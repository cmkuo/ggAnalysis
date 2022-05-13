// -*- C++ -*-
//
// Package:    ggAnalysis/ggNtuplizer
// Class:      ggNtuplizer
//
/**\class ggNtuplizer ggNtuplizer.cc ggAnalysis/ggNtuplizer/plugins/ggNtuplizer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Prasant Kumar Rout
//         Created:  Thu, 12 May 2022 19:26:52 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "../interface/ggNtuplizer.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

using namespace std;
using namespace edm;

void setbit(UShort_t& x, UShort_t bit) {
  UShort_t a = 1;
  x |= (a << bit);
}


ggNtuplizer::ggNtuplizer(const edm::ParameterSet& iConfig) :
  ecalClusterToolsESGetTokens_{consumesCollector()}
{
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  setupDataToken_ = esConsumes<SetupData, SetupRecord>();
#endif

  doGenParticles_            =                                          iConfig.getParameter<bool>("doGenParticles");
  runOnParticleGun_          =                                          iConfig.getParameter<bool>("runOnParticleGun");
  runOnSherpa_               =                                          iConfig.getParameter<bool>("runOnSherpa");
  dumpCrystalinfo_           =                                          iConfig.getParameter<bool>("dumpCrystalinfo"); 
  dumpPDFSystWeight_         =                                          iConfig.getParameter<bool>("dumpPDFSystWeight");
  dumpGenScaleSystWeights_   =                                          iConfig.getParameter<bool>("dumpGenScaleSystWeights");
  newparticles_              =                                          iConfig.getParameter< vector<int > >("newParticles");

  vtxLabel_                  = consumes<reco::VertexCollection>        (iConfig.getParameter<InputTag>("VtxLabel"));
  //vtxBSLabel_                = consumes<reco::VertexCollection>        (iConfig.getParameter<InputTag>("VtxBSLabel"));                     
  
  generatorLabel_            = consumes<GenEventInfoProduct>           (iConfig.getParameter<InputTag>("generatorLabel"));
  lheEventLabel_             = consumes<LHEEventProduct>               (iConfig.getParameter<InputTag>("LHEEventLabel"));
  puCollection_              = consumes<vector<PileupSummaryInfo> >    (iConfig.getParameter<InputTag>("pileupCollection"));
  genParticlesCollection_    = consumes<vector<reco::GenParticle> >    (iConfig.getParameter<InputTag>("genParticleSrc"));
  
  rhoLabel_                  = consumes<double>                        (iConfig.getParameter<InputTag>("rhoLabel"));
  rhoCentralLabel_           = consumes<double>                        (iConfig.getParameter<InputTag>("rhoCentralLabel"));
  //trgEventLabel_             = consumes<trigger::TriggerEvent>         (iConfig.getParameter<InputTag>("triggerEvent"));                    
  //triggerObjectsLabel_       = consumes<pat::TriggerObjectStandAloneCollection>(iConfig.getParameter<edm::InputTag>("triggerEvent"));      
  
  trgResultsLabel_           = consumes<edm::TriggerResults>           (iConfig.getParameter<InputTag>("triggerResults"));
  //patTrgResultsLabel_        = consumes<edm::TriggerResults>           (iConfig.getParameter<InputTag>("patTriggerResults"));              
  
  trgResultsProcess_         =                                          iConfig.getParameter<InputTag>("triggerResults").process();

  electronCollection_         = consumes<View<pat::Electron> > (iConfig.getParameter<InputTag>("electronSrc"));
  photonCollection_           = consumes<View<pat::Photon> >   (iConfig.getParameter<InputTag>("photonSrc"));
  barrelClusterToken_         = consumes<reco::BasicClusterCollection> (iConfig.getParameter<InputTag>("barrelClusterCollection"));
  endcapClusterToken_         = consumes<reco::BasicClusterCollection> (iConfig.getParameter<InputTag>("endcapClusterCollection"));
  ebReducedRecHitCollection_  = consumes<EcalRecHitCollection> (iConfig.getParameter<InputTag>("ebReducedRecHitCollection"));
  eeReducedRecHitCollection_  = consumes<EcalRecHitCollection> (iConfig.getParameter<InputTag>("eeReducedRecHitCollection"));
  esReducedRecHitCollection_  = consumes<EcalRecHitCollection> (iConfig.getParameter<InputTag>("esReducedRecHitCollection"));

  edm::Service<TFileService> fs;
  tree = fs->make<TTree>("EventTree", "EventInfo");

  branchesGlobalEvent(tree);
  if (doGenParticles_) {
    branchesGenInfo(tree, fs);
    branchesGenPart(tree);
  }

  branchesPhotons(tree);
  branchesElectrons(tree);
  branchesRechits(tree);

  //now do what ever initialization is needed
}

ggNtuplizer::~ggNtuplizer() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called for each event  ------------
void ggNtuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  //for (const auto& track : iEvent.get(tracksToken_)) {
    // do something with track parameters, e.g, plot the charge.
    // int charge = track.charge();

  fillGlobalEvent(iEvent, iSetup);
  if (!iEvent.isRealData()) {
    fillGenInfo(iEvent);
    if (doGenParticles_)
      fillGenPart(iEvent);
  }

  edm::Handle<reco::VertexCollection> vtxHandle;
  iEvent.getByToken(vtxLabel_, vtxHandle);

  reco::Vertex vtx;

  // best-known primary vertex coordinates                                                                                                   
  
  math::XYZPoint pv(0, 0, 0);
  for (vector<reco::Vertex>::const_iterator v = vtxHandle->begin(); v != vtxHandle->end(); ++v) {
    // replace isFake() for miniAOD since it requires tracks while miniAOD vertices don't have tracks:      
                                   
    // Vertex.h: bool isFake() const {return (chi2_==0 && ndof_==0 && tracks_.empty());}                                                     
  
    bool isFake = (v->chi2() == 0 && v->ndof() == 0);

    if (!isFake) {
      pv.SetXYZ(v->x(), v->y(), v->z());
      vtx = *v;
      break;
    }
  }


  fillPhotons(iEvent, iSetup);
  fillElectrons(iEvent, iSetup, pv);
  fillRechits(iEvent, iSetup);

  
  tree->Fill();


#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  // if the SetupData is always needed
  auto setup = iSetup.getData(setupToken_);
  // if need the ESHandle to check if the SetupData was there or not
  auto pSetup = iSetup.getHandle(setupToken_);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
/*void ggNtuplizer::beginJob() {
  // please remove this method if not needed
}

// ------------ method called once each job just after ending the event loop  ------------
void ggNtuplizer::endJob() {
  // please remove this method if not needed
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void ggNtuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(ggNtuplizer);
