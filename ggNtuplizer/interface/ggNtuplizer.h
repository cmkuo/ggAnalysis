// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/Common/interface/TriggerResults.h"
//#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"

#include "CondFormats/EcalObjects/interface/EcalIntercalibConstants.h"
#include "CondFormats/EcalObjects/interface/EcalADCToGeVConstant.h"
#include "CalibCalorimetry/EcalLaserCorrection/interface/EcalLaserDbService.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgoRcd.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"

#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateFwd.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDetId/interface/ESDetId.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "Geometry/Records/interface/CaloTopologyRecord.h"
#include "CondFormats/DataRecord/interface/EcalIntercalibConstantsRcd.h"
#include "CondFormats/DataRecord/interface/EcalADCToGeVConstantRcd.h"
#include "CalibCalorimetry/EcalLaserCorrection/interface/EcalLaserDbRecord.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgoRcd.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "CondFormats/EcalObjects/interface/EcalPFRecHitThresholds.h"
#include "CondFormats/DataRecord/interface/EcalPFRecHitThresholdsRcd.h"
#include "RecoEcal/EgammaCoreTools/interface/EgammaLocalCovParamDefaults.h"

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"           
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "TTree.h"
#include <cmath>
#include "TVector3.h"
#include "TLorentzVector.h"
#include "iterator"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

using namespace std;
using namespace edm;
using reco::TrackCollection;

void setbit(UShort_t& x, UShort_t bit);

class ggNtuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
 public:
  explicit ggNtuplizer(const edm::ParameterSet&);
  ~ggNtuplizer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  TTree   *tree; 

  bool doGenParticles_;
  bool dumpPDFSystWeight_;
  bool dumpGenScaleSystWeights_;
  bool runOnParticleGun_;
  bool runOnSherpa_;
  bool dumpCrystalinfo_;

  vector<int> newparticles_;

 private:
  //void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&); // override;
  //Void- Endjob() override;

  void branchesGlobalEvent(TTree*);
  void branchesGenInfo    (TTree*, edm::Service<TFileService>&);
  void branchesGenPart    (TTree*);
  void branchesPhotons    (TTree*);
  void branchesElectrons  (TTree*);
  void branchesRechits    (TTree*);

  void fillGlobalEvent(const edm::Event&, const edm::EventSetup&);
  void fillGenInfo    (const edm::Event&);
  void fillGenPart    (const edm::Event&);
  void fillPhotons    (const edm::Event&, const edm::EventSetup&);
  void fillElectrons  (const edm::Event&, const edm::EventSetup&, math::XYZPoint&);
  void fillRechits    (const edm::Event&, const edm::EventSetup&);

  // ----------member data ---------------------------
  //edm::EDGetTokenT<TrackCollection> tracksToken_;  //used to select what tracks to read from configuration file
  //edm::EDGetTokenT<edm::View<reco::GsfElectron> > electronCollection_;

  edm::EDGetTokenT<reco::VertexCollection>         vtxLabel_;
  //edm::EDGetTokenT<reco::VertexCollection>         vtxBSLabel_;
  edm::EDGetTokenT<GenEventInfoProduct>            generatorLabel_;
  edm::EDGetTokenT<LHEEventProduct>                lheEventLabel_;
  edm::EDGetTokenT<vector<PileupSummaryInfo> >     puCollection_;
  edm::EDGetTokenT<vector<reco::GenParticle> >     genParticlesCollection_;

  edm::EDGetTokenT<double>                         rhoLabel_;
  edm::EDGetTokenT<double>                         rhoCentralLabel_;
  edm::EDGetTokenT<trigger::TriggerEvent>          trgEventLabel_;
  edm::EDGetTokenT<pat::TriggerObjectStandAloneCollection> triggerObjectsLabel_;
  edm::EDGetTokenT<edm::TriggerResults>            trgResultsLabel_;
  string                                           trgResultsProcess_;
  edm::EDGetTokenT<edm::TriggerResults>            patTrgResultsLabel_;
  //HLTPrescaleProvider                              hltPrescaleProvider_;

  edm::EDGetTokenT<edm::View<pat::Electron> >      electronCollection_;
  edm::EDGetTokenT<edm::View<pat::Photon> >        photonCollection_;
  edm::EDGetToken barrelClusterToken_;
  edm::EDGetToken endcapClusterToken_;

  edm::ESHandle<CaloGeometry>                      pG_;
  edm::EDGetTokenT<EcalRecHitCollection>           ebReducedRecHitCollection_;
  edm::EDGetTokenT<EcalRecHitCollection>           eeReducedRecHitCollection_;
  edm::EDGetTokenT<EcalRecHitCollection>           esReducedRecHitCollection_;
  const EcalClusterLazyTools::ESGetTokens ecalClusterToolsESGetTokens_;

    
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  edm::ESGetToken<SetupData, SetupRecord> setupToken_;
#endif
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
