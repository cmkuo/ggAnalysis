#include "ggAnalysis/ggNtuplizer/interface/ggNtuplizer.h"

using namespace std;
using namespace edm;

void setbit(UShort_t& x, UShort_t bit) {
  UShort_t a = 1;
  x |= (a << bit);
}

ggNtuplizer::ggNtuplizer(const edm::ParameterSet& ps) {

  development_               = ps.getParameter<bool>("development");
  addFilterInfoAOD_          = ps.getParameter<bool>("addFilterInfoAOD");
  addFilterInfoMINIAOD_      = ps.getParameter<bool>("addFilterInfoMINIAOD");
  doNoHFMET_                 = ps.getParameter<bool>("doNoHFMET");

  doGenParticles_            = ps.getParameter<bool>("doGenParticles");
  runOnParticleGun_          = ps.getParameter<bool>("runOnParticleGun");
  dumpPhotons_               = ps.getParameter<bool>("dumpPhotons");
  dumpJets_                  = ps.getParameter<bool>("dumpJets");
  dumpSubJets_               = ps.getParameter<bool>("dumpSubJets");
  dumpSoftDrop_               = ps.getParameter<bool>("dumpSoftDrop");
  dumpTaus_                  = ps.getParameter<bool>("dumpTaus");
  isAOD_                     = ps.getParameter<bool>("isAOD");

  runphoIDVID_               = ps.getParameter<bool>("runphoIDVID");
  runeleIDVID_               = ps.getParameter<bool>("runeleIDVID");
  runeleMVAID_               = ps.getParameter<bool>("runeleMVAID");
  runphoMVAID_               = ps.getParameter<bool>("runphoMVAID");

  runHFElectrons_            = ps.getParameter<bool>("runHFElectrons");

  trgFilterDeltaPtCut_       = ps.getParameter<double>("trgFilterDeltaPtCut");
  trgFilterDeltaRCut_        = ps.getParameter<double>("trgFilterDeltaRCut");

  vtxLabel_                  = consumes<reco::VertexCollection>     (ps.getParameter<InputTag>("VtxLabel"));
  vtxBSLabel_                = consumes<reco::VertexCollection>     (ps.getParameter<InputTag>("VtxBSLabel"));
  rhoLabel_                  = consumes<double>                     (ps.getParameter<InputTag>("rhoLabel"));
  rhoCentralLabel_           = consumes<double>                     (ps.getParameter<InputTag>("rhoCentralLabel"));
  trgEventLabel_             = consumes<trigger::TriggerEvent>      (ps.getParameter<InputTag>("triggerEvent"));
  triggerObjectsLabel_       = consumes<pat::TriggerObjectStandAloneCollection>(ps.getParameter<edm::InputTag>("triggerEvent"));
  trgResultsLabel_           = consumes<edm::TriggerResults>        (ps.getParameter<InputTag>("triggerResults"));
  patTrgResultsLabel_        = consumes<edm::TriggerResults>        (ps.getParameter<InputTag>("patTriggerResults"));
  trgResultsProcess_         =                                       ps.getParameter<InputTag>("triggerResults").process();
  generatorLabel_            = consumes<GenEventInfoProduct>        (ps.getParameter<InputTag>("generatorLabel"));
  puCollection_              = consumes<vector<PileupSummaryInfo> > (ps.getParameter<InputTag>("pileupCollection"));
  genParticlesCollection_    = consumes<vector<reco::GenParticle> > (ps.getParameter<InputTag>("genParticleSrc"));
  pfMETlabel_                = consumes<View<pat::MET> >            (ps.getParameter<InputTag>("pfMETLabel"));
  electronCollection_        = consumes<View<pat::Electron> >       (ps.getParameter<InputTag>("electronSrc"));
  gsfTracks_                 = consumes<View<reco::GsfTrack>>       (ps.getParameter<InputTag>("gsfTrackSrc"));

  photonCollection_          = consumes<View<pat::Photon> >         (ps.getParameter<InputTag>("photonSrc"));
  muonCollection_            = consumes<View<pat::Muon> >           (ps.getParameter<InputTag>("muonSrc"));
  ebReducedRecHitCollection_ = consumes<EcalRecHitCollection>       (ps.getParameter<InputTag>("ebReducedRecHitCollection"));
  eeReducedRecHitCollection_ = consumes<EcalRecHitCollection>       (ps.getParameter<InputTag>("eeReducedRecHitCollection"));
  esReducedRecHitCollection_ = consumes<EcalRecHitCollection>       (ps.getParameter<InputTag>("esReducedRecHitCollection")); 
  recophotonCollection_      = consumes<reco::PhotonCollection>     (ps.getParameter<InputTag>("recoPhotonSrc"));
  tracklabel_                = consumes<reco::TrackCollection>      (ps.getParameter<InputTag>("TrackLabel"));
  gsfElectronlabel_          = consumes<reco::GsfElectronCollection>(ps.getParameter<InputTag>("gsfElectronLabel"));
  tauCollection_             = consumes<vector<pat::Tau> >          (ps.getParameter<InputTag>("tauSrc"));
  pfAllParticles_            = consumes<reco::PFCandidateCollection>(ps.getParameter<InputTag>("PFAllCandidates"));
  jetsAK4Label_ = consumes<View<pat::Jet> > (ps.getParameter<InputTag>("ak4JetSrc"));
  jetsAK8Label_ = consumes<View<pat::Jet> > (ps.getParameter<InputTag>("ak8JetSrc"));
  newparticles_              = ps.getParameter< vector<int > >("newParticles");

  //pfLooseId_                 = ps.getParameter<ParameterSet>("pfLooseId");
  pckPFCdsLabel_ = ps.getParameter<edm::InputTag>("packedPFCands");

  cicPhotonId_ = new CiCPhotonID(ps);

  // electron ID 
  eleVetoIdMapToken_    = consumes<edm::ValueMap<bool> >(ps.getParameter<edm::InputTag>("eleVetoIdMap"));
  eleLooseIdMapToken_   = consumes<edm::ValueMap<bool> >(ps.getParameter<edm::InputTag>("eleLooseIdMap"));
  eleMediumIdMapToken_  = consumes<edm::ValueMap<bool> >(ps.getParameter<edm::InputTag>("eleMediumIdMap"));
  eleTightIdMapToken_   = consumes<edm::ValueMap<bool> >(ps.getParameter<edm::InputTag>("eleTightIdMap"));
  eleHEEPIdMapToken_    = consumes<edm::ValueMap<bool> >(ps.getParameter<edm::InputTag>("eleHEEPIdMap"));
  eleNonTrgMVAValuesMapToken_ = consumes<edm::ValueMap<float> >(ps.getParameter<edm::InputTag>("eleNonTrgMVAValuesMap"));
  eleTrgMVAValuesMapToken_    = consumes<edm::ValueMap<float> >(ps.getParameter<edm::InputTag>("eleTrgMVAValuesMap"));

  // Photon ID in VID framwork 
  phoLooseIdMapToken_             = consumes<edm::ValueMap<bool> >(ps.getParameter<edm::InputTag>("phoLooseIdMap"));
  phoMediumIdMapToken_            = consumes<edm::ValueMap<bool> >(ps.getParameter<edm::InputTag>("phoMediumIdMap"));
  phoTightIdMapToken_             = consumes<edm::ValueMap<bool> >(ps.getParameter<edm::InputTag>("phoTightIdMap"));
  phoMVAValuesMapToken_           = consumes<edm::ValueMap<float> >(ps.getParameter<edm::InputTag>("phoMVAValuesMap")); 
  phoChargedIsolationToken_       = consumes <edm::ValueMap<float> >(ps.getParameter<edm::InputTag>("phoChargedIsolation"));
  phoNeutralHadronIsolationToken_ = consumes <edm::ValueMap<float> >(ps.getParameter<edm::InputTag>("phoNeutralHadronIsolation"));
  phoPhotonIsolationToken_        = consumes <edm::ValueMap<float> >(ps.getParameter<edm::InputTag>("phoPhotonIsolation"));
  phoWorstChargedIsolationToken_  = consumes <edm::ValueMap<float> >(ps.getParameter<edm::InputTag>("phoWorstChargedIsolation"));


  Service<TFileService> fs;
  tree_    = fs->make<TTree>("EventTree", "Event data");
  hEvents_ = fs->make<TH1F>("hEvents",    "total processed and skimmed events",   2,  0,   2);

  branchesGlobalEvent(tree_);

  if (doGenParticles_) {
    branchesGenInfo(tree_, fs);
    branchesGenPart(tree_);
  }

  branchesMET(tree_);
  if (dumpPhotons_) branchesPhotons(tree_);
  branchesElectrons(tree_);
  if (isAOD_ && runHFElectrons_) branchesHFElectrons(tree_);
  branchesMuons(tree_);
  if (dumpTaus_) branchesTaus(tree_);
  if (dumpJets_) branchesJets(tree_);
}

ggNtuplizer::~ggNtuplizer() {
  cleanupPhotons();
  delete cicPhotonId_;
}

void ggNtuplizer::analyze(const edm::Event& e, const edm::EventSetup& es) {

  hEvents_->Fill(0.5);

  edm::Handle<reco::VertexCollection> vtxHandle;
  e.getByToken(vtxLabel_, vtxHandle);

  reco::Vertex vtx;

  // best-known primary vertex coordinates
  math::XYZPoint pv(0, 0, 0);
  for (vector<reco::Vertex>::const_iterator v = vtxHandle->begin(); v != vtxHandle->end(); ++v) {
    // replace isFake() for miniAOD since it requires tracks while miniAOD vertices don't have tracks:
    // Vertex.h: bool isFake() const {return (chi2_==0 && ndof_==0 && tracks_.empty());}
    bool isFake = isAOD_ ? v->isFake() : (v->chi2() == 0 && v->ndof() == 0);

    if (!isFake) {
      pv.SetXYZ(v->x(), v->y(), v->z());
      vtx = *v;
      break;
    }
  }

  initTriggerFilters(e);

  fillGlobalEvent(e, es);
  if (!e.isRealData()) {
    fillGenInfo(e);
    if (doGenParticles_)
      fillGenPart(e);
  }

  fillMET(e, es);
  fillPhotons(e, es); // FIXME: photons have different vertex (not pv)
  fillElectrons(e, es, pv);
  if (isAOD_ && runHFElectrons_ ) fillHFElectrons(e);
  fillMuons(e, pv, vtx);
  if (dumpTaus_) fillTaus(e);
  if (dumpJets_) fillJets(e,es);

  hEvents_->Fill(1.5);
  tree_->Fill();

}


// void ggNtuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions)
// {
//   //The following says we do not know what parameters are allowed so do no validation
//   // Please change this to state exactly what you do use, even if it is no parameters
//   edm::ParameterSetDescription desc;
//   desc.setUnknown();
//   descriptions.addDefault(desc);
// }

DEFINE_FWK_MODULE(ggNtuplizer);
