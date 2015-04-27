#include "ggAnalysis/ggNtuplizer/interface/ggNtuplizer.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "PhysicsTools/SelectorUtils/interface/JetIDSelectionFunctor.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"

#include <TMVA/Reader.h>

ggNtuplizer::ggNtuplizer(const edm::ParameterSet& ps) {

  doGenParticles_            = ps.getParameter<bool>("doGenParticles");
  runOnParticleGun_          = ps.getParameter<bool>("runOnParticleGun");
  vtxLabel_                  = ps.getParameter<InputTag>("VtxLabel");
  vtxBSLabel_                = ps.getParameter<InputTag>("VtxBSLabel");
  rhoLabel_                  = ps.getParameter<InputTag>("rhoLabel");
  generatorLabel_            = ps.getParameter<InputTag>("generatorLabel");
  puCollection_              = ps.getParameter<InputTag>("pileupCollection");
  genParticlesCollection_ = consumes<vector<reco::GenParticle> >   (ps.getParameter<edm::InputTag>("genParticleSrc"));
  pfMETlabel_                = ps.getParameter<InputTag>("pfMETLabel");
  electronCollection_        = ps.getParameter<InputTag>("electronSrc");
  photonCollection_          = ps.getParameter<InputTag>("photonSrc");
  muonCollection_            = ps.getParameter<InputTag>("muonSrc");
  ebReducedRecHitCollection_ = consumes<EcalRecHitCollection>(ps.getParameter<InputTag>("ebReducedRecHitCollection"));
  eeReducedRecHitCollection_ = consumes<EcalRecHitCollection>(ps.getParameter<InputTag>("eeReducedRecHitCollection"));
  esReducedRecHitCollection_ = consumes<EcalRecHitCollection>(ps.getParameter<InputTag>("esReducedRecHitCollection"));
  recophotonCollection_      = ps.getParameter<InputTag>("recoPhotonSrc");
  tracklabel_                = ps.getParameter<InputTag>("TrackLabel");
  gsfElectronlabel_          = ps.getParameter<InputTag>("gsfElectronLabel");
  pfAllParticles_            = ps.getParameter<InputTag>("PFAllCandidates");
  jetsCHSLabel_              = edm::InputTag("selectedPatJetsCA8PFCHS");
  jetCollection_             = ps.getParameter<InputTag>("jetSrc");
  dumpJets_                  = ps.getParameter<bool>("dumpJets");
  dumpSubJets_               = ps.getParameter<bool>("dumpSubJets");
  dumpTaus_                  = ps.getParameter<bool>("dumpTaus");
  pfLooseId_                 = ps.getParameter<edm::ParameterSet>("pfLooseId");
  tauCollection_             = ps.getParameter<InputTag>("tauSrc");

  cicPhotonId_ = new CiCPhotonID(ps);

  Service<TFileService> fs;
  hEvents_ = fs->make<TH1F>("hEvents",    "total processed and skimmed events",   2,  0,   2);
  hPU_     = fs->make<TH1F>("hPU",        "number of pileup",                   200,  0, 200);
  hPUTrue_ = fs->make<TH1F>("hPUTrue",    "number of true pilepu"             , 1000, 0, 200);
  tree_    = fs->make<TTree>("EventTree", "Event data");

  tree_->Branch("run", &run_, "run/I");
  tree_->Branch("event", &event_, "event/L");
  tree_->Branch("lumis", &lumis_, "lumis/I");
  tree_->Branch("isData", &isData_, "isData/O");
  tree_->Branch("nVtx", &nVtx_, "nVtx/I");
  tree_->Branch("nTrks", &nTrks_, "nTrks/I");
  tree_->Branch("rho", &rho_, "rho/F");
  if (doGenParticles_) {
    tree_->Branch("pdf", &pdf_);
    tree_->Branch("pthat", &pthat_, "pthat/F");
    tree_->Branch("processID", &processID_, "processID/F");

    makeBranchesGenParticles(tree_);

    // Gen MET
    //tree_->Branch("genMET", &genMET_, "genMET/F");
    //tree_->Branch("genMETPhi", &genMETPhi_, "genMETPhi/F");
    // PU Info
    tree_->Branch("nPUInfo", &nPUInfo_, "nPUInfo/I");
    tree_->Branch("nPU", &nPU_);
    tree_->Branch("puBX", &puBX_);
    tree_->Branch("puTrue", &puTrue_);
  }
  // pfMET Type1  
  tree_->Branch("pfMET",       &pfMET_,       "pfMET/F");
  tree_->Branch("pfMETPhi",    &pfMETPhi_,    "pfMETPhi/F");
  tree_->Branch("pfMETsumEt",  &pfMETsumEt_,  "pfMETsumEt/F");
  tree_->Branch("pfMETmEtSig", &pfMETmEtSig_, "pfMETmEtSig/F");
  tree_->Branch("pfMETSig",    &pfMETSig_,    "pfMETSig/F");
  // Electron
  tree_->Branch("nEle",                  &nEle_, "nEle/I");
  tree_->Branch("eleCharge",             &eleCharge_);
  tree_->Branch("eleChargeConsistent",   &eleChargeConsistent_);
  tree_->Branch("eleEn",                 &eleEn_);
  tree_->Branch("eleSCEn",               &eleSCEn_);
  tree_->Branch("eleESEn",               &eleESEn_);
  tree_->Branch("eleD0",                 &eleD0_);
  tree_->Branch("eleDz",                 &eleDz_);
  tree_->Branch("elePt",                 &elePt_);
  tree_->Branch("eleEta",                &eleEta_);

  tree_->Branch("eleTheta",                &eleTheta_);

  tree_->Branch("elePhi",                &elePhi_);
  tree_->Branch("eleSCEta",              &eleSCEta_);
  tree_->Branch("eleSCPhi",              &eleSCPhi_);
  tree_->Branch("eleSCRawEn",            &eleSCRawEn_);
  tree_->Branch("eleSCEtaWidth",         &eleSCEtaWidth_);
  tree_->Branch("eleSCPhiWidth",         &eleSCPhiWidth_);
  tree_->Branch("eleHoverE",             &eleHoverE_);
  tree_->Branch("eleEoverP",             &eleEoverP_);
  tree_->Branch("eleEoverPInv",          &eleEoverPInv_);
  tree_->Branch("eleBrem",               &eleBrem_);
  tree_->Branch("eledEtaAtVtx",          &eledEtaAtVtx_);
  tree_->Branch("eledPhiAtVtx",          &eledPhiAtVtx_);
  tree_->Branch("eleSigmaIEtaIEta",      &eleSigmaIEtaIEta_);
  tree_->Branch("eleSigmaIEtaIPhi",      &eleSigmaIEtaIPhi_);
  tree_->Branch("eleSigmaIPhiIPhi",      &eleSigmaIPhiIPhi_);
  tree_->Branch("eleSigmaIEtaIEta_2012", &eleSigmaIEtaIEta_2012_);
  tree_->Branch("eleConvVeto",           &eleConvVeto_);
  tree_->Branch("eleMissHits",           &eleMissHits_);
  tree_->Branch("eleESEffSigmaRR",       &eleESEffSigmaRR_);
  tree_->Branch("elePFChIso",            &elePFChIso_);
  tree_->Branch("elePFPhoIso",           &elePFPhoIso_);
  tree_->Branch("elePFNeuIso",           &elePFNeuIso_);
  tree_->Branch("elePFPUIso",            &elePFPUIso_);
  tree_->Branch("eleBC1E",               &eleBC1E_);
  tree_->Branch("eleBC1Eta",             &eleBC1Eta_);
  tree_->Branch("eleBC2E",               &eleBC2E_);
  tree_->Branch("eleBC2Eta",             &eleBC2Eta_);
  tree_->Branch("eleIDMVA",              &eleIDMVA_);

  tree_->Branch("eledEtaseedAtVtx",           &eledEtaseedAtVtx_);

  tree_->Branch("eleE1x5",           &eleE1x5_);
  tree_->Branch("eleE2x5",           &eleE2x5_);
  tree_->Branch("eleE5x5",           &eleE5x5_);
  tree_->Branch("eleRelIsoWithDBeta",            &eleRelIsoWithDBeta_);

  tree_->Branch("eleE1x5_2012",           &eleE1x5_2012_);
  tree_->Branch("eleE2x5_2012",           &eleE2x5_2012_);
  tree_->Branch("eleE5x5_2012",           &eleE5x5_2012_);

  tree_->Branch("eleEcalDrivenSeed",           &eleEcalDrivenSeed_);
  tree_->Branch("eleDr03EcalRecHitSumEt",           &eleDr03EcalRecHitSumEt_);
  tree_->Branch("eleDr03HcalDepth1TowerSumEt",           &eleDr03HcalDepth1TowerSumEt_);
  tree_->Branch("eleDr03HcalDepth2TowerSumEt",           &eleDr03HcalDepth2TowerSumEt_);
  tree_->Branch("eleDr03HcalTowerSumEt",           &eleDr03HcalTowerSumEt_);
  tree_->Branch("eleDr03TkSumPt",           &eleDr03TkSumPt_);  
  tree_->Branch("elecaloEnergy",           &elecaloEnergy_);

  tree_->Branch("eleTrkdxy",           &eleTrkdxy_);
 
  // Photon
  tree_->Branch("nPho",                  &nPho_, "nPho/I");
  tree_->Branch("phoE",                  &phoE_);
  tree_->Branch("phoEt",                 &phoEt_);
  tree_->Branch("phoEta",                &phoEta_);
  tree_->Branch("phoPhi",                &phoPhi_);
  tree_->Branch("phoSCE",                &phoSCE_);
  tree_->Branch("phoSCRawE",             &phoSCRawE_);
  tree_->Branch("phoESEn",               &phoESEn_);
  tree_->Branch("phoSCEta",              &phoSCEta_);
  tree_->Branch("phoSCPhi",              &phoSCPhi_);
  tree_->Branch("phoSCEtaWidth",         &phoSCEtaWidth_);
  tree_->Branch("phoSCPhiWidth",         &phoSCPhiWidth_);
  tree_->Branch("phoSCBrem",             &phoSCBrem_);
  tree_->Branch("phohasPixelSeed",       &phohasPixelSeed_);
  tree_->Branch("phoEleVeto",            &phoEleVeto_);
  tree_->Branch("phoR9",                 &phoR9_);
  tree_->Branch("phoHoverE",             &phoHoverE_);
  tree_->Branch("phoSigmaIEtaIEta",      &phoSigmaIEtaIEta_);
  tree_->Branch("phoSigmaIEtaIPhi",      &phoSigmaIEtaIPhi_);
  tree_->Branch("phoSigmaIPhiIPhi",      &phoSigmaIPhiIPhi_);
  tree_->Branch("phoE1x3",               &phoE1x3_);
  tree_->Branch("phoE2x2",               &phoE2x2_);
  tree_->Branch("phoE2x5Max",            &phoE2x5Max_);
  tree_->Branch("phoE5x5",               &phoE5x5_);
  tree_->Branch("phoESEffSigmaRR",       &phoESEffSigmaRR_);
  tree_->Branch("phoSigmaIEtaIEta_2012", &phoSigmaIEtaIEta_2012_);
  tree_->Branch("phoSigmaIEtaIPhi_2012", &phoSigmaIEtaIPhi_2012_);
  tree_->Branch("phoSigmaIPhiIPhi_2012", &phoSigmaIPhiIPhi_2012_);
  tree_->Branch("phoE1x3_2012",          &phoE1x3_2012_);
  tree_->Branch("phoE2x2_2012",          &phoE2x2_2012_);
  tree_->Branch("phoE2x5Max_2012",       &phoE2x5Max_2012_);
  tree_->Branch("phoE5x5_2012",          &phoE5x5_2012_);
  tree_->Branch("phoPFChIso",            &phoPFChIso_);
  tree_->Branch("phoPFPhoIso",           &phoPFPhoIso_);
  tree_->Branch("phoPFNeuIso",           &phoPFNeuIso_);
  tree_->Branch("phoPFChWorstIso",       &phoPFChWorstIso_);
  tree_->Branch("phoPFChIsoFrix1",       &phoPFChIsoFrix1_);
  tree_->Branch("phoPFChIsoFrix2",       &phoPFChIsoFrix2_);
  tree_->Branch("phoPFChIsoFrix3",       &phoPFChIsoFrix3_);
  tree_->Branch("phoPFChIsoFrix4",       &phoPFChIsoFrix4_);
  tree_->Branch("phoPFChIsoFrix5",       &phoPFChIsoFrix5_);
  tree_->Branch("phoPFChIsoFrix6",       &phoPFChIsoFrix6_);
  tree_->Branch("phoPFChIsoFrix7",       &phoPFChIsoFrix7_);
  tree_->Branch("phoPFChIsoFrix8",       &phoPFChIsoFrix8_);
  tree_->Branch("phoPFPhoIsoFrix1",      &phoPFPhoIsoFrix1_);
  tree_->Branch("phoPFPhoIsoFrix2",      &phoPFPhoIsoFrix2_);
  tree_->Branch("phoPFPhoIsoFrix3",      &phoPFPhoIsoFrix3_);
  tree_->Branch("phoPFPhoIsoFrix4",      &phoPFPhoIsoFrix4_);
  tree_->Branch("phoPFPhoIsoFrix5",      &phoPFPhoIsoFrix5_);
  tree_->Branch("phoPFPhoIsoFrix6",      &phoPFPhoIsoFrix6_);
  tree_->Branch("phoPFPhoIsoFrix7",      &phoPFPhoIsoFrix7_);
  tree_->Branch("phoPFPhoIsoFrix8",      &phoPFPhoIsoFrix8_);
  tree_->Branch("phoPFNeuIsoFrix1",      &phoPFNeuIsoFrix1_);
  tree_->Branch("phoPFNeuIsoFrix2",      &phoPFNeuIsoFrix2_);
  tree_->Branch("phoPFNeuIsoFrix3",      &phoPFNeuIsoFrix3_);
  tree_->Branch("phoPFNeuIsoFrix4",      &phoPFNeuIsoFrix4_);
  tree_->Branch("phoPFNeuIsoFrix5",      &phoPFNeuIsoFrix5_);
  tree_->Branch("phoPFNeuIsoFrix6",      &phoPFNeuIsoFrix6_);
  tree_->Branch("phoPFNeuIsoFrix7",      &phoPFNeuIsoFrix7_);
  tree_->Branch("phoPFNeuIsoFrix8",      &phoPFNeuIsoFrix8_);
  tree_->Branch("phoBC1E",               &phoBC1E_);
  tree_->Branch("phoBC1Eta",             &phoBC1Eta_);
  tree_->Branch("phoBC2E",               &phoBC2E_);
  tree_->Branch("phoBC2Eta",             &phoBC2Eta_);

  tree_->Branch("phoIDMVA",                  &phoIDMVA_);
  tree_->Branch("phoEcalRecHitSumEtConeDR03",                  &phoEcalRecHitSumEtConeDR03_);
  tree_->Branch("phohcalDepth1TowerSumEtConeDR03",                  &phohcalDepth1TowerSumEtConeDR03_);
  tree_->Branch("phohcalDepth2TowerSumEtConeDR03",                  &phohcalDepth2TowerSumEtConeDR03_);
  tree_->Branch("phohcalTowerSumEtConeDR03",                  &phohcalTowerSumEtConeDR03_);

  tree_->Branch("photrkSumPtHollowConeDR03",                  &photrkSumPtHollowConeDR03_);

  // muon
  tree_->Branch("nMu",                   &nMu_, "nMu/I");
  tree_->Branch("muPt",                  &muPt_);
  tree_->Branch("muEta",                 &muEta_);
  tree_->Branch("muPhi",                 &muPhi_);
  tree_->Branch("muCharge",              &muCharge_);
  tree_->Branch("muType",                &muType_);
  tree_->Branch("muIsGood",              &muIsGood_);
  //tree_->Branch("muID",                  &muID_);
  tree_->Branch("muD0",                  &muD0_);
  tree_->Branch("muDz",                  &muDz_);
  tree_->Branch("muChi2NDF",             &muChi2NDF_);
  tree_->Branch("muInnerD0",             &muInnerD0_);
  tree_->Branch("muInnerDz",             &muInnerDz_);
  tree_->Branch("muTrkLayers",           &muTrkLayers_);
  tree_->Branch("muPixelLayers",         &muPixelLayers_);
  tree_->Branch("muPixelHits",           &muPixelHits_);
  tree_->Branch("muMuonHits",            &muMuonHits_);
  tree_->Branch("muStations",            &muStations_);
  tree_->Branch("muTrkQuality",          &muTrkQuality_);
  tree_->Branch("muIsoTrk",              &muIsoTrk_);
  tree_->Branch("muPFChIso",             &muPFChIso_);
  tree_->Branch("muPFPhoIso",            &muPFPhoIso_);
  tree_->Branch("muPFNeuIso",            &muPFNeuIso_);
  tree_->Branch("muPFPUIso",             &muPFPUIso_);

  ///SJ
  tree_->Branch("muInnervalidFraction",         &muInnervalidFraction_);
  tree_->Branch("musegmentCompatibility",                  &musegmentCompatibility_);
  tree_->Branch("muchi2LocalPosition",                  &muchi2LocalPosition_);
  tree_->Branch("mutrkKink",                  &mutrkKink_);
  tree_->Branch("muBestTrkPtError",                  &muBestTrkPtError_);
  tree_->Branch("muBestTrkPt",                  &muBestTrkPt_);

  //tau : Lvdp 
  if(dumpTaus_){
    tree_->Branch("nTau", &nTau_,"nTau/I");
    tree_->Branch("tauByLooseElectronRejection", &tauByLooseElectronRejection_);
    tree_->Branch("tauByMediumElectronRejection", &tauByMediumElectronRejection_);
    tree_->Branch("tauByTightElectronRejection", &tauByTightElectronRejection_);
    tree_->Branch("tauByMVA5LooseElectronRejection", &tauByMVA5LooseElectronRejection_);
    tree_->Branch("tauByMVA5MediumElectronRejection", &tauByMVA5MediumElectronRejection_);
    tree_->Branch("tauByMVA5TightElectronRejection", &tauByMVA5TightElectronRejection_);
    tree_->Branch("tauByMVA5VTightElectronRejection", &tauByMVA5VTightElectronRejection_);
    tree_->Branch("tauByLooseMuonRejection", &tauByLooseMuonRejection_);
    tree_->Branch("tauByMediumMuonRejection", &tauByMediumMuonRejection_);
    tree_->Branch("tauByTightMuonRejection", &tauByTightMuonRejection_);
    tree_->Branch("tauByLooseMuonRejection3", &tauByLooseMuonRejection3_);
    tree_->Branch("tauByTightMuonRejection3", &tauByTightMuonRejection3_);
    tree_->Branch("tauByMVALooseMuonRejection", &tauByMVALooseMuonRejection_);
    tree_->Branch("tauByMVAMediumMuonRejection", &tauByMVAMediumMuonRejection_);
    tree_->Branch("tauByMVATightMuonRejection", &tauByMVATightMuonRejection_);
    tree_->Branch("tauByMVArawMuonRejection", &tauByMVArawMuonRejection_);
    //    tree_->Branch("pfTausDiscriminationByDecayModeFinding", &pfTausDiscriminationByDecayModeFinding_);
    tree_->Branch("tauByVLooseIsolation", &tauByVLooseIsolation_);
    tree_->Branch("tauByVLooseCombinedIsolationDBSumPtCorr", &tauByVLooseCombinedIsolationDBSumPtCorr_);
    tree_->Branch("tauByLooseCombinedIsolationDBSumPtCorr", &tauByLooseCombinedIsolationDBSumPtCorr_);
    tree_->Branch("tauByMediumCombinedIsolationDBSumPtCorr", &tauByMediumCombinedIsolationDBSumPtCorr_);
    tree_->Branch("tauByTightCombinedIsolationDBSumPtCorr", &tauByTightCombinedIsolationDBSumPtCorr_);
    tree_->Branch("tauByLooseCombinedIsolationDBSumPtCorr3Hits", &tauByLooseCombinedIsolationDBSumPtCorr3Hits_);
    tree_->Branch("tauByMediumCombinedIsolationDBSumPtCorr3Hits", &tauByMediumCombinedIsolationDBSumPtCorr3Hits_);
    tree_->Branch("tauByTightCombinedIsolationDBSumPtCorr3Hits", &tauByTightCombinedIsolationDBSumPtCorr3Hits_);
    tree_->Branch("tauByVLooseIsolationMVA3newDMwoLT", &tauByVLooseIsolationMVA3newDMwoLT_);
    tree_->Branch("tauByLooseIsolationMVA3newDMwoLT", &tauByLooseIsolationMVA3newDMwoLT_);
    tree_->Branch("tauByMediumIsolationMVA3newDMwoLT", &tauByMediumIsolationMVA3newDMwoLT_);
    tree_->Branch("tauByTightIsolationMVA3newDMwoLT", &tauByTightIsolationMVA3newDMwoLT_);
    tree_->Branch("tauByVTightIsolationMVA3newDMwoLT", &tauByVTightIsolationMVA3newDMwoLT_);
    tree_->Branch("tauByVVTightIsolationMVA3newDMwoLT", &tauByVVTightIsolationMVA3newDMwoLT_);
    tree_->Branch("tauByIsolationMVA3newDMwoLTraw", &tauByIsolationMVA3newDMwoLTraw_);
    tree_->Branch("tauByVLooseIsolationMVA3oldDMwLT", &tauByVLooseIsolationMVA3oldDMwLT_);
    tree_->Branch("tauByLooseIsolationMVA3oldDMwLT", &tauByLooseIsolationMVA3oldDMwLT_);
    tree_->Branch("tauByMediumIsolationMVA3oldDMwLT", &tauByMediumIsolationMVA3oldDMwLT_);
    tree_->Branch("tauByTightIsolationMVA3oldDMwLT", &tauByTightIsolationMVA3oldDMwLT_);
    tree_->Branch("tauByVTightIsolationMVA3oldDMwLT", &tauByVTightIsolationMVA3oldDMwLT_);
    tree_->Branch("tauByVVTightIsolationMVA3oldDMwLT", &tauByVVTightIsolationMVA3oldDMwLT_);
    tree_->Branch("tauByIsolationMVA3oldDMwLTraw", &tauByIsolationMVA3oldDMwLTraw_);
    tree_->Branch("tauByVLooseIsolationMVA3oldDMwoLT", &tauByVLooseIsolationMVA3oldDMwoLT_);
    tree_->Branch("tauByLooseIsolationMVA3oldDMwoLT", &tauByLooseIsolationMVA3oldDMwoLT_);
    tree_->Branch("tauByTightIsolationMVA3oldDMwoLT", &tauByTightIsolationMVA3oldDMwoLT_);
    tree_->Branch("tauByVTightIsolationMVA3oldDMwoLT", &tauByVTightIsolationMVA3oldDMwoLT_);
    tree_->Branch("tauByVVTightIsolationMVA3oldDMwoLT", &tauByVVTightIsolationMVA3oldDMwoLT_);
    tree_->Branch("tauByIsolationMVA3newDMwoLTraw", &tauByIsolationMVA3newDMwoLTraw_);
    tree_->Branch("tauByLooseIsolationMVA3newDMwLT", &tauByLooseIsolationMVA3newDMwLT_);
    tree_->Branch("tauByVLooseIsolationMVA3newDMwLT", &tauByVLooseIsolationMVA3newDMwLT_);
    tree_->Branch("tauByMediumIsolationMVA3newDMwLT", &tauByMediumIsolationMVA3newDMwLT_);
    tree_->Branch("tauByTightIsolationMVA3newDMwLT", &tauByTightIsolationMVA3newDMwLT_);
    tree_->Branch("tauByVTightIsolationMVA3newDMwLT", &tauByVTightIsolationMVA3newDMwLT_);
    tree_->Branch("tauByVVTightIsolationMVA3newDMwLT", &tauByVVTightIsolationMVA3newDMwLT_);
    tree_->Branch("tauByIsolationMVA3newDMwLTraw", &tauByIsolationMVA3newDMwLTraw_);

    tree_->Branch("tauEta"  ,&tauEta_);
    tree_->Branch("tauPhi"  ,&tauPhi_);
    tree_->Branch("tauPt"  ,&tauPt_);
    tree_->Branch("tauEt"  ,&tauEt_);
    tree_->Branch("tauCharge"  ,&tauCharge_);
    tree_->Branch("tauDecayMode"  ,&tauDecayMode_);
    tree_->Branch("tauEMFraction"  ,&tauEMFraction_);
    tree_->Branch("tauHCAL3x3OverPLead"  ,&tauHCAL3x3OverPLead_);
    tree_->Branch("tauHCALMaxOverPLead"  ,&tauHCALMaxOverPLead_);
    tree_->Branch("tauHCALTotOverPLead"  ,&tauHCALTotOverPLead_);
    tree_->Branch("tauIsolationPFChargedHadrCandsPtSum"  ,&tauIsolationPFChargedHadrCandsPtSum_);
    tree_->Branch("tauIsolationPFGammaCandsEtSum"  ,&tauIsolationPFGammaCandsEtSum_);
    tree_->Branch("tauLeadPFChargedHadrCandsignedSipt"  ,&tauLeadPFChargedHadrCandsignedSipt_);
    tree_->Branch("tauLeadChargedHadronExists"  ,&tauLeadChargedHadronExists_);
    tree_->Branch("tauLeadChargedHadronEta"  ,&tauLeadChargedHadronEta_);
    tree_->Branch("tauLeadChargedHadronPhi"  ,&tauLeadChargedHadronPhi_);
    tree_->Branch("tauLeadChargedHadronPt"  ,&tauLeadChargedHadronPt_);
  }  
  
  
  
  if (dumpJets_) {
    tree_->Branch("nJet", &nJet_, "nJet/I");
    tree_->Branch("jetPt", &jetPt_);
    tree_->Branch("jetEta", &jetEta_);
    tree_->Branch("jetPhi", &jetPhi_);
    tree_->Branch("jetCHF", &jetCHF_);
    tree_->Branch("jetNHF", &jetNHF_);
    tree_->Branch("jetCEF", &jetCEF_);
    tree_->Branch("jetNEF", &jetNEF_);
    tree_->Branch("jetNCH", &jetNCH_);
    tree_->Branch("jetHFHAE", &jetHFHAE_);
    tree_->Branch("jetHFEME", &jetHFEME_);
    tree_->Branch("jetNConstituents", &jetNConstituents_);
    tree_->Branch("jetCombinedSecondaryVtxBJetTags", &jetCombinedSecondaryVtxBJetTags_);
    tree_->Branch("jetJetProbabilityBJetTags", &jetJetProbabilityBJetTags_);
    tree_->Branch("jetJetBProbabilityBJetTags", &jetJetBProbabilityBJetTags_);
    tree_->Branch("jetTrackCountingHighPurBJetTags", &jetTrackCountingHighPurBJetTags_);
    tree_->Branch("jetTrackCountingHighEffBJetTags", &jetTrackCountingHighEffBJetTags_);
    tree_->Branch("jetSimpleSecondaryVertexHighEffBJetTags", &jetSimpleSecondaryVertexHighEffBJetTags_);
    tree_->Branch("jetSimpleSecondaryVertexHighPurBJetTags", &jetSimpleSecondaryVertexHighPurBJetTags_);
if (doGenParticles_) tree_->Branch("jetPartonID", &jetPartonID_);
    tree_->Branch("jetPFLooseId", &jetPFLooseId_);

}

  // SubJet
  if (dumpSubJets_) {
    tree_->Branch("nCA8Jet",&nCA8Jet_, "nCA8Jet/I");
    tree_->Branch("CA8JetPt", &CA8JetPt_);
    tree_->Branch("CA8JetEta", &CA8JetEta_);
    tree_->Branch("CA8JetPhi", &CA8JetPhi_);
    tree_->Branch("CA8JetMass", &CA8JetMass_);
    tree_->Branch("CA8JetArea", &CA8JetArea_);
    tree_->Branch("CA8Jet_tau1", &CA8Jet_tau1_);
    tree_->Branch("CA8Jet_tau2", &CA8Jet_tau2_);
    tree_->Branch("CA8Jet_tau3", &CA8Jet_tau3_);
    tree_->Branch("CA8JetCHF", &CA8JetCHF_);
    tree_->Branch("CA8JetNHF",&CA8JetNHF_);
    tree_->Branch("CA8JetCEF",&CA8JetCEF_);
    tree_->Branch("CA8JetNEF",&CA8JetNEF_);
    tree_->Branch("CA8JetNCH",&CA8JetNCH_);
    tree_->Branch("CA8Jetnconstituents",&CA8Jetnconstituents_);
    tree_->Branch("CA8prunedJetMass", &CA8prunedJetMass_);
  }


  ////////////Prepare for photon ID MVA for Run II///////////
  //
  // Create and configure barrel MVA
  //
  tmvaReader_[0] = new TMVA::Reader( "!Color:!Silent:Error" );  
  tmvaReader_[0]->SetVerbose(kFALSE);
  // Add all the vars, we take the string with variable name from the weights file (the Expression field)
  tmvaReader_[0]->AddVariable("recoPhi"   , &varPhi_);
  tmvaReader_[0]->AddVariable("r9"        , &varR9_);
  tmvaReader_[0]->AddVariable("sieie_2012", &varSieie_);
  tmvaReader_[0]->AddVariable("sieip_2012", &varSieip_);
  tmvaReader_[0]->AddVariable("e1x3_2012/e5x5_2012"        , &varE1x3overE5x5_);
  tmvaReader_[0]->AddVariable("e2x2_2012/e5x5_2012"        , &varE2x2overE5x5_);
  tmvaReader_[0]->AddVariable("e2x5_2012/e5x5_2012"        , &varE2x5overE5x5_);
  tmvaReader_[0]->AddVariable("recoSCEta" , &varSCEta_);
  tmvaReader_[0]->AddVariable("rawE"      , &varRawE_);
  tmvaReader_[0]->AddVariable("scEtaWidth", &varSCEtaWidth_);
  tmvaReader_[0]->AddVariable("scPhiWidth", &varSCPhiWidth_);
  tmvaReader_[0]->AddVariable("rho"       , &varRho_);
  tmvaReader_[0]->AddVariable("phoIsoRaw" , &varPhoIsoRaw_);
  tmvaReader_[0]->AddVariable("chIsoRaw"  , &varChIsoRaw_);
  tmvaReader_[0]->AddVariable("chWorstRaw", &varWorstChRaw_);
  // Add spectators
  tmvaReader_[0]->AddSpectator("recoPt" , &varPt_);
  tmvaReader_[0]->AddSpectator("recoEta", &varEta_);

  //
  // Create and configure endcap MVA
  //
  tmvaReader_[1] = new TMVA::Reader( "!Color:!Silent:Error" );  
  tmvaReader_[1]->SetVerbose(kFALSE);
  // Add all the vars, we take the string with variable name from the weights file (the Expression field)
  tmvaReader_[1]->AddVariable("recoPhi"   , &varPhi_);
  tmvaReader_[1]->AddVariable("r9"        , &varR9_);
  tmvaReader_[1]->AddVariable("sieie_2012", &varSieie_);
  tmvaReader_[1]->AddVariable("sieip_2012", &varSieip_);
  tmvaReader_[1]->AddVariable("e1x3_2012/e5x5_2012"        , &varE1x3overE5x5_);
  tmvaReader_[1]->AddVariable("e2x2_2012/e5x5_2012"        , &varE2x2overE5x5_);
  tmvaReader_[1]->AddVariable("e2x5_2012/e5x5_2012"        , &varE2x5overE5x5_);
  tmvaReader_[1]->AddVariable("recoSCEta" , &varSCEta_);
  tmvaReader_[1]->AddVariable("rawE"      , &varRawE_);
  tmvaReader_[1]->AddVariable("scEtaWidth", &varSCEtaWidth_);
  tmvaReader_[1]->AddVariable("scPhiWidth", &varSCPhiWidth_);
  tmvaReader_[1]->AddVariable("esEn/rawE" , &varESEnOverRawE_);
  tmvaReader_[1]->AddVariable("esRR"      , &varESEffSigmaRR_);
  tmvaReader_[1]->AddVariable("rho"       , &varRho_);
  tmvaReader_[1]->AddVariable("phoIsoRaw" , &varPhoIsoRaw_);
  tmvaReader_[1]->AddVariable("chIsoRaw"  , &varChIsoRaw_);
  tmvaReader_[1]->AddVariable("chWorstRaw", &varWorstChRaw_);
  // Add spectators
  tmvaReader_[1]->AddSpectator("recoPt" , &varPt_);
  tmvaReader_[1]->AddSpectator("recoEta", &varEta_);

  //
  // Book the MVA method for each category
  //
  std::string cmssw_base_src = getenv("CMSSW_BASE");
  cmssw_base_src += "/src/";
  //
  TString localFileName1 = "EgammaAnalysis/PhotonTools/data/PHYS14/photon_general_MVA_phys14_pu20bx25_EB_V1.weights.xml";
  TString weightsFileName1 = TString(cmssw_base_src) + localFileName1;
  //methodName_[0] = "BDT photons barrel";
  methodName_[0] = "BDT";
  tmvaReader_[0]->BookMVA(methodName_[0], weightsFileName1);
  //
  TString localFileName2 = "EgammaAnalysis/PhotonTools/data/PHYS14/photon_general_MVA_phys14_pu20bx25_EE_V1.weights.xml";
  TString weightsFileName2 = TString(cmssw_base_src) + localFileName2;
  //methodName_[1] = "BDT photons endcap";
  methodName_[1] = "BDT";
  tmvaReader_[1]->BookMVA(methodName_[1], weightsFileName2);



}

ggNtuplizer::~ggNtuplizer() {
  delete cicPhotonId_;

  
  delete tmvaReader_[0];
  delete tmvaReader_[1];

}

void ggNtuplizer::getHandles(const edm::Event & event,
			     edm::Handle<VertexCollection>                & recVtxs,
			     edm::Handle<VertexCollection>                & recVtxsBS,
			     edm::Handle<double>                          & rhoHandle,
			     edm::Handle<edm::View<pat::MET> >            & pfMETHandle,
			     edm::Handle<edm::View<pat::Electron> >       & electronHandle,
			     edm::Handle<edm::View<pat::Photon> >         & photonHandle,
			     edm::Handle<edm::View<pat::Muon> >           & muonHandle,
			     edm::Handle<std::vector<pat::Tau> >            & tauHandle,
			     edm::Handle<EcalRecHitCollection>            & EBReducedRecHits,
			     edm::Handle<EcalRecHitCollection>            & EEReducedRecHits,
			     edm::Handle<EcalRecHitCollection>            & ESReducedRecHits,
			     edm::Handle<reco::PhotonCollection>          & recoPhotonHandle,
			     edm::Handle<TrackCollection>                 & tracksHandle,
			     edm::Handle<GsfElectronCollection>           & gsfElectronHandle,
			     edm::Handle<PFCandidateCollection>           & pfAllCandidates,
                             edm::Handle<edm::View<pat::Jet> > &        jetHandle
			     ) { 

  event.getByLabel(vtxLabel_,                  recVtxs);
  event.getByLabel(vtxBSLabel_,                recVtxsBS);
  event.getByLabel(rhoLabel_,                  rhoHandle);
  event.getByLabel(pfMETlabel_,                pfMETHandle);
  event.getByLabel(electronCollection_,        electronHandle);
  event.getByLabel(photonCollection_,          photonHandle);
  event.getByLabel(muonCollection_,            muonHandle);
  event.getByLabel(tauCollection_,             tauHandle);
  event.getByToken(ebReducedRecHitCollection_, EBReducedRecHits);
  event.getByToken(eeReducedRecHitCollection_, EEReducedRecHits);
  event.getByToken(esReducedRecHitCollection_, ESReducedRecHits);
  event.getByLabel(recophotonCollection_,      recoPhotonHandle);
  event.getByLabel(tracklabel_,                tracksHandle);
  event.getByLabel(gsfElectronlabel_,          gsfElectronHandle);
  event.getByLabel(pfAllParticles_,            pfAllCandidates);
  event.getByLabel(jetCollection_         , jetHandle);

}

void ggNtuplizer::analyze(const edm::Event& e, const edm::EventSetup& es) {

  this->getHandles(e, recVtxs_, recVtxsBS_, rhoHandle_, pfMETHandle_, electronHandle_,
		   photonHandle_, muonHandle_, tauHandle_, EBReducedRecHits_, EEReducedRecHits_, ESReducedRecHits_, 
		   recoPhotonHandle_, tracksHandle_, gsfElectronHandle_, pfAllCandidates_,
                    jetHandle_); 
  clearVectors();
  hEvents_->Fill(0.5);

  GEDPhoIDTools *GEDIdTool = new GEDPhoIDTools(e);

  lazyTool = new EcalClusterLazyTools(e, es, ebReducedRecHitCollection_, eeReducedRecHitCollection_, esReducedRecHitCollection_);
  lazyToolnoZS = new noZS::EcalClusterLazyTools(e, es, ebReducedRecHitCollection_, eeReducedRecHitCollection_, esReducedRecHitCollection_);

  run_    = e.id().run();
  event_  = e.id().event();
  lumis_  = e.luminosityBlock();
  isData_ = e.isRealData();
  rho_    = *(rhoHandle_.product());

  cicPhotonId_->configure(recVtxsBS_, tracksHandle_, gsfElectronHandle_, pfAllCandidates_, rho_); 

  // vertex
  nVtx_ = 0;
  math::XYZPoint pv(0, 0, 0);
  if (recVtxs_.isValid()) {
    for (size_t i=0; i<recVtxs_->size(); ++i) {
      if (!((*recVtxs_)[i].isFake())) {
	if (nVtx_ == 0) {
	  pv.SetXYZ((*recVtxs_)[i].x(), (*recVtxs_)[i].y(), (*recVtxs_)[i].z());
	  nTrks_ = (*recVtxs_)[i].tracksSize();
	}
	nVtx_++;
      }
    }
  }

  // PDF information, pthat and processID
  nPUInfo_ = 0; 
  if (!isData_) {
    
    pthat_ = -99;
    processID_ = -99;

    Handle<GenEventInfoProduct> pdfInfoHandle;
    if (e.getByLabel(generatorLabel_, pdfInfoHandle)) {
      if (pdfInfoHandle->pdf()) {
	pdf_.push_back(pdfInfoHandle->pdf()->id.first);    // PDG ID of incoming parton #1
	pdf_.push_back(pdfInfoHandle->pdf()->id.second);   // PDG ID of incoming parton #2
	pdf_.push_back(pdfInfoHandle->pdf()->x.first);     // x value of parton #1
	pdf_.push_back(pdfInfoHandle->pdf()->x.second);    // x value of parton #2
	pdf_.push_back(pdfInfoHandle->pdf()->xPDF.first);  // PDF weight for parton #1
	pdf_.push_back(pdfInfoHandle->pdf()->xPDF.second); // PDF weight for parton #2
	pdf_.push_back(pdfInfoHandle->pdf()->scalePDF);    // scale of the hard interaction
      }
    }

    Handle<GenEventInfoProduct> genEventScale;
    if (e.getByLabel(generatorLabel_, genEventScale)) {
      if (genEventScale->hasBinningValues()) pthat_ = genEventScale->binningValues()[0];
      processID_ = genEventScale->signalProcessID();
    }

    Handle< vector<PileupSummaryInfo> > pileupHandle;
    if (e.getByLabel(puCollection_, pileupHandle)) {
      vector<PileupSummaryInfo>::const_iterator PVI;
      for (PVI = pileupHandle->begin(); PVI != pileupHandle->end(); ++PVI) {
	if (PVI->getBunchCrossing() == 0) {
	  hPU_->Fill(PVI->getPU_NumInteractions());
	  hPUTrue_->Fill(PVI->getTrueNumInteractions());
	}
	nPU_.push_back(PVI->getPU_NumInteractions());
	puTrue_.push_back(PVI->getTrueNumInteractions());
	puBX_.push_back(PVI->getBunchCrossing());
	nPUInfo_++;
      }
    }
  }
  
  if (!isData_ && doGenParticles_)
    fillGenParticles(e);

  // Gen, pfMET  
  if (pfMETHandle_.isValid()) {
    const pat::MET *pfMET = 0;
    pfMET = &(pfMETHandle_->front());

    pfMET_       = pfMET->et();
    pfMETPhi_    = pfMET->phi();
    pfMETsumEt_  = pfMET->sumEt();
    pfMETmEtSig_ = (pfMET->mEtSig() < 1.e10) ? pfMET->mEtSig() : 0;
    pfMETSig_    = (pfMET->significance() < 1.e10) ? pfMET->significance() : 0;;

    if (!isData_) {
      genMET_    = pfMET->genMET()->et();
      genMETPhi_ = pfMET->genMET()->phi();
    }
  }
  
  // electrons
  nEle_ = 0;
  if (electronHandle_.isValid()) {
    for (View<pat::Electron>::const_iterator iEle = electronHandle_->begin(); iEle != electronHandle_->end(); ++iEle) {

      eleCharge_          .push_back(iEle->charge());
      eleChargeConsistent_.push_back((Int_t)iEle->isGsfCtfScPixChargeConsistent());
      eleEn_              .push_back(iEle->energy());
      eleD0_              .push_back(iEle->gsfTrack()->dxy(pv));
      eleDz_              .push_back(iEle->gsfTrack()->dz(pv));
      elePt_              .push_back(iEle->pt());
      eleEta_             .push_back(iEle->eta());
      eleTheta_             .push_back(iEle->theta());

      elePhi_             .push_back(iEle->phi());
      eleSCEn_            .push_back(iEle->superCluster()->energy());
      eleESEn_            .push_back(iEle->superCluster()->preshowerEnergy());
      eleSCEta_           .push_back(iEle->superCluster()->eta());
      eleSCPhi_           .push_back(iEle->superCluster()->phi());
      eleSCRawEn_         .push_back(iEle->superCluster()->rawEnergy());
      eleSCEtaWidth_      .push_back(iEle->superCluster()->etaWidth());
      eleSCPhiWidth_      .push_back(iEle->superCluster()->phiWidth());

      ///SJ - 16th April
      ///https://cmssdt.cern.ch/SDT/doxygen/CMSSW_7_2_2/doc/html/d8/dac/GsfElectron_8h_source.html
      //eleHoverE_          .push_back(iEle->hcalOverEcalBc());
      eleHoverE_          .push_back(iEle->hcalOverEcal());

      ///https://cmssdt.cern.ch/SDT/doxygen/CMSSW_7_2_2/doc/html/d8/dac/GsfElectron_8h_source.html
      eleEoverP_          .push_back(iEle->eSuperClusterOverP());
      eleEoverPInv_       .push_back(fabs(1./iEle->ecalEnergy()-1./iEle->trackMomentumAtVtx().R()));
      eleBrem_            .push_back(iEle->fbrem());
      eledEtaAtVtx_       .push_back(iEle->deltaEtaSuperClusterTrackAtVtx());
      eledPhiAtVtx_       .push_back(iEle->deltaPhiSuperClusterTrackAtVtx());
      eleSigmaIEtaIEta_   .push_back(iEle->sigmaIetaIeta()); ///new sigmaietaieta
      eleSigmaIEtaIPhi_   .push_back(iEle->sigmaIetaIphi());
      eleSigmaIPhiIPhi_   .push_back(iEle->sigmaIphiIphi());
      eleConvVeto_        .push_back((Int_t)iEle->passConversionVeto()); // ConvVtxFit || missHit == 0
      eleMissHits_        .push_back(iEle->gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS));
      eleESEffSigmaRR_    .push_back(lazyTool->eseffsirir(*((*iEle).superCluster())));
      //elePFChIso_         .push_back(iEle->chargedHadronIso());
      //elePFPhoIso_        .push_back(iEle->photonIso());
      //elePFNeuIso_        .push_back(iEle->neutralHadronIso());

      ///HEEP ID
      double eledEtaseedAtVtx = iEle->superCluster().isNonnull() && iEle->superCluster()->seed().isNonnull() ? 
	iEle->deltaEtaSuperClusterTrackAtVtx() - iEle->superCluster()->eta() + iEle->superCluster()->seed()->eta() : std::numeric_limits<float>::max();
      
      eledEtaseedAtVtx_   .push_back(eledEtaseedAtVtx);
      
      eleE1x5_            .push_back(iEle->e1x5());
      eleE2x5_            .push_back(iEle->e2x5Max());
      eleE5x5_            .push_back(iEle->e5x5());

      GsfElectron::PflowIsolationVariables pfIso = iEle->pfIsolationVariables();
      elePFChIso_         .push_back(pfIso.sumChargedHadronPt);
      elePFPhoIso_        .push_back(pfIso.sumPhotonEt);
      elePFNeuIso_        .push_back(pfIso.sumNeutralHadronEt);
      elePFPUIso_         .push_back(pfIso.sumPUPt);

      elecaloEnergy_      .push_back(iEle->caloEnergy());
      ///SJ - https://github.com/lgray/cmssw/blob/common_isolation_selection_70X/TestElectronID/ElectronIDAnalyzer/plugins/ElectronIDAnalyzer.cc
      
      float absiso = elePFChIso_[nEle_]
	+ std::max(0.0 , elePFNeuIso_[nEle_] + elePFPhoIso_[nEle_] - 0.5 * elePFPUIso_[nEle_]);
      eleRelIsoWithDBeta_ .push_back(absiso/elePt_[nEle_]);


      /////quantities which were used for Run1 - these do not 
      ///calculated through PF (meaning no energy is subtracted 
      ///using PF)
      ///https://cmssdt.cern.ch/SDT/doxygen/CMSSW_7_2_2/doc/html/d9/d44/ElectronIDValueMapProducer_8cc_source.html
      ///line 120

      std::vector<float> vCovEle = lazyToolnoZS->localCovariances( *((*iEle).superCluster()->seed()) );
      const float seeEle = (isnan(vCovEle[0]) ? 0. : sqrt(vCovEle[0]));

      eleSigmaIEtaIEta_2012_ .push_back(seeEle); //sigmaietaieta of run1
      
      const reco::CaloClusterPtr seed = (*iEle).superCluster()->seed();
      eleE1x5_2012_            .push_back(lazyToolnoZS->e1x5(*seed));
      eleE2x5_2012_            .push_back(lazyToolnoZS->e2x5Max(*seed));
      eleE5x5_2012_            .push_back(lazyToolnoZS->e5x5(*seed));

      ///For HEEP ID
      eleEcalDrivenSeed_                       .push_back(iEle->ecalDrivenSeed());
      eleDr03EcalRecHitSumEt_                  .push_back(iEle->dr03EcalRecHitSumEt());
      eleDr03HcalDepth1TowerSumEt_             .push_back(iEle->dr03HcalDepth1TowerSumEt());
      eleDr03HcalDepth2TowerSumEt_             .push_back(iEle->dr03HcalDepth2TowerSumEt());
      eleDr03HcalTowerSumEt_                   .push_back(iEle->dr03HcalTowerSumEt());
      eleDr03TkSumPt_                          .push_back(iEle->dr03TkSumPt());


      reco::GsfTrackRef trackref = iEle->gsfTrack();
      if(iEle->gsfTrack().isNonnull())
	{
	  if(recVtxs_->size() > 0 ){
	    eleTrkdxy_        .push_back(trackref->dxy(recVtxs_->front().position()));
	  }//if(recVtxs->size() > 0 )
	  
	}//if(iEle->gsfTrack().isNonnull())
      
      
      eleBC1E_               .push_back((*iEle).superCluster()->seed()->energy());
      eleBC1Eta_             .push_back((*iEle).superCluster()->seed()->eta());

      Int_t nBCEle = 0;
      for (CaloCluster_iterator itbc = iEle->superCluster()->clustersBegin(); itbc != iEle->superCluster()->clustersEnd(); ++itbc) {
	if (nBCEle == 1) {
	  eleBC2E_  .push_back((*itbc)->energy());
	  eleBC2Eta_.push_back((*itbc)->eta());
	}
	nBCEle++;
      }
      if (nBCEle == 1) {
	eleBC2E_  .push_back(-99.);
	eleBC2Eta_.push_back(-99.);
      }
      
      ////////9th April, 2015 - SHILPI JAIN
      ///MVA for electrons in PHYS14
      
      eleIDMVA_   .push_back(iEle->electronID("trigMVAid"));

      nEle_++;
    }
  } 

  // photons
  nPho_ = 0;
  if (photonHandle_.isValid()) {
    for (View<pat::Photon>::const_iterator iPho = photonHandle_->begin(); iPho != photonHandle_->end(); ++iPho) {

      phoE_             .push_back(iPho->energy());
      phoEt_            .push_back(iPho->et());
      phoEta_           .push_back(iPho->eta());
      phoPhi_           .push_back(iPho->phi());
      phoSCE_           .push_back((*iPho).superCluster()->energy());
      phoSCRawE_        .push_back((*iPho).superCluster()->rawEnergy());
      phoESEn_          .push_back((*iPho).superCluster()->preshowerEnergy());
      phoSCEta_         .push_back((*iPho).superCluster()->eta());
      phoSCPhi_         .push_back((*iPho).superCluster()->phi());
      phoSCEtaWidth_    .push_back((*iPho).superCluster()->etaWidth());
      phoSCPhiWidth_    .push_back((*iPho).superCluster()->phiWidth());
      phoSCBrem_        .push_back((*iPho).superCluster()->phiWidth()/(*iPho).superCluster()->etaWidth());
      phohasPixelSeed_  .push_back((Int_t)iPho->hasPixelSeed());
      phoEleVeto_       .push_back((Int_t)iPho->passElectronVeto());
      phoR9_            .push_back(iPho->r9());
      phoHoverE_        .push_back(iPho->hadTowOverEm());
      phoSigmaIEtaIEta_ .push_back(iPho->see());
      phoSigmaIEtaIPhi_ .push_back(iPho->sep());
      phoSigmaIPhiIPhi_ .push_back(iPho->spp());
      phoE1x3_          .push_back(lazyTool->e1x3(*((*iPho).superCluster()->seed())));
      phoE2x2_          .push_back(lazyTool->e2x2(*((*iPho).superCluster()->seed())));
      phoE2x5Max_       .push_back(lazyTool->e2x5Max(*((*iPho).superCluster()->seed())));
      phoE5x5_          .push_back(lazyTool->e5x5(*((*iPho).superCluster()->seed())));
      phoESEffSigmaRR_  .push_back(lazyTool->eseffsirir(*((*iPho).superCluster())));
      //phoPFChIso_       .push_back(iPho->chargedHadronIso());
      //phoPFPhoIso_      .push_back(iPho->photonIso());
      //phoPFNeuIso_      .push_back(iPho->neutralHadronIso());

      std::vector<float> vCov = lazyToolnoZS->localCovariances( *((*iPho).superCluster()->seed()) );
      const float see = (isnan(vCov[0]) ? 0. : sqrt(vCov[0]));
      const float spp = (isnan(vCov[2]) ? 0. : sqrt(vCov[2]));
      const float sep = vCov[1];

      phoSigmaIEtaIEta_2012_ .push_back(see);
      phoSigmaIEtaIPhi_2012_ .push_back(sep);
      phoSigmaIPhiIPhi_2012_ .push_back(spp);
      phoE1x3_2012_          .push_back(lazyToolnoZS->e1x3(*((*iPho).superCluster()->seed())));
      phoE2x2_2012_          .push_back(lazyToolnoZS->e2x2(*((*iPho).superCluster()->seed())));
      phoE2x5Max_2012_       .push_back(lazyToolnoZS->e2x5Max(*((*iPho).superCluster()->seed())));
      phoE5x5_2012_          .push_back(lazyToolnoZS->e5x5(*((*iPho).superCluster()->seed())));

      size_t rightRecoPho = -1;
      for (size_t iv = 0; iv < recoPhotonHandle_->size(); ++iv) {
	reco::PhotonRef recophoRef2(recoPhotonHandle_, iv);
	if (deltaR(iPho->eta(), iPho->phi(), recophoRef2->eta(), recophoRef2->phi()) < 0.01) rightRecoPho = iv;
      }
      reco::PhotonRef recophoRef(recoPhotonHandle_, rightRecoPho);
      reco::Vertex pv = recVtxs_->at(0);
      GEDIdTool->setPhotonP4(recophoRef, pv);

      phoPFChIso_       .push_back(GEDIdTool->SolidConeIso(0.3, reco::PFCandidate::h));
      phoPFPhoIso_      .push_back(GEDIdTool->SolidConeIso(0.3, reco::PFCandidate::gamma));
      phoPFNeuIso_      .push_back(GEDIdTool->SolidConeIso(0.3, reco::PFCandidate::h0));

      std::vector<double> IsoRings;
      GEDIdTool->FrixioneIso(0.1, 8, reco::PFCandidate::h, IsoRings);
      phoPFChIsoFrix1_.push_back(IsoRings[0]);
      phoPFChIsoFrix2_.push_back(IsoRings[1]);
      phoPFChIsoFrix3_.push_back(IsoRings[2]);
      phoPFChIsoFrix4_.push_back(IsoRings[3]);
      phoPFChIsoFrix5_.push_back(IsoRings[4]);
      phoPFChIsoFrix6_.push_back(IsoRings[5]);
      phoPFChIsoFrix7_.push_back(IsoRings[6]);
      phoPFChIsoFrix8_.push_back(IsoRings[7]);
      IsoRings.resize(0);

      GEDIdTool->FrixioneIso(0.1, 8, reco::PFCandidate::gamma, IsoRings);
      phoPFPhoIsoFrix1_.push_back(IsoRings[0]);
      phoPFPhoIsoFrix2_.push_back(IsoRings[1]);
      phoPFPhoIsoFrix3_.push_back(IsoRings[2]);
      phoPFPhoIsoFrix4_.push_back(IsoRings[3]);
      phoPFPhoIsoFrix5_.push_back(IsoRings[4]);
      phoPFPhoIsoFrix6_.push_back(IsoRings[5]);
      phoPFPhoIsoFrix7_.push_back(IsoRings[6]);
      phoPFPhoIsoFrix8_.push_back(IsoRings[7]);
      IsoRings.resize(0);

      GEDIdTool->FrixioneIso(0.1, 8, reco::PFCandidate::h0, IsoRings);
      phoPFNeuIsoFrix1_.push_back(IsoRings[0]);
      phoPFNeuIsoFrix2_.push_back(IsoRings[1]);
      phoPFNeuIsoFrix3_.push_back(IsoRings[2]);
      phoPFNeuIsoFrix4_.push_back(IsoRings[3]);
      phoPFNeuIsoFrix5_.push_back(IsoRings[4]);
      phoPFNeuIsoFrix6_.push_back(IsoRings[5]);
      phoPFNeuIsoFrix7_.push_back(IsoRings[6]);
      phoPFNeuIsoFrix8_.push_back(IsoRings[7]);

      std::vector<float> vtxIsolations03 = cicPhotonId_->pfTkIsoWithVertex(recophoRef, 0.3, 0.0, 0.0, 0.0, 0.2, 0.1, reco::PFCandidate::h);
      phoPFChWorstIso_  .push_back(*max_element(vtxIsolations03.begin(), vtxIsolations03.end()));

      phoBC1E_          .push_back((*iPho).superCluster()->seed()->energy());
      phoBC1Eta_        .push_back((*iPho).superCluster()->seed()->eta());

      Int_t nBCPho = 0;
      for (CaloCluster_iterator itbc = iPho->superCluster()->clustersBegin(); itbc != iPho->superCluster()->clustersEnd(); ++itbc) {
	if (nBCPho == 1) {
	  phoBC2E_  .push_back((*itbc)->energy());
	  phoBC2Eta_.push_back((*itbc)->eta());
	}
	nBCPho++;
      }
      if (nBCPho == 1) {
	phoBC2E_  .push_back(-99.);
	phoBC2Eta_.push_back(-99.);
      }

      ///SJ - isolation variables
      phoEcalRecHitSumEtConeDR03_                   .push_back(iPho->ecalRecHitSumEtConeDR03());
      phohcalDepth1TowerSumEtConeDR03_              .push_back(iPho->hcalDepth1TowerSumEtConeDR03());
      phohcalDepth2TowerSumEtConeDR03_              .push_back(iPho->hcalDepth2TowerSumEtConeDR03());
      phohcalTowerSumEtConeDR03_                    .push_back(iPho->hcalTowerSumEtConeDR03());
      photrkSumPtHollowConeDR03_                    .push_back(iPho->trkSumPtHollowConeDR03());
    
      ////////9th April, 2015 - SJ
      ///MVA for photons in PHYS14
      // set MVA variables
      
      varPhi_          = phoPhi_[nPho_];
      varR9_           = phoR9_[nPho_];
      varSieie_        = phoSigmaIEtaIEta_2012_[nPho_];
      varSieip_        = phoSigmaIEtaIPhi_2012_[nPho_];
      varE1x3overE5x5_ = phoE1x3_2012_[nPho_]/phoE5x5_2012_[nPho_];
      varE2x2overE5x5_ = phoE2x2_2012_[nPho_]/phoE5x5_2012_[nPho_];
      varSCEta_        = phoSCEta_[nPho_];
      varRawE_         = phoSCRawE_[nPho_];
      varSCEtaWidth_   = phoSCEtaWidth_[nPho_];
      varSCPhiWidth_   = phoSCPhiWidth_[nPho_];
      varESEnOverRawE_ = phoESEn_[nPho_]/phoSCRawE_[nPho_];
      varESEffSigmaRR_ = phoESEffSigmaRR_[nPho_];
      varRho_          = rho_;
      varPhoIsoRaw_    = phoPFPhoIso_[nPho_];
      varChIsoRaw_     = phoPFChIso_[nPho_];
      varWorstChRaw_   = phoPFNeuIso_[nPho_];
      // Declare spectator vars 
      varPt_           = phoEt_[nPho_];
      varEta_          = phoEta_[nPho_];

      // 0=ECAL barrel or 1=ECAL endcaps
      int iBE = (fabs(phoSCEta_[nPho_]) < 1.479) ? 0 : 1;

      phoIDMVA_                   .push_back(tmvaReader_[iBE]->EvaluateMVA("BDT"));
      


      nPho_++;
    }
  }

  // muons
  nMu_ = 0;
  if ( muonHandle_.isValid() ) {
    for (View<pat::Muon>::const_iterator iMu = muonHandle_->begin(); iMu != muonHandle_->end(); ++iMu) {

      if (iMu->pt() < 5) continue;
      if (! (iMu->isPFMuon() || iMu->isGlobalMuon() || iMu->isTrackerMuon())) continue; 

      muPt_    .push_back(iMu->pt());
      muEta_   .push_back(iMu->eta());
      muPhi_   .push_back(iMu->phi());
      muCharge_.push_back(iMu->charge());
      muType_  .push_back(iMu->type());
      muIsGood_.push_back((int) iMu->isGood("TMOneStationTight"));
      muD0_    .push_back(iMu->muonBestTrack()->dxy(pv));
      muDz_    .push_back(iMu->muonBestTrack()->dz(pv));

      ///SJ
      muBestTrkPtError_  .push_back(iMu->muonBestTrack()->ptError());
      muBestTrkPt_       .push_back(iMu->muonBestTrack()->pt());

      ///SJ
      musegmentCompatibility_  .push_back(iMu->segmentCompatibility());
      muchi2LocalPosition_     .push_back(iMu->combinedQuality().chi2LocalPosition);

      mutrkKink_               .push_back(iMu->combinedQuality().trkKink);
      
      
      const reco::TrackRef glbmu = iMu->globalTrack();
      const reco::TrackRef innmu = iMu->innerTrack();

      if (glbmu.isNull()) {
	muChi2NDF_ .push_back(-99.);
	muMuonHits_.push_back(-99);
      } else {
	muChi2NDF_.push_back(glbmu->normalizedChi2());
        muMuonHits_.push_back(glbmu->hitPattern().numberOfValidMuonHits());
      }

      if (innmu.isNull()) {
	muInnerD0_     .push_back(-99.);
	muInnerDz_     .push_back(-99.);
	muTrkLayers_   .push_back(-99);
	muPixelLayers_ .push_back(-99);
	muPixelHits_   .push_back(-99);
	muTrkQuality_  .push_back(-99);

	///SJ
	muInnervalidFraction_ .push_back(-99);
	
      } else {
	muInnerD0_     .push_back(innmu->dxy(pv));
        muInnerDz_     .push_back(innmu->dz(pv));
	muTrkLayers_   .push_back(innmu->hitPattern().trackerLayersWithMeasurement());
	muPixelLayers_ .push_back(innmu->hitPattern().pixelLayersWithMeasurement());
        muPixelHits_   .push_back(innmu->hitPattern().numberOfValidPixelHits());
	muTrkQuality_  .push_back(innmu->quality(reco::TrackBase::highPurity));

	muInnervalidFraction_ .push_back(innmu->validFraction()); 
      }

      muStations_  .push_back(iMu->numberOfMatchedStations());
      muIsoTrk_    .push_back(iMu->trackIso());
      muPFChIso_   .push_back(iMu->pfIsolationR04().sumChargedHadronPt);
      muPFPhoIso_  .push_back(iMu->pfIsolationR04().sumPhotonEt);
      muPFNeuIso_  .push_back(iMu->pfIsolationR04().sumNeutralHadronEt);
      muPFPUIso_   .push_back(iMu->pfIsolationR04().sumPUPt);

      nMu_++;
    }
  }
  
  //start jets Lvdp
if(dumpJets_){
  nJet_ = 0;
  if (jetHandle_.isValid())
    for (View<pat::Jet>::const_iterator iJet = jetHandle_->begin(); iJet != jetHandle_->end(); ++iJet) {
      jetPt_.push_back(    iJet->pt());
      jetEta_.push_back(   iJet->eta());
      jetPhi_.push_back(   iJet->phi());
      jetCEF_.push_back(   iJet->chargedEmEnergyFraction());
      jetNEF_.push_back(   iJet->neutralEmEnergyFraction());
      jetCHF_.push_back(   iJet->chargedHadronEnergyFraction());
      jetNHF_.push_back(   iJet->neutralHadronEnergyFraction());
      jetHFHAE_.push_back( iJet->HFHadronEnergy());
      jetHFEME_.push_back( iJet->HFEMEnergy());
      jetNCH_.push_back(   iJet->chargedMultiplicity());
      jetNConstituents_.push_back(iJet->getPFConstituents().size());
//b-tagging
      jetCombinedSecondaryVtxBJetTags_.push_back(iJet->bDiscriminator("combinedSecondaryVertexBJetTags"));
      jetJetProbabilityBJetTags_.push_back(iJet->bDiscriminator("jetProbabilityBJetTags")); 
      jetJetBProbabilityBJetTags_.push_back(iJet->bDiscriminator("jetBProbabilityBJetTags"));
      jetTrackCountingHighPurBJetTags_.push_back(iJet->bDiscriminator("trackCountingHighPurBJetTags"));
      jetTrackCountingHighEffBJetTags_.push_back(iJet->bDiscriminator("trackCountingHighEffBJetTags"));
      jetSimpleSecondaryVertexHighEffBJetTags_.push_back(iJet->bDiscriminator("simpleSecondaryVertexHighEffBJetTags"));
      jetSimpleSecondaryVertexHighPurBJetTags_.push_back(iJet->bDiscriminator("simpleSecondaryVertexHighPurBJetTags"));
//parton id
      jetPartonID_.push_back(iJet->partonFlavour());
//jet PF Loose ID
      pat::strbitset retjet = pfLooseId_.getBitTemplate();
      jetPFLooseId_.push_back(pfLooseId_(*iJet, retjet));
      nJet_++;
}
}
  if(dumpSubJets_){
    nCA8Jet_ = 0;
    //jet substructure
    
    edm::Handle<edm::View<pat::Jet> > jetsCHS;
    e.getByLabel(jetsCHSLabel_, jetsCHS);
    
    View<pat::Jet>::const_iterator beginCHS = jetsCHS->begin();
    View<pat::Jet>::const_iterator endCHS = jetsCHS->end();
    View<pat::Jet>::const_iterator ijetCHS = beginCHS;
    
    
    // Loop over the "hard" jets
    for(ijetCHS = beginCHS; ijetCHS != endCHS; ++ijetCHS ){
      if( ijetCHS->pt() < 30.0 ) continue;
      nCA8Jet_++;
      CA8JetPt_.push_back( ijetCHS->pt() );
      CA8JetEta_.push_back( ijetCHS->eta() );
      CA8JetPhi_.push_back( ijetCHS->phi() );
      CA8JetMass_.push_back( ijetCHS->mass() );
      CA8Jet_tau1_.push_back( ijetCHS->userFloat("NjettinessCA8:tau1") );
      CA8Jet_tau2_.push_back( ijetCHS->userFloat("NjettinessCA8:tau2") );
      CA8Jet_tau3_.push_back( ijetCHS->userFloat("NjettinessCA8:tau3") );
      
      CA8JetCHF_.push_back( ijetCHS->chargedHadronEnergyFraction()); // 0.0
      CA8JetNHF_.push_back( ( ijetCHS->neutralHadronEnergy() + ijetCHS->HFHadronEnergy() ) / ijetCHS->energy()); //0.99
      CA8JetCEF_.push_back( ijetCHS->chargedEmEnergyFraction()); //0.99
      CA8JetNEF_.push_back( ijetCHS->neutralEmEnergyFraction()); //0.99
      CA8JetNCH_.push_back( ijetCHS->chargedMultiplicity()); //0
      CA8Jetnconstituents_.push_back( ijetCHS->numberOfDaughters()); //1
      CA8prunedJetMass_.push_back(ijetCHS->userFloat("ca8PFJetsCHSPrunedLinks"));
} 
 }//endjets Lvdp
  
//startTaus Lvdp
 nTau_ = 0;
if ( tauHandle_.isValid() ) {
for(vector<pat::Tau>::const_iterator itau = tauHandle_->begin(); itau != tauHandle_->end(); ++itau) {

  tauByLooseElectronRejection_.push_back(itau->tauID("byLooseElectronRejection"));
  tauByMediumElectronRejection_.push_back(itau->tauID("byMediumElectronRejection"));
  tauByTightElectronRejection_.push_back(itau->tauID("byTightElectronRejection"));
  tauByMVA5LooseElectronRejection_.push_back(itau->tauID("byMVA5LooseElectronRejection"));
  tauByMVA5MediumElectronRejection_.push_back(itau->tauID("byMVA5MediumElectronRejection"));
  tauByMVA5TightElectronRejection_.push_back(itau->tauID("byMVA5TightElectronRejection"));
  tauByMVA5VTightElectronRejection_.push_back(itau->tauID("byMVA5VTightElectronRejection"));
  tauByLooseMuonRejection_.push_back(itau->tauID("byLooseMuonRejection"));
  tauByMediumMuonRejection_.push_back(itau->tauID("byMediumMuonRejection"));
  tauByTightMuonRejection_.push_back(itau->tauID("byTightMuonRejection"));
  tauByLooseMuonRejection3_.push_back(itau->tauID("byLooseMuonRejection3"));
  tauByTightMuonRejection3_.push_back(itau->tauID("byTightMuonRejection3"));
  tauByMVALooseMuonRejection_.push_back(itau->tauID("byMVALooseMuonRejection"));
  tauByMVAMediumMuonRejection_.push_back(itau->tauID("byMVAMediumMuonRejection"));
  tauByMVATightMuonRejection_.push_back(itau->tauID("byMVATightMuonRejection"));
  tauByMVArawMuonRejection_.push_back(itau->tauID("byMVArawMuonRejection"));
  //  pfTausDiscriminationByDecayModeFinding_.push_back(itau->tauID("pfTausDiscriminationByDecayModeFinding"));
  tauByVLooseIsolation_.push_back(itau->tauID("byVLooseIsolation"));
  tauByVLooseCombinedIsolationDBSumPtCorr_.push_back(itau->tauID("byVLooseCombinedIsolationDBSumPtCorr"));
  tauByLooseCombinedIsolationDBSumPtCorr_.push_back(itau->tauID("byLooseCombinedIsolationDBSumPtCorr"));
  tauByMediumCombinedIsolationDBSumPtCorr_.push_back(itau->tauID("byMediumCombinedIsolationDBSumPtCorr"));
  tauByTightCombinedIsolationDBSumPtCorr_.push_back(itau->tauID("byTightCombinedIsolationDBSumPtCorr"));
  tauByLooseCombinedIsolationDBSumPtCorr3Hits_.push_back(itau->tauID("byLooseCombinedIsolationDBSumPtCorr3Hits"));
  tauByMediumCombinedIsolationDBSumPtCorr3Hits_.push_back(itau->tauID("byMediumCombinedIsolationDBSumPtCorr3Hits"));
  tauByTightCombinedIsolationDBSumPtCorr3Hits_.push_back(itau->tauID("byTightCombinedIsolationDBSumPtCorr3Hits"));
  tauByVLooseIsolationMVA3newDMwoLT_.push_back(itau->tauID("byVLooseIsolationMVA3newDMwoLT"));
  tauByLooseIsolationMVA3newDMwoLT_.push_back(itau->tauID("byLooseIsolationMVA3newDMwoLT"));
  tauByMediumIsolationMVA3newDMwoLT_.push_back(itau->tauID("byMediumIsolationMVA3newDMwoLT"));
  tauByTightIsolationMVA3newDMwoLT_.push_back(itau->tauID("byTightIsolationMVA3newDMwoLT"));
  tauByVTightIsolationMVA3newDMwoLT_.push_back(itau->tauID("byVTightIsolationMVA3newDMwoLT"));
  tauByVVTightIsolationMVA3newDMwoLT_.push_back(itau->tauID("byVVTightIsolationMVA3newDMwoLT"));
  tauByIsolationMVA3newDMwoLTraw_.push_back(itau->tauID("byIsolationMVA3newDMwoLTraw"));
  tauByVLooseIsolationMVA3oldDMwLT_.push_back(itau->tauID("byVLooseIsolationMVA3oldDMwLT"));
  tauByLooseIsolationMVA3oldDMwLT_.push_back(itau->tauID("byLooseIsolationMVA3oldDMwLT"));
  tauByMediumIsolationMVA3oldDMwLT_.push_back(itau->tauID("byMediumIsolationMVA3oldDMwLT"));
  tauByTightIsolationMVA3oldDMwLT_.push_back(itau->tauID("byTightIsolationMVA3oldDMwLT"));
  tauByVTightIsolationMVA3oldDMwLT_.push_back(itau->tauID("byVTightIsolationMVA3oldDMwLT"));
  tauByVVTightIsolationMVA3oldDMwLT_.push_back(itau->tauID("byVVTightIsolationMVA3oldDMwLT"));
  tauByIsolationMVA3oldDMwLTraw_.push_back(itau->tauID("byIsolationMVA3oldDMwLTraw"));
  tauByVLooseIsolationMVA3oldDMwoLT_.push_back(itau->tauID("byVLooseIsolationMVA3oldDMwoLT"));
  tauByLooseIsolationMVA3oldDMwoLT_.push_back(itau->tauID("byLooseIsolationMVA3oldDMwoLT"));
  tauByTightIsolationMVA3oldDMwoLT_.push_back(itau->tauID("byTightIsolationMVA3oldDMwoLT"));
  tauByVTightIsolationMVA3oldDMwoLT_.push_back(itau->tauID("byVTightIsolationMVA3oldDMwoLT"));
  tauByVVTightIsolationMVA3oldDMwoLT_.push_back(itau->tauID("byVVTightIsolationMVA3oldDMwoLT"));
  tauByIsolationMVA3newDMwoLTraw_.push_back(itau->tauID("byIsolationMVA3newDMwoLTraw"));
  tauByLooseIsolationMVA3newDMwLT_.push_back(itau->tauID("byLooseIsolationMVA3newDMwLT"));
  tauByVLooseIsolationMVA3newDMwLT_.push_back(itau->tauID("byVLooseIsolationMVA3newDMwLT"));
  tauByMediumIsolationMVA3newDMwLT_.push_back(itau->tauID("byMediumIsolationMVA3newDMwLT"));
  tauByTightIsolationMVA3newDMwLT_.push_back(itau->tauID("byTightIsolationMVA3newDMwLT"));
  tauByVTightIsolationMVA3newDMwLT_.push_back(itau->tauID("byVTightIsolationMVA3newDMwLT"));
  tauByVVTightIsolationMVA3newDMwLT_.push_back(itau->tauID("byVVTightIsolationMVA3newDMwLT"));
  tauByIsolationMVA3newDMwLTraw_.push_back(itau->tauID("byIsolationMVA3newDMwLTraw"));






  tauEta_.push_back(itau->eta());
  tauPhi_.push_back(itau->phi());
  tauPt_ .push_back(itau->pt());
  tauEt_ .push_back(itau->et());
  tauCharge_.push_back(itau->charge());
  tauDecayMode_.push_back(itau->decayMode());
  tauEMFraction_.push_back(itau->emFraction());
  tauHCAL3x3OverPLead_.push_back(itau->hcal3x3OverPLead());
  tauHCALMaxOverPLead_.push_back(itau->hcalMaxOverPLead());
  tauHCALTotOverPLead_.push_back(itau->hcalTotOverPLead());
  tauIsolationPFChargedHadrCandsPtSum_.push_back(itau->isolationPFChargedHadrCandsPtSum());
  tauIsolationPFGammaCandsEtSum_.push_back(itau->isolationPFGammaCandsEtSum());
  tauLeadPFChargedHadrCandsignedSipt_.push_back(itau->leadPFChargedHadrCandsignedSipt());
  const reco::PFCandidatePtr& leadPFChargedHadrCand_Ref = itau->leadPFChargedHadrCand();
  if(leadPFChargedHadrCand_Ref.isNonnull()){// this check is needed in case hpsTau fails decayModeFinding
    tauLeadChargedHadronExists_.push_back(true);
    tauLeadChargedHadronEta_.push_back(leadPFChargedHadrCand_Ref->eta());
    tauLeadChargedHadronPhi_.push_back(leadPFChargedHadrCand_Ref->phi());
    tauLeadChargedHadronPt_.push_back(leadPFChargedHadrCand_Ref->pt());
  } else {
    tauLeadChargedHadronExists_.push_back(false);
    tauLeadChargedHadronEta_.push_back(0);
    tauLeadChargedHadronPhi_.push_back(0);
tauLeadChargedHadronPt_.push_back(0);
  }
  ++nTau_;
 } // loop over tau candidates
 } // tau handle is valid
 



  hEvents_->Fill(1.5);
  tree_->Fill();
  
  delete lazyTool;
  delete lazyToolnoZS;
  delete GEDIdTool;
}

void ggNtuplizer::beginJob() {
}

void ggNtuplizer::endJob() {
}

void ggNtuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

void ggNtuplizer::clearVectors() {

  pdf_.clear();
  nPU_.clear();
  puBX_.clear();
  puTrue_.clear();

  eleCharge_.clear();
  eleChargeConsistent_.clear();
  eleEn_.clear();
  eleSCEn_.clear();
  eleESEn_.clear();
  eleD0_.clear();
  eleDz_.clear();
  elePt_.clear();
  eleEta_.clear();
  eleTheta_.clear();
  elePhi_.clear();
  eleSCEta_.clear();
  eleSCPhi_.clear();
  eleSCRawEn_.clear();
  eleSCEtaWidth_.clear();
  eleSCPhiWidth_.clear();
  eleHoverE_.clear();
  eleEoverP_.clear();
  eleEoverPInv_.clear();
  eleBrem_.clear();
  eledEtaAtVtx_.clear();
  eledPhiAtVtx_.clear();
  eleSigmaIEtaIEta_.clear();
  eleSigmaIEtaIPhi_.clear();
  eleSigmaIPhiIPhi_.clear();
  eleSigmaIEtaIEta_2012_.clear();
  eleConvVeto_.clear();
  eleMissHits_.clear();
  eleESEffSigmaRR_.clear();
  elePFChIso_.clear();
  elePFPhoIso_.clear();
  elePFNeuIso_.clear();
  elePFPUIso_.clear();
  eleBC1E_.clear();
  eleBC1Eta_.clear();
  eleBC2E_.clear();
  eleBC2Eta_.clear();

  eleIDMVA_.clear();
  eledEtaseedAtVtx_.clear();
  eleE1x5_.clear();
  eleE2x5_.clear();
  eleE5x5_.clear();
  eleRelIsoWithDBeta_.clear();
  eleEcalDrivenSeed_.clear();
  eleDr03EcalRecHitSumEt_.clear();
  eleDr03HcalDepth1TowerSumEt_.clear();
  eleDr03HcalDepth2TowerSumEt_.clear();
  eleDr03HcalTowerSumEt_.clear();
  eleDr03TkSumPt_.clear();

  eleE1x5_2012_.clear();
  eleE2x5_2012_.clear();
  eleE5x5_2012_.clear();

  elecaloEnergy_.clear();
  eleTrkdxy_.clear();
  
  phoE_.clear();
  phoEt_.clear();
  phoEta_.clear();
  phoPhi_.clear();
  phoSCE_.clear();
  phoSCRawE_.clear();
  phoESEn_.clear();
  phoSCEta_.clear();
  phoSCPhi_.clear();
  phoSCEtaWidth_.clear();
  phoSCPhiWidth_.clear();
  phoSCBrem_.clear();
  phohasPixelSeed_.clear();
  phoEleVeto_.clear();
  phoR9_.clear();
  phoHoverE_.clear();
  phoSigmaIEtaIEta_.clear();
  phoSigmaIEtaIPhi_.clear();
  phoSigmaIPhiIPhi_.clear();
  phoE1x3_.clear();
  phoE2x2_.clear();
  phoE2x5Max_.clear();
  phoE5x5_.clear();
  phoESEffSigmaRR_.clear();
  phoSigmaIEtaIEta_2012_.clear();
  phoSigmaIEtaIPhi_2012_.clear();
  phoSigmaIPhiIPhi_2012_.clear();
  phoE1x3_2012_.clear();
  phoE2x2_2012_.clear();
  phoE2x5Max_2012_.clear();
  phoE5x5_2012_.clear();
  phoPFChIso_.clear();
  phoPFPhoIso_.clear();
  phoPFNeuIso_.clear();
  phoPFChWorstIso_.clear();
  phoPFChIsoFrix1_.clear();
  phoPFChIsoFrix2_.clear();
  phoPFChIsoFrix3_.clear();
  phoPFChIsoFrix4_.clear();
  phoPFChIsoFrix5_.clear();
  phoPFChIsoFrix6_.clear();
  phoPFChIsoFrix7_.clear();
  phoPFChIsoFrix8_.clear();
  phoPFPhoIsoFrix1_.clear();
  phoPFPhoIsoFrix2_.clear();
  phoPFPhoIsoFrix3_.clear();
  phoPFPhoIsoFrix4_.clear();
  phoPFPhoIsoFrix5_.clear();
  phoPFPhoIsoFrix6_.clear();
  phoPFPhoIsoFrix7_.clear();
  phoPFPhoIsoFrix8_.clear();
  phoPFNeuIsoFrix1_.clear();
  phoPFNeuIsoFrix2_.clear();
  phoPFNeuIsoFrix3_.clear();
  phoPFNeuIsoFrix4_.clear();
  phoPFNeuIsoFrix5_.clear();
  phoPFNeuIsoFrix6_.clear();
  phoPFNeuIsoFrix7_.clear();
  phoPFNeuIsoFrix8_.clear();
  phoBC1E_.clear();
  phoBC1Eta_.clear();
  phoBC2E_.clear();
  phoBC2Eta_.clear();

  phoIDMVA_.clear();

  phoEcalRecHitSumEtConeDR03_.clear();
  phohcalDepth1TowerSumEtConeDR03_.clear();
  phohcalDepth2TowerSumEtConeDR03_.clear();
  phohcalTowerSumEtConeDR03_.clear();
  photrkSumPtHollowConeDR03_.clear();

  muPt_.clear();
  muEta_.clear();
  muPhi_.clear();
  muCharge_.clear();
  muType_.clear();
  muIsGood_.clear();
  //muID_.clear();
  muD0_.clear();
  muDz_.clear();
  muChi2NDF_.clear();
  muInnerD0_.clear();
  muInnerDz_.clear();
  muTrkLayers_.clear();
  muPixelLayers_.clear();
  muPixelHits_.clear();
  muMuonHits_.clear();
  muStations_.clear();
  muTrkQuality_.clear();
  muIsoTrk_.clear();
  muPFChIso_.clear();
  muPFPhoIso_.clear();
  muPFNeuIso_.clear();
  muPFPUIso_.clear();
  ///SJ
  muInnervalidFraction_.clear();
  musegmentCompatibility_.clear();
  muchi2LocalPosition_.clear();
  mutrkKink_.clear();
  muBestTrkPtError_.clear();
  muBestTrkPt_.clear();

  //Lovedeep
  //@Lovedeep: 09/14 (Updated TauId Info)
  tauByLooseElectronRejection_.clear();
  tauByMediumElectronRejection_.clear();
  tauByTightElectronRejection_.clear();
  tauByMVA5LooseElectronRejection_.clear();
  tauByMVA5MediumElectronRejection_.clear();
  tauByMVA5TightElectronRejection_.clear();
  tauByMVA5VTightElectronRejection_.clear();
  tauByLooseMuonRejection_.clear();
  tauByMediumMuonRejection_.clear();
  tauByTightMuonRejection_.clear();
  tauByLooseMuonRejection3_.clear();
  tauByTightMuonRejection3_.clear();
  tauByMVALooseMuonRejection_.clear();
  tauByMVAMediumMuonRejection_.clear();
  tauByMVATightMuonRejection_.clear();
  tauByMVArawMuonRejection_.clear();
  pfTausDiscriminationByDecayModeFinding_.clear();
  tauByVLooseIsolation_.clear();
  tauByVLooseCombinedIsolationDBSumPtCorr_.clear();
  tauByLooseCombinedIsolationDBSumPtCorr_.clear();
  tauByMediumCombinedIsolationDBSumPtCorr_.clear();
  tauByTightCombinedIsolationDBSumPtCorr_.clear();
  tauByLooseCombinedIsolationDBSumPtCorr3Hits_.clear();
  tauByMediumCombinedIsolationDBSumPtCorr3Hits_.clear();
  tauByTightCombinedIsolationDBSumPtCorr3Hits_.clear();
  tauByVLooseIsolationMVA3newDMwoLT_.clear();
  tauByLooseIsolationMVA3newDMwoLT_.clear();
  tauByMediumIsolationMVA3newDMwoLT_.clear();
  tauByTightIsolationMVA3newDMwoLT_.clear();
  tauByVTightIsolationMVA3newDMwoLT_.clear();
  tauByVVTightIsolationMVA3newDMwoLT_.clear();
  tauByIsolationMVA3newDMwoLTraw_.clear();
  tauByVLooseIsolationMVA3oldDMwLT_.clear();
  tauByLooseIsolationMVA3oldDMwLT_.clear();
  tauByMediumIsolationMVA3oldDMwLT_.clear();
  tauByTightIsolationMVA3oldDMwLT_.clear();
  tauByVTightIsolationMVA3oldDMwLT_.clear();
  tauByVVTightIsolationMVA3oldDMwLT_.clear();
  tauByIsolationMVA3oldDMwLTraw_.clear();
  tauByVLooseIsolationMVA3oldDMwoLT_.clear();
  tauByLooseIsolationMVA3oldDMwoLT_.clear();
  tauByTightIsolationMVA3oldDMwoLT_.clear();
  tauByVTightIsolationMVA3oldDMwoLT_.clear();
  tauByVVTightIsolationMVA3oldDMwoLT_.clear();
  tauByIsolationMVA3newDMwoLTraw_.clear();
  tauByLooseIsolationMVA3newDMwLT_.clear();
  tauByVLooseIsolationMVA3newDMwLT_.clear();
  tauByMediumIsolationMVA3newDMwLT_.clear();
  tauByTightIsolationMVA3newDMwLT_.clear();
  tauByVTightIsolationMVA3newDMwLT_.clear();
  tauByVVTightIsolationMVA3newDMwLT_.clear();
  tauByIsolationMVA3newDMwLTraw_.clear();

  tauEta_.clear();
  tauPhi_.clear();
  tauPt_.clear();
  tauEt_.clear();
  tauCharge_.clear();
  tauDecayMode_.clear();
  tauEMFraction_.clear();
  tauHCAL3x3OverPLead_.clear();
  tauHCALMaxOverPLead_.clear();
  tauHCALTotOverPLead_.clear();
  tauIsolationPFChargedHadrCandsPtSum_.clear();
  tauIsolationPFGammaCandsEtSum_.clear();
  tauLeadPFChargedHadrCandsignedSipt_.clear();
  tauLeadChargedHadronExists_.clear();
  tauLeadChargedHadronEta_.clear();
  tauLeadChargedHadronPhi_.clear();
  tauLeadChargedHadronPt_.clear();






  // jets
  jetPt_.clear();
  jetEta_.clear();
  jetPhi_.clear();
  jetCHF_.clear();
  jetNHF_.clear();
  jetCEF_.clear();
  jetNEF_.clear();
  jetNCH_.clear();
  jetHFHAE_.clear();
  jetHFEME_.clear();
  jetNConstituents_.clear();
  jetCombinedSecondaryVtxBJetTags_.clear();
  jetJetProbabilityBJetTags_.clear();
  jetJetBProbabilityBJetTags_.clear();
  jetTrackCountingHighPurBJetTags_.clear();
  jetTrackCountingHighEffBJetTags_.clear();
  jetSimpleSecondaryVertexHighEffBJetTags_.clear();
  jetSimpleSecondaryVertexHighPurBJetTags_.clear();
  jetPartonID_.clear();
  jetPFLooseId_.clear();

// SubJet
  CA8JetPt_.clear();
  CA8JetEta_.clear();
  CA8JetPhi_.clear();
  CA8JetMass_.clear();
  CA8JetArea_.clear();
  CA8Jet_tau1_.clear();
  CA8Jet_tau2_.clear();
  CA8Jet_tau3_.clear();
  CA8JetCHF_.clear();
  CA8JetNHF_.clear();
  CA8JetCEF_.clear();
  CA8JetNEF_.clear();
  CA8JetNCH_.clear();
  CA8Jetnconstituents_.clear();
  CA8prunedJetMass_.clear();
}

DEFINE_FWK_MODULE(ggNtuplizer);
