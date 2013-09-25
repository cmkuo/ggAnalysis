#include "ggAnalysis/ggNtuplizer/interface/ggNtuplizer.h"
#include "ggAnalysis/ggNtuplizer/interface/GenParticleParentage.h"

#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "CondFormats/DataRecord/interface/GBRWrapperRcd.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDetId/interface/ESDetId.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HepMCCandidate/interface/PdfInfo.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/METReco/interface/BeamHaloSummary.h"
#include "DataFormats/METReco/interface/HcalNoiseSummary.h"
#include "DataFormats/MuonReco/interface/MuonCocktails.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/Scalers/interface/DcsStatus.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "EGamma/EGammaAnalysisTools/src/PFIsolationEstimator.cc"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/CaloGeometry/interface/TruncatedPyramid.h"
#include "Geometry/EcalAlgo/interface/EcalPreshowerGeometry.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "PFIsolation/SuperClusterFootprintRemoval/interface/SuperClusterFootprintRemoval.h"
#include "PhysicsTools/CandUtils/interface/AddFourMomenta.h"
#include "RecoCaloTools/Navigation/interface/EcalPreshowerNavigator.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "RecoEcal/EgammaCoreTools/interface/Mustache.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"
#include "RecoEgamma/EgammaIsolationAlgos/plugins/EgammaTowerIsolationProducer.h"
#include "RecoEgamma/EgammaTools/interface/ConversionFinder.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "RecoEgamma/EgammaTools/interface/ggPFPhotons.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "ggAnalysis/ggNtuplizer/interface/fCorrs.h"
#include "ggAnalysis/ggNtuplizer/interface/ggPFIsolation.h"
#include <iostream>
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"

bool compareMass(TLorentzVector l1, TLorentzVector l2){
  return (l1.M()>l2.M());
}

ggNtuplizer::ggNtuplizer(const edm::ParameterSet& ps) : verbosity_(0) {

  trgResults_        = ps.getParameter<InputTag>("triggerResults");
  trgEvent_          = ps.getParameter<InputTag>("triggerEvent"); 
  getBlocks_         = ps.getParameter<bool>("getBlocks");
  useAllPF_          = ps.getParameter<bool>("useAllPF");
  doGenParticles_    = ps.getParameter<bool>("doGenParticles");
  doSkim_            = ps.getParameter<bool>("doSkim");
  dumpESHits_        = ps.getParameter<bool>("dumpESHits");
  dumpESClusterInfo_ = ps.getParameter<bool>("dumpESClusterInfo");
  dumpTrks_          = ps.getParameter<bool>("dumpTrks"); 
  dumpJets_          = ps.getParameter<bool>("dumpJets");
  dumpSubJets_       = ps.getParameter<bool>("dumpSubJets");
  runOnParticleGun_  = ps.getParameter<bool>("runOnParticleGun");
  develop_           = ps.getParameter<bool>("development");
  doCentrality_      = ps.getParameter<bool>("doCentrality");
  
  vtxlabel_          = ps.getParameter<InputTag>("VtxLabel");
  tracklabel_        = ps.getParameter<InputTag>("TrackLabel");
  gsfElectronlabel_  = ps.getParameter<InputTag>("gsfElectronLabel");
  pfMETlabel_        = ps.getParameter<InputTag>("pfMETLabel");
  pfType01METlabel_  = ps.getParameter<InputTag>("pfType01METLabel");
  recoPfMETlabel_    = ps.getParameter<InputTag>("recoPfMETLabel");

  genParticlesCollection_    = ps.getParameter<InputTag>("genParticleSrc");
  METCollection_             = ps.getParameter<InputTag>("METSrc");
  electronCollection_        = ps.getParameter<InputTag>("electronSrc");
  photonCollection_          = ps.getParameter<InputTag>("photonSrc");
  recophotonCollection_      = ps.getParameter<InputTag>("Photons");
  muonCollection_            = ps.getParameter<InputTag>("muonSrc");
  tauCollection_             = ps.getParameter<InputTag>("tauSrc");
  jetCollection_             = ps.getParameter<InputTag>("jetSrc");
  ebReducedRecHitCollection_ = ps.getParameter<InputTag>("ebReducedRecHitCollection");
  eeReducedRecHitCollection_ = ps.getParameter<InputTag>("eeReducedRecHitCollection");
  esRecHitCollection_        = ps.getParameter<InputTag>("esRecHitCollection");
  towerCollection_           = ps.getParameter<InputTag>("towerCollection");
  beamSpotCollection_        = ps.getParameter<InputTag>("BeamSpotCollection");
  puCollection_              = ps.getParameter<InputTag>("pileupCollection");
  
  rhoLepPFisoCollection_     = ps.getParameter<InputTag>("rhoCollectionLepPFiso");
  rhoCollection25_           = ps.getParameter<InputTag>("rhoCollection25");
  rhoCollection25_neu_       = ps.getParameter<InputTag>("rhoCollection25_neu");
  rhoCollection25_eleLabel_  = edm::InputTag("kt6PFJetsCentralChargedPileUp","rho");  
  rhoCollection44_           = ps.getParameter<InputTag>("rhoCollection44");
  rho2011Label_		     = ps.getParameter<InputTag>("rho2011Label");
  rho2012Label_              = edm::InputTag("kt6PFJets","rho");
  pfAllParticles_            = ps.getParameter<InputTag>("PFAllCandidates");
  pfParticles_               = ps.getParameter<InputTag>("PFCandidates");
  pfPhotonParticles_         = ps.getParameter<InputTag>("PFCandidatePhotons");
  allConversionsColl_        = ps.getParameter<InputTag>("ConvertedPhotonColl");
  pfPhotonCollection_        = ps.getParameter<InputTag>("PFPhotons");
  generatorLabel_            = edm::InputTag("generator");
  HBHENoiseFilterLabel_      = edm::InputTag("HBHENoiseFilterResultProducer","HBHENoiseFilterResult");
  HcalLaserFilterLabel_      = edm::InputTag("hcalLaserEventFilter");
  EcalDeadCellFilterLabel_   = edm::InputTag("EcalDeadCellTriggerPrimitiveFilter");
  TrackingFailureFilterLabel_ = edm::InputTag("trackingFailureFilter");
  EEBadScFilterLabel_        = edm::InputTag("eeBadScFilter");
  EcalLaserFilterLabel_      = edm::InputTag("ecalLaserCorrFilter");
  Manystripclus53XLabel_     = edm::InputTag("manystripclus53X");
  Toomanystripclus53XLabel_  = edm::InputTag("toomanystripclus53X");
  LogErrorTooManyClustersLabel_ = edm::InputTag("logErrorTooManyClusters");
  muonNoCutsLabel_           = edm::InputTag("muons");
  eleNoCutsLabel_            = edm::InputTag("gsfElectrons");
  recVtxsBSLabel_            = edm::InputTag("offlinePrimaryVerticesWithBS","");

  QGTagsHandleMLPLabel_      = edm::InputTag("QGTagger","qgMLP");
  QGTagsHandleLikelihoodLabel_ = edm::InputTag("QGTagger","qgLikelihood");
  QGtagjetLabel_             = edm::InputTag("selectedPatJetsAK5PF");
  jetsCHSprunedLabel_        = edm::InputTag("selectedPatJetsCA8CHSpruned");
  jetsCHSLabel_              = edm::InputTag("selectedPatJetsCA8CHSwithNsub");

  jetMVAAlgos_               = ps.getUntrackedParameter<std::vector<edm::ParameterSet> >("puJetIDAlgos");
  pfLooseId_                 = ps.getParameter<edm::ParameterSet>("pfLooseId");
  
  // jetMVAs_.resize(jetMVAAlgos_.size());
  // jetWPLevels_.resize(jetMVAAlgos_.size());
  // jetMVAsExt_.resize(jetMVAAlgos_.size());
  // jetWPLevelsExt_.resize(jetMVAAlgos_.size());
  pujetIDalgos_.resize(jetMVAAlgos_.size());
  for(unsigned int imva=0; imva<jetMVAAlgos_.size(); ++imva){
    pujetIDalgos_[imva] = new PileupJetIdAlgo((jetMVAAlgos_.at(imva)));
  }

  inputTagPhotonIsoDeposits_ = ps.getParameter<std::vector<edm::InputTag> >("PhotonIsoDeposits");

  // PF isolation
  inputTagIsoDepElectrons_     = ps.getParameter< std::vector<edm::InputTag> >("IsoDepElectron");
  inputTagIsoDepPhotons_       = ps.getParameter< std::vector<edm::InputTag> >("IsoDepPhoton");
  inputTagIsoValElectronsPFId_ = ps.getParameter< std::vector<edm::InputTag> >("IsoValElectronPF");
  inputTagIsoValPhotonsPFId_   = ps.getParameter< std::vector<edm::InputTag> >("IsoValPhoton");

  // SC footprint remover parameters
  scRemover03Pset_ = ps.getParameter<edm::ParameterSet>("scRemover03");
  scRemover04Pset_ = ps.getParameter<edm::ParameterSet>("scRemover04");

  cicPhotonId_ = new CiCPhotonID(ps);
  trackMET_ = new trackMET(ps);

  Service<TFileService> fs;
  hEvents_ = fs->make<TH1F>("hEvents", "total processed and skimmed events",   2, 0,   2);
  hPU_     = fs->make<TH1F>("hPU",     "number of pileup",                   100, 0, 100);
  hPUTrue_ = fs->make<TH1F>("hPUTrue", "number of true pilepu"             , 500, 0, 100);
  tree_    = fs->make<TTree>("EventTree", "Event data");

  tree_->Branch("run", &run_, "run/I");
  tree_->Branch("event", &event_, "event/L");
  tree_->Branch("lumis", &lumis_, "lumis/I");
  tree_->Branch("isData", &isData_, "isData/O");
  tree_->Branch("nHLT", &nHLT_, "nHLT/I");
  tree_->Branch("HLT", HLT_, "HLT[nHLT]/I");
  tree_->Branch("HLTIndex", HLTIndex_, "HLTIndex[70]/I");
  tree_->Branch("bspotPos", bspotPos_, "bspotPos[3]/F");
  tree_->Branch("nVtx", &nVtx_, "nVtx/I");
  tree_->Branch("vtx_x", &vtx_x_);
  tree_->Branch("vtx_y", &vtx_y_);
  tree_->Branch("vtx_z", &vtx_z_);
  tree_->Branch("IsVtxGood", &IsVtxGood_, "IsVtxGood/I");
  tree_->Branch("nGoodVtx", &nGoodVtx_, "nGoodVtx/I");
  if (doCentrality_) tree_->Branch("centrality", centrality_, "centrality[5]/F");
  tree_->Branch("nVtxBS", &nVtxBS_, "nVtxBS/I");
  tree_->Branch("vtxbs_x", &vtxbs_x_);
  tree_->Branch("vtxbs_y", &vtxbs_y_);
  tree_->Branch("vtxbs_z", &vtxbs_z_);
  if (dumpTrks_) {
    tree_->Branch("vtxbsPtMod", &vtxbsPtMod_);
    tree_->Branch("vtxbsSumPt2", &vtxbsSumPt2_);
    tree_->Branch("vtxbsTkIndex", &vtxbsTkIndex_);
    tree_->Branch("vtxbsTkWeight", &vtxbsTkWeight_);
    tree_->Branch("nTrk", &nTrk_, "nTrk/I");
    tree_->Branch("trkP_x", &trkP_x_);
    tree_->Branch("trkP_y", &trkP_y_);
    tree_->Branch("trkP_z", &trkP_z_);
    tree_->Branch("trkVtx_x", &trkVtx_x_);
    tree_->Branch("trkVtx_y", &trkVtx_y_);
    tree_->Branch("trkVtx_z", &trkVtx_z_);
    tree_->Branch("trkd0", &trkd0_);
    tree_->Branch("trkd0Err", &trkd0Err_);
    tree_->Branch("trkdz", &trkdz_);
    tree_->Branch("trkdzErr", &trkdzErr_);
    tree_->Branch("trkPtErr", &trkPtErr_);
    tree_->Branch("trkQuality", &trkQuality_);
    tree_->Branch("nGoodTrk", &nGoodTrk_, "nGoodTrk/I");
    tree_->Branch("IsTracksGood", &IsTracksGood_, "IsTracksGood/I");
  }
  if (doGenParticles_) {
    tree_->Branch("pdf", pdf_, "pdf[7]/F");
    tree_->Branch("pthat", &pthat_, "pthat/F");
    tree_->Branch("processID", &processID_, "processID/F");
    // genParticle
    tree_->Branch("nMC", &nMC_, "nMC/I");
    tree_->Branch("mcPID", &mcPID);
    tree_->Branch("mcVtx_x", &mcVtx_x);
    tree_->Branch("mcVtx_y", &mcVtx_y);
    tree_->Branch("mcVtx_z", &mcVtx_z);
    tree_->Branch("mcPt", &mcPt);
    tree_->Branch("mcMass", &mcMass);
    tree_->Branch("mcEta", &mcEta);
    tree_->Branch("mcPhi", &mcPhi);
    tree_->Branch("mcE", &mcE);
    tree_->Branch("mcEt", &mcEt);
    tree_->Branch("mcGMomPID", &mcGMomPID);
    tree_->Branch("mcMomPID", &mcMomPID);
    tree_->Branch("mcMomPt", &mcMomPt);
    tree_->Branch("mcMomMass", &mcMomMass);
    tree_->Branch("mcMomEta", &mcMomEta);
    tree_->Branch("mcMomPhi", &mcMomPhi);
    tree_->Branch("mcIndex", &mcIndex);
    tree_->Branch("mcDecayType", &mcDecayType); //-999:non W or Z, 1:hardronic, 2:e, 3:mu, 4:tau
    tree_->Branch("mcParentage", &mcParentage); // 16*lepton + 8*boson + 4*non-prompt + 2*qcd + exotics
    tree_->Branch("mcStatus", &mcStatus); // status of the particle
    // Gen MET
    tree_->Branch("genMET", &genMET_, "genMET/F");
    tree_->Branch("genMETPhi", &genMETPhi_, "genMETPhi/F");
    // PU Info
    tree_->Branch("nPUInfo", &nPUInfo_, "nPUInfo/I");
    tree_->Branch("nPU", &nPU_);
    tree_->Branch("puBX", &puBX_);
    tree_->Branch("puTrue", &puTrue_);
  }
  // pfMET Type1
  tree_->Branch("pfMET", &pfMET_, "pfMET/F");
  tree_->Branch("pfMETPhi", &pfMETPhi_, "pfMETPhi/F");
  tree_->Branch("pfMETsumEt", &pfMETsumEt_, "pfMETsumEt/F");
  tree_->Branch("pfMETmEtSig", &pfMETmEtSig_, "pfMETmEtSig/F");
  tree_->Branch("pfMETSig", &pfMETSig_, "pfMETSig/F");
  // pfMET Type0+1
  tree_->Branch("pfType01MET",       &pfType01MET_,       "pfType01MET/F");
  tree_->Branch("pfType01METPhi",    &pfType01METPhi_,    "pfType01METPhi/F");
  tree_->Branch("pfType01METsumEt",  &pfType01METsumEt_,  "pfType01METsumEt/F");
  tree_->Branch("pfType01METmEtSig", &pfType01METmEtSig_, "pfType01METmEtSig/F");
  tree_->Branch("pfType01METSig",    &pfType01METSig_,    "pfType01METSig/F");
  // reco pfMET
  tree_->Branch("recoPfMET", &recoPfMET_, "recoPfMET/F");
  tree_->Branch("recoPfMETPhi", &recoPfMETPhi_, "recoPfMETPhi/F");
  tree_->Branch("recoPfMETsumEt", &recoPfMETsumEt_, "recoPfMETsumEt/F");
  tree_->Branch("recoPfMETmEtSig", &recoPfMETmEtSig_, "recoPfMETmEtSig/F");
  tree_->Branch("recoPfMETSig", &recoPfMETSig_, "recoPfMETSig/F");
  // track MET
  tree_->Branch("trkMETxPV",   &trkMETxPV_,   "trkMETxPV/F");
  tree_->Branch("trkMETyPV",   &trkMETyPV_,   "trkMETyPV/F");
  tree_->Branch("trkMETPhiPV", &trkMETPhiPV_, "trkMETPhiPV/F");
  tree_->Branch("trkMETPV",    &trkMETPV_,    "trkMETPV/F");
  tree_->Branch("trkMETx",     &trkMETx_);
  tree_->Branch("trkMETy",     &trkMETy_);
  tree_->Branch("trkMETPhi",   &trkMETPhi_);
  tree_->Branch("trkMET",      &trkMET_);
  // MET filters
  tree_->Branch("metFilters", metFilters_, "metFilters[10]/I");
  // Electron
  tree_->Branch("nEle", &nEle_, "nEle/I");
  tree_->Branch("eleTrg", &eleTrg_);
  tree_->Branch("eleClass", &eleClass_);
  tree_->Branch("eleIsEcalDriven", &eleIsEcalDriven_);
  tree_->Branch("eleCharge", &eleCharge_);
  tree_->Branch("eleChargeConsistent", &eleChargeConsistent_);
  tree_->Branch("eleEn", &eleEn_);
  tree_->Branch("eleEcalEn", &eleEcalEn_);
  tree_->Branch("eleSCRawEn", &eleSCRawEn_);
  tree_->Branch("eleSCEn", &eleSCEn_);
  tree_->Branch("eleESEn", &eleESEn_);
  tree_->Branch("elePt", &elePt_);
  tree_->Branch("eleEta", &eleEta_);
  tree_->Branch("elePhi", &elePhi_);
  tree_->Branch("eleR9", &eleR9_);
  tree_->Branch("eleEtaVtx", &eleEtaVtx_);
  tree_->Branch("elePhiVtx", &elePhiVtx_);
  tree_->Branch("eleEtVtx", &eleEtVtx_);
  tree_->Branch("eleSCEta", &eleSCEta_);
  tree_->Branch("eleSCPhi", &eleSCPhi_);
  tree_->Branch("eleSCEtaWidth", &eleSCEtaWidth_);
  tree_->Branch("eleSCPhiWidth", &eleSCPhiWidth_);
  tree_->Branch("eleVtx_x", &eleVtx_x_);
  tree_->Branch("eleVtx_y", &eleVtx_y_);
  tree_->Branch("eleVtx_z", &eleVtx_z_);
  tree_->Branch("eleD0", &eleD0_);
  tree_->Branch("eleDz", &eleDz_);
  tree_->Branch("eleD0GV",   &eleD0GV_);
  tree_->Branch("eleDzGV",   &eleDzGV_);
  tree_->Branch("eleD0Vtx",  &eleD0Vtx_);
  tree_->Branch("eleDzVtx",  &eleDzVtx_);
  tree_->Branch("eleHoverE", &eleHoverE_);
  tree_->Branch("eleHoverE12", &eleHoverE12_);
  tree_->Branch("eleEoverP", &eleEoverP_);
  tree_->Branch("elePin", &elePin_);
  tree_->Branch("elePout", &elePout_);
  tree_->Branch("eleTrkMomErr", &eleTrkMomErr_);
  tree_->Branch("eleBrem", &eleBrem_);
  tree_->Branch("eledEtaAtVtx", &eledEtaAtVtx_);
  tree_->Branch("eledPhiAtVtx", &eledPhiAtVtx_);
  tree_->Branch("eleSigmaIEtaIEta", &eleSigmaIEtaIEta_);
  tree_->Branch("eleSigmaIEtaIPhi", &eleSigmaIEtaIPhi_);
  tree_->Branch("eleSigmaIPhiIPhi", &eleSigmaIPhiIPhi_);
  tree_->Branch("eleEmax", &eleEmax_);
  tree_->Branch("eleE2ndMax", &eleE2ndMax_);
  tree_->Branch("eleETop", &eleETop_);
  tree_->Branch("eleEBottom", &eleEBottom_);
  tree_->Branch("eleELeft", &eleELeft_);
  tree_->Branch("eleERight", &eleERight_);  
  tree_->Branch("eleE1x5", &eleE1x5_);
  tree_->Branch("eleE3x3", &eleE3x3_);
  tree_->Branch("eleE5x5", &eleE5x5_);
  tree_->Branch("eleE2x5Max", &eleE2x5Max_);
  tree_->Branch("eleE2x5Top", &eleE2x5Top_);
  tree_->Branch("eleE2x5Bottom", &eleE2x5Bottom_);
  tree_->Branch("eleE2x5Left", &eleE2x5Left_);
  tree_->Branch("eleE2x5Right", &eleE2x5Right_);
  tree_->Branch("eleSeedEta", &eleSeedEta_);
  tree_->Branch("eleSeedE", &eleSeedE_);
  tree_->Branch("eleSeedPhi", &eleSeedPhi_);
  tree_->Branch("eleCrysEta", &eleCrysEta_);
  tree_->Branch("eleCrysPhi", &eleCrysPhi_);
  tree_->Branch("eleCrysIEta", &eleCrysIEta_);
  tree_->Branch("eleCrysIPhi", &eleCrysIPhi_);
  tree_->Branch("eleRegrE", &eleRegrE_);
  tree_->Branch("eleRegrEerr", &eleRegrEerr_);
  tree_->Branch("elePhoRegrE", &elePhoRegrE_);
  tree_->Branch("elePhoRegrEerr", &elePhoRegrEerr_);
  tree_->Branch("eleSeedTime", &eleSeedTime_);
  // If Flag == 2, it means that rechit is out of time
  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/EcalFirstBeam09Anomalous#Spike_identification_in_collision
  tree_->Branch("eleRecoFlag", &eleRecoFlag_);
  tree_->Branch("elePos", &elePos_);
  if (doGenParticles_) {
    tree_->Branch("eleGenIndex", &eleGenIndex_);
    tree_->Branch("eleGenGMomPID", &eleGenGMomPID_);
    tree_->Branch("eleGenMomPID", &eleGenMomPID_);
    tree_->Branch("eleGenMomPt", &eleGenMomPt_);
  }
  tree_->Branch("eleIsoTrkDR03", &eleIsoTrkDR03_);
  tree_->Branch("eleIsoEcalDR03", &eleIsoEcalDR03_);
  tree_->Branch("eleIsoHcalDR03", &eleIsoHcalDR03_);
  tree_->Branch("eleIsoHcalDR0312", &eleIsoHcalDR0312_);
  tree_->Branch("eleIsoTrkDR04", &eleIsoTrkDR04_);
  tree_->Branch("eleIsoEcalDR04", &eleIsoEcalDR04_);
  tree_->Branch("eleIsoHcalDR04", &eleIsoHcalDR04_);
  tree_->Branch("eleIsoHcalDR0412", &eleIsoHcalDR0412_);
  tree_->Branch("eleModIsoTrk", &eleModIsoTrk_);
  tree_->Branch("eleModIsoEcal", &eleModIsoEcal_);
  tree_->Branch("eleModIsoHcal", &eleModIsoHcal_);
  tree_->Branch("eleMissHits", &eleMissHits_);
  tree_->Branch("eleConvDist", &eleConvDist_);
  tree_->Branch("eleConvDcot", &eleConvDcot_);
  tree_->Branch("eleConvVtxFit", &eleConvVtxFit_);
  tree_->Branch("eleIP3D", &eleIP3D_);
  tree_->Branch("eleIP3DErr", &eleIP3DErr_);
  tree_->Branch("eleIDMVANonTrig", &eleIDMVANonTrig_);
  tree_->Branch("eleIDMVATrig", &eleIDMVATrig_);
  tree_->Branch("elePFChIso03", &elePFChIso03_);
  tree_->Branch("elePFPhoIso03", &elePFPhoIso03_);
  tree_->Branch("elePFNeuIso03", &elePFNeuIso03_);
  tree_->Branch("elePFChIso04", &elePFChIso04_);
  tree_->Branch("elePFPhoIso04", &elePFPhoIso04_);
  tree_->Branch("elePFNeuIso04", &elePFNeuIso04_);
  tree_->Branch("eleESEffSigmaRR_x", &eleESEffSigmaRR_x_);
  tree_->Branch("eleESEffSigmaRR_y", &eleESEffSigmaRR_y_);
  tree_->Branch("eleESEffSigmaRR_z", &eleESEffSigmaRR_z_);
  if (dumpESHits_) {
    tree_->Branch("eleESDetId", eleESDetId_, "eleESDetId_[nEle][2]/I" );
    tree_->Branch("eleESHits",  eleESHits_,  "eleESHits_[nEle][3][62]/F");
  }
  if (dumpESClusterInfo_) {
    tree_->Branch("eleESE1",  eleESE1_,  "eleESE1_[nEle][2]/F");
    tree_->Branch("eleESE3",  eleESE3_,  "eleESE3_[nEle][2]/F");
    tree_->Branch("eleESE5",  eleESE5_,  "eleESE5_[nEle][2]/F");
    tree_->Branch("eleESE7",  eleESE7_,  "eleESE7_[nEle][2]/F");
    tree_->Branch("eleESE11", eleESE11_, "eleESE11_[nEle][2]/F");
    tree_->Branch("eleESE21", eleESE21_, "eleESE21_[nEle][2]/F");
  }
  if (develop_) {
    tree_->Branch("eleNBC", &eleNBC_);
    tree_->Branch("eleBrLinear"  , &eleBrLinear_  );  
    tree_->Branch("eleCetaCorrE" , &eleCetaCorrE_ );
    tree_->Branch("eleCetaCorrEt", &eleCetaCorrEt_);
    tree_->Branch("eleBremCorrE" , &eleBremCorrE_ );
    tree_->Branch("eleBremCorrEt", &eleBremCorrEt_);
    tree_->Branch("eleFullCorrE" , &eleFullCorrE_ );
    tree_->Branch("eleFullCorrEt", &eleFullCorrEt_);  
  }
  // Photon
  tree_->Branch("nPho", &nPho_, "nPho/I");
  tree_->Branch("phoTrg", &phoTrg_);
  tree_->Branch("phoTrgFilter", &phoTrgFilter_);
  tree_->Branch("phoIsPhoton", &phoIsPhoton_);
  tree_->Branch("phoSCPos_x", &phoSCPos_x_);
  tree_->Branch("phoSCPos_y", &phoSCPos_y_);
  tree_->Branch("phoSCPos_z", &phoSCPos_z_);
  tree_->Branch("phoCaloPos_x", &phoCaloPos_x_);
  tree_->Branch("phoCaloPos_y", &phoCaloPos_y_);
  tree_->Branch("phoCaloPos_z", &phoCaloPos_z_);

  tree_->Branch("phoE", &phoE_);
  tree_->Branch("phoEt", &phoEt_);
  tree_->Branch("phoEta", &phoEta_);
  tree_->Branch("phoVtx_x", &phoVtx_x_);
  tree_->Branch("phoVtx_y", &phoVtx_y_);
  tree_->Branch("phoVtx_z", &phoVtx_z_);
  tree_->Branch("phoPhi", &phoPhi_);
  tree_->Branch("phoEtVtx", &phoEtVtx_);
  tree_->Branch("phoEtaVtx", &phoEtaVtx_);
  tree_->Branch("phoPhiVtx", &phoPhiVtx_);
  tree_->Branch("phoR9", &phoR9_);
  tree_->Branch("phoNClus", &phoNClus_);
  if (develop_) {
    tree_->Branch("phoCetaCorrE" , &phoCetaCorrE_ );
    tree_->Branch("phoCetaCorrEt", &phoCetaCorrEt_);
    tree_->Branch("phoBremCorrE" , &phoBremCorrE_ );
    tree_->Branch("phoBremCorrEt", &phoBremCorrEt_);
    tree_->Branch("phoFullCorrE" , &phoFullCorrE_ );
    tree_->Branch("phoFullCorrEt", &phoFullCorrEt_);
  }

  //tree_->Branch("phoTrkIsoSolidDR03", &phoTrkIsoSolidDR03_);
  tree_->Branch("phoTrkIsoHollowDR03", &phoTrkIsoHollowDR03_);
  tree_->Branch("phoEcalIsoDR03", &phoEcalIsoDR03_);
  tree_->Branch("phoHcalIsoDR03", &phoHcalIsoDR03_);
  tree_->Branch("phoHcalIsoDR0312", &phoHcalIsoDR0312_);
  //tree_->Branch("phoHcalIsoSolidDR03", &phoHcalIsoSolidDR03_);
  //tree_->Branch("phoTrkIsoSolidDR04", &phoTrkIsoSolidDR04_);
  tree_->Branch("phoTrkIsoHollowDR04", &phoTrkIsoHollowDR04_);
  //tree_->Branch("phoCiCTrkIsoDR03", &phoCiCTrkIsoDR03_);
  //tree_->Branch("phoCiCTrkIsoDR04", &phoCiCTrkIsoDR04_);
  tree_->Branch("phoCiCdRtoTrk", &phoCiCdRtoTrk_);
  tree_->Branch("phoEcalIsoDR04", &phoEcalIsoDR04_);
  tree_->Branch("phoHcalIsoDR04", &phoHcalIsoDR04_);
  tree_->Branch("phoHcalIsoDR0412", &phoHcalIsoDR0412_);
  //tree_->Branch("phoHcalIsoSolidDR04", &phoHcalIsoSolidDR04_);
  tree_->Branch("phoHoverE", &phoHoverE_);
  tree_->Branch("phoHoverE12", &phoHoverE12_);
  tree_->Branch("phoEleVeto", &phoEleVeto_);
  tree_->Branch("phoSigmaIEtaIEta", &phoSigmaIEtaIEta_);
  tree_->Branch("phoSigmaIEtaIPhi", &phoSigmaIEtaIPhi_);
  tree_->Branch("phoSigmaIPhiIPhi", &phoSigmaIPhiIPhi_);
  if (dumpESClusterInfo_) {
    tree_->Branch("phoCiCPF4phopfIso005", &phoCiCPF4phopfIso005_);
    tree_->Branch("phoCiCPF4phopfIso01", &phoCiCPF4phopfIso01_);
    tree_->Branch("phoCiCPF4phopfIso02", &phoCiCPF4phopfIso02_);
  }
  tree_->Branch("phoCiCPF4phopfIso03", &phoCiCPF4phopfIso03_);
  tree_->Branch("phoCiCPF4phopfIso04", &phoCiCPF4phopfIso04_);
  if (dumpESClusterInfo_) {
    tree_->Branch("phoCiCPF4phopfIso05", &phoCiCPF4phopfIso05_);
    tree_->Branch("phoCiCPF4phopfIso06", &phoCiCPF4phopfIso06_);
    tree_->Branch("phoCiCPF4phopfIso07", &phoCiCPF4phopfIso07_);
    tree_->Branch("phoCiCPF4phopfIso08", &phoCiCPF4phopfIso08_);
    tree_->Branch("phoCiCPF4chgpfIso005", &phoCiCPF4chgpfIso005_);
    tree_->Branch("phoCiCPF4chgpfIso01", &phoCiCPF4chgpfIso01_);
  }
  tree_->Branch("phoCiCPF4chgpfIso02", &phoCiCPF4chgpfIso02_);
  tree_->Branch("phoCiCPF4chgpfIso03", &phoCiCPF4chgpfIso03_);
  tree_->Branch("phoCiCPF4chgpfIso04", &phoCiCPF4chgpfIso04_);
  if (dumpESClusterInfo_) {
    tree_->Branch("phoCiCPF4chgpfIso05", &phoCiCPF4chgpfIso05_);
    tree_->Branch("phoCiCPF4chgpfIso06", &phoCiCPF4chgpfIso06_);
    tree_->Branch("phoCiCPF4chgpfIso07", &phoCiCPF4chgpfIso07_);
    tree_->Branch("phoCiCPF4chgpfIso08", &phoCiCPF4chgpfIso08_);
    tree_->Branch("phoCiCPF4phopfIsoNoVETO005", &phoCiCPF4phopfIsoNoVETO005_);
    tree_->Branch("phoCiCPF4phopfIsoNoVETO01", &phoCiCPF4phopfIsoNoVETO01_);
    tree_->Branch("phoCiCPF4phopfIsoNoVETO02", &phoCiCPF4phopfIsoNoVETO02_);
    tree_->Branch("phoCiCPF4phopfIsoNoVETO03", &phoCiCPF4phopfIsoNoVETO03_);
    tree_->Branch("phoCiCPF4phopfIsoNoVETO04", &phoCiCPF4phopfIsoNoVETO04_);
    tree_->Branch("phoCiCPF4phopfIsoNoVETO05", &phoCiCPF4phopfIsoNoVETO05_);
    tree_->Branch("phoCiCPF4phopfIsoNoVETO06", &phoCiCPF4phopfIsoNoVETO06_);
    tree_->Branch("phoCiCPF4phopfIsoNoVETO07", &phoCiCPF4phopfIsoNoVETO07_);
    tree_->Branch("phoCiCPF4phopfIsoNoVETO08", &phoCiCPF4phopfIsoNoVETO08_);
    tree_->Branch("phoCiCPF4chgpfIsoNoVETO005", &phoCiCPF4chgpfIsoNoVETO005_);
    tree_->Branch("phoCiCPF4chgpfIsoNoVETO01", &phoCiCPF4chgpfIsoNoVETO01_);
    tree_->Branch("phoCiCPF4chgpfIsoNoVETO02", &phoCiCPF4chgpfIsoNoVETO02_);
    tree_->Branch("phoCiCPF4chgpfIsoNoVETO03", &phoCiCPF4chgpfIsoNoVETO03_);
    tree_->Branch("phoCiCPF4chgpfIsoNoVETO04", &phoCiCPF4chgpfIsoNoVETO04_);
    tree_->Branch("phoCiCPF4chgpfIsoNoVETO05", &phoCiCPF4chgpfIsoNoVETO05_);
    tree_->Branch("phoCiCPF4chgpfIsoNoVETO06", &phoCiCPF4chgpfIsoNoVETO06_);
    tree_->Branch("phoCiCPF4chgpfIsoNoVETO07", &phoCiCPF4chgpfIsoNoVETO07_);
    tree_->Branch("phoCiCPF4chgpfIsoNoVETO08", &phoCiCPF4chgpfIsoNoVETO08_);
  }
  tree_->Branch("phoEmax", &phoEmax_);
  tree_->Branch("phoETop", &phoETop_);
  tree_->Branch("phoEBottom", &phoEBottom_);
  tree_->Branch("phoELeft", &phoELeft_);
  tree_->Branch("phoERight", &phoERight_);
  tree_->Branch("phoE2ndMax", &phoE2ndMax_);
  tree_->Branch("phoE3x3", &phoE3x3_);
  tree_->Branch("phoE3x1", &phoE3x1_);
  tree_->Branch("phoE1x3", &phoE1x3_);
  tree_->Branch("phoE5x5", &phoE5x5_);
  tree_->Branch("phoE1x5", &phoE1x5_);
  tree_->Branch("phoE2x2", &phoE2x2_);
  tree_->Branch("phoE2x5Max", &phoE2x5Max_);
  tree_->Branch("phoE2x5Top", &phoE2x5Top_);
  tree_->Branch("phoE2x5Bottom", &phoE2x5Bottom_);
  tree_->Branch("phoE2x5Left", &phoE2x5Left_);
  tree_->Branch("phoE2x5Right", &phoE2x5Right_);
  tree_->Branch("phoSeedE", &phoSeedE_);
  tree_->Branch("phoSeedEta", &phoSeedEta_);
  tree_->Branch("phoSeedPhi", &phoSeedPhi_);
  tree_->Branch("phoCrysEta", &phoCrysEta_);
  tree_->Branch("phoCrysPhi", &phoCrysPhi_);
  tree_->Branch("phoCrysEta", &phoCrysEta_);
  tree_->Branch("phoCrysIEta", &phoCrysIEta_);
  tree_->Branch("phoCrysIPhi", &phoCrysIPhi_);
  tree_->Branch("phoPFChIso", &phoPFChIso_);
  tree_->Branch("phoPFPhoIso", &phoPFPhoIso_);
  tree_->Branch("phoPFNeuIso", &phoPFNeuIso_);
  tree_->Branch("phoSCRChIso", &phoSCRChIso_);
  tree_->Branch("phoSCRPhoIso", &phoSCRPhoIso_);
  tree_->Branch("phoSCRNeuIso", &phoSCRNeuIso_);
  tree_->Branch("phoSCRChIso04", &phoSCRChIso04_);
  tree_->Branch("phoSCRPhoIso04", &phoSCRPhoIso04_);
  tree_->Branch("phoSCRNeuIso04", &phoSCRNeuIso04_);
  tree_->Branch("phoRandConeChIso", &phoRandConeChIso_);
  tree_->Branch("phoRandConePhoIso", &phoRandConePhoIso_);
  tree_->Branch("phoRandConeNeuIso", &phoRandConeNeuIso_);
  tree_->Branch("phoRandConeChIso04", &phoRandConeChIso04_);
  tree_->Branch("phoRandConePhoIso04", &phoRandConePhoIso04_);
  tree_->Branch("phoRandConeNeuIso04", &phoRandConeNeuIso04_);
  tree_->Branch("phoRegrE", &phoRegrE_);
  tree_->Branch("phoRegrEerr", &phoRegrEerr_);
  tree_->Branch("phoSeedTime", &phoSeedTime_);
  tree_->Branch("phoSeedDetId1", &phoSeedDetId1_);
  tree_->Branch("phoSeedDetId2", &phoSeedDetId2_);
  tree_->Branch("phoLICTD", &phoLICTD_);
  // If Flag == 2, it means that rechit is out of time
  tree_->Branch("phoRecoFlag", &phoRecoFlag_);
  tree_->Branch("phoPos", &phoPos_);
  if (doGenParticles_) {
    tree_->Branch("phoGenIndex", &phoGenIndex_);
    tree_->Branch("phoGenGMomPID", &phoGenGMomPID_);
    tree_->Branch("phoGenMomPID", &phoGenMomPID_);
    tree_->Branch("phoGenMomPt", &phoGenMomPt_);
  }
  tree_->Branch("phoSCE", &phoSCE_);
  tree_->Branch("phoSCRawE", &phoSCRawE_);
  tree_->Branch("phoESEn", &phoESEn_);
  tree_->Branch("phoSCEt", &phoSCEt_);
  tree_->Branch("phoSCEta", &phoSCEta_);
  tree_->Branch("phoSCPhi", &phoSCPhi_);
  tree_->Branch("phoSCEtaWidth", &phoSCEtaWidth_);
  tree_->Branch("phoSCPhiWidth", &phoSCPhiWidth_);
  tree_->Branch("phoSCBrem", &phoSCBrem_);
  tree_->Branch("phoOverlap", &phoOverlap_);
  tree_->Branch("phohasPixelSeed", &phohasPixelSeed_);
  tree_->Branch("pho_hasConvPf", &pho_hasConvPf_);
  tree_->Branch("pho_hasSLConvPf", &pho_hasSLConvPf_);
  if (develop_) {
    tree_->Branch("MustacheEin", &MustacheEin_);
    tree_->Branch("MustacheEOut", &MustacheEOut_);
    tree_->Branch("MustacheEtOut", &MustacheEtOut_);
    tree_->Branch("PFRecoMatch", &PFRecoMatch_);
    tree_->Branch("PFEleMatch", &PFEleMatch_);
    tree_->Branch("PFEleVeto", &PFEleVeto_);
    tree_->Branch("PFLowestClustE", &PFLowestClustE_);
    tree_->Branch("PFClustdEta", &PFClustdEta_);
    tree_->Branch("PFClustdPhi", &PFClustdPhi_);
    tree_->Branch("PFClustRMSPhi", &PFClustRMSPhi_);
    tree_->Branch("PFClustRMSPhiMust", &PFClustRMSPhiMust_);
    tree_->Branch("PFPreShowerE1", &PFPreShowerE1_);
    tree_->Branch("PFPreShowerE2", &PFPreShowerE2_);
  }
  tree_->Branch("pho_pfconvVtxZ", &pho_pfconvVtxZ_);
  tree_->Branch("pho_pfconvVtxZErr", &pho_pfconvVtxZErr_);
  tree_->Branch("pho_nSLConv", &pho_nSLConv_);
  tree_->Branch("pho_pfSLConvPos_x", &pho_pfSLConvPos_x_);
  tree_->Branch("pho_pfSLConvPos_y", &pho_pfSLConvPos_y_);
  tree_->Branch("pho_pfSLConvPos_z", &pho_pfSLConvPos_z_);
  tree_->Branch("pho_pfSLConvVtxZ", &pho_pfSLConvVtxZ_);
  //Conversion Bran ches
  tree_->Branch("phoIsConv", &phoIsConv_);
  tree_->Branch("phoNConv", &phoNConv_);
  tree_->Branch("phoConvInvMass", &phoConvInvMass_);
  tree_->Branch("phoConvCotTheta", &phoConvCotTheta_);
  tree_->Branch("phoConvEoverP", &phoConvEoverP_);
  tree_->Branch("phoConvZofPVfromTrks", &phoConvZofPVfromTrks_);
  tree_->Branch("phoConvMinDist", &phoConvMinDist_);
  tree_->Branch("phoConvdPhiAtVtx", &phoConvdPhiAtVtx_);
  tree_->Branch("phoConvdPhiAtCalo", &phoConvdPhiAtCalo_);
  tree_->Branch("phoConvdEtaAtCalo", &phoConvdEtaAtCalo_);
  tree_->Branch("phoConvTrkd0_x", &phoConvTrkd0_x_);
  tree_->Branch("phoConvTrkd0_y", &phoConvTrkd0_y_);
  tree_->Branch("phoConvTrkPin_x", &phoConvTrkPin_x_);
  tree_->Branch("phoConvTrkPin_y", &phoConvTrkPin_y_);
  tree_->Branch("phoConvTrkPout_x", &phoConvTrkPout_x_);
  tree_->Branch("phoConvTrkPout_y", &phoConvTrkPout_y_);
  tree_->Branch("phoConvTrkdz_x", &phoConvTrkdz_x_);
  tree_->Branch("phoConvTrkdz_y", &phoConvTrkdz_y_);
  tree_->Branch("phoConvTrkdzErr_x", &phoConvTrkdzErr_x_);
  tree_->Branch("phoConvTrkdzErr_y", &phoConvTrkdzErr_y_);
  tree_->Branch("phoConvChi2", &phoConvChi2_);
  tree_->Branch("phoConvChi2Prob", &phoConvChi2Prob_);
  tree_->Branch("phoConvNTrks", &phoConvNTrks_);
  tree_->Branch("phoConvCharge1", &phoConvCharge1_);
  tree_->Branch("phoConvCharge2", &phoConvCharge2_);
  tree_->Branch("phoConvValidVtx", &phoConvValidVtx_);
  tree_->Branch("phoConvLikeLihood", &phoConvLikeLihood_);
  tree_->Branch("phoConvP4_0", &phoConvP4_0_);
  tree_->Branch("phoConvP4_1", &phoConvP4_1_);
  tree_->Branch("phoConvP4_2", &phoConvP4_2_);
  tree_->Branch("phoConvP4_3", &phoConvP4_3_);
  tree_->Branch("phoConvVtx_x", &phoConvVtx_x_);
  tree_->Branch("phoConvVtx_y", &phoConvVtx_y_);
  tree_->Branch("phoConvVtx_z", &phoConvVtx_z_);
  tree_->Branch("phoConvVtxErr_x", &phoConvVtxErr_x_);
  tree_->Branch("phoConvVtxErr_y", &phoConvVtxErr_y_);
  tree_->Branch("phoConvVtxErr_z", &phoConvVtxErr_z_);
  tree_->Branch("phoConvPairMomentum_x", &phoConvPairMomentum_x_);
  tree_->Branch("phoConvPairMomentum_y", &phoConvPairMomentum_y_);
  tree_->Branch("phoConvPairMomentum_z", &phoConvPairMomentum_z_);
  tree_->Branch("phoConvRefittedMomentum_x", &phoConvRefittedMomentum_x_);
  tree_->Branch("phoConvRefittedMomentum_y", &phoConvRefittedMomentum_y_);
  tree_->Branch("phoConvRefittedMomentum_z", &phoConvRefittedMomentum_z_);
  tree_->Branch("SingleLegConv", &SingleLegConv_);
  tree_->Branch("phoPFConvVtx_x", &phoPFConvVtx_x_);
  tree_->Branch("phoPFConvVtx_y", &phoPFConvVtx_y_);
  tree_->Branch("phoPFConvVtx_z", &phoPFConvVtx_z_);
  tree_->Branch("phoPFConvMom_x", &phoPFConvMom_x_);
  tree_->Branch("phoPFConvMom_y", &phoPFConvMom_y_);
  tree_->Branch("phoPFConvMom_z", &phoPFConvMom_z_);
  tree_->Branch("phoESEffSigmaRR_x", &phoESEffSigmaRR_x_);
  tree_->Branch("phoESEffSigmaRR_y", &phoESEffSigmaRR_y_);
  tree_->Branch("phoESEffSigmaRR_z", &phoESEffSigmaRR_z_);
  if (dumpESHits_) {
    tree_->Branch("phoESDetId", phoESDetId_, "phoESDetId_[nPho][2]/I" );
    tree_->Branch("phoESHits",  phoESHits_,  "phoESHits_[nPho][3][62]/F");
  }
  if (dumpESClusterInfo_) {
    tree_->Branch("phoESE1",  phoESE1_,  "phoESE1_[nPho][2]/F");
    tree_->Branch("phoESE3",  phoESE3_,  "phoESE3_[nPho][2]/F");
    tree_->Branch("phoESE5",  phoESE5_,  "phoESE5_[nPho][2]/F");
    tree_->Branch("phoESE7",  phoESE7_,  "phoESE7_[nPho][2]/F");
    tree_->Branch("phoESE11", phoESE11_, "phoESE11_[nPho][2]/F");
    tree_->Branch("phoESE21", phoESE21_, "phoESE21_[nPho][2]/F");
  }
  // Muon
  tree_->Branch("nMu", &nMu_, "nMu/I");
  tree_->Branch("muTrg", &muTrg_);
  tree_->Branch("muEta", &muEta_);
  tree_->Branch("muPhi", &muPhi_);
  tree_->Branch("muCharge", &muCharge_);
  tree_->Branch("muPt", &muPt_);
  tree_->Branch("muPz", &muPz_);
  tree_->Branch("muVtx_x", &muVtx_x_);
  tree_->Branch("muVtx_y", &muVtx_y_);
  tree_->Branch("muVtx_z", &muVtx_z_);
  tree_->Branch("muVtxGlb_x", &muVtxGlb_x_);
  tree_->Branch("muVtxGlb_y", &muVtxGlb_y_);
  tree_->Branch("muVtxGlb_z", &muVtxGlb_z_);
  if (doGenParticles_)
    tree_->Branch("muGenIndex", &muGenIndex_);
  tree_->Branch("mucktPt",    &mucktPt_);
  tree_->Branch("mucktPtErr", &mucktPtErr_);
  tree_->Branch("mucktEta",   &mucktEta_);
  tree_->Branch("mucktPhi",   &mucktPhi_);
  tree_->Branch("mucktdxy",   &mucktdxy_);
  tree_->Branch("mucktdz",    &mucktdz_);

  tree_->Branch("muIsoTrk",  &muIsoTrk_);
  tree_->Branch("muIsoCalo", &muIsoCalo_);
  tree_->Branch("muIsoEcal", &muIsoEcal_);
  tree_->Branch("muIsoHcal", &muIsoHcal_);

  tree_->Branch("muChi2NDF", &muChi2NDF_);
  tree_->Branch("muInnerChi2NDF", &muInnerChi2NDF_);
  tree_->Branch("muPFIsoR04_CH", &muPFIsoR04_CH_);
  tree_->Branch("muPFIsoR04_NH", &muPFIsoR04_NH_);
  tree_->Branch("muPFIsoR04_Pho", &muPFIsoR04_Pho_);
  tree_->Branch("muPFIsoR04_PU", &muPFIsoR04_PU_);
  tree_->Branch("muPFIsoR04_CPart", &muPFIsoR04_CPart_);
  tree_->Branch("muPFIsoR04_NHHT", &muPFIsoR04_NHHT_);
  tree_->Branch("muPFIsoR04_PhoHT", &muPFIsoR04_PhoHT_);
  tree_->Branch("muPFIsoR03_CH", &muPFIsoR03_CH_);
  tree_->Branch("muPFIsoR03_NH", &muPFIsoR03_NH_);
  tree_->Branch("muPFIsoR03_Pho", &muPFIsoR03_Pho_);
  tree_->Branch("muPFIsoR03_PU", &muPFIsoR03_PU_);
  tree_->Branch("muPFIsoR03_CPart", &muPFIsoR03_CPart_);
  tree_->Branch("muPFIsoR03_NHHT", &muPFIsoR03_NHHT_);
  tree_->Branch("muPFIsoR03_PhoHT", &muPFIsoR03_PhoHT_);
  tree_->Branch("muType", &muType_);
  tree_->Branch("muD0", &muD0_);
  tree_->Branch("muDz", &muDz_);
  tree_->Branch("muD0GV", &muD0GV_);
  tree_->Branch("muDzGV", &muDzGV_);
  tree_->Branch("muD0Vtx", &muD0Vtx_);
  tree_->Branch("muDzVtx", &muDzVtx_);
  tree_->Branch("muInnerD0", &muInnerD0_);
  tree_->Branch("muInnerDz", &muInnerDz_);
  tree_->Branch("muInnerD0GV", &muInnerD0GV_);
  tree_->Branch("muInnerDzGV", &muInnerDzGV_);
  tree_->Branch("muInnerPt", &muInnerPt_);
  tree_->Branch("muInnerPtErr", &muInnerPtErr_);
  tree_->Branch("muNumberOfValidTrkLayers", &muNumberOfValidTrkLayers_);
  tree_->Branch("muNumberOfValidTrkHits", &muNumberOfValidTrkHits_);
  tree_->Branch("muNumberOfValidPixelLayers", &muNumberOfValidPixelLayers_);
  tree_->Branch("muNumberOfValidPixelHits", &muNumberOfValidPixelHits_);
  tree_->Branch("muNumberOfValidMuonHits", &muNumberOfValidMuonHits_);
  tree_->Branch("muStations", &muStations_);
  tree_->Branch("muChambers", &muChambers_);
  tree_->Branch("muIP3D", &muIP3D_);
  tree_->Branch("muIP3DErr", &muIP3DErr_);

  // taus
  tree_->Branch("nTau", &nTau_, "nTau/I");
  tree_->Branch("tauDecayModeFinding", &tauDecayModeFinding_);
  tree_->Branch("tauAgainstElectronLooseMVA3", &tauAgainstElectronLooseMVA3_);
  tree_->Branch("tauAgainstElectronMediumMVA3", &tauAgainstElectronMediumMVA3_);
  tree_->Branch("tauAgainstElectronTightMVA3", &tauAgainstElectronTightMVA3_);
  tree_->Branch("tauAgainstElectronVTightMVA3", &tauAgainstElectronVTightMVA3_);
  tree_->Branch("tauAgainstElectronDeadECAL", &tauAgainstElectronDeadECAL_);
  tree_->Branch("tauAgainstMuonLoose2", &tauAgainstMuonLoose2_);
  tree_->Branch("tauAgainstMuonMedium2", &tauAgainstMuonMedium2_);
  tree_->Branch("tauAgainstMuonTight2", &tauAgainstMuonTight2_);
  tree_->Branch("tauCombinedIsolationDeltaBetaCorrRaw3Hits", &tauCombinedIsolationDeltaBetaCorrRaw3Hits_);
  tree_->Branch("tauLooseCombinedIsolationDeltaBetaCorr3Hits", &tauLooseCombinedIsolationDeltaBetaCorr3Hits_);
  tree_->Branch("tauMediumCombinedIsolationDeltaBetaCorr3Hits", &tauMediumCombinedIsolationDeltaBetaCorr3Hits_);
  tree_->Branch("tauTightCombinedIsolationDeltaBetaCorr3Hits", &tauTightCombinedIsolationDeltaBetaCorr3Hits_);
  tree_->Branch("tauEta", &tauEta_);
  tree_->Branch("tauPhi", &tauPhi_);
  tree_->Branch("tauPt", &tauPt_);
  tree_->Branch("tauEt", &tauEt_);
  tree_->Branch("tauCharge", &tauCharge_);
  tree_->Branch("tauDecayMode", &tauDecayMode_);
  tree_->Branch("tauEMFraction", &tauEMFraction_);
  tree_->Branch("tauHCAL3x3OverPLead", &tauHCAL3x3OverPLead_);
  tree_->Branch("tauHCALMaxOverPLead", &tauHCALMaxOverPLead_);
  tree_->Branch("tauHCALTotOverPLead", &tauHCALTotOverPLead_);
  tree_->Branch("tauIsolationPFChargedHadrCandsPtSum", &tauIsolationPFChargedHadrCandsPtSum_);
  tree_->Branch("tauIsolationPFGammaCandsEtSum", &tauIsolationPFGammaCandsEtSum_);
  tree_->Branch("tauLeadPFChargedHadrCandsignedSipt", &tauLeadPFChargedHadrCandsignedSipt_);
  tree_->Branch("tauLeadChargedHadronExists", &tauLeadChargedHadronExists_);
  tree_->Branch("tauLeadChargedHadronEta", &tauLeadChargedHadronEta_);
  tree_->Branch("tauLeadChargedHadronPhi", &tauLeadChargedHadronPhi_);
  tree_->Branch("tauLeadChargedHadronPt", &tauLeadChargedHadronPt_);
  if (develop_) {
    tree_->Branch("nPFPho", &nPFPho_, "nPFPho_/I");    
    tree_->Branch("PFPhoE", &PFPhoE_);
    tree_->Branch("PFPhoEt", &PFPhoEt_);
    tree_->Branch("PFPhoEta", &PFPhoEta_);
    tree_->Branch("PFPhoPhi", &PFPhoPhi_);
    tree_->Branch("PFPhoType", &PFPhoType_);
    tree_->Branch("PFPhoIso", &PFPhoIso_);
    tree_->Branch("PFPhoMaxEtaWidth_", PFPhoMaxEtaWidth_, "PFPhoMaxEtaWidth_[nPFPho_]/F");
    tree_->Branch("PFPhoMaxPhiWidth_", PFPhoMaxPhiWidth_, "PFPhoMaxPhiWidth_[nPFPho_]/F");
    tree_->Branch("nPFPhoClust_", nPFPhoClust_, "nPFPhoClust_[nPFPho_]/I");
    tree_->Branch("nPFPhoCrys_", nPFPhoCrys_, "nPFPhoCrys_[nPFPho_][20]/I");
    tree_->Branch("PFPho_clusteta_", PFPho_clusteta_, "PFPho_clusteta_[nPFPho_][20]/F");
    tree_->Branch("PFPho_clustphi_", PFPho_clustphi_, "PFPho_clustphi_[nPFPho_][20]/F");
    tree_->Branch("PFPho_clustE_", PFPho_clustE_, "PFPho_clustE_[nPFPho_][20]/F");   
    tree_->Branch("PFPho_clustSCEfrac_", PFPho_clustSCEfrac_, "PFPho_clustSCEfrac_[nPFPho_][20]/F");
    tree_->Branch("PFPho_clustEt_", PFPho_clustEt_, "PFPho_clustEt_[nPFPho_][20]/F");
    tree_->Branch("PFPho_clustEseed_", PFPho_clustEseed_, "PFPho_clustEseed_[nPFPho_][20]/F");
    tree_->Branch("PFPho_clustEtop_", PFPho_clustEtop_, "PFPho_clustEtop_[nPFPho_][20]/F");
    tree_->Branch("PFPho_clustEbottom_", PFPho_clustEbottom_, "PFPho_clustEbottom_[nPFPho_][20]/F");
    tree_->Branch("PFPho_clustEleft_", PFPho_clustEleft_, "PFPho_clustEleft_[nPFPho_][20]/F");
    tree_->Branch("PFPho_clustEright_", PFPho_clustEright_, "PFPho_clustEright_[nPFPho_][20]/F");
    tree_->Branch("PFPho_clustE3x3_", PFPho_clustE3x3_, "PFPho_clustE3x3_[nPFPho_][20]/F");
    tree_->Branch("PFPho_clustE1x3_", PFPho_clustE1x3_, "PFPho_clustE1x3_[nPFPho_][20]/F");
    tree_->Branch("PFPho_clustE3x1_", PFPho_clustE3x1_, "PFPho_clustE3x1_[nPFPho_][20]/F");
    tree_->Branch("PFPho_clustE1x5_", PFPho_clustE1x5_, "PFPho_clustE1x5_[nPFPho_][20]/F");
    tree_->Branch("PFPho_clustE5x5_", PFPho_clustE5x5_, "PFPho_clustE5x5_[nPFPho_][20]/F");
    tree_->Branch("PFPho_clustE2x5Max_", PFPho_clustE2x5Max_, "PFPho_clustE2x5Max_[nPFPho_][20]/F");
    tree_->Branch("PFPho_clustE2x5Top_", PFPho_clustE2x5Top_, "PFPho_clustE2x5Top_[nPFPho_][20]/F");
    tree_->Branch("PFPho_clustE2x5Bottom_", PFPho_clustE2x5Bottom_, "PFPho_clustE2x5Bottom_[nPFPho_][20]/F");
    tree_->Branch("PFPho_clustE2x5Left_", PFPho_clustE2x5Left_, "PFPho_clustE2x5Left_[nPFPho_][20]/F");
    tree_->Branch("PFPho_clustE2x5Right_", PFPho_clustE2x5Right_, "PFPho_clustE2x5Right_[nPFPho_][20]/F");
    tree_->Branch("PFPho_crysIeta_", PFPho_crysIeta_, "PFPho_crysIeta_[nPFPho_][20]/I");
    tree_->Branch("PFPho_crysIphi_", PFPho_crysIphi_, "PFPho_crysIphi_[nPFPho_][20]/I");
    tree_->Branch("PFPho_crysIX_", PFPho_crysIX_, "PFPho_crysIX_[nPFPho_][20]/I");
    tree_->Branch("PFPho_crysIY_", PFPho_crysIY_, "PFPho_crysIY_[nPFPho_][20]/I");
    tree_->Branch("PFPho_crysetafix_", PFPho_crysetafix_, "PFPho_crysetafix_[nPFPho_][20]/F"); 	 
    tree_->Branch("PFPho_crysphifix_", PFPho_crysphifix_, "PFPho_crysphifix_[nPFPho_][20]/F"); 	 
    tree_->Branch("PFPho_crysxfix_", PFPho_crysxfix_, "PFPho_crysxfix_[nPFPho_][20]");
    tree_->Branch("PFPho_crysyfix_", PFPho_crysyfix_, "PFPho_crysyfix_[nPFPho_][20]");    
    tree_->Branch("hasGSF_", hasGSF_, "hasGSF_[nPFPho_]/I");
    tree_->Branch("PFPhoGSFPin_", PFPhoGSFPin_, "PFPhoGSFPin_[nPFPho_][3]/F");
    tree_->Branch("PFPhoGSFPout_", PFPhoGSFPout_, "PFPhoGSFPout_[nPFPho_][3]/F"); 
    tree_->Branch("PFPhoGsf_In_", PFPhoGsf_In_, "PFPhoGsf_In_[nPFPho_][3]/F"); 
    tree_->Branch("PFPhoGsf_Theta_", PFPhoGsf_Theta_, "PFPhoGsf_Theta_[nPFPho_]/F"); 
    tree_->Branch("PFPhoGsf_ThetaErr_", PFPhoGsf_ThetaErr_, "PFPhoGsf_ThetaErr_[nPFPho_]/F");   
    tree_->Branch("PFPhoGsfeta_In_", PFPhoGsfeta_In_, "PFPhoGsfeta_In_[nPFPho_]/F");  
    tree_->Branch("PFPhoGsfeta_Out_", PFPhoGsfeta_Out_, "PFPhoGsfeta_Out_[nPFPho_]/F");  
    tree_->Branch("PFPhoGsfphi_In_", PFPhoGsfphi_In_, "PFPhoGsfphi_In_[nPFPho_]/F");  
    tree_->Branch("PFPhoGsfphi_Out_", PFPhoGsfphi_Out_, "PFPhoGsfphi_Out_[nPFPho_]/F");
    tree_->Branch("nPFPhoTks_", nPFPhoTks_, "nPFPhoTks_[nPFPho_]/I");
    tree_->Branch("PFPho_tkPos_", PFPho_tkPos_, "PFPho_tkPos_[nPFPho_][20][3]/F");
    tree_->Branch("PFPho_tkpt_", PFPho_tkpt_, "PFPho_tkpt_[nPFPho_][20]/F");	
    tree_->Branch("PFPho_tkR_", PFPho_tkR_, "PFPho_tkR_[nPFPho_][20]/F");
    tree_->Branch("PFPho_tkTheta_", PFPho_tkTheta_, "PFPho_tkTheta_[nPFPho_][20]");
    tree_->Branch("PFPho_tkerreta_", PFPho_tkerreta_, "PFPho_tkerreta_[nPFPho_][20]/F");
    tree_->Branch("PFPho_tkerrphi_", PFPho_tkerrphi_, "PFPho_tkerrphi_[nPFPho_][20]/F");
    tree_->Branch("PFPho_tkerrthet_", PFPho_tkerrthet_, "PFPho_tkerrthet_[nPFPho_][20]/F");
    tree_->Branch("PFPho_sigIetaIeta_", PFPho_sigIetaIeta_, "PFPho_sigIetaIeta_[nPFPho_][20]/F");
    tree_->Branch("PFPho_sigIphiIphi_", PFPho_sigIphiIphi_, "PFPho_sigIphiIphi_[nPFPho_][20]/F");
    tree_->Branch("PFPho_ES1Energy_", PFPho_ES1Energy_, "PFPho_ES1Energy_[nPFPho_][20]/F");
    tree_->Branch("PFPho_ES2Energy_", PFPho_ES2Energy_, "PFPho_ES2Energy_[nPFPho_][20]/F");
    tree_->Branch("PFPhoMustEtOut_", PFPhoMustEtOut_, "PFPhoMustEtOut_[nPFPho_]/F");
    tree_->Branch("PFPhoMustExcl_", PFPhoMustExcl_, "PFPhoMustExcl_[nPFPho_]/F");
    tree_->Branch("PFPhoMustEin_", PFPhoMustEin_, "PFPhoMustEin_[nPFPho_]/F");
    tree_->Branch("PFPhoMustEout_", PFPhoMustEout_, "PFPhoMustEout_[nPFPho_]/F");
  }
  if (dumpESClusterInfo_) {
    tree_->Branch("nPFPhoES1Clust_", nPFPhoES1Clust_, "nPFPhoES1Clust_[nPFPho_]/I");
    tree_->Branch("PFPho_ES1clusteta_", PFPho_ES1clusteta_, "PFPho_ES1clusteta_[nPFPho_][30]/F");
    tree_->Branch("PFPho_ES1clustphi_", PFPho_ES1clustphi_, "PFPho_ES1clustphi_[nPFPho_][30]/F");
    tree_->Branch("PFPho_ES1clustz_", PFPho_ES1clustz_, "PFPho_ES1clustz_[nPFPho_][30]/F");
    tree_->Branch("PFPho_ES1clustE_", PFPho_ES1clustE_, "PFPho_ES1clustE_[nPFPho_][30]/F");
    tree_->Branch("nPFPho_ES1Clust_perEEclust_",nPFPho_ES1Clust_perEEclust_,"nPFPho_ES1Clust_perEEclust_[nPFPho_][20]/I");
    tree_->Branch("PFPho_ES1size_perEEclust_",PFPho_ES1size_perEEclust_,"PFPho_ES1size_perEEclust_[nPFPho_][20]/I");
    tree_->Branch("PFPho_ES1weightedX_perEEclust_",PFPho_ES1weightedX_perEEclust_,"PFPho_ES1weightedX_perEEclust_[nPFPho_][20]/F");
    tree_->Branch("PFPho_ES1linkD_perEEclust_",PFPho_ES1linkD_perEEclust_,"PFPho_ES1linkD_perEEclust_[nPFPho_][20]/F");
    tree_->Branch("PFPho_ES1clustx_", PFPho_ES1clustx_, "PFPho_ES1clustx_[nPFPho_][30]/F");
    tree_->Branch("PFPho_ES1clusty_", PFPho_ES1clusty_, "PFPho_ES1clusty_[nPFPho_][30]/F");
    tree_->Branch("PFPho_ES1clustLinkDist_",PFPho_ES1clustLinkDist_,"PFPho_ES1clustLinkDist_[nPFPho_][30]/F");
    tree_->Branch("PFPho_ES1_stripsDetId_",PFPho_ES1_stripsDetId_,"PFPho_ES1_stripsDetId_[nPFPho_][30][30]/I");
    tree_->Branch("PFPho_ES1_stripsX_",PFPho_ES1_stripsX_,"PFPho_ES1_stripsX_[nPFPho_][30][30]/F");
    tree_->Branch("PFPho_ES1_stripsY_",PFPho_ES1_stripsY_,"PFPho_ES1_stripsY_[nPFPho_][30][30]/F");
    tree_->Branch("PFPho_ES1_stripsZ_",PFPho_ES1_stripsZ_,"PFPho_ES1_stripsZ_[nPFPho_][30][30]/F");
    tree_->Branch("PFPho_ES1_stripsEta_",PFPho_ES1_stripsEta_,"PFPho_ES1_stripsEta_[nPFPho_][30][30]/F");
    tree_->Branch("PFPho_ES1_stripsPhi_",PFPho_ES1_stripsPhi_,"PFPho_ES1_stripsPhi_[nPFPho_][30][30]/F");
    tree_->Branch("PFPho_ES1_stripsFrac_",PFPho_ES1_stripsFrac_,"PFPho_ES1_stripsFrac_[nPFPho_][30][30]/F");
    tree_->Branch("PFPho_ES1_stripsE_",PFPho_ES1_stripsE_,"PFPho_ES1_stripsE_[nPFPho_][30][30]/F");
    tree_->Branch("PFPho_ES1size_",PFPho_ES1size_,"PFPho_ES1size_[nPFPho_][30]/I");
    tree_->Branch("PFPho_weightedXst_perES1clust_",PFPho_weightedXst_perES1clust_,"PFPho_weightedXst_perES1clust_[nPFPho_][30]/F");
    tree_->Branch("nPFPhoES2Clust_", nPFPhoES2Clust_, "nPFPhoES2Clust_[nPFPho_]/I");
    tree_->Branch("PFPho_ES2clusteta_", PFPho_ES2clusteta_, "PFPho_ES2clusteta_[nPFPho_][30]/F");
    tree_->Branch("PFPho_ES2clustphi_", PFPho_ES2clustphi_, "PFPho_ES2clustphi_[nPFPho_][30]/F");
    tree_->Branch("PFPho_ES2clustz_", PFPho_ES2clustz_, "PFPho_ES2clustz_[nPFPho_][30]/F");
    tree_->Branch("PFPho_ES2clustE_", PFPho_ES2clustE_, "PFPho_ES2clustE_[nPFPho_][30]/F");
    tree_->Branch("nPFPho_ES2Clust_perEEclust_",nPFPho_ES2Clust_perEEclust_,"nPFPho_ES2Clust_perEEclust_[nPFPho_][20]/I");
    tree_->Branch("PFPho_ES2size_perEEclust_",PFPho_ES2size_perEEclust_,"PFPho_ES2size_perEEclust_[nPFPho_][20]/I");
    tree_->Branch("PFPho_ES2weightedY_perEEclust_",PFPho_ES2weightedY_perEEclust_,"PFPho_ES2weightedY_perEEclust_[nPFPho_][20]/F");
    tree_->Branch("PFPho_ES_siphiiphi_",PFPho_ES_siphiiphi_,"PFPho_ES_siphiiphi_[nPFPho_][20]/F");
    tree_->Branch("PFPho_ES_sietaieta_",PFPho_ES_sietaieta_,"PFPho_ES_sietaieta_[nPFPho_][20]/F");
    tree_->Branch("PFPho_ES2linkD_perEEclust_",PFPho_ES2linkD_perEEclust_,"PFPho_ES2linkD_perEEclust_[nPFPho_][20]/F");
    tree_->Branch("PFPho_ES2clustx_", PFPho_ES2clustx_, "PFPho_ES2clustx_[nPFPho_][30]/F");
    tree_->Branch("PFPho_ES2clusty_", PFPho_ES2clusty_, "PFPho_ES2clusty_[nPFPho_][30]/F");
    tree_->Branch("PFPho_ES2clustLinkDist_",PFPho_ES2clustLinkDist_,"PFPho_ES2clustLinkDist_[nPFPho_][30]/F");
    tree_->Branch("PFPho_ES2_stripsDetId_",PFPho_ES2_stripsDetId_,"PFPho_ES2_stripsDetId_[nPFPho_][30][30]/I");
    tree_->Branch("PFPho_ES2_stripsX_",PFPho_ES2_stripsX_,"PFPho_ES2_stripsX_[nPFPho_][30][30]/F");
    tree_->Branch("PFPho_ES2_stripsY_",PFPho_ES2_stripsY_,"PFPho_ES2_stripsY_[nPFPho_][30][30]/F");
    tree_->Branch("PFPho_ES2_stripsZ_",PFPho_ES2_stripsZ_,"PFPho_ES2_stripsZ_[nPFPho_][30][30]/F");
    tree_->Branch("PFPho_ES2_stripsEta_",PFPho_ES2_stripsEta_,"PFPho_ES2_stripsEta_[nPFPho_][30][30]/F");
    tree_->Branch("PFPho_ES2_stripsPhi_",PFPho_ES2_stripsPhi_,"PFPho_ES2_stripsPhi_[nPFPho_][30][30]/F");
    tree_->Branch("PFPho_ES2_stripsFrac_",PFPho_ES2_stripsFrac_,"PFPho_ES2_stripsFrac_[nPFPho_][30][30]/F");
    tree_->Branch("PFPho_ES2_stripsE_",PFPho_ES2_stripsE_,"PFPho_ES2_stripsE_[nPFPho_][30][30]/F");
    tree_->Branch("PFPho_ES2size_",PFPho_ES2size_,"PFPho_ES2size_[nPFPho_][30]/I");
    tree_->Branch("PFPho_weightedYst_perES2clust_",PFPho_weightedYst_perES2clust_,"PFPho_weightedYst_perES2clust_[nPFPho_][30]/F");
  }
  //PFElectrons  
  if (develop_) {
    tree_->Branch("nPFEle_", &nPFEle_, "nPFEle_/I");  
    tree_->Branch("PFElePt_", PFElePt_,"PFElePt_[nPFEle_]/F");  
    tree_->Branch("PFEleEta_", PFEleEta_,"PFEleEta_[nPFEle_]/F");  
    tree_->Branch("PFElePhi_", PFElePhi_,"PFElePhi_[nPFEle_]/F");   
    tree_->Branch("PFEleEn_", PFEleEn_, "PFEleEn_[nPFEle_]/F");  
    tree_->Branch("PFEleCharge", PFEleCharge_, "PFEleCharge_[nPFEle_]/I");
    tree_->Branch("PFEleMaxEtaWidth_", PFEleMaxEtaWidth_, "PFEleMaxEtaWidth_[nPFEle_]/F");
    tree_->Branch("PFEleMaxPhiWidth_", PFEleMaxPhiWidth_, "PFEleMaxPhiWidth_[nPFEle_]/F");
    tree_->Branch("nPFEleClust_", nPFEleClust_, "nPFEleClust_[nPFEle_]/I");
    tree_->Branch("PFEle_clusteta_", PFEle_clusteta_,"PFEle_clusteta_[nPFEle_][20]/F");
    tree_->Branch("PFEle_clustphi_", PFEle_clustphi_,"PFEle_clustphi_[nPFEle_][20]/F");
    tree_->Branch("PFEle_clustE_", PFEle_clustE_,"PFEle_clustE_[nPFEle_][20]/F");
    tree_->Branch("PFEle_clustSCEfrac_", PFEle_clustSCEfrac_, "PFEle_clustSCEfrac_[nPFEle_][20]/F");
    tree_->Branch("PFEle_clustEt_", PFEle_clustEt_,"PFEle_clustEt_[nPFEle_][20]/F");
    tree_->Branch("PFEle_crysIeta_", PFEle_crysIeta_, "PFEle_crysIeta_[nPFEle_][20]/I");
    tree_->Branch("PFEle_crysIphi_", PFEle_crysIphi_, "PFEle_crysIphi_[nPFEle_][20]/I");
    tree_->Branch("PFEle_crysIX_", PFEle_crysIX_, "PFEle_crysIX_[nPFEle_][20]/I");
    tree_->Branch("PFEle_crysIY_", PFEle_crysIY_, "PFEle_crysIY_[nPFEle_][20]/I");
    tree_->Branch("PFEle_crysetafix_", PFEle_crysetafix_, "PFEle_crysetafix_[nPFEle_][20]/F"); 	 
    tree_->Branch("PFEle_crysphifix_", PFEle_crysphifix_, "PFEle_crysphifix_[nPFEle_][20]/F"); 	 
    tree_->Branch("PFEle_modetafix_", PFEle_modetafix_, "PFEle_modetafix_[nPFEle_][20]/F"); 	 
    tree_->Branch("PFEle_modphifix_", PFEle_modphifix_, "PFEle_modphifix_[nPFEle_][20]/F"); 	 
    tree_->Branch("PFEle_Smodetafix_", PFEle_Smodetafix_, "PFEle_Smodetafix_[nPFEle_][20]/F"); 	 
    tree_->Branch("PFEle_Smodphifix_", PFEle_Smodphifix_, "PFEle_Smodphifix_[nPFEle_][20]/F");
    tree_->Branch("PFEle_crysxfix_", PFEle_crysxfix_, "PFEle_crysxfix_[nPFEle_][20]");
    tree_->Branch("PFEle_crysyfix_", PFEle_crysyfix_, "PFEle_crysyfix_[nPFEle_][20]");
    tree_->Branch("PFEle_modxfix_", PFEle_modxfix_, "PFEle_modxfix_[nPFEle_][20]");
    tree_->Branch("PFEle_modyfix_", PFEle_modyfix_, "PFEle_modyfix_[nPFEle_][20]");
    tree_->Branch("PFEle_Smodxfix_", PFEle_Smodxfix_, "PFEle_Smodxfix_[nPFEle_][20]");
    tree_->Branch("PFEle_Smodyfix_", PFEle_Smodyfix_, "PFEle_Smodyfix_[nPFEle_][20]");
    tree_->Branch("PFEle_clustEseed_", PFEle_clustEseed_, "PFEle_clustEseed_[nPFEle_][20]/F");
    tree_->Branch("PFEle_clustEtop_", PFEle_clustEtop_, "PFEle_clustEtop_[nPFEle_][20]/F");
    tree_->Branch("PFEle_clustEbottom_", PFEle_clustEbottom_, "PFEle_clustEbottom_[nPFEle_][20]/F");
    tree_->Branch("PFEle_clustEleft_", PFEle_clustEleft_, "PFEle_clustEleft_[nPFEle_][20]/F");
    tree_->Branch("PFEle_clustEright_", PFEle_clustEright_, "PFEle_clustEright_[nPFEle_][20]/F");
    tree_->Branch("PFEle_clustE3x3_", PFEle_clustE3x3_, "PFEle_clustE3x3_[nPFEle_][20]/F");
    tree_->Branch("PFEle_clustE1x3_", PFEle_clustE1x3_, "PFEle_clustE1x3_[nPFEle_][20]/F");
    tree_->Branch("PFEle_clustE3x1_", PFEle_clustE3x1_, "PFEle_clustE3x1_[nPFEle_][20]/F");
    tree_->Branch("PFEle_clustE1x5_", PFEle_clustE1x5_, "PFEle_clustE1x5_[nPFEle_][20]/F");
    tree_->Branch("PFEle_clustE5x5_", PFEle_clustE5x5_, "PFEle_clustE5x5_[nPFEle_][20]/F");
    tree_->Branch("PFEle_clustE2x5Max_", PFEle_clustE2x5Max_, "PFEle_clustE2x5Max_[nPFEle_][20]/F");
    tree_->Branch("PFEle_clustE2x5Top_", PFEle_clustE2x5Top_, "PFEle_clustE2x5Top_[nPFEle_][20]/F");
    tree_->Branch("PFEle_clustE2x5Bottom_", PFEle_clustE2x5Bottom_, "PFEle_clustE2x5Bottom_[nPFEle_][20]/F");
    tree_->Branch("PFEle_clustE2x5Left_", PFEle_clustE2x5Left_, "PFEle_clustE2x5Left_[nPFEle_][20]/F");
    tree_->Branch("PFEle_clustE2x5Right_", PFEle_clustE2x5Right_, "PFEle_clustE2x5Right_[nPFEle_][20]/F");
    tree_->Branch("PFElePin_", PFElePin_,"PFElePin_[nPFEle_][3]/F");  
    tree_->Branch("PFElePout_", PFElePout_,"PFElePout_[nPFEle_]/F"); 
    tree_->Branch("PFEleChi2NDF_", PFEleChi2NDF_, "PFEleChi2NDF_[nPFEle_]/F");
    tree_->Branch("PFGsf_In_", PFGsf_In_, "PFGsf_In_[nPFEle_][3]/F");
    tree_->Branch("PFGsf_Theta_", PFGsf_Theta_, "PFGsf_Theta_[nPFEle_]/F");
    tree_->Branch("PFGsf_ThetaErr_", PFGsf_ThetaErr_, "PFGsf_ThetaErr_[nPFEle_]/F");
    tree_->Branch("PFGsfeta_In_", PFGsfeta_In_, "PFGsfeta_In_[nPFEle_]/F");
    tree_->Branch("PFGsfeta_Out_", PFGsfeta_Out_, "PFGsfeta_Out_[nPFEle_]/F");
    tree_->Branch("PFGsfphi_In_", PFGsfphi_In_, "PFGsfphi_In_[nPFEle_]/F");
    tree_->Branch("PFGsfphi_Out_", PFGsfphi_Out_, "PFGsfphi_Out_[nPFEle_]/F");
    tree_->Branch("PFEleMustEtOut_", PFEleMustEtOut_, "PFEleMustEtOut_[nPFEle_]/F");
    tree_->Branch("PFEleMustExcl_", PFEleMustExcl_, "PFEleMustExcl_[nPFEle_]/F");
    tree_->Branch("PFEleMustEin_", PFEleMustEin_, "PFEleMustEin_[nPFEle_]/F");
    tree_->Branch("PFEleMustEout_", PFEleMustEout_, "PFEleMustEout_[nPFEle_]/F");
  }
  if (useAllPF_) {
    //PFCharged Hadrons:  
    tree_->Branch("nPFchad_",&nPFchad_, "nPFchad_/I");  
    tree_->Branch("PFhad_charge_",PFhad_charge_, "PFhad_charge_[nPFchad_]/I");  
    tree_->Branch("PFchad_Eta_", PFchad_Eta_, "PFchad_Eta_[nPFchad_]/F");  
    tree_->Branch("PFchad_Phi_", PFchad_Phi_, "PFchad_Phi_[nPFchad_]/F");    
    tree_->Branch("PFchad_Pt_", PFchad_Pt_, "PFchad_Pt_[nPFchad_]/F");    
    tree_->Branch("PFchad_P_", PFchad_P_, "PFchad_P_[nPFchad_]/F");  
    tree_->Branch("PFchad_E_", PFchad_E_, "PFchad_E_[nPFchad_]/F");  
    tree_->Branch("PFchad_Ecal_", PFchad_Ecal_, "PFchad_Ecal_[nPFchad_]/F");  
    tree_->Branch("PFchad_Hcal_", PFchad_Hcal_, "PFchad_Hcal_[nPFchad_]/F");  
    //PFNeutral Hadrons  
    tree_->Branch("nPFnhad_",&nPFnhad_, "nPFnhad_/I");  
    tree_->Branch("PFnhad_Eta_", PFnhad_Eta_, "PFnhad_Eta_[nPFnhad_]/F");  
    tree_->Branch("PFnhad_Phi_", PFnhad_Phi_, "PFnhad_Phi_[nPFnhad_]/F");    
    tree_->Branch("PFnhad_Pt_", PFnhad_Pt_, "PFnhad_Pt_[nPFnhad_]/F");    
    tree_->Branch("PFnhad_P_", PFnhad_P_, "PFnhad_P_[nPFnhad_]/F");  
    tree_->Branch("PFnhad_E_", PFnhad_E_, "PFnhad_E_[nPFnhad_]/F");  
    tree_->Branch("PFnhad_Hcal_", PFnhad_Hcal_, "PFnhad_Hcal_[nPFnhad_]/F");
  }
  // rho
  tree_->Branch("rho25", &rho25_, "rho25/F");
  tree_->Branch("rho25_neu", &rho25_neu_, "rho25_neu/F");
  tree_->Branch("rho25_muPFiso", &rho25_muPFiso_, "rho25_muPFiso/F");
  tree_->Branch("rho25_elePFiso", &rho25_elePFiso_, "rho25_elePFiso/F");
  tree_->Branch("rho2011", &rho2011_, "rho2011/F");
  tree_->Branch("rho2012", &rho2012_, "rho2012/F");

  //QGTag
  tree_->Branch("QGTag_MLP", &QGTag_MLP_, "QGTag_MLP/F");
  tree_->Branch("QGTag_likelihood", &QGTag_likelihood_, "QGTag_likelihood/F");
  
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
    tree_->Branch("CA8prunedJetMass", &CA8prunedJetMass_);
    tree_->Branch("CA8prunedJet_nSubJets", &CA8prunedJet_nSubJets_) ;
    tree_->Branch("CA8prunedJet_SubjetPt", &CA8prunedJet_SubjetPt_);
    tree_->Branch("CA8prunedJet_SubjetEta", &CA8prunedJet_SubjetEta_);
    tree_->Branch("CA8prunedJet_SubjetPhi", &CA8prunedJet_SubjetPhi_);
    tree_->Branch("CA8prunedJet_SubjetMass", &CA8prunedJet_SubjetMass_);
  }

  // Jet
  if (dumpJets_) {
    tree_->Branch("nJet", &nJet_, "nJet/I");
    tree_->Branch("jetTrg", &jetTrg_);
    tree_->Branch("jetEn", &jetEn_);
    tree_->Branch("jetPt", &jetPt_);
    tree_->Branch("jetEta", &jetEta_);
    tree_->Branch("jetPhi", &jetPhi_);
    tree_->Branch("jetCharge", &jetCharge_);
    tree_->Branch("jetEt", &jetEt_);
    tree_->Branch("jetRawPt", &jetRawPt_);
    tree_->Branch("jetRawEn", &jetRawEn_);
    tree_->Branch("jetArea", &jetArea_);
    tree_->Branch("jetCHF", &jetCHF_);
    tree_->Branch("jetNHF", &jetNHF_);
    tree_->Branch("jetCEF", &jetCEF_);
    tree_->Branch("jetNEF", &jetNEF_);
    tree_->Branch("jetNCH", &jetNCH_);
    tree_->Branch("jetHFHAE", &jetHFHAE_);
    tree_->Branch("jetHFEME", &jetHFEME_);
    tree_->Branch("jetNConstituents", &jetNConstituents_);
    tree_->Branch("jetCombinedSecondaryVtxBJetTags", &jetCombinedSecondaryVtxBJetTags_);
    tree_->Branch("jetCombinedSecondaryVtxMVABJetTags", &jetCombinedSecondaryVtxMVABJetTags_);
    tree_->Branch("jetJetProbabilityBJetTags", &jetJetProbabilityBJetTags_);
    tree_->Branch("jetJetBProbabilityBJetTags", &jetJetBProbabilityBJetTags_);
    tree_->Branch("jetBetaStar", &jetBetaStar_);
    // CMG Jet Id Variables
    tree_->Branch("jetPFLooseId", &jetPFLooseId_);
    tree_->Branch("jetDRMean", &jetDRMean_);
    tree_->Branch("jetDR2Mean", &jetDR2Mean_);
    tree_->Branch("jetDZ", &jetDZ_);
    tree_->Branch("jetFrac01", &jetFrac01_);
    tree_->Branch("jetFrac02", &jetFrac02_);
    tree_->Branch("jetFrac03", &jetFrac03_);
    tree_->Branch("jetFrac04", &jetFrac04_);
    tree_->Branch("jetFrac05", &jetFrac05_);
    tree_->Branch("jetFrac06", &jetFrac06_);
    tree_->Branch("jetFrac07", &jetFrac07_);
    tree_->Branch("jetBeta", &jetBeta_);
    tree_->Branch("jetBetaStarCMG", &jetBetaStarCMG_);
    tree_->Branch("jetBetaStarClassic", &jetBetaStarClassic_);
    tree_->Branch("jetBetaExt", &jetBetaExt_);
    tree_->Branch("jetBetaStarCMGExt", &jetBetaStarCMGExt_);
    tree_->Branch("jetBetaStarClassicExt", &jetBetaStarClassicExt_);
    tree_->Branch("jetNNeutrals", &jetNNeutrals_);
    tree_->Branch("jetNCharged", &jetNCharged_);
    tree_->Branch("jetMVAs", &jetMVAs_); // 0: simple, 1: full, 2: cut based, 3: philv1
    tree_->Branch("jetWPLevels", &jetWPLevels_); // 0: simple, 1: full, 2: cut based, 3: philv1
//     tree_->Branch("jetMVAsExt", jetMVAsExt_, "jetMVAsExt[nJet][4][100]/F");             // 0: simple, 1: full, 2: cut based, 3: philv1; vtxbs
//     tree_->Branch("jetWPLevelsExt", jetWPLevelsExt_, "jetWPLevelsExt[nJet][4][100]/I"); // 0: simple, 1: full, 2: cut based, 3: philv1; vtxbs
    tree_->Branch("jetMVAsExt_simple", &jetMVAsExt_simple_);
    tree_->Branch("jetWPLevelsExt_simple", &jetWPLevelsExt_simple_);
    tree_->Branch("jetMVAsExt_full", &jetMVAsExt_full_);
    tree_->Branch("jetWPLevelsExt_full", &jetWPLevelsExt_full_);
    tree_->Branch("jetMVAsExt_cutBased", &jetMVAsExt_cutBased_);
    tree_->Branch("jetWPLevelsExt_cutBased", &jetWPLevelsExt_cutBased_);
    tree_->Branch("jetMVAsExt_philv1", &jetMVAsExt_philv1_);
    tree_->Branch("jetWPLevelsExt_philv1", &jetWPLevelsExt_philv1_);

    //b-jet regression variables
    tree_->Branch("jetMt", &jetMt_);
    tree_->Branch("jetJECUnc", &jetJECUnc_);
    tree_->Branch("jetLeadTrackPt", &jetLeadTrackPt_);
    tree_->Branch("jetVtxPt", &jetVtxPt_);
    tree_->Branch("jetVtxMass", &jetVtxMass_);
    tree_->Branch("jetVtx3dL", &jetVtx3dL_);
    tree_->Branch("jetVtx3deL", &jetVtx3deL_);
    tree_->Branch("jetSoftLeptPt", &jetSoftLeptPt_);
    tree_->Branch("jetSoftLeptPtRel", &jetSoftLeptPtRel_);
    tree_->Branch("jetSoftLeptdR", &jetSoftLeptdR_);
    tree_->Branch("jetSoftLeptIdlooseMu", &jetSoftLeptIdlooseMu_);
    tree_->Branch("jetSoftLeptIdEle95", &jetSoftLeptIdEle95_);
    tree_->Branch("jetDPhiMETJet", &jetDPhiMETJet_);
    tree_->Branch("jetPuJetIdL", &jetPuJetIdL_);
    tree_->Branch("jetPuJetIdM", &jetPuJetIdM_);
    tree_->Branch("jetPuJetIdT", &jetPuJetIdT_);
    if (doGenParticles_) {
      tree_->Branch("jetPartonID", &jetPartonID_);
      tree_->Branch("jetGenJetIndex", &jetGenJetIndex_);
      tree_->Branch("jetGenJetEn", &jetGenJetEn_);
      tree_->Branch("jetGenJetPt", &jetGenJetPt_);
      tree_->Branch("jetGenJetEta", &jetGenJetEta_);
      tree_->Branch("jetGenJetPhi", &jetGenJetPhi_);
      tree_->Branch("jetGenPartonID", &jetGenPartonID_);
      tree_->Branch("jetGenEn", &jetGenEn_);
      tree_->Branch("jetGenPt", &jetGenPt_);
      tree_->Branch("jetGenEta", &jetGenEta_);
      tree_->Branch("jetGenPhi", &jetGenPhi_);
      tree_->Branch("jetGenPartonMomID", &jetGenPartonMomID_);
    }
    // Low Pt Jets
    if (dumpTrks_) {
      tree_->Branch("nLowPtJet", &nLowPtJet_, "nLowPtJet/I");
      tree_->Branch("jetLowPtEn", &jetLowPtEn_);
      tree_->Branch("jetLowPtPt", &jetLowPtPt_);
      tree_->Branch("jetLowPtEta", &jetLowPtEta_);
      tree_->Branch("jetLowPtPhi", &jetLowPtPhi_);
      tree_->Branch("jetLowPtCharge", &jetLowPtCharge_);
      tree_->Branch("jetLowPtEt", &jetLowPtEt_);
      tree_->Branch("jetLowPtRawPt", &jetLowPtRawPt_);
      tree_->Branch("jetLowPtRawEn", &jetLowPtRawEn_);
      tree_->Branch("jetLowPtArea", &jetLowPtArea_);
      if (doGenParticles_) {
	tree_->Branch("jetLowPtPartonID", &jetLowPtPartonID_);
	tree_->Branch("jetLowPtGenJetEn", &jetLowPtGenJetEn_);
	tree_->Branch("jetLowPtGenJetPt", &jetLowPtGenJetPt_);
	tree_->Branch("jetLowPtGenJetEta", &jetLowPtGenJetEta_);
	tree_->Branch("jetLowPtGenJetPhi", &jetLowPtGenJetPhi_);
	tree_->Branch("jetLowPtGenPartonID", &jetLowPtGenPartonID_);
	tree_->Branch("jetLowPtGenEn", &jetLowPtGenEn_);
	tree_->Branch("jetLowPtGenPt", &jetLowPtGenPt_);
	tree_->Branch("jetLowPtGenEta", &jetLowPtGenEta_);
	tree_->Branch("jetLowPtGenPhi", &jetLowPtGenPhi_);
      }
    }
  }
  // Converted Photon Collection
  tree_->Branch("nConv", &nConv_, "nConv/I");
  tree_->Branch("convP4_x", &convP4_x_);
  tree_->Branch("convP4_y", &convP4_y_);
  tree_->Branch("convP4_z", &convP4_z_);
  tree_->Branch("convP4_E", &convP4_E_);
  tree_->Branch("convVtx_x", &convVtx_x_);
  tree_->Branch("convVtx_y", &convVtx_y_);
  tree_->Branch("convVtx_z", &convVtx_z_);
  tree_->Branch("convVtxErr_x", &convVtxErr_x_);
  tree_->Branch("convVtxErr_y", &convVtxErr_y_);
  tree_->Branch("convVtxErr_z", &convVtxErr_z_);
  tree_->Branch("convPairMomentum_x", &convPairMomentum_x_);
  tree_->Branch("convPairMomentum_y", &convPairMomentum_y_);
  tree_->Branch("convPairMomentum_z", &convPairMomentum_z_);
  tree_->Branch("convRefittedMomentum_x", &convRefittedMomentum_x_);
  tree_->Branch("convRefittedMomentum_y", &convRefittedMomentum_y_);
  tree_->Branch("convRefittedMomentum_z", &convRefittedMomentum_z_);
  tree_->Branch("convNTracks", &convNTracks_);
  tree_->Branch("convPairInvMass", &convPairInvMass_);
  tree_->Branch("convPairCotThetaSep", &convPairCotThetaSep_);
  tree_->Branch("convEoverP", &convEoverP_);
  tree_->Branch("convDistOfMinApproach", &convDistOfMinApproach_);
  tree_->Branch("convDPhiTrksAtVtx", &convDPhiTrksAtVtx_);
  tree_->Branch("convDPhiTrksAtEcal", &convDPhiTrksAtEcal_);
  tree_->Branch("convDEtaTrksAtEcal", &convDEtaTrksAtEcal_);
  tree_->Branch("convDxy", &convDxy_);
  tree_->Branch("convDz", &convDz_);
  tree_->Branch("convLxy", &convLxy_);
  tree_->Branch("convLz", &convLz_);
  tree_->Branch("convZofPrimVtxFromTrks", &convZofPrimVtxFromTrks_);
  tree_->Branch("convNHitsBeforeVtx_0", &convNHitsBeforeVtx_0_);
  tree_->Branch("convNHitsBeforeVtx_1", &convNHitsBeforeVtx_1_);
  tree_->Branch("convNSharedHits", &convNSharedHits_);
  tree_->Branch("convValidVtx", &convValidVtx_);
  tree_->Branch("convMVALikelihood", &convMVALikelihood_);
  tree_->Branch("convChi2", &convChi2_);
  tree_->Branch("convChi2Probability", &convChi2Probability_);
  // per track quantities
  tree_->Branch("convTk1Dz", &convTk1Dz_);
  tree_->Branch("convTk2Dz", &convTk2Dz_);
  tree_->Branch("convTk1DzErr", &convTk1DzErr_);
  tree_->Branch("convTk2DzErr", &convTk2DzErr_);
  tree_->Branch("convCh1Ch2", &convCh1Ch2_);
  tree_->Branch("convTk1D0", &convTk1D0_);
  tree_->Branch("convTk1Pout", &convTk1Pout_);
  tree_->Branch("convTk1Pin", &convTk1Pin_);
  tree_->Branch("convTk2D0", &convTk2D0_);
  tree_->Branch("convTk2Pout", &convTk2Pout_);
  tree_->Branch("convTk2Pin", &convTk2Pin_);

  vtxbsTkIndex_  = new std::vector<std::vector<Int_t> >;   vtxbsTkIndex_->clear();
  vtxbsTkWeight_ = new std::vector<std::vector<Float_t> >; vtxbsTkWeight_->clear();

}

ggNtuplizer::~ggNtuplizer() {

  delete cicPhotonId_;
  delete trackMET_;
  delete vtxbsTkIndex_;
  delete vtxbsTkWeight_;

}

void ggNtuplizer::clearVectors() {
  vtx_x_.clear();
  vtx_y_.clear();
  vtx_z_.clear();
  vtxbs_x_.clear();
  vtxbs_y_.clear();
  vtxbs_z_.clear();
  vtxbsTkIndex_->clear();
  vtxbsTkWeight_->clear();
  trkP_x_.clear();
  trkP_y_.clear();
  trkP_z_.clear();
  trkVtx_x_.clear();
  trkVtx_y_.clear();
  trkVtx_z_.clear();
  trkd0_.clear();
  trkd0Err_.clear();
  trkdz_.clear();
  trkdzErr_.clear();
  trkPtErr_.clear();
  trkQuality_.clear();

  mcPID.clear();
  mcVtx_x.clear();
  mcVtx_y.clear();
  mcVtx_z.clear();
  mcPt.clear();
  mcMass.clear();
  mcEta.clear();
  mcPhi.clear();
  mcE.clear();
  mcEt.clear();
  mcGMomPID.clear();
  mcMomPID.clear();
  mcMomPt.clear();
  mcMomMass.clear();
  mcMomEta.clear();
  mcMomPhi.clear();
  mcIndex.clear();
  mcDecayType.clear();
  mcParentage.clear();
  mcStatus.clear();
  
  nPU_.clear();
  puBX_.clear();
  puTrue_.clear();

  trkMETx_.clear();
  trkMETy_.clear();
  trkMETPhi_.clear();
  trkMET_.clear();

  eleTrg_.clear();
  eleClass_.clear();
  eleIsEcalDriven_.clear();
  eleCharge_.clear();
  eleChargeConsistent_.clear();
  eleEn_.clear();
  eleEcalEn_.clear();
  eleSCEn_.clear();
  eleESEn_.clear();
  eleVtx_x_.clear();
  eleVtx_y_.clear();
  eleVtx_z_.clear();
  eleD0_.clear();
  eleDz_.clear();
  eleD0GV_.clear();
  eleDzGV_.clear();
  eleD0Vtx_.clear();
  eleDzVtx_.clear();
  elePt_.clear();
  eleEta_.clear();
  elePhi_.clear();
  eleSCEta_.clear();
  eleSCPhi_.clear();
  eleEtVtx_.clear();
  eleEtaVtx_.clear();
  elePhiVtx_.clear();
  eleSCRawEn_.clear();
  eleSCEtaWidth_.clear();
  eleSCPhiWidth_.clear();
  eleHoverE_.clear();
  eleHoverE12_.clear(); 
  eleEoverP_.clear();
  elePin_.clear();
  elePout_.clear();
  eleTrkMomErr_.clear();
  eleBrem_.clear();
  eledEtaAtVtx_.clear();
  eledPhiAtVtx_.clear();
  eleSigmaIEtaIEta_.clear();
  eleSigmaIEtaIPhi_.clear();
  eleSigmaIPhiIPhi_.clear();
  eleEmax_.clear();
  eleE1x5_.clear();
  eleE3x3_.clear();
  eleE5x5_.clear();
  eleE2x5Max_.clear();
  eleRegrE_.clear();
  eleRegrEerr_.clear();
  elePhoRegrE_.clear();
  elePhoRegrEerr_.clear();
  eleSeedTime_.clear();
  eleRecoFlag_.clear();
  elePos_.clear();
  eleGenIndex_.clear();
  eleGenGMomPID_.clear();
  eleGenMomPID_.clear();
  eleGenMomPt_.clear();
  eleIsoTrkDR03_.clear();
  eleIsoEcalDR03_.clear();
  eleIsoHcalDR03_.clear();
  eleIsoHcalDR0312_.clear();
  eleIsoTrkDR04_.clear();
  eleIsoEcalDR04_.clear();
  eleIsoHcalDR04_.clear();
  eleIsoHcalDR0412_.clear();
  eleModIsoTrk_.clear();
  eleModIsoEcal_.clear();
  eleModIsoHcal_.clear();
  eleChi2NDF_.clear();
  eleMissHits_.clear();
  eleConvDist_.clear();
  eleConvDcot_.clear();
  eleConvVtxFit_.clear();
  eleIP3D_.clear();
  eleIP3DErr_.clear();
  eleIDMVANonTrig_.clear();
  eleIDMVATrig_.clear();
  eleID2012_0_.clear();
  eleID2012_1_.clear();
  eleID2012_2_.clear();
  eleID2012_3_.clear();
  eleESEffSigmaRR_x_.clear();
  eleESEffSigmaRR_y_.clear();
  eleESEffSigmaRR_z_.clear();
  eleNBC_.clear();
  eleBrLinear_.clear();
  eleCetaCorrE_.clear();
  eleCetaCorrEt_.clear();
  eleBremCorrE_.clear();
  eleBremCorrEt_.clear();
  eleFullCorrE_.clear();
  eleFullCorrEt_.clear();
  elePFChIso03_.clear();
  elePFPhoIso03_.clear();
  elePFNeuIso03_.clear();
  elePFChIso04_.clear();
  elePFPhoIso04_.clear();
  elePFNeuIso04_.clear();
  
  // photon vectors
  phoTrg_.clear();
  phoTrgFilter_.clear();
  phoIsPhoton_.clear();
  phoE_.clear();
  phoEt_.clear();
  phoEta_.clear();
  phoPhi_.clear();
  phoEtVtx_.clear();
  phoEtaVtx_.clear();
  phoPhiVtx_.clear();
  phoVtx_x_.clear();
  phoVtx_y_.clear();
  phoVtx_z_.clear();
  phoSCPos_x_.clear();
  phoSCPos_y_.clear();
  phoSCPos_z_.clear();
  phoCaloPos_x_.clear(); // EXTRA
  phoCaloPos_y_.clear(); // EXTRA
  phoCaloPos_z_.clear(); // EXTRA
  phoR9_.clear();
  phoCetaCorrE_.clear();
  phoCetaCorrEt_.clear();
  phoBremCorrE_.clear();
  phoBremCorrEt_.clear();
  phoFullCorrE_.clear();
  phoFullCorrEt_.clear();
  phoTrkIsoSolidDR03_.clear();
  phoTrkIsoHollowDR03_.clear();
  phoEcalIsoDR03_.clear();
  phoHcalIsoSolidDR03_.clear();
  phoHcalIsoDR03_.clear();
  phoHcalIsoDR0312_.clear();
  phoTrkIsoSolidDR04_.clear();
  phoTrkIsoHollowDR04_.clear();
  phoEcalIsoDR04_.clear();
  phoHcalIsoDR04_.clear();
  phoHcalIsoDR0412_.clear();
  phoHcalIsoSolidDR04_.clear();
  
  phoCiCTrkIsoDR03_.clear();
  phoCiCTrkIsoDR04_.clear();
  phoCiCdRtoTrk_.clear();
  phoHoverE_.clear();
  phoHoverEBCdepth1_.clear();
  phoHoverEBCdepth2_.clear();
  phoHoverE12_.clear();
  phoSigmaIEtaIEta_.clear();
  phoSigmaIEtaIPhi_.clear();
  phoSigmaIPhiIPhi_.clear();
  phoEmax_.clear();
  phoE3x3_.clear();
  phoE5x5_.clear();
  phoE2x5Max_.clear();
  phoE5x1_.clear();
  phoE1x5_.clear();
  phoE3x1_.clear();
  phoE1x3_.clear();
  phoE2x2_.clear();
  phoRegrE_.clear();
  phoRegrEerr_.clear();
  phoPFChIso_.clear();
  phoPFPhoIso_.clear();
  phoPFNeuIso_.clear();
  phoSCRChIso_.clear();
  phoSCRPhoIso_.clear();
  phoSCRNeuIso_.clear();
  phoSCRChIso04_.clear();
  phoSCRPhoIso04_.clear();
  phoSCRNeuIso04_.clear();
  phoRandConeChIso_.clear();
  phoRandConePhoIso_.clear();
  phoRandConeNeuIso_.clear();
  phoRandConeChIso04_.clear();
  phoRandConePhoIso04_.clear();
  phoRandConeNeuIso04_.clear();
  phoSeedTime_.clear();
  phoLICTD_.clear();
  phoCiCPF4phopfIso005_.clear();
  phoCiCPF4phopfIso01_.clear();
  phoCiCPF4phopfIso02_.clear();
  phoCiCPF4phopfIso03_.clear();
  phoCiCPF4phopfIso04_.clear();
  phoCiCPF4phopfIso05_.clear();
  phoCiCPF4phopfIso06_.clear();
  phoCiCPF4phopfIso07_.clear();
  phoCiCPF4phopfIso08_.clear();
  phoCiCPF4chgpfIso005_.clear();
  phoCiCPF4chgpfIso01_.clear();
  phoCiCPF4chgpfIso02_.clear();
  phoCiCPF4chgpfIso03_.clear();
  phoCiCPF4chgpfIso04_.clear();
  phoCiCPF4chgpfIso05_.clear();
  phoCiCPF4chgpfIso06_.clear();
  phoCiCPF4chgpfIso07_.clear();
  phoCiCPF4chgpfIso08_.clear();
  
  phoCiCPF4phopfIsoNoVETO005_.clear();
  phoCiCPF4phopfIsoNoVETO01_.clear();
  phoCiCPF4phopfIsoNoVETO02_.clear();
  phoCiCPF4phopfIsoNoVETO03_.clear();
  phoCiCPF4phopfIsoNoVETO04_.clear();
  phoCiCPF4phopfIsoNoVETO05_.clear();
  phoCiCPF4phopfIsoNoVETO06_.clear();
  phoCiCPF4phopfIsoNoVETO07_.clear();
  phoCiCPF4phopfIsoNoVETO08_.clear();
  phoCiCPF4chgpfIsoNoVETO005_.clear();
  phoCiCPF4chgpfIsoNoVETO01_.clear();
  phoCiCPF4chgpfIsoNoVETO02_.clear();
  phoCiCPF4chgpfIsoNoVETO03_.clear();
  phoCiCPF4chgpfIsoNoVETO04_.clear();
  phoCiCPF4chgpfIsoNoVETO05_.clear();
  phoCiCPF4chgpfIsoNoVETO06_.clear();
  phoCiCPF4chgpfIsoNoVETO07_.clear();
  phoCiCPF4chgpfIsoNoVETO08_.clear();
  
  phoCiCPF4chgpfIso02AppVeto_.clear();
  phoCiCPF4chgpfIso03AppVeto_.clear();
  phoCiCPF4chgpfIso04AppVeto_.clear();
  phoCiCPF4phopfIso03AppVeto_.clear();
  phoCiCPF4phopfIso04AppVeto_.clear();
  phoCiCPF4phopfIso03AppVetoClean_.clear();
  phoCiCPF4phopfIso04AppVetoClean_.clear();
  phoCiCPF4chgpfIso03Mod_.clear();
  phoCiCPF4chgpfIso04Mod_.clear();
  
  phoSeedDetId1_.clear();
  phoSeedDetId2_.clear();
  phoRecoFlag_.clear();
  phoPos_.clear();
  phoGenIndex_.clear();
  phoGenGMomPID_.clear();
  phoGenMomPID_.clear();
  phoGenMomPt_.clear();
  phoSCE_.clear();
  phoSCRawE_.clear();
  phoESEn_.clear();
  phoSCEt_.clear();
  phoSCEta_.clear();
  phoSCPhi_.clear();
  phoSCEtaWidth_.clear();
  phoSCPhiWidth_.clear();
  phoSCBrem_.clear();
  phoOverlap_.clear();
  phohasPixelSeed_.clear();
  phoIsConv_.clear();
  phoEleVeto_.clear();
  phoESEffSigmaRR_x_.clear();
  phoESEffSigmaRR_y_.clear();
  phoESEffSigmaRR_z_.clear();
  phoNConv_.clear();
  phoConvNTrks_.clear();
  phoConvInvMass_.clear();
  phoConvCotTheta_.clear();
  phoConvEoverP_.clear();
  phoConvZofPVfromTrks_.clear();
  phoConvMinDist_.clear();
  phoConvdPhiAtVtx_.clear();
  phoConvdPhiAtCalo_.clear();
  phoConvdEtaAtCalo_.clear();
  phoConvTrkd0_x_.clear();
  phoConvTrkd0_y_.clear();
  phoConvTrkPin_x_.clear();
  phoConvTrkPin_y_.clear();
  phoConvTrkPout_x_.clear();
  phoConvTrkPout_y_.clear();
  phoConvTrkdz_x_.clear();
  phoConvTrkdz_y_.clear();
  phoConvTrkdzErr_x_.clear();
  phoConvTrkdzErr_y_.clear();
  phoConvChi2_.clear();
  phoConvChi2Prob_.clear();
  phoConvCharge1_.clear();
  phoConvCharge2_.clear();
  phoConvValidVtx_.clear();
  phoConvLikeLihood_.clear();
  phoConvP4_0_.clear();
  phoConvP4_1_.clear();
  phoConvP4_2_.clear();
  phoConvP4_3_.clear();
  phoConvVtx_x_.clear();
  phoConvVtx_y_.clear();
  phoConvVtx_z_.clear();
  phoConvVtxErr_x_.clear();
  phoConvVtxErr_y_.clear();
  phoConvVtxErr_z_.clear();
  phoConvPairMomentum_x_.clear();
  phoConvPairMomentum_y_.clear();
  phoConvPairMomentum_z_.clear();
  phoConvRefittedMomentum_x_.clear();
  phoConvRefittedMomentum_y_.clear();
  phoConvRefittedMomentum_z_.clear();
  SingleLegConv_.clear();
  phoPFConvRefitMom_x_.clear();
  phoPFConvRefitMom_y_.clear();
  phoPFConvRefitMom_z_.clear();
  phoPFConvMom_x_.clear();
  phoPFConvMom_y_.clear();
  phoPFConvMom_z_.clear();
  phoPFConvVtx_x_.clear();
  phoPFConvVtx_y_.clear();
  phoPFConvVtx_z_.clear();
  
  //new PF Variables stored when Matched to Reco:
  PFRecoMatch_.clear();
  PFEleMatch_.clear();
  PFEleVeto_.clear();    
  PFPreShowerE1_.clear();
  PFPreShowerE2_.clear();
  MustacheEin_.clear();
  MustacheEOut_.clear();
  MustacheEtOut_.clear();
  PFLowestClustE_.clear();
  PFClustdEta_.clear();
  PFClustdPhi_.clear();
  PFClustRMSPhi_.clear();
  PFClustRMSPhiMust_.clear();
  PFClustEneCorr_.clear();
  pho_hasSLConvPf_.clear();
  pho_hasConvPf_.clear();
  pho_pfconvVtxZ_.clear();
  pho_pfconvVtxZErr_.clear();
  pho_nSLConv_.clear();
  pho_pfSLConvPos_x_.clear();
  pho_pfSLConvPos_y_.clear();
  pho_pfSLConvPos_z_.clear();
  pho_pfSLConvVtxZ_.clear();

  muTrg_.clear();
  muEta_.clear();
  muPhi_.clear();
  muCharge_.clear();
  muPt_.clear();
  muPz_.clear();
  muVtx_x_.clear();
  muVtx_y_.clear();
  muVtx_z_.clear();
  muVtxGlb_x_.clear();
  muVtxGlb_y_.clear();
  muVtxGlb_z_.clear();
  muGenIndex_.clear();
  mucktPt_.clear();
  mucktPtErr_.clear();
  mucktEta_.clear();
  mucktPhi_.clear();
  mucktdxy_.clear();
  mucktdz_.clear();
  muIsoTrk_.clear();
  muIsoCalo_.clear();
  muIsoEcal_.clear();
  muIsoHcal_.clear();
  muChi2NDF_.clear();
  muInnerChi2NDF_.clear();
  muPFIsoR04_CH_.clear();
  muPFIsoR04_NH_.clear();
  muPFIsoR04_Pho_.clear();
  muPFIsoR04_PU_.clear();
  muPFIsoR04_CPart_.clear();
  muPFIsoR04_NHHT_.clear();
  muPFIsoR04_PhoHT_.clear();
  muPFIsoR03_CH_.clear();
  muPFIsoR03_NH_.clear();
  muPFIsoR03_Pho_.clear();
  muPFIsoR03_PU_.clear();
  muPFIsoR03_CPart_.clear();
  muPFIsoR03_NHHT_.clear();
  muPFIsoR03_PhoHT_.clear();
  muType_.clear();
  muD0_.clear();
  muDz_.clear();
  muD0GV_.clear();
  muDzGV_.clear();
  muD0Vtx_.clear();
  muDzVtx_.clear();
  muInnerD0_.clear();
  muInnerDz_.clear();
  muInnerD0GV_.clear();
  muInnerDzGV_.clear();
  muInnerPt_.clear();
  muInnerPtErr_.clear();
  muNumberOfValidTrkLayers_.clear();
  muNumberOfValidTrkHits_.clear();
  muNumberOfValidPixelLayers_.clear();
  muNumberOfValidPixelHits_.clear();
  muNumberOfValidMuonHits_.clear();
  muStations_.clear();
  muChambers_.clear();
  muIP3D_.clear();
  muIP3DErr_.clear();

  // taus
  tauDecayModeFinding_.clear();
  tauAgainstElectronLooseMVA3_.clear();
  tauAgainstElectronMediumMVA3_.clear();
  tauAgainstElectronTightMVA3_.clear();
  tauAgainstElectronVTightMVA3_.clear();
  tauAgainstElectronDeadECAL_.clear();
  tauAgainstMuonLoose2_.clear();
  tauAgainstMuonMedium2_.clear();
  tauAgainstMuonTight2_.clear();
  tauCombinedIsolationDeltaBetaCorrRaw3Hits_.clear();
  tauLooseCombinedIsolationDeltaBetaCorr3Hits_.clear();
  tauMediumCombinedIsolationDeltaBetaCorr3Hits_.clear();
  tauTightCombinedIsolationDeltaBetaCorr3Hits_.clear();
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


  // pfphotons
  PFPhoE_.clear();
  PFPhoEt_.clear();
  PFPhoEta_.clear();
  PFPhoPhi_.clear();  
  PFPhoType_.clear();
  PFPhoIso_.clear();

  // SubJet
  CA8JetPt_.clear();
  CA8JetEta_.clear();
  CA8JetPhi_.clear();
  CA8JetMass_.clear();
  CA8JetArea_.clear();
  CA8Jet_tau1_.clear();
  CA8Jet_tau2_.clear();
  CA8Jet_tau3_.clear();
  CA8prunedJetMass_.clear();
  CA8prunedJet_nSubJets_.clear();
  CA8prunedJet_SubjetPt_.clear();
  CA8prunedJet_SubjetEta_.clear();
  CA8prunedJet_SubjetPhi_.clear();
  CA8prunedJet_SubjetMass_.clear();

  // jets
  jetTrg_.clear();
  jetAlgo_.clear();
  jetEn_.clear();
  jetPt_.clear();
  jetEta_.clear();
  jetPhi_.clear();
  jetEt_.clear();
  jetRawPt_.clear();
  jetRawEn_.clear();
  jetCharge_.clear();
  jetArea_.clear();
  jetCHF_.clear();
  jetNHF_.clear();
  jetCEF_.clear();
  jetNEF_.clear();
  jetNCH_.clear();
  jetHFHAE_.clear();
  jetHFEME_.clear();
  jetPartonID_.clear();
  jetNConstituents_.clear();
  jetCombinedSecondaryVtxBJetTags_.clear();
  jetCombinedSecondaryVtxMVABJetTags_.clear();
  jetJetProbabilityBJetTags_.clear();
  jetJetBProbabilityBJetTags_.clear();
  jetBetaStar_.clear();
  jetGenJetIndex_.clear();
  jetGenJetEn_.clear();
  jetGenJetPt_.clear();
  jetGenJetEta_.clear();
  jetGenJetPhi_.clear();
  jetGenPartonID_.clear();
  jetGenEn_.clear();
  jetGenPt_.clear();
  jetGenEta_.clear();
  jetGenPhi_.clear();
  jetGenPartonMomID_.clear();
  //JetIDMVAvariable.clear();
  jetMVAs_.clear();
  jetWPLevels_.clear();

  //pujetIDalgos_.clear();//EXTRA
  jetDRMean_.clear();
  jetDR2Mean_.clear();
  jetDZ_.clear();
  jetFrac01_.clear();
  jetFrac02_.clear();
  jetFrac03_.clear();
  jetFrac04_.clear();
  jetFrac05_.clear();
  jetFrac06_.clear();
  jetFrac07_.clear();
  jetBeta_.clear();
  jetBetaStarCMG_.clear();
  jetBetaStarClassic_.clear();
  jetBetaExt_.clear();
  jetBetaStarCMGExt_.clear();
  jetBetaStarClassicExt_.clear();
  jetNNeutrals_.clear();
  jetNCharged_.clear();
  jetPFLooseId_.clear();
  //b-jetregressionvariable.clear();
  jetMt_.clear();
  jetJECUnc_.clear();
  jetLeadTrackPt_.clear();
  jetVtxPt_.clear();
  jetVtxMass_.clear();
  jetVtx3dL_.clear();
  jetVtx3deL_.clear();
  jetSoftLeptPt_.clear();
  jetSoftLeptPtRel_.clear();
  jetSoftLeptdR_.clear();
  jetSoftLeptIdlooseMu_.clear();
  jetSoftLeptIdEle95_.clear();
  jetDPhiMETJet_.clear();
  jetPuJetIdL_.clear();
  jetPuJetIdM_.clear();
  jetPuJetIdT_.clear();

  jetLowPtEn_.clear();
  jetLowPtPt_.clear();
  jetLowPtEta_.clear();
  jetLowPtPhi_.clear();
  jetLowPtCharge_.clear();
  jetLowPtEt_.clear();
  jetLowPtRawPt_.clear();
  jetLowPtRawEn_.clear();
  jetLowPtArea_.clear();
  jetLowPtGenJetEn_.clear();
  jetLowPtGenJetPt_.clear();
  jetLowPtGenJetEta_.clear();
  jetLowPtGenJetPhi_.clear();
  jetLowPtPartonID_.clear();
  jetLowPtGenPartonID_.clear();
  jetLowPtGenEn_.clear();
  jetLowPtGenPt_.clear();
  jetLowPtGenEta_.clear();
  jetLowPtGenPhi_.clear();

  // conversions
  convBarrel_.clear();
  convScInd_.clear();
  convP4_x_.clear();
  convP4_y_.clear();
  convP4_z_.clear();
  convP4_E_.clear();
  convVtx_x_.clear();
  convVtx_y_.clear();
  convVtx_z_.clear();
  convVtxErr_x_.clear();
  convVtxErr_y_.clear();
  convVtxErr_z_.clear();
  convPairMomentum_x_.clear();
  convPairMomentum_y_.clear();
  convPairMomentum_z_.clear();
  convRefittedMomentum_x_.clear();
  convRefittedMomentum_y_.clear();
  convRefittedMomentum_z_.clear();
  convNTracks_.clear();
  convPairInvMass_.clear();
  convPairCotThetaSep_.clear();
  convEoverP_.clear();
  convDistOfMinApproach_.clear();
  convDPhiTrksAtVtx_.clear();
  convDPhiTrksAtEcal_.clear();
  convDEtaTrksAtEcal_.clear();
  convDxy_.clear();
  convDz_.clear();
  convLxy_.clear();
  convLz_.clear();
  convZofPrimVtxFromTrks_.clear();
  convNHitsBeforeVtx_0_.clear();
  convNHitsBeforeVtx_1_.clear();
  convNSharedHits_.clear();
  convValidVtx_.clear();
  convMVALikelihood_.clear();
  convChi2_.clear();
  convChi2Probability_.clear();
  convTk1Dz_.clear();
  convTk2Dz_.clear();
  convTk1DzErr_.clear();
  convTk2DzErr_.clear();
  convCh1Ch2_.clear();
  convTk1D0_.clear();
  convTk1Pout_.clear();
  convTk1Pin_.clear();
  convTk2D0_.clear();
  convTk2Pout_.clear();
  convTk2Pin_.clear();

}

void ggNtuplizer::getHandles(edm::Event & event,
			     edm::Handle<std::vector<reco::GenParticle> > & genParticlesHandle,
			     edm::Handle<VertexCollection> &            recVtxs,
			     edm::Handle<VertexCollection> &            recVtxsBS,
			     edm::Handle<TriggerResults> &              trgResultsHandle,
			     edm::Handle<TriggerEvent> &                triggerEvent,
			     edm::Handle<TrackCollection> &             tracksHandle,
			     edm::Handle<GsfElectronCollection> &       gsfElectronHandle,
			     edm::Handle<edm::View<pat::MET> > &        METHandle,
			     edm::Handle<edm::View<pat::MET> > &        pfMETHandle,
			     edm::Handle<edm::View<pat::MET> > &        pfType01METHandle,
			     edm::Handle<reco::PFMETCollection> &       recoPfMETHandle,
			     edm::Handle<edm::View<pat::Electron> > &   electronHandle,
			     edm::Handle<edm::View<pat::Photon> > &     photonHandle,
			     edm::Handle<reco::PhotonCollection> &      recoPhotonHandle,
			     edm::Handle<edm::View<pat::Muon> > &       muonHandle,
			     edm::Handle<vector<pat::Tau> > &           tauHandle,
			     edm::Handle<edm::View<pat::Jet> > &        jetHandle,
			     edm::Handle<BeamSpot> &                    beamSpotHandle,
			     edm::Handle<EcalRecHitCollection> &        EBReducedRecHits,
			     edm::Handle<EcalRecHitCollection> &        EEReducedRecHits,
			     edm::Handle<EcalRecHitCollection> &        ESRecHits, 
			     edm::Handle<CaloTowerCollection> &         towers,
                             edm::Handle<PFCandidateCollection>&        pfAllCandidates,
                             edm::Handle<PFCandidateCollection>&        pfCandidates,
			     edm::Handle<PFCandidateCollection>&        pfCandidatePhotons,
                             edm::Handle<reco::ConversionCollection>&   convH,
			     edm::Handle<reco::PhotonCollection>& pfPhoTranslator
			     // edm::Handle<reco::SuperClusterCollection>& pfPhoSClusters,
			     //  edm::Handle<reco::SuperClusterCollection>& pfEleSClusters,
			     // edm::Handle<reco::PhotonCoreCollection>& pfPhotonCore,
			     //     edm::Handle<reco::GsfElectronCoreCollection>& pfElectronCore_

			     ) {

  if (doGenParticles_) event.getByLabel(genParticlesCollection_, genParticlesHandle);
  event.getByLabel(vtxlabel_              , recVtxs);
  event.getByLabel(recVtxsBSLabel_        , recVtxsBS);
  event.getByLabel(trgResults_            , trgResultsHandle_);
  event.getByLabel(trgEvent_              , triggerEvent);
  event.getByLabel(tracklabel_            , tracksHandle);
  event.getByLabel(gsfElectronlabel_      , gsfElectronHandle);
  event.getByLabel(METCollection_         , METHandle);
  event.getByLabel(pfMETlabel_            , pfMETHandle);
  event.getByLabel(pfType01METlabel_      , pfType01METHandle);
  event.getByLabel(recoPfMETlabel_        , recoPfMETHandle);
  event.getByLabel(electronCollection_    , electronHandle);
  event.getByLabel(photonCollection_      , photonHandle);
  event.getByLabel(recophotonCollection_  , recoPhotonHandle);
  event.getByLabel(muonCollection_        , muonHandle);
  event.getByLabel(tauCollection_         , tauHandle);
  event.getByLabel(jetCollection_         , jetHandle);
  event.getByLabel(beamSpotCollection_    , beamSpotHandle);
  event.getByLabel(ebReducedRecHitCollection_, EBReducedRecHits);
  event.getByLabel(eeReducedRecHitCollection_, EEReducedRecHits);
  event.getByLabel(esRecHitCollection_       , ESRecHits);
  event.getByLabel(towerCollection_          , towers);
  event.getByLabel(pfAllParticles_           , pfAllCandidates);
  event.getByLabel(pfParticles_              , pfCandidates);
  event.getByLabel(pfPhotonCollection_       , pfPhoTranslator); 
  event.getByLabel(pfParticles_              , pfCandidatePhotons);
  event.getByLabel(allConversionsColl_       , convH);
}

void ggNtuplizer::produce(edm::Event & e, const edm::EventSetup & es) {

  //PhotonFix::initialiseGeometry(es);
  hEvents_->Fill(0.5);

  //hcalHelper->readEvent(const_cast<edm::Event &>(e));
  // hcalHelperPflow->readEvent(const_cast<edm::Event &>(e));

  this->getHandles(e,
		   genParticlesHandle_,
		   recVtxs_,
		   recVtxsBS_,
		   trgResultsHandle_,
		   triggerEvent_,
		   tracksHandle_,
		   gsfElectronHandle_,
		   METHandle_,
		   pfMETHandle_,
		   pfType01METHandle_,
		   recoPfMETHandle_,
		   electronHandle_,
		   photonHandle_,
		   recoPhotonHandle_,
		   muonHandle_,
		   tauHandle_,
		   jetHandle_,
		   beamSpotHandle_,
		   EBReducedRecHits_,
		   EEReducedRecHits_,
		   ESRecHits_,
		   towers_,
                   pfAllCandidates,
                   pfCandidates,
		   pfCandidatePhotons,
                   convH_,
		   pfPhoTranslator_
		   //pfPhoSClusters_,
		   // pfEleSClusters_,
		   //  pfPhotonCore_,
		   //   pfElectronCore_
		   );

  // clear vectors
  clearVectors();

  // this is for trigger results, ugly, I know.
  ULong_t powerOfTwo[42];
  for(unsigned int i = 0; i != 42; ++i) 
    powerOfTwo[i] = pow(2, i);

  // get the iso deposits. 3 (charged hadrons, photons, neutral hadrons)                                                                   
  unsigned nTypes=3;
  IsoDepositMaps electronIsoDep(nTypes);
  for (size_t j = 0; j<inputTagIsoDepElectrons_.size(); ++j) {
    e.getByLabel(inputTagIsoDepElectrons_[j], electronIsoDep[j]);
  }

  IsoDepositMaps photonIsoDep(nTypes);
  for (size_t j = 0; j<inputTagIsoDepPhotons_.size(); ++j) {
    e.getByLabel(inputTagIsoDepPhotons_[j], photonIsoDep[j]);
  }

  IsoDepositVals electronIsoValPFId03(nTypes);
  IsoDepositVals electronIsoValPFId04(nTypes);
  IsoDepositVals photonIsoValPFId(nTypes);
  const IsoDepositVals * electronIsoVals03 = &electronIsoValPFId03;
  const IsoDepositVals * electronIsoVals04 = &electronIsoValPFId04;
  const IsoDepositVals * photonIsoVals = &photonIsoValPFId;

  for (size_t j = 0; j<inputTagIsoValElectronsPFId_.size(); ++j) {
    if (j < 3) e.getByLabel(inputTagIsoValElectronsPFId_[j], electronIsoValPFId03[j]);
    else e.getByLabel(inputTagIsoValElectronsPFId_[j], electronIsoValPFId04[j-3]);
  }

  for (size_t j = 0; j<inputTagIsoValPhotonsPFId_.size(); ++j) {
    e.getByLabel(inputTagIsoValPhotonsPFId_[j], photonIsoValPFId[j]);
  }

  // supercluster footprint removal
  SuperClusterFootprintRemoval scRemover03(e, es, scRemover03Pset_);
  SuperClusterFootprintRemoval scRemover04(e, es, scRemover04Pset_);

  // rho collection
  Handle<double> rhoHandle25;
  e.getByLabel(rhoCollection25_, rhoHandle25);
  rho25_ = *(rhoHandle25.product());

  Handle<double> rhoHandle25_neu;
  e.getByLabel(rhoCollection25_neu_, rhoHandle25_neu);
  rho25_neu_ = *(rhoHandle25_neu.product()) ;

  Handle<double> rhoHandle_elePFiso;
  e.getByLabel(rhoCollection25_eleLabel_, rhoHandle_elePFiso);
  rho25_elePFiso_ = *(rhoHandle_elePFiso.product());

  Handle<double> rhoHandle_muPFiso;
  e.getByLabel(rhoLepPFisoCollection_ ,rhoHandle_muPFiso);
  rho25_muPFiso_  = *(rhoHandle_muPFiso.product());

  edm::Handle<double> rho2011Handle;
  e.getByLabel(rho2011Label_, rho2011Handle);
  rho2011_ = *(rho2011Handle.product());

  Handle<double> rhoHandle_enReg;
  e.getByLabel(rho2012Label_, rhoHandle_enReg); 
  rho2012_ = *(rhoHandle_enReg.product());

  cicPhotonId_->configure(recVtxsBS_, tracksHandle_, gsfElectronHandle_, pfAllCandidates, rho2012_); 
  trackMET_->configure(recVtxsBS_, pfAllCandidates);

  lazyTool = new EcalClusterLazyTools(e, es, ebReducedRecHitCollection_, eeReducedRecHitCollection_);

  //const CaloTowerCollection* towers = towers_.product();
  //EgammaTowerIsolation hcaliso03(0.3, 0., 0, -1, towers);
  //EgammaTowerIsolation hcaliso04(0.4, 0., 0, -1, towers);

  // ES geometry
  ESHandle<CaloGeometry> geoHandle;
  es.get<CaloGeometryRecord>().get(geoHandle);
  const CaloSubdetectorGeometry *geometry = geoHandle->getSubdetectorGeometry(DetId::Ecal, EcalPreshower);
  const CaloSubdetectorGeometry *& geometry_p = geometry;

  CaloSubdetectorTopology *topology_p = 0;
  if (geometry) topology_p = new EcalPreshowerTopology(geoHandle);

  // ECAL Geometry
  edm::ESHandle<CaloGeometry> pG;
  es.get<CaloGeometryRecord>().get(pG);
  const CaloSubdetectorGeometry *geomBar = pG->getSubdetectorGeometry(DetId::Ecal,1);
  const CaloSubdetectorGeometry *geomEnd = pG->getSubdetectorGeometry(DetId::Ecal,2);

  // make the map of rechits
  rechits_map_.clear();
  if (ESRecHits_.isValid()) {
    EcalRecHitCollection::const_iterator it;
    for (it = ESRecHits_->begin(); it != ESRecHits_->end(); ++it) {
      // remove bad ES rechits
      if (it->recoFlag()==1 || it->recoFlag()==14 || (it->recoFlag()<=10 && it->recoFlag()>=5)) continue;
      //Make the map of DetID, EcalRecHit pairs
      rechits_map_.insert(std::make_pair(it->id(), *it));   
    }
  }

  edm::ESHandle<TransientTrackBuilder> builder;
  es.get<TransientTrackRecord>().get("TransientTrackBuilder", builder);
  TransientTrackBuilder thebuilder = *(builder.product());

  const TriggerMatchHelper matchHelper;

  if (doCentrality_) {
    if (!centProvider_) centProvider_ = new CentralityProvider(es);
    centProvider_->newEvent(e, es); 
    const reco::Centrality* centrality = centProvider_->raw();

    centrality_[0] = centrality->EtHFtowerSum();
    centrality_[1] = centrality->EtHFtowerSumPlus();
    centrality_[2] = centrality->EtHFtowerSumMinus();
    centrality_[3] = centrality->EtHFtruncatedPlus();
    centrality_[4] = centrality->EtHFtruncatedMinus();
  }

  run_    = e.id().run();
  event_  = e.id().event();
  lumis_  = e.luminosityBlock();
  isData_ = e.isRealData();

  // beam spot position 
  bspotPos_[0] = beamSpotHandle_->position().x();
  bspotPos_[1] = beamSpotHandle_->position().y();
  bspotPos_[2] = beamSpotHandle_->position().z();

  // vertex
  nVtx_ = 0;
  IsVtxGood_ = -1;
  Int_t firstGoodVtx = -1;
  nGoodVtx_ = 0;
  if (recVtxs_.isValid()) {
    for (size_t i=0; i<recVtxs_->size(); ++i) {
      if (!((*recVtxs_)[i].isFake())) {
	vtx_x_.push_back((*recVtxs_)[i].x());
	vtx_y_.push_back((*recVtxs_)[i].y());
	vtx_z_.push_back((*recVtxs_)[i].z());
	vtxNTrk_.push_back((*recVtxs_)[i].tracksSize());
	vtxNDF_.push_back((*recVtxs_)[i].ndof());
	vtxD0_.push_back((*recVtxs_)[i].position().rho());
	
	if ((*recVtxs_)[i].ndof() > 4 && 
	    fabs((*recVtxs_)[i].z()) <= 24  && 
	    fabs((*recVtxs_)[i].position().rho()) <= 2) {
	  nGoodVtx_++;
	  if (nGoodVtx_ == 1) firstGoodVtx = i;
	}
	nVtx_++;
      }
    }
  }
  if (nGoodVtx_ > 0) IsVtxGood_ = firstGoodVtx;

  // Set PV and first good vertex
  math::XYZPoint pv(vtx_x_[0], vtx_y_[0], vtx_z_[0]);
  if (firstGoodVtx < 0) firstGoodVtx = 0; 
  math::XYZPoint gv(vtx_x_[firstGoodVtx], vtx_y_[firstGoodVtx], vtx_z_[firstGoodVtx]);

  // vertex with beamspot
  vtxbsTkIndex_->clear();
  vtxbsTkWeight_->clear();

  nVtxBS_ = 0;
  if (recVtxsBS_.isValid()) {
    for (size_t i=0; i<recVtxsBS_->size(); ++i) {
      
      reco::VertexRef vtxRef(recVtxsBS_, i);
      
      float vtxbsSumPt2 = 0.;
      
      vtxbs_x_.push_back((*recVtxsBS_)[i].x());
      vtxbs_y_.push_back((*recVtxsBS_)[i].y());
      vtxbs_z_.push_back((*recVtxsBS_)[i].z());
      
      // porting H2ggolbe vertex info
      int nVtxBSTrk = 0;
      TVector3 vtxbsP_(0.,0.,0.);
      std::vector<Int_t>   tmp_itrk;
      std::vector<Float_t> tmp_weight;

      std::vector<reco::TrackBaseRef>::const_iterator tk;
      for (tk=(*recVtxsBS_)[i].tracks_begin(); tk!=(*recVtxsBS_)[i].tracks_end(); ++tk) {
	int itrk = 0;
	for (reco::TrackCollection::size_type j = 0; j<tracksHandle_->size(); ++j) {
	  reco::TrackRef track(tracksHandle_, j);
	  if (&(**tk) == &(*track)) {
	    const TVector3 tkPVec(track->px(),track->py(),track->pz());
	    TVector2 tkPtVec = tkPVec.XYvector();
	    vtxbsSumPt2 += tkPtVec.Mod2();

            tmp_itrk.push_back(static_cast<Int_t>(itrk));
            tmp_weight.push_back(static_cast<Float_t>(vtxRef->trackWeight(track)));

	    vtxbsP_ += tkPVec;
            
	    nVtxBSTrk++;
	    break;
	  }
	  itrk++;
	}
      }
      vtxbsTkIndex_->push_back(tmp_itrk);
      vtxbsTkWeight_->push_back(tmp_weight);

      vtxbsPtMod_.push_back(vtxbsP_.XYvector().Mod());
      vtxbsSumPt2_.push_back(vtxbsSumPt2);
      nVtxBS_++; 
    }
  }

  // track quality
  TrackBase::TrackQuality trkPurity_;
  if (tracksHandle_.isValid()) {

    const reco::TrackCollection *track = tracksHandle_.product();

    trkPurity_ = reco::TrackBase::qualityByName("highPurity");

    nTrk_ = 0;
    nGoodTrk_ = 0;
    for (reco::TrackCollection::const_iterator aTrk = track->begin(); aTrk != track->end(); ++aTrk) {

      if (dumpTrks_) {
	trkP_x_.push_back(aTrk->px()); 
	trkP_y_.push_back(aTrk->py()); 
	trkP_z_.push_back(aTrk->pz()); 

	trkVtx_x_.push_back(aTrk->vx());
	trkVtx_y_.push_back(aTrk->vy());
	trkVtx_z_.push_back(aTrk->vz());

	trkd0_.push_back(aTrk->d0());
	trkd0Err_.push_back(aTrk->d0Error());
	trkdz_.push_back(aTrk->dz());
	trkdzErr_.push_back(aTrk->dzError());
	trkPtErr_.push_back(aTrk->ptError());
	trkQuality_.push_back(aTrk->qualityMask());
      }

      if (aTrk->quality(trkPurity_)) nGoodTrk_++;
      nTrk_++;
    }
  }
  IsTracksGood_ = 0;
  if (nTrk_ > 10) {
    if (((float)nGoodTrk_/(float)nTrk_) > 0.25) IsTracksGood_ = 1;
  } else {
    IsTracksGood_ = 1;
  }

  // PDF information, pthat and processID
  nPUInfo_ = 0; 
  if (!isData_) {

    for (unsigned int i=0; i<7; ++i) pdf_[i] = -99;
    pthat_ = -99;
    processID_ = -99;

    Handle<GenEventInfoProduct> pdfInfoHandle;
    if (e.getByLabel(generatorLabel_, pdfInfoHandle)) {
      if (pdfInfoHandle->pdf()) {
	pdf_[0] = pdfInfoHandle->pdf()->id.first;    // PDG ID of incoming parton #1
	pdf_[1] = pdfInfoHandle->pdf()->id.second;   // PDG ID of incoming parton #2
	pdf_[2] = pdfInfoHandle->pdf()->x.first;     // x value of parton #1
	pdf_[3] = pdfInfoHandle->pdf()->x.second;    // x value of parton #2
	pdf_[4] = pdfInfoHandle->pdf()->xPDF.first;  // PDF weight for parton #1
	pdf_[5] = pdfInfoHandle->pdf()->xPDF.second; // PDF weight for parton #2
	pdf_[6] = pdfInfoHandle->pdf()->scalePDF;    // scale of the hard interaction
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

  // GenParticle
  if (!isData_ && genParticlesHandle_.isValid() ) {
    
    nMC_ = 0;
    int genIndex = 0;
    
    for (vector<GenParticle>::const_iterator ip = genParticlesHandle_->begin(); 
	 ip != genParticlesHandle_->end(); ++ip) {
      genIndex++;

      int status = ip->status() - 10*(ip->status()/10);
      
      bool stableFinalStateParticle = status == 1 && ip->pt() > 5.0;
      
      // keep all the photons with pT > 5.0 and all leptons;
      bool photonOrLepton = 
	(status == 1 && ip->pdgId() == 22 && ip->pt() > 5.0 ) ||
	(status == 1 && ( abs(ip->pdgId()) >= 11 && abs(ip->pdgId()) <= 16 ))  ||
	(status < 10 && abs(ip->pdgId()) == 15 );
      
      // select also Z, W, H, and top
      bool heavyParticle =  
	(ip->pdgId() == 23 || abs(ip->pdgId()) == 24 || ip->pdgId() == 25 || 
	 abs(ip->pdgId()) == 6 || abs(ip->pdgId()) == 5);
      
      if ( stableFinalStateParticle || heavyParticle || photonOrLepton ) {
	const Candidate *p = (const Candidate*)&(*ip);
	if (!runOnParticleGun_ && !p->mother()) continue;
	
	reco::GenParticleRef partRef = reco::GenParticleRef(genParticlesHandle_,
							    ip-genParticlesHandle_->begin());
	genpartparentage::GenParticleParentage particleHistory(partRef);
	
	mcPID    .push_back(p->pdgId());
	mcVtx_x  .push_back(p->vx());
	mcVtx_y  .push_back(p->vy());
	mcVtx_z  .push_back(p->vz());
	mcPt     .push_back(p->pt());
	mcMass   .push_back(p->mass());
	mcEta    .push_back(p->eta());
	mcPhi    .push_back(p->phi());
	mcE      .push_back(p->energy());
	mcEt     .push_back(p->et());
	mcParentage.push_back( 
			      particleHistory.hasLeptonParent()*16   + 
			      particleHistory.hasBosonParent()*8     + 
			      particleHistory.hasNonPromptParent()*4 +
			      particleHistory.hasQCDParent()*2       +
			      particleHistory.hasExoticParent());
	mcStatus.push_back(p->status());

	int mcDecayType_ = -999;
	// if genParticle is W or Z, check its decay type
	if ( ip->pdgId() == 23 || abs(ip->pdgId()) == 24 ) {
	  for (size_t k=0; k < p->numberOfDaughters(); ++k) {
	    const Candidate *dp = p->daughter(k);
	    if (abs(dp->pdgId())<=6)
	      mcDecayType_ = 1;
	    else if (abs(dp->pdgId())==11 || abs(dp->pdgId())==12)
	      mcDecayType_ = 2;
	    else if (abs(dp->pdgId())==13 || abs(dp->pdgId())==14)
	      mcDecayType_ = 3;
	    else if (abs(dp->pdgId())==15 || abs(dp->pdgId())==16)
	      mcDecayType_ = 4;
	  }
	}
	mcDecayType.push_back(mcDecayType_);

	int mcGMomPID_ = -999;
	int mcMomPID_  = -999;
	float mcMomPt_    = -999.;
	float mcMomMass_  = -999.;
	float mcMomEta_   = -999.;
	float mcMomPhi_   = -999.;
	if ( particleHistory.hasRealParent() ) {
	  reco::GenParticleRef momRef = particleHistory.parent();
	  if ( momRef.isNonnull() && momRef.isAvailable() ) {
	    mcMomPID_  = momRef->pdgId();
	    mcMomPt_   = momRef->pt();
	    mcMomMass_ = momRef->mass();
	    mcMomEta_  = momRef->eta();
	    mcMomPhi_  = momRef->phi();
	    
	    // get Granny
	    genpartparentage::GenParticleParentage motherParticle(momRef);
	    if ( motherParticle.hasRealParent() ) {
	      reco::GenParticleRef granny = motherParticle.parent();
	      mcGMomPID_ = granny->pdgId();
	    }
	  }
	}
	mcGMomPID.push_back(mcGMomPID_);
	mcMomPID.push_back(mcMomPID_);
	mcMomPt.push_back(mcMomPt_);
	mcMomMass.push_back(mcMomMass_);
	mcMomEta.push_back(mcMomEta_);
	mcMomPhi.push_back(mcMomPhi_);

	mcIndex.push_back(genIndex-1);
	nMC_++;
      } // save info on particles of interest
    } // loop over gen-level particles
  }

  // HLT
  bool changed (true);
  hltConfigProvider_.init(e.getRun(), es, trgResults_.process(), changed);
 
  for (int i=0; i<70; ++i) HLTIndex_[i] = -1;
 
  nHLT_ = 0;
  const TriggerNames &trgNames = e.triggerNames(*trgResultsHandle_);
  vector<string> hlNames = trgNames.triggerNames();
  nHLT_ = trgNames.size();

  for (size_t i=0; i<trgNames.size(); ++i) {
    HLT_[i] = (trgResultsHandle_->accept(i) == true) ? hltConfigProvider_.prescaleValue(e, es, hlNames[i]) : 0;
    if (hlNames[i].find("HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_R9Id85_OR_CaloId10_Iso50_Mass60_v") != string::npos) HLTIndex_[0] = i;
    else if (hlNames[i].find("HLT_Photon26_Photon18_v") != string::npos) HLTIndex_[1] = i;
    else if (hlNames[i].find("HLT_Photon26_CaloId10_Iso50_Photon18_CaloId10_Iso50_Mass60_v") != string::npos) HLTIndex_[2] = i;
    else if (hlNames[i].find("HLT_Photon26_CaloId10_Iso50_Photon18_R9Id85_Mass60_v") != string::npos) HLTIndex_[3] = i;
    else if (hlNames[i].find("HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_v") != string::npos) HLTIndex_[4] = i;
    else if (hlNames[i].find("HLT_Photon26_R9Id85_Photon18_CaloId10_Iso50_Mass60_v") != string::npos) HLTIndex_[5] = i;
    else if (hlNames[i].find("HLT_Photon26_R9Id85_Photon18_R9Id85_Mass60_v") != string::npos) HLTIndex_[6] = i;
    else if (hlNames[i].find("HLT_Photon26_R9Id85_OR_CaloId10_Iso50_Photon18_R9Id85_OR_CaloId10_Iso50_Mass70_v") != string::npos) HLTIndex_[7] = i;
    else if (hlNames[i].find("HLT_Photon36_R9Id85_OR_CaloId10_Iso50_Photon22_R9Id85_OR_CaloId10_Iso50_v") != string::npos) HLTIndex_[8] = i;
    else if (hlNames[i].find("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v") != string::npos) HLTIndex_[9] = i;
    else if (hlNames[i].find("HLT_Ele27_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele15_CaloIdT_CaloIsoVL_trackless_v") != string::npos) HLTIndex_[10] = i;
    else if (hlNames[i].find("HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_v") != string::npos) HLTIndex_[11] = i;
    else if (hlNames[i].find("HLT_Mu22_TkMu22_v") != string::npos) HLTIndex_[12] = i;
    else if (hlNames[i].find("HLT_Mu17_TkMu8_v") != string::npos) HLTIndex_[13] = i;
    else if (hlNames[i].find("HLT_Mu17_Mu8_v") != string::npos) HLTIndex_[14] = i;
    else if (hlNames[i].find("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v") != string::npos) HLTIndex_[15] = i;
    else if (hlNames[i].find("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v") != string::npos) HLTIndex_[16] = i;
    else if (hlNames[i].find("HLT_Ele27_WP80_v") != string::npos) HLTIndex_[17] = i;
    else if (hlNames[i].find("HLT_IsoMu24_eta2p1_v") != string::npos) HLTIndex_[18] = i;
    else if (hlNames[i].find("HLT_IsoMu24_v") != string::npos) HLTIndex_[19] = i;
    else if (hlNames[i].find("HLT_Mu22_TkMu8_v") != string::npos) HLTIndex_[20] = i;
    else if (hlNames[i].find("HLT_Mu22_Mu8_v") != string::npos) HLTIndex_[21] = i;
    else if (hlNames[i].find("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v") != string::npos) HLTIndex_[22] = i;
    else if (hlNames[i].find("HLT_Photon36_Photon22_v") != string::npos) HLTIndex_[23] = i;
    else if (hlNames[i].find("HLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50_v") != string::npos) HLTIndex_[24] = i;
    else if (hlNames[i].find("HLT_Photon36_R9Id85_OR_CaloId10_Iso50_Photon10_R9Id85_OR_CaloId10_Iso50_Mass80_v") != string::npos) HLTIndex_[25] = i;
    else if (hlNames[i].find("HLT_PAL1DoubleMuOpen_v") != string::npos) HLTIndex_[26] = i;
    else if (hlNames[i].find("HLT_DiJet35_MJJ650_AllJets_DEta3p5_VBF_v") != string::npos) HLTIndex_[27] = i;
    else if (hlNames[i].find("HLT_DiJet35_MJJ700_AllJets_DEta3p5_VBF_v") != string::npos) HLTIndex_[28] = i;
    else if (hlNames[i].find("HLT_PAMu12_v") != string::npos) HLTIndex_[29] = i;
    else if (hlNames[i].find("HLT_PAPhoton15_TightCaloIdVL_v") != string::npos) HLTIndex_[30] = i;
    else if (hlNames[i].find("HLT_Photon30_R9Id90_CaloId_HE10_Iso40_EBOnly_v") != string::npos) HLTIndex_[31] = i;
    else if (hlNames[i].find("HLT_Photon30_R9Id90_CaloId_HE10_Iso40_EBOnly_Met25_HBHENoiseCleaned_v") != string::npos) HLTIndex_[32] = i;
    else if (hlNames[i].find("HLT_Photon30_v") != string::npos) HLTIndex_[33] = i;
    else if (hlNames[i].find("HLT_PFJet40_v") != string::npos) HLTIndex_[34] = i;
    else if (hlNames[i].find("HLT_PFJet80_v") != string::npos) HLTIndex_[35] = i; 
    else if (hlNames[i].find("HLT_PFJet140_v") != string::npos) HLTIndex_[36] = i;
    else if (hlNames[i].find("HLT_PFJet200_v") != string::npos) HLTIndex_[37] = i; 
    else if (hlNames[i].find("HLT_PFJet260_v") != string::npos) HLTIndex_[38] = i; 
    else if (hlNames[i].find("HLT_Photon30_CaloIdVL_v") != string::npos) HLTIndex_[39] = i;
    else if (hlNames[i].find("HLT_Photon50_CaloIdVL_IsoL_v") != string::npos) HLTIndex_[40] = i;
    else if (hlNames[i].find("HLT_Photon50_CaloIdVL_v") != string::npos) HLTIndex_[41] = i;
    else if (hlNames[i].find("HLT_Photon75_CaloIdVL_v") != string::npos) HLTIndex_[42] = i; 
    else if (hlNames[i].find("HLT_Photon90_CaloIdVL_v") != string::npos) HLTIndex_[43] = i; 
    else if (hlNames[i].find("HLT_Photon135_v") != string::npos) HLTIndex_[44] = i;
    else if (hlNames[i].find("HLT_Photon150_v") != string::npos) HLTIndex_[45] = i; 
    else if (hlNames[i].find("HLT_Photon160_v") != string::npos) HLTIndex_[46] = i; 
    else if (hlNames[i].find("HLT_Mu22_Photon22_CaloIdL_v") != string::npos) HLTIndex_[47] = i; 
  }

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

  if (pfType01METHandle_.isValid()) {
    const pat::MET *pfType01MET = 0;
    pfType01MET = &(pfType01METHandle_->front());

    pfType01MET_       = pfType01MET->et();
    pfType01METPhi_    = pfType01MET->phi();
    pfType01METsumEt_  = pfType01MET->sumEt();
    pfType01METmEtSig_ = (pfType01MET->mEtSig() < 1.e10) ? pfType01MET->mEtSig() : 0;
    pfType01METSig_    = (pfType01MET->significance() < 1.e10) ? pfType01MET->significance() : 0;;
  }

  if (recoPfMETHandle_.isValid()) {
    recoPfMET_       = recoPfMETHandle_->begin()->et();
    recoPfMETPhi_    = recoPfMETHandle_->begin()->phi();
    recoPfMETsumEt_  = recoPfMETHandle_->begin()->sumEt();
    recoPfMETmEtSig_ = (recoPfMETHandle_->begin()->mEtSig() < 1.e10) ? recoPfMETHandle_->begin()->mEtSig() : 0;
    recoPfMETSig_    = (recoPfMETHandle_->begin()->significance() < 1.e10) ? recoPfMETHandle_->begin()->significance() : 0;;
  }

  // track MET
  std::vector<float> trkMETPV = trackMET_->trackMETWithPV(0, 0.2, 10000.);
  trkMETxPV_   = trkMETPV[0];
  trkMETyPV_   = trkMETPV[1];
  trkMETPhiPV_ = trkMETPV[2];
  trkMETPV_    = trkMETPV[3];

  std::vector<float> trkMET = trackMET_->trackMETWithVertex(0, 0.2, 10000.);
  for (unsigned int v=0; v<recVtxsBS_->size(); ++v) {
    trkMETx_.push_back(trkMET[v*4]);
    trkMETy_.push_back(trkMET[v*4+1]);
    trkMETPhi_.push_back(trkMET[v*4+2]);
    trkMET_.push_back(trkMET[v*4+3]);
  }

  for (unsigned int i=0; i<10; ++i) metFilters_[i] = -1;
  if (isData_) {
    // MET filters
    // [0]CSC Beam Halo [1]HBHE Noise [2]HCAL laser [3]ECAL dead cell [4]Tracking failure
    // [5]Bad EE SC [6]ECAL large laser calibration [7]Many strips [8]too many strips [9]log error too many clusters
    
    //Handle<BeamHaloSummary> theBeamHaloSummary;
    //e.getByLabel("BeamHaloSummary", theBeamHaloSummary);
    //const BeamHaloSummary theBHSummary = (*theBeamHaloSummary.product());
    //metFilters_[0] = (int) (theBHSummary.CSCTightHaloId());       
    metFilters_[0] = -99;
    
    Handle<bool> HBHENoiseFilter;
    if (e.getByLabel(HBHENoiseFilterLabel_, HBHENoiseFilter)) metFilters_[1] = (int) (*HBHENoiseFilter);
    
    Handle<bool> HcalLaserFilter;
    if (e.getByLabel(HcalLaserFilterLabel_, HcalLaserFilter)) metFilters_[2] = (int) (*HcalLaserFilter);
    
    Handle<bool> EcalDeadCellFilter;
    if (e.getByLabel(EcalDeadCellFilterLabel_, EcalDeadCellFilter)) metFilters_[3] = (int) (*EcalDeadCellFilter);
    
    Handle<bool> TrackingFailureFilter;
    if (e.getByLabel(TrackingFailureFilterLabel_, TrackingFailureFilter)) metFilters_[4] = (int) (*TrackingFailureFilter);
    
    Handle<bool> EEBadScFilter;
    if (e.getByLabel(EEBadScFilterLabel_, EEBadScFilter)) metFilters_[5] = (int) (*EEBadScFilter);
    
    Handle<bool> EcalLaserFilter;
    if (e.getByLabel(EcalLaserFilterLabel_, EcalLaserFilter)) metFilters_[6] = (int) (*EcalLaserFilter);
    
    Handle<bool> Manystripclus53X;
    if (e.getByLabel(Manystripclus53XLabel_, Manystripclus53X)) metFilters_[7] = (int) (*Manystripclus53X);
    
    Handle<bool> Toomanystripclus53X;
    if (e.getByLabel(Toomanystripclus53XLabel_, Toomanystripclus53X)) metFilters_[8] = (int) (*Toomanystripclus53X);
    
    Handle<bool> LogErrorTooManyClusters;
    if (e.getByLabel(LogErrorTooManyClustersLabel_, LogErrorTooManyClusters)) metFilters_[9] = (int) (*LogErrorTooManyClusters);
  } 

  // PF Candidates
  const PFCandidateCollection thePfColl = *(pfCandidates.product());
  
  PFIsolationEstimator isolator;
  isolator.initializePhotonIsolation(kTRUE);
  isolator.setConeSize(0.3);
  VertexRef myVtxRef(recVtxs_, 0);

  // Electron (note that trigger match 2,3,4 are useless for 8TeV for the time being
  const TriggerObjectMatch *eleTriggerMatch1(triggerEvent_->triggerObjectMatchResult("electronTriggerMatchHLTEle27WP80"));
  const TriggerObjectMatch *eleTriggerMatch2(triggerEvent_->triggerObjectMatchResult("electronTriggerMatchHLTPAPhoton15TightCaloIdVL"));  
  const TriggerObjectMatch *eleTriggerMatch3(triggerEvent_->triggerObjectMatchResult("electronTriggerMatchHLTEle17CaloIdLCaloIsoVLEle8CaloIdLCaloIsoVL"));
  const TriggerObjectMatch *eleTriggerMatch4(triggerEvent_->triggerObjectMatchResult("electronTriggerMatchHLTEle27WP80PFMT50"));
  const TriggerObjectMatch *eleTriggerMatch5(triggerEvent_->triggerObjectMatchResult("electronTriggerMatchEle17CaloIdTCaloIsoVLTrkIdVLTrkIsoVLEle8CaloIdTCaloIsoVLTrkIdVLTrkIsoVL"));
  const TriggerObjectMatch *eleTriggerMatch6(triggerEvent_->triggerObjectMatchResult("electronTriggerMatchEle17CaloIdTTrkIdVLCaloIsoVLTrkIsoVLEle8CaloIdTTrkIdVLCaloIsoVLTrkIsoVL"));
  const TriggerObjectMatch *eleTriggerMatch7(triggerEvent_->triggerObjectMatchResult("electronTriggerMatchHLTEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoFilter")); // matched to Ele17 lag
  const TriggerObjectMatch *eleTriggerMatch8(triggerEvent_->triggerObjectMatchResult("electronTriggerMatchHLTEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDoubleFilter")); // matched to Ele17 and Ele8
  const TriggerObjectMatch *eleTriggerMatch9(triggerEvent_->triggerObjectMatchResult("electronTriggerMatchHLTEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDZ")); // matched to Ele17 and Ele8 and pass DZ cut
  const TriggerObjectMatch *eleTriggerMatch10(triggerEvent_->triggerObjectMatchResult("electronTriggerMatchHLTEle20CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC4Mass50"));
  const TriggerObjectMatch *eleTriggerMatch11(triggerEvent_->triggerObjectMatchResult("electronFilterMatchEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVT"));
  const TriggerObjectMatch *eleTriggerMatch12(triggerEvent_->triggerObjectMatchResult("electronFilterMatchEle8Mass50"));
  const TriggerObjectMatch *eleTriggerMatch13(triggerEvent_->triggerObjectMatchResult("electronFilterMatchEle20CaloIdVTCaloIsoVTTrkIdTTrkIsoVT"));
  const TriggerObjectMatch *eleTriggerMatch14(triggerEvent_->triggerObjectMatchResult("electronFilterMatchSC4Mass50"));
  const TriggerObjectMatch *eleTriggerMatch15(triggerEvent_->triggerObjectMatchResult("electronFilterMatchEle32CaloIdTCaloIsoTTrkIdTTrkIsoT"));
  const TriggerObjectMatch *eleTriggerMatch16(triggerEvent_->triggerObjectMatchResult("electronFilterMatchSC17Mass50"));

  const TriggerObjectMatch *eleTriggerFMatch1(triggerEvent_->triggerObjectMatchResult("eleTrgFMatchEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVT"));
  const TriggerObjectMatch *eleTriggerFMatch2(triggerEvent_->triggerObjectMatchResult("eleTrgFMatchEle8Mass50"));
  const TriggerObjectMatch *eleTriggerFMatch3(triggerEvent_->triggerObjectMatchResult("eleTrgFMatchEle20CaloIdVTCaloIsoVTTrkIdTTrkIsoVT"));
  const TriggerObjectMatch *eleTriggerFMatch4(triggerEvent_->triggerObjectMatchResult("eleTrgFMatchSC4Mass50"));
  const TriggerObjectMatch *eleTriggerFMatch5(triggerEvent_->triggerObjectMatchResult("eleTrgFMatchEle32CaloIdTCaloIsoTTrkIdTTrkIsoT"));
  const TriggerObjectMatch *eleTriggerFMatch6(triggerEvent_->triggerObjectMatchResult("eleTrgFMatchSC17Mass50"));

  nEle_ = 0;
  const Candidate *elemom = 0;
  if ( electronHandle_.isValid() )
    for (View<pat::Electron>::const_iterator iEle = electronHandle_->begin(); iEle != electronHandle_->end(); ++iEle) {
      
      edm::RefToBase<pat::Electron> eleRef = electronHandle_->refAt(nEle_);
      reco::CandidateBaseRef eleBaseRef(eleRef);
      const TriggerObjectRef eleTrigRef1( matchHelper.triggerMatchObject( eleBaseRef, eleTriggerMatch1, e, *triggerEvent_ ) );
      const TriggerObjectRef eleTrigRef2( matchHelper.triggerMatchObject( eleBaseRef, eleTriggerMatch2, e, *triggerEvent_ ) );
      const TriggerObjectRef eleTrigRef3( matchHelper.triggerMatchObject( eleBaseRef, eleTriggerMatch3, e, *triggerEvent_ ) );
      const TriggerObjectRef eleTrigRef4( matchHelper.triggerMatchObject( eleBaseRef, eleTriggerMatch4, e, *triggerEvent_ ) );
      const TriggerObjectRef eleTrigRef5( matchHelper.triggerMatchObject( eleBaseRef, eleTriggerMatch5, e, *triggerEvent_ ) );
      const TriggerObjectRef eleTrigRef6( matchHelper.triggerMatchObject( eleBaseRef, eleTriggerMatch6, e, *triggerEvent_ ) );
      const TriggerObjectRef eleTrigRef7( matchHelper.triggerMatchObject( eleBaseRef, eleTriggerMatch7, e, *triggerEvent_ ) );
      const TriggerObjectRef eleTrigRef8( matchHelper.triggerMatchObject( eleBaseRef, eleTriggerMatch8, e, *triggerEvent_ ) );
      const TriggerObjectRef eleTrigRef9( matchHelper.triggerMatchObject( eleBaseRef, eleTriggerMatch9, e, *triggerEvent_ ) );
      const TriggerObjectRef eleTrigRef10( matchHelper.triggerMatchObject( eleBaseRef, eleTriggerMatch10, e, *triggerEvent_ ) );
      const TriggerObjectRef eleTrigRef11( matchHelper.triggerMatchObject( eleBaseRef, eleTriggerMatch11, e, *triggerEvent_ ) );
      const TriggerObjectRef eleTrigRef12( matchHelper.triggerMatchObject( eleBaseRef, eleTriggerMatch12, e, *triggerEvent_ ) );
      const TriggerObjectRef eleTrigRef13( matchHelper.triggerMatchObject( eleBaseRef, eleTriggerMatch13, e, *triggerEvent_ ) );
      const TriggerObjectRef eleTrigRef14( matchHelper.triggerMatchObject( eleBaseRef, eleTriggerMatch14, e, *triggerEvent_ ) );
      const TriggerObjectRef eleTrigRef15( matchHelper.triggerMatchObject( eleBaseRef, eleTriggerMatch15, e, *triggerEvent_ ) );
      const TriggerObjectRef eleTrigRef16( matchHelper.triggerMatchObject( eleBaseRef, eleTriggerMatch16, e, *triggerEvent_ ) );
      ULong_t results = 0;
      results +=(eleTrigRef1.isAvailable())  ? powerOfTwo[0]  : 0; // 2^0 - set 0th bit
      results +=(eleTrigRef2.isAvailable())  ? powerOfTwo[1]  : 0; // 2^1 - set 1st bit
      results +=(eleTrigRef3.isAvailable())  ? powerOfTwo[2]  : 0;
      results +=(eleTrigRef4.isAvailable())  ? powerOfTwo[3]  : 0;
      results +=(eleTrigRef5.isAvailable())  ? powerOfTwo[4]  : 0;
      results +=(eleTrigRef6.isAvailable())  ? powerOfTwo[5]  : 0;
      results +=(eleTrigRef7.isAvailable())  ? powerOfTwo[6]  : 0;
      results +=(eleTrigRef8.isAvailable())  ? powerOfTwo[7]  : 0;
      results +=(eleTrigRef9.isAvailable())  ? powerOfTwo[8]  : 0;
      results +=(eleTrigRef10.isAvailable()) ? powerOfTwo[9]  : 0;
      results +=(eleTrigRef11.isAvailable()) ? powerOfTwo[10] : 0;
      results +=(eleTrigRef12.isAvailable()) ? powerOfTwo[11] : 0;
      results +=(eleTrigRef13.isAvailable()) ? powerOfTwo[12] : 0;
      results +=(eleTrigRef14.isAvailable()) ? powerOfTwo[13] : 0;
      results +=(eleTrigRef15.isAvailable()) ? powerOfTwo[14] : 0;
      results +=(eleTrigRef16.isAvailable()) ? powerOfTwo[15] : 0;
      eleTrg_.push_back(results);

      edm::Ptr<reco::Candidate> recoEleRef = iEle->originalObjectRef();                                                                                  
      const reco::GsfElectron *recoElectron = dynamic_cast<const reco::GsfElectron *>(recoEleRef.get());        

      eleClass_            .push_back(iEle->classification());
      eleIsEcalDriven_     .push_back(iEle->ecalDrivenSeed());
      eleCharge_           .push_back(iEle->charge());
      eleChargeConsistent_ .push_back((Int_t)iEle->isGsfCtfScPixChargeConsistent());

      eleEn_           .push_back(recoElectron->energy()); // iEle->energy();
      eleEcalEn_       .push_back(iEle->ecalEnergy());
      elePt_           .push_back(recoElectron->pt()); // iEle->pt();
      eleEta_          .push_back(iEle->eta());
      elePhi_          .push_back(iEle->phi());
      eleR9_           .push_back(iEle->r9());
      eleHoverE_       .push_back(iEle->hadronicOverEm());
      eleHoverE12_     .push_back(iEle->hcalOverEcalBc());
      eleEoverP_       .push_back(iEle->eSuperClusterOverP());

      elePin_       .push_back(iEle->trackMomentumAtVtx().R());
      elePout_      .push_back(iEle->trackMomentumOut().R());
      eleTrkMomErr_ .push_back(iEle->trackMomentumError());
      eleBrem_      .push_back(iEle->fbrem());
      
      eleD0_   .push_back(iEle->gsfTrack()->dxy(pv));
      eleDz_   .push_back(iEle->gsfTrack()->dz(pv));
      eleD0GV_ .push_back(iEle->gsfTrack()->dxy(gv));
      eleDzGV_ .push_back(iEle->gsfTrack()->dz(gv));

      eledEtaAtVtx_ .push_back(iEle->deltaEtaSuperClusterTrackAtVtx());
      eledPhiAtVtx_ .push_back(iEle->deltaPhiSuperClusterTrackAtVtx());
      
      eleMissHits_   .push_back(iEle->gsfTrack()->trackerExpectedHitsInner().numberOfHits());
      eleConvDist_   .push_back(iEle->convDist());
      eleConvDcot_   .push_back(iEle->convDcot());

      eleConvVtxFit_ .push_back((int) ConversionTools::hasMatchedConversion(*recoElectron, convH_, beamSpotHandle_->position(), true, 2.0, 1e-6, 0));      

      // IP3D
      float eleIP3D = -999.;
      float eleIP3DErr = -99.;
      if (iEle->gsfTrack().isNonnull()) {
	const double gsfsign = ((- iEle->gsfTrack()->dxy(pv)) >= 0) ? 1. : -1.;
	const reco::TransientTrack tt_ele = thebuilder.build(iEle->gsfTrack());
	const std::pair<bool,Measurement1D> ip3dpv_ele = IPTools::absoluteImpactParameter3D(tt_ele, *recVtxs_->begin());
	if (ip3dpv_ele.first) {
	  eleIP3D    = gsfsign*ip3dpv_ele.second.value();
	  eleIP3DErr = ip3dpv_ele.second.error();
	}
      }
      eleIP3D_.push_back(eleIP3D);
      eleIP3DErr_.push_back(eleIP3DErr);

      // 2012 MVA eID
      eleIDMVANonTrig_ .push_back(iEle->electronID("mvaNonTrigV0"));
      eleIDMVATrig_    .push_back(iEle->electronID("mvaTrigV0"));

      // Access PF isolation
      elePFChIso03_  .push_back((*(*electronIsoVals03)[0])[recoEleRef]);
      elePFPhoIso03_ .push_back((*(*electronIsoVals03)[1])[recoEleRef]);
      elePFNeuIso03_ .push_back((*(*electronIsoVals03)[2])[recoEleRef]);

      elePFChIso04_  .push_back((*(*electronIsoVals04)[0])[recoEleRef]);
      elePFPhoIso04_ .push_back((*(*electronIsoVals04)[1])[recoEleRef]);
      elePFNeuIso04_ .push_back((*(*electronIsoVals04)[2])[recoEleRef]);

      // Access super cluster
      eleSCEta_      .push_back(iEle->superCluster()->eta());
      eleSCPhi_      .push_back(iEle->superCluster()->phi());
      eleSCRawEn_    .push_back(iEle->superCluster()->rawEnergy());
      eleSCEn_       .push_back(iEle->superCluster()->energy());
      eleESEn_       .push_back(iEle->superCluster()->preshowerEnergy());
      eleSCEtaWidth_ .push_back(iEle->superCluster()->etaWidth());
      eleSCPhiWidth_ .push_back(iEle->superCluster()->phiWidth());
      
      eleVtx_x_.push_back(iEle->trackPositionAtVtx().x());
      eleVtx_y_.push_back(iEle->trackPositionAtVtx().y());
      eleVtx_z_.push_back(iEle->trackPositionAtVtx().z());
      
      // Gen Particle
      int eleGenIndex = -1;
      int EleGenIndex_count = 0;
      int eleGenMomPID = 999;
      int eleGenGMomPID = 999;
      float eleGenMomPt  = -1.;

      if (!isData_) {
        if ((*iEle).genLepton() && genParticlesHandle_.isValid() ) {
	  
          for (vector<GenParticle>::const_iterator iGen = genParticlesHandle_->begin(); 
	       iGen != genParticlesHandle_->end(); ++iGen) {
	    
            if (iGen->p4() == (*iEle).genLepton()->p4() && 
		iGen->pdgId() == (*iEle).genLepton()->pdgId() && 
		iGen->status() == (*iEle).genLepton()->status()) {
	      
	      eleGenIndex = EleGenIndex_count;
	      
	      const Candidate *elep = (const Candidate*)&(*iGen);
	      
	      for (size_t j=0; j<elep->numberOfMothers(); ++j) {
		
		elemom = elep->mother(j);
		eleGenMomPID = elemom->pdgId();
		eleGenMomPt  = elemom->pt();
		if (elemom->mother()) eleGenGMomPID = elemom->mother()->pdgId();
	      }
            }
            EleGenIndex_count++;
          }
        }
      }
      eleGenIndex_.push_back(eleGenIndex);
      eleGenMomPID_.push_back(eleGenMomPID);
      eleGenGMomPID_.push_back(eleGenGMomPID);
      eleGenMomPt_.push_back(eleGenMomPt);
      
      eleIsoTrkDR03_       .push_back(iEle->dr03TkSumPt());
      eleIsoEcalDR03_      .push_back(iEle->dr03EcalRecHitSumEt());
      eleIsoHcalDR03_      .push_back(iEle->dr03HcalTowerSumEt());
      eleIsoHcalDR0312_    .push_back(iEle->dr03HcalDepth1TowerSumEt());
      eleIsoTrkDR04_       .push_back(iEle->dr04TkSumPt());
      eleIsoEcalDR04_      .push_back(iEle->dr04EcalRecHitSumEt());
      eleIsoHcalDR04_      .push_back(iEle->dr04HcalTowerSumEt());
      eleIsoHcalDR0412_    .push_back(iEle->dr04HcalDepth1TowerSumEt());

      eleModIsoTrk_        .push_back(iEle->userIso(0));
      eleModIsoEcal_       .push_back(iEle->userIso(1));
      eleModIsoHcal_       .push_back(iEle->userIso(2));
      
      const reco::CaloClusterPtr eleSeed = (*iEle).superCluster()->seed();
      eleNClus_.push_back((*iEle).superCluster()->clustersSize());
      vector<float> eleCov;
      eleCov = lazyTool->localCovariances(*eleSeed);
      eleSigmaIEtaIEta_.push_back(iEle->sigmaIetaIeta());
      eleSigmaIEtaIPhi_.push_back(eleCov[1]);
      eleSigmaIPhiIPhi_.push_back(eleCov[2]);
      
      eleEmax_     .push_back(lazyTool->eMax(*eleSeed));
      eleE2ndMax_  .push_back(lazyTool->e2nd(*eleSeed));
      eleETop_     .push_back(lazyTool->eTop(*eleSeed));
      eleEBottom_  .push_back(lazyTool->eBottom(*eleSeed));
      eleELeft_    .push_back(lazyTool->eLeft(*eleSeed));
      eleERight_   .push_back(lazyTool->eRight(*eleSeed));
      eleSeedE_  .push_back(eleSeed->energy());
      eleSeedEta_  .push_back(eleSeed->eta());
      eleSeedPhi_  .push_back(eleSeed->phi());
      eleE1x5_     .push_back(lazyTool->e1x5(*eleSeed));
      eleE3x3_     .push_back(lazyTool->e3x3(*eleSeed));
      eleE5x5_     .push_back(iEle->e5x5());
      eleE2x5Max_  .push_back(lazyTool->e2x5Max(*eleSeed));
      eleE2x5Left_  .push_back(lazyTool->e2x5Left(*eleSeed));
      eleE2x5Right_  .push_back(lazyTool->e2x5Right(*eleSeed));
      eleE2x5Top_  .push_back(lazyTool->e2x5Top(*eleSeed));     
      eleE2x5Bottom_  .push_back(lazyTool->e2x5Bottom(*eleSeed));   
      if (iEle->isEB()){
       float betacry, bphicry, bthetatilt, bphitilt;
       int bieta, biphi;
         _ecalLocal.localCoordsEB(*eleSeed,es,betacry,bphicry,bieta,biphi,bthetatilt,bphitilt);
         eleCrysEta_.push_back(betacry);
         eleCrysPhi_.push_back(bphicry);
         eleCrysIEta_.push_back(bieta);
         eleCrysIPhi_.push_back(biphi);
      }
       else{
         phoCrysEta_.push_back(-99);
         phoCrysPhi_.push_back(-99);
         phoCrysIEta_.push_back(-9999);
         phoCrysIPhi_.push_back(-9999);
      }
       
      // Energy Regression Correction
      if (!egCorrEle_.IsInitialized()) {   
	//std::string filenameEle = "http://homepages.spa.umn.edu/~rekovic/cms/regweights52xV3/gbrv3ele_52x.root";
	std::string filenameEle = "gbrv3ele_52x.root";
	egCorrEle_.Initialize(es, filenameEle);
      }

      eleRegrE_    .push_back(iEle->p4(reco::GsfElectron::P4_COMBINATION).t());
      eleRegrEerr_ .push_back(iEle->p4Error(reco::GsfElectron::P4_COMBINATION));

      std::pair<double,double> regrCorEle = egCorrEle_.CorrectedEnergyWithErrorV3(*eleRef, *(recVtxsBS_.product()), *rhoHandle_enReg, *lazyTool, es);
      elePhoRegrE_   .push_back(regrCorEle.first);
      elePhoRegrEerr_.push_back(regrCorEle.second);
      
      float eleSeedTime = -999.;
      int   eleRecoFlag = -999;

      // where is electron ? (0: EB, 1: EE, 2: EBGap, 3: EEGap, 4: EBEEGap)
      int elePos = -1;
      if (iEle->isEB() == true) elePos = 0;
      if (iEle->isEE() == true) elePos = 1;
      if (iEle->isEBGap() == true) elePos = 2;
      if (iEle->isEEGap() == true) elePos = 3;
      if (iEle->isEBEEGap() == true) elePos = 4;
      elePos_.push_back(elePos);

      DetId eleSeedDetId = lazyTool->getMaximum(*eleSeed).first;
      
      if ( iEle->isEB() && EBReducedRecHits_.isValid() ) {
        EcalRecHitCollection::const_iterator eleebrhit = EBReducedRecHits_->find(eleSeedDetId);
        if ( eleebrhit != EBReducedRecHits_->end() ) { 
	  eleSeedTime = eleebrhit->time(); 
	  eleRecoFlag = eleebrhit->recoFlag();
	}
      } else if ( EEReducedRecHits_.isValid() ) {
        EcalRecHitCollection::const_iterator eleeerhit = EEReducedRecHits_->find(eleSeedDetId);
        if ( eleeerhit != EEReducedRecHits_->end() ) { 
	  eleSeedTime = eleeerhit->time(); 
	  eleRecoFlag = eleeerhit->recoFlag();
	}
      }
      eleSeedTime_.push_back(eleSeedTime);
      eleRecoFlag_.push_back(eleRecoFlag);
      
      for (unsigned int i=0; i<2; ++i) eleESDetId_[nEle_][i] = 0;
      for (unsigned int i=0; i<3; ++i) 
	for (unsigned int j=0; j<62; ++j)
	  eleESHits_[nEle_][i][j] = 0;
      float eleESEffSigmaRRx = 0.;
      float eleESEffSigmaRRy = 0.;
      float eleESEffSigmaRRz = 0.;
      if (dumpESClusterInfo_) {
	for (unsigned int i=0; i<2; ++i) {
	  eleESE1_[nEle_][i] = 0.;
	  eleESE3_[nEle_][i] = 0.;
	  eleESE5_[nEle_][i] = 0.;
	  eleESE7_[nEle_][i] = 0.;
	  eleESE11_[nEle_][i] = 0.;
	  eleESE21_[nEle_][i] = 0.;
	}
      }

      if (ESRecHits_.isValid() && (fabs(eleSCEta_[nEle_]) > 1.6 && fabs(eleSCEta_[nEle_]) < 3)) {
        const GlobalPoint scPoint((*iEle).superCluster()->x(), (*iEle).superCluster()->y(), (*iEle).superCluster()->z());
        DetId extrapolatedESId1 = (dynamic_cast<const EcalPreshowerGeometry*>(geometry_p))->getClosestCellInPlane(scPoint, 1);
        DetId extrapolatedESId2 = (dynamic_cast<const EcalPreshowerGeometry*>(geometry_p))->getClosestCellInPlane(scPoint, 2);
        eleESDetId_[nEle_][0] = extrapolatedESId1.rawId();
        eleESDetId_[nEle_][1] = extrapolatedESId2.rawId();
	
        vector<float> eleESHits0 = getESHits((*iEle).superCluster()->x(), (*iEle).superCluster()->y(), (*iEle).superCluster()->z(), rechits_map_, geometry_p, topology_p, 0);
        vector<float> eleESHits1 = getESHits((*iEle).superCluster()->x(), (*iEle).superCluster()->y(), (*iEle).superCluster()->z(), rechits_map_, geometry_p, topology_p, 1);
        vector<float> eleESHits2 = getESHits((*iEle).superCluster()->x(), (*iEle).superCluster()->y(), (*iEle).superCluster()->z(), rechits_map_, geometry_p, topology_p, -1);
        for (unsigned int i=0; i<62; ++i) {
          eleESHits_[nEle_][0][i] = eleESHits0[i];
          eleESHits_[nEle_][1][i] = eleESHits1[i];
          eleESHits_[nEle_][2][i] = eleESHits2[i];
        }

        vector<float> eleESShape = getESEffSigmaRR(eleESHits0);
        eleESEffSigmaRRx = eleESShape[0];
        eleESEffSigmaRRy = eleESShape[1];
        eleESEffSigmaRRz = eleESShape[2];
	
	if (dumpESClusterInfo_) {
	  vector<float> eleESEn = getESEn(eleESHits0);
	  eleESE1_[nEle_][0]  = eleESEn[0];   eleESE1_[nEle_][1]  = eleESEn[1];
	  eleESE3_[nEle_][0]  = eleESEn[2];   eleESE3_[nEle_][1]  = eleESEn[3];
	  eleESE5_[nEle_][0]  = eleESEn[4];   eleESE5_[nEle_][1]  = eleESEn[5];
	  eleESE7_[nEle_][0]  = eleESEn[6];   eleESE7_[nEle_][1]  = eleESEn[7];
	  eleESE11_[nEle_][0] = eleESEn[8];   eleESE11_[nEle_][1] = eleESEn[9];
	  eleESE21_[nEle_][0] = eleESEn[10];  eleESE21_[nEle_][1] = eleESEn[11];
	}
      }
      eleESEffSigmaRR_x_.push_back(eleESEffSigmaRRx);
      eleESEffSigmaRR_y_.push_back(eleESEffSigmaRRy);
      eleESEffSigmaRR_z_.push_back(eleESEffSigmaRRz);
      
      //------------ Applying corrections ----------
      if (develop_) {
	float phiwidth = iEle->superCluster()->phiWidth();
	float etawidth = iEle->superCluster()->etaWidth();
	float rawE     = iEle->superCluster()->rawEnergy();
	float scEta    = iEle->superCluster()->eta();

	eleBrLinear_.push_back(phiwidth/etawidth);
	eleNBC_.push_back(iEle->superCluster()->clustersSize());
	
	float absEta = fabs(scEta);
	float energyWithEtaCorrection;
	float pTWithEtaCorrection;
	
	// C(eta) for EB only
	if (iEle->isEB()){
	  energyWithEtaCorrection = rawE/fcorrs::f5x5((int)(absEta*(5/0.087)));
	  pTWithEtaCorrection     = energyWithEtaCorrection / cosh(scEta);
	  
	  eleCetaCorrE_.push_back(energyWithEtaCorrection);
	  eleCetaCorrEt_.push_back(pTWithEtaCorrection);
	  
	  // f(brem)-corrected energies
	  float eleBremCorrE = fcorrs::fBrem(phiwidth/etawidth, energyWithEtaCorrection);
	  eleBremCorrE_.push_back(eleBremCorrE);
	  eleBremCorrEt_[nEle_] = eleBremCorrE / cosh(scEta);
	  
	  // fully-corrected energies
	  float eleFullCorrEt = fcorrs::fullCorr(eleBremCorrE/cosh(scEta), absEta);
	  eleFullCorrEt_.push_back(eleFullCorrEt);
	  eleFullCorrE_.push_back(eleFullCorrEt * cosh(scEta));
	}
	else{
	  energyWithEtaCorrection = rawE + (*iEle).superCluster()->preshowerEnergy();
	  pTWithEtaCorrection     = energyWithEtaCorrection / cosh(scEta);
	  
	  eleCetaCorrE_.push_back(energyWithEtaCorrection);
	  eleCetaCorrEt_.push_back(pTWithEtaCorrection);
	  
	  // f(brem)-corrected energies
	  float eleBremCorrE = fcorrs::fBrem_ee(phiwidth/etawidth, energyWithEtaCorrection);
	  eleBremCorrE_.push_back(eleBremCorrE);
	  eleBremCorrEt_.push_back(eleBremCorrE / cosh(scEta));
	  
	  // fully-corrected energies
	  float eleFullCorrEt = fcorrs::fullCorr_ee(eleBremCorrE/cosh(scEta), absEta);
	  eleFullCorrEt_.push_back(eleFullCorrEt);
	  eleFullCorrE_.push_back(eleFullCorrEt * cosh(scEta));
	}
      }

      vector<math::XYZPoint> vtxPoint;
      vector<float> eleEtaVtx;
      vector<float> elePhiVtx;
      vector<float> eleEtVtx;
      vector<float> eleD0Vtx;
      vector<float> eleDzVtx;
      for (size_t iv=0; iv<recVtxsBS_->size(); ++iv) { 
	vtxPoint.push_back(math::XYZPoint(vtxbs_x_[iv], vtxbs_y_[iv], vtxbs_z_[iv]));
	TVector3 *v3 = new TVector3(iEle->caloPosition().x()-(*recVtxsBS_)[iv].x(), 
				    iEle->caloPosition().y()-(*recVtxsBS_)[iv].y(), 
				    iEle->caloPosition().z()-(*recVtxsBS_)[iv].z());
	eleEtaVtx.push_back(v3->Eta());
	elePhiVtx.push_back(v3->Phi());
	eleEtVtx.push_back(iEle->energy() * TMath::Sin(2*TMath::ATan(TMath::Exp( - v3->Eta() ))));
	eleD0Vtx.push_back(iEle->gsfTrack()->dxy(vtxPoint[iv]));
	eleDzVtx.push_back(iEle->gsfTrack()->dz(vtxPoint[iv]));
	delete v3;
      }
      eleEtaVtx_.push_back(eleEtaVtx);
      elePhiVtx_.push_back(elePhiVtx);
      eleEtVtx_.push_back(eleEtVtx);
      eleD0Vtx_.push_back(eleD0Vtx);
      eleDzVtx_.push_back(eleDzVtx);
      nEle_++;	    
    }

  // Photon
  const Candidate *phomom = 0;
  nPho_ = 0;
  Int_t nGoodPho = 0;  

  const TriggerObjectMatch *phoTriggerMatch1(triggerEvent_->triggerObjectMatchResult("photonTriggerMatchHLTPhoton26Photon18"));
  const TriggerObjectMatch *phoTriggerMatch2(triggerEvent_->triggerObjectMatchResult("photonTriggerMatchHLTPhoton36Photon22"));
  const TriggerObjectMatch *phoTriggerMatch3(triggerEvent_->triggerObjectMatchResult("photonTriggerMatchHLTPhoton26R9Id85ORCaloId10Iso50Photon18R9Id85ORCaloId10Iso50Mass60"));
  const TriggerObjectMatch *phoTriggerMatch4(triggerEvent_->triggerObjectMatchResult("photonTriggerMatchHLTPhoton26R9Id85ORCaloId10Iso50Photon18R9Id85ORCaloId10Iso50Mass70"));
  const TriggerObjectMatch *phoTriggerMatch5(triggerEvent_->triggerObjectMatchResult("photonTriggerMatchHLTPhoton36R9Id85ORCaloId10Iso50Photon22R9Id85ORCaloId10Iso50"));
  const TriggerObjectMatch *phoTriggerMatch6(triggerEvent_->triggerObjectMatchResult("photonTriggerMatchHLTDoublePhoton60")); // useless now
  const TriggerObjectMatch *phoTriggerMatch7(triggerEvent_->triggerObjectMatchResult("photonTriggerMatchHLTDoublePhoton80")); // useless now

  const TriggerObjectMatch *phoTriggerFMatch1(triggerEvent_->triggerObjectMatchResult("phoTrgFMatch26IdIso18IdIso18"));
  const TriggerObjectMatch *phoTriggerFMatch2(triggerEvent_->triggerObjectMatchResult("phoTrgFMatch26IdIso18IdIso26"));
  const TriggerObjectMatch *phoTriggerFMatch3(triggerEvent_->triggerObjectMatchResult("phoTrgFMatch26R918R918"));
  const TriggerObjectMatch *phoTriggerFMatch4(triggerEvent_->triggerObjectMatchResult("phoTrgFMatch26R918R926"));
  const TriggerObjectMatch *phoTriggerFMatch5(triggerEvent_->triggerObjectMatchResult("phoTrgFMatch26R918IdIso18"));
  const TriggerObjectMatch *phoTriggerFMatch6(triggerEvent_->triggerObjectMatchResult("phoTrgFMatch26R918IdIso26"));
  const TriggerObjectMatch *phoTriggerFMatch7(triggerEvent_->triggerObjectMatchResult("phoTrgFMatch26IdIso18R918"));
  const TriggerObjectMatch *phoTriggerFMatch8(triggerEvent_->triggerObjectMatchResult("phoTrgFMatch26IdIso18R926"));
  const TriggerObjectMatch *phoTriggerFMatch9(triggerEvent_->triggerObjectMatchResult("phoTrgFMatch36IdIso22IdIso22"));
  const TriggerObjectMatch *phoTriggerFMatch10(triggerEvent_->triggerObjectMatchResult("phoTrgFMatch36IdIso22IdIso36"));
  const TriggerObjectMatch *phoTriggerFMatch11(triggerEvent_->triggerObjectMatchResult("phoTrgFMatch36R922R922"));
  const TriggerObjectMatch *phoTriggerFMatch12(triggerEvent_->triggerObjectMatchResult("phoTrgFMatch36R922R936"));
  const TriggerObjectMatch *phoTriggerFMatch13(triggerEvent_->triggerObjectMatchResult("phoTrgFMatch36IdIso22R922"));
  const TriggerObjectMatch *phoTriggerFMatch14(triggerEvent_->triggerObjectMatchResult("phoTrgFMatch36IdIso22R936"));
  const TriggerObjectMatch *phoTriggerFMatch15(triggerEvent_->triggerObjectMatchResult("phoTrgFMatch26IdIsoXL18IdIsoXL18"));
  const TriggerObjectMatch *phoTriggerFMatch16(triggerEvent_->triggerObjectMatchResult("phoTrgFMatch26IdIsoXL18IdIsoXL26"));
  const TriggerObjectMatch *phoTriggerFMatch17(triggerEvent_->triggerObjectMatchResult("phoTrgFMatch26R9T18R9T18"));
  const TriggerObjectMatch *phoTriggerFMatch18(triggerEvent_->triggerObjectMatchResult("phoTrgFMatch26R9T18R9T26"));
  const TriggerObjectMatch *phoTriggerFMatch19(triggerEvent_->triggerObjectMatchResult("phoTrgFMatch26R918IdIsoXL18"));
  const TriggerObjectMatch *phoTriggerFMatch20(triggerEvent_->triggerObjectMatchResult("phoTrgFMatch26R918IdIsoXL26"));
  const TriggerObjectMatch *phoTriggerFMatch21(triggerEvent_->triggerObjectMatchResult("phoTrgFMatch26IdIsoXL18R918"));
  const TriggerObjectMatch *phoTriggerFMatch22(triggerEvent_->triggerObjectMatchResult("phoTrgFMatch26IdIsoXL18R926")); // 7 TeV
  const TriggerObjectMatch *phoTriggerFMatch23(triggerEvent_->triggerObjectMatchResult("phoTrgFMatch26Id10Iso50HcalIso"));
  const TriggerObjectMatch *phoTriggerFMatch24(triggerEvent_->triggerObjectMatchResult("phoTrgFMatch26R9Id85"));
  const TriggerObjectMatch *phoTriggerFMatch25(triggerEvent_->triggerObjectMatchResult("phoTrgFMatch18Id10Iso50TrkIsolDouble"));
  const TriggerObjectMatch *phoTriggerFMatch26(triggerEvent_->triggerObjectMatchResult("phoTrgFMatch18R9Id85"));
  const TriggerObjectMatch *phoTriggerFMatch27(triggerEvent_->triggerObjectMatchResult("phoTrgFMatch18Id10Iso50TrkIso"));
  const TriggerObjectMatch *phoTriggerFMatch28(triggerEvent_->triggerObjectMatchResult("phoTrgFMatch18Id10Iso50TrkIsoDouble"));
  const TriggerObjectMatch *phoTriggerFMatch29(triggerEvent_->triggerObjectMatchResult("phoTrgFMatch26Id10Isol50HcalIsol"));
  const TriggerObjectMatch *phoTriggerFMatch30(triggerEvent_->triggerObjectMatchResult("phoTrgFMatch18Id10Iso50TrkIsol"));
  const TriggerObjectMatch *phoTriggerFMatch31(triggerEvent_->triggerObjectMatchResult("phoTrgFMatch30EBOnly"));

  if ( photonHandle_.isValid() )
    for (View<pat::Photon>::const_iterator iPho = photonHandle_->begin(); iPho != photonHandle_->end(); ++iPho) {
      
      edm::RefToBase<pat::Photon> phoRef = photonHandle_->refAt(nPho_);
      reco::CandidateBaseRef phoBaseRef(phoRef);
      const TriggerObjectRef phoTrigRef1( matchHelper.triggerMatchObject( phoBaseRef, phoTriggerMatch1, e, *triggerEvent_ ) );
      const TriggerObjectRef phoTrigRef2( matchHelper.triggerMatchObject( phoBaseRef, phoTriggerMatch2, e, *triggerEvent_ ) );
      const TriggerObjectRef phoTrigRef3( matchHelper.triggerMatchObject( phoBaseRef, phoTriggerMatch3, e, *triggerEvent_ ) );
      const TriggerObjectRef phoTrigRef4( matchHelper.triggerMatchObject( phoBaseRef, phoTriggerMatch4, e, *triggerEvent_ ) );
      const TriggerObjectRef phoTrigRef5( matchHelper.triggerMatchObject( phoBaseRef, phoTriggerMatch5, e, *triggerEvent_ ) );
      const TriggerObjectRef phoTrigRef6( matchHelper.triggerMatchObject( phoBaseRef, phoTriggerMatch6, e, *triggerEvent_ ) );
      const TriggerObjectRef phoTrigRef7( matchHelper.triggerMatchObject( phoBaseRef, phoTriggerMatch7, e, *triggerEvent_ ) );
      ULong_t results = 0;
      results += phoTrigRef1.isAvailable() ? powerOfTwo[0] : 0;
      results += phoTrigRef2.isAvailable() ? powerOfTwo[1] : 0;
      results += phoTrigRef3.isAvailable() ? powerOfTwo[2] : 0;
      results += phoTrigRef4.isAvailable() ? powerOfTwo[3] : 0;
      results += phoTrigRef5.isAvailable() ? powerOfTwo[4] : 0;
      results += phoTrigRef6.isAvailable() ? powerOfTwo[5] : 0;
      results += phoTrigRef7.isAvailable() ? powerOfTwo[6] : 0;
      phoTrg_.push_back(results);

      const TriggerObjectRef eleTrigFRef1( matchHelper.triggerMatchObject( phoBaseRef, eleTriggerFMatch1, e, *triggerEvent_ ) );
      const TriggerObjectRef eleTrigFRef2( matchHelper.triggerMatchObject( phoBaseRef, eleTriggerFMatch2, e, *triggerEvent_ ) );
      const TriggerObjectRef eleTrigFRef3( matchHelper.triggerMatchObject( phoBaseRef, eleTriggerFMatch3, e, *triggerEvent_ ) );
      const TriggerObjectRef eleTrigFRef4( matchHelper.triggerMatchObject( phoBaseRef, eleTriggerFMatch4, e, *triggerEvent_ ) );
      const TriggerObjectRef eleTrigFRef5( matchHelper.triggerMatchObject( phoBaseRef, eleTriggerFMatch5, e, *triggerEvent_ ) );
      const TriggerObjectRef eleTrigFRef6( matchHelper.triggerMatchObject( phoBaseRef, eleTriggerFMatch6, e, *triggerEvent_ ) );

      const TriggerObjectRef phoTrigFRef1( matchHelper.triggerMatchObject( phoBaseRef, phoTriggerFMatch1, e, *triggerEvent_ ) );
      const TriggerObjectRef phoTrigFRef2( matchHelper.triggerMatchObject( phoBaseRef, phoTriggerFMatch2, e, *triggerEvent_ ) );
      const TriggerObjectRef phoTrigFRef3( matchHelper.triggerMatchObject( phoBaseRef, phoTriggerFMatch3, e, *triggerEvent_ ) );
      const TriggerObjectRef phoTrigFRef4( matchHelper.triggerMatchObject( phoBaseRef, phoTriggerFMatch4, e, *triggerEvent_ ) );
      const TriggerObjectRef phoTrigFRef5( matchHelper.triggerMatchObject( phoBaseRef, phoTriggerFMatch5, e, *triggerEvent_ ) );
      const TriggerObjectRef phoTrigFRef6( matchHelper.triggerMatchObject( phoBaseRef, phoTriggerFMatch6, e, *triggerEvent_ ) );
      const TriggerObjectRef phoTrigFRef7( matchHelper.triggerMatchObject( phoBaseRef, phoTriggerFMatch7, e, *triggerEvent_ ) );
      const TriggerObjectRef phoTrigFRef8( matchHelper.triggerMatchObject( phoBaseRef, phoTriggerFMatch8, e, *triggerEvent_ ) );
      const TriggerObjectRef phoTrigFRef9( matchHelper.triggerMatchObject( phoBaseRef, phoTriggerFMatch9, e, *triggerEvent_ ) );
      const TriggerObjectRef phoTrigFRef10( matchHelper.triggerMatchObject( phoBaseRef, phoTriggerFMatch10, e, *triggerEvent_ ) );
      const TriggerObjectRef phoTrigFRef11( matchHelper.triggerMatchObject( phoBaseRef, phoTriggerFMatch11, e, *triggerEvent_ ) );
      const TriggerObjectRef phoTrigFRef12( matchHelper.triggerMatchObject( phoBaseRef, phoTriggerFMatch12, e, *triggerEvent_ ) );
      const TriggerObjectRef phoTrigFRef13( matchHelper.triggerMatchObject( phoBaseRef, phoTriggerFMatch13, e, *triggerEvent_ ) );
      const TriggerObjectRef phoTrigFRef14( matchHelper.triggerMatchObject( phoBaseRef, phoTriggerFMatch14, e, *triggerEvent_ ) );
      const TriggerObjectRef phoTrigFRef15( matchHelper.triggerMatchObject( phoBaseRef, phoTriggerFMatch15, e, *triggerEvent_ ) );
      const TriggerObjectRef phoTrigFRef16( matchHelper.triggerMatchObject( phoBaseRef, phoTriggerFMatch16, e, *triggerEvent_ ) );
      const TriggerObjectRef phoTrigFRef17( matchHelper.triggerMatchObject( phoBaseRef, phoTriggerFMatch17, e, *triggerEvent_ ) );
      const TriggerObjectRef phoTrigFRef18( matchHelper.triggerMatchObject( phoBaseRef, phoTriggerFMatch18, e, *triggerEvent_ ) );
      const TriggerObjectRef phoTrigFRef19( matchHelper.triggerMatchObject( phoBaseRef, phoTriggerFMatch19, e, *triggerEvent_ ) );
      const TriggerObjectRef phoTrigFRef20( matchHelper.triggerMatchObject( phoBaseRef, phoTriggerFMatch20, e, *triggerEvent_ ) );
      const TriggerObjectRef phoTrigFRef21( matchHelper.triggerMatchObject( phoBaseRef, phoTriggerFMatch21, e, *triggerEvent_ ) );
      const TriggerObjectRef phoTrigFRef22( matchHelper.triggerMatchObject( phoBaseRef, phoTriggerFMatch22, e, *triggerEvent_ ) );
      const TriggerObjectRef phoTrigFRef23( matchHelper.triggerMatchObject( phoBaseRef, phoTriggerFMatch23, e, *triggerEvent_ ) );
      const TriggerObjectRef phoTrigFRef24( matchHelper.triggerMatchObject( phoBaseRef, phoTriggerFMatch24, e, *triggerEvent_ ) );
      const TriggerObjectRef phoTrigFRef25( matchHelper.triggerMatchObject( phoBaseRef, phoTriggerFMatch25, e, *triggerEvent_ ) );
      const TriggerObjectRef phoTrigFRef26( matchHelper.triggerMatchObject( phoBaseRef, phoTriggerFMatch26, e, *triggerEvent_ ) );
      const TriggerObjectRef phoTrigFRef27( matchHelper.triggerMatchObject( phoBaseRef, phoTriggerFMatch27, e, *triggerEvent_ ) );
      const TriggerObjectRef phoTrigFRef28( matchHelper.triggerMatchObject( phoBaseRef, phoTriggerFMatch28, e, *triggerEvent_ ) );
      const TriggerObjectRef phoTrigFRef29( matchHelper.triggerMatchObject( phoBaseRef, phoTriggerFMatch29, e, *triggerEvent_ ) );
      const TriggerObjectRef phoTrigFRef30( matchHelper.triggerMatchObject( phoBaseRef, phoTriggerFMatch30, e, *triggerEvent_ ) );
      const TriggerObjectRef phoTrigFRef31( matchHelper.triggerMatchObject( phoBaseRef, phoTriggerFMatch31, e, *triggerEvent_ ) );
      
      results = 0;
      results += eleTrigFRef1.isAvailable() ? powerOfTwo[0] : 0;
      results += eleTrigFRef2.isAvailable() ? powerOfTwo[1] : 0;
      results += eleTrigFRef3.isAvailable() ? powerOfTwo[2] : 0;
      results += eleTrigFRef4.isAvailable() ? powerOfTwo[3] : 0;
      results += eleTrigFRef5.isAvailable() ? powerOfTwo[4] : 0;
      results += eleTrigFRef6.isAvailable() ? powerOfTwo[5] : 0;
      //phoTrgFilter.push_back(-99); // 7th bit  -  64 
      //phoTrgFilter.push_back(-99); // 8th bit  - 128
      //phoTrgFilter.push_back(-99); // 9th bit  - 256
      //phoTrgFilter.push_back(-99); // 10th bit - 512
      results += phoTrigFRef1.isAvailable() ? powerOfTwo[10] : 0;
      results += phoTrigFRef2.isAvailable() ? powerOfTwo[11] : 0;
      results += phoTrigFRef3.isAvailable() ? powerOfTwo[12] : 0;
      results += phoTrigFRef4.isAvailable() ? powerOfTwo[13] : 0;
      results += phoTrigFRef5.isAvailable() ? powerOfTwo[14] : 0;
      results += phoTrigFRef6.isAvailable() ? powerOfTwo[15] : 0;
      results += phoTrigFRef7.isAvailable() ? powerOfTwo[16] : 0;
      results += phoTrigFRef8.isAvailable() ? powerOfTwo[17] : 0;
      results += phoTrigFRef9.isAvailable() ? powerOfTwo[18] : 0;
      results += phoTrigFRef10.isAvailable() ? powerOfTwo[19] : 0;
      results += phoTrigFRef11.isAvailable() ? powerOfTwo[20] : 0;
      results += phoTrigFRef12.isAvailable() ? powerOfTwo[21] : 0;
      results += phoTrigFRef13.isAvailable() ? powerOfTwo[22] : 0;
      results += phoTrigFRef14.isAvailable() ? powerOfTwo[23] : 0;
      results += phoTrigFRef15.isAvailable() ? powerOfTwo[24] : 0;
      results += phoTrigFRef16.isAvailable() ? powerOfTwo[25] : 0;
      results += phoTrigFRef17.isAvailable() ? powerOfTwo[26] : 0;
      results += phoTrigFRef18.isAvailable() ? powerOfTwo[27] : 0;
      results += phoTrigFRef19.isAvailable() ? powerOfTwo[28] : 0;
      results += phoTrigFRef20.isAvailable() ? powerOfTwo[29] : 0;
      results += phoTrigFRef21.isAvailable() ? powerOfTwo[30] : 0;
      results += phoTrigFRef22.isAvailable() ? powerOfTwo[31] : 0;
      results += phoTrigFRef23.isAvailable() ? powerOfTwo[32] : 0;
      results += phoTrigFRef24.isAvailable() ? powerOfTwo[33] : 0;
      results += phoTrigFRef25.isAvailable() ? powerOfTwo[34] : 0;
      results += phoTrigFRef26.isAvailable() ? powerOfTwo[35] : 0;
      results += phoTrigFRef27.isAvailable() ? powerOfTwo[36] : 0;
      results += phoTrigFRef28.isAvailable() ? powerOfTwo[37] : 0;
      results += phoTrigFRef29.isAvailable() ? powerOfTwo[38] : 0;
      results += phoTrigFRef30.isAvailable() ? powerOfTwo[39] : 0;
      results += phoTrigFRef31.isAvailable() ? powerOfTwo[40] : 0;
      phoTrgFilter_.push_back(results);

      math::XYZPoint photonXYZ(iPho->caloPosition().x(),iPho->caloPosition().y(),iPho->caloPosition().z());
     
      edm::Ptr<reco::Candidate> recoPhoRef = iPho->originalObjectRef();
      const reco::Photon *recoPhoton = dynamic_cast<const reco::Photon *>(recoPhoRef.get());
      //isolation.mvaID(pfCandidates.product(),recoPhoton,recVtxs_);

      // PF isolation from Alternate code
      isolator.fGetIsolation(recoPhoton, &thePfColl, myVtxRef, recVtxs_);
      phoPFChIso_ .push_back(isolator.getIsolationCharged());
      phoPFNeuIso_.push_back(isolator.getIsolationNeutral());
      phoPFPhoIso_.push_back(isolator.getIsolationPhoton());

      phoIsPhoton_.push_back(iPho->isPhoton());
      phoE_       .push_back(iPho->energy());
      phoEt_      .push_back(iPho->et());
      phoEta_.push_back(iPho->eta());
      phoPhi_.push_back(iPho->phi());

      phoVtx_x_.push_back(iPho->vx());
      phoVtx_y_.push_back(iPho->vy());
      phoVtx_z_.push_back(iPho->vz());

      phoSCPos_x_.push_back((*iPho).superCluster()->x());
      phoSCPos_y_.push_back((*iPho).superCluster()->y());
      phoSCPos_z_.push_back((*iPho).superCluster()->z());
      phoCaloPos_x_.push_back(iPho->caloPosition().x());
      phoCaloPos_y_.push_back(iPho->caloPosition().y());
      phoCaloPos_z_.push_back(iPho->caloPosition().z());

      phoR9_.push_back(iPho->r9());
      phoNClus_.push_back((*iPho).superCluster()->clustersSize());
      phoTrkIsoHollowDR03_.push_back(iPho->trkSumPtHollowConeDR03());
      phoEcalIsoDR03_     .push_back(iPho->ecalRecHitSumEtConeDR03());
      phoHcalIsoDR03_     .push_back(iPho->hcalTowerSumEtConeDR03());
      phoHcalIsoDR0312_   .push_back(iPho->hcalTowerSumEtConeDR03() + 
				     (iPho->hadronicOverEm() - iPho->hadTowOverEm())*(*iPho).superCluster()->energy()/cosh((*iPho).superCluster()->eta()));
      
      phoTrkIsoHollowDR04_.push_back(iPho->trkSumPtHollowConeDR04());
      phoEcalIsoDR04_    .push_back(iPho->ecalRecHitSumEtConeDR04());
      phoHcalIsoDR04_    .push_back(iPho->hcalTowerSumEtConeDR04());
      phoHcalIsoDR0412_  .push_back(iPho->hcalTowerSumEtConeDR04() + (iPho->hadronicOverEm() - iPho->hadTowOverEm())*(*iPho).superCluster()->energy()/cosh((*iPho).superCluster()->eta()));

      phoHoverE_ .push_back(iPho->hadronicOverEm());     
      phoHoverE12_.push_back(iPho->hadTowOverEm());

      phoOverlap_.push_back(iPho->hasOverlaps("electrons"));
      phohasPixelSeed_.push_back(iPho->hasPixelSeed());
      phoEleVeto_.push_back((int)ConversionTools::hasMatchedPromptElectron(recoPhoton->superCluster(), gsfElectronHandle_, convH_, beamSpotHandle_->position()));

      size_t rightRecoPho = -1;
      for (size_t iv = 0; iv < recoPhotonHandle_->size(); ++iv) {
	reco::PhotonRef recophoRef2(recoPhotonHandle_, iv);
	if (deltaR(iPho->eta(), iPho->phi(), recophoRef2->eta(), recophoRef2->phi()) < 0.01) rightRecoPho = iv;
      }
      reco::PhotonRef recophoRef(recoPhotonHandle_, rightRecoPho);
      ggPFIsolation PFIso;
      
      std::vector<reco::PFCandidate::ParticleType> temp;
      temp.clear();
      temp.push_back(reco::PFCandidate::gamma);
      if (dumpESClusterInfo_) {
	float val_pfiso_photon005 = cicPhotonId_->pfEcalIso(recophoRef, 0.05, 0.0, 0.070, 0.015, 0.0, 0.0, 0.0, reco::PFCandidate::gamma);
	float val_pfiso_photon01 = cicPhotonId_->pfEcalIso(recophoRef, 0.1, 0.0, 0.070, 0.015, 0.0, 0.0, 0.0, reco::PFCandidate::gamma);
	float val_pfiso_photon02 = cicPhotonId_->pfEcalIso(recophoRef, 0.2, 0.0, 0.070, 0.015, 0.0, 0.0, 0.0, reco::PFCandidate::gamma);
	phoCiCPF4phopfIso005_.push_back( val_pfiso_photon005);
	phoCiCPF4phopfIso01_.push_back(  val_pfiso_photon01);
	phoCiCPF4phopfIso02_.push_back(  val_pfiso_photon02);
      }
      float val_pfiso_photon03 = cicPhotonId_->pfEcalIso(recophoRef, 0.3, 0.0, 0.070, 0.015, 0.0, 0.0, 0.0, reco::PFCandidate::gamma);
      float val_pfiso_photon04 = cicPhotonId_->pfEcalIso(recophoRef, 0.4, 0.0, 0.070, 0.015, 0.0, 0.0, 0.0, reco::PFCandidate::gamma);
      phoCiCPF4phopfIso03_.push_back( val_pfiso_photon03);
      phoCiCPF4phopfIso04_.push_back( val_pfiso_photon04);
      if (dumpESClusterInfo_) {
	float val_pfiso_photon05 = cicPhotonId_->pfEcalIso(recophoRef, 0.5, 0.0, 0.070, 0.015, 0.0, 0.0, 0.0, reco::PFCandidate::gamma);
	float val_pfiso_photon06 = cicPhotonId_->pfEcalIso(recophoRef, 0.6, 0.0, 0.070, 0.015, 0.0, 0.0, 0.0, reco::PFCandidate::gamma);
	float val_pfiso_photon07 = cicPhotonId_->pfEcalIso(recophoRef, 0.7, 0.0, 0.070, 0.015, 0.0, 0.0, 0.0, reco::PFCandidate::gamma);
	float val_pfiso_photon08 = cicPhotonId_->pfEcalIso(recophoRef, 0.8, 0.0, 0.070, 0.015, 0.0, 0.0, 0.0, reco::PFCandidate::gamma);
	phoCiCPF4phopfIso05_.push_back(val_pfiso_photon05);
	phoCiCPF4phopfIso06_.push_back(val_pfiso_photon06);
	phoCiCPF4phopfIso07_.push_back(val_pfiso_photon07);
	phoCiCPF4phopfIso08_.push_back(val_pfiso_photon08);
	//no vetoes
	float val_pfiso_photon005NoV = cicPhotonId_->pfEcalIso(recophoRef, 0.05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, reco::PFCandidate::gamma);
	float val_pfiso_photon01NoV = cicPhotonId_->pfEcalIso(recophoRef, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, reco::PFCandidate::gamma);
	float val_pfiso_photon02NoV = cicPhotonId_->pfEcalIso(recophoRef, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, reco::PFCandidate::gamma);
	float val_pfiso_photon04NoV = cicPhotonId_->pfEcalIso(recophoRef, 0.4, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, reco::PFCandidate::gamma);
	float val_pfiso_photon03NoV = cicPhotonId_->pfEcalIso(recophoRef, 0.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, reco::PFCandidate::gamma);
	float val_pfiso_photon05NoV = cicPhotonId_->pfEcalIso(recophoRef, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, reco::PFCandidate::gamma);
	float val_pfiso_photon06NoV = cicPhotonId_->pfEcalIso(recophoRef, 0.6, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, reco::PFCandidate::gamma);
	float val_pfiso_photon07NoV = cicPhotonId_->pfEcalIso(recophoRef, 0.7, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, reco::PFCandidate::gamma);
	float val_pfiso_photon08NoV = cicPhotonId_->pfEcalIso(recophoRef, 0.8, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, reco::PFCandidate::gamma);
	phoCiCPF4phopfIsoNoVETO005_.push_back(val_pfiso_photon005NoV);
	phoCiCPF4phopfIsoNoVETO01_.push_back( val_pfiso_photon01NoV);
	phoCiCPF4phopfIsoNoVETO02_.push_back( val_pfiso_photon02NoV);
	phoCiCPF4phopfIsoNoVETO03_.push_back( val_pfiso_photon03NoV);
	phoCiCPF4phopfIsoNoVETO04_.push_back( val_pfiso_photon04NoV);
	phoCiCPF4phopfIsoNoVETO05_.push_back( val_pfiso_photon05NoV);
	phoCiCPF4phopfIsoNoVETO06_.push_back( val_pfiso_photon06NoV);
	phoCiCPF4phopfIsoNoVETO07_.push_back( val_pfiso_photon07NoV);
	phoCiCPF4phopfIsoNoVETO08_.push_back( val_pfiso_photon08NoV);
      }
      temp.clear();
      temp.push_back(reco::PFCandidate::h); 
      if (dumpESClusterInfo_) {
	vector<float> vtxIsolations005 = 
	  cicPhotonId_->pfTkIsoWithVertex(recophoRef, 0.05, 0.02, 0.02, 0.0, 0.2, 0.1, reco::PFCandidate::h);
	vector<float> vtxIsolations01  = 
	  cicPhotonId_->pfTkIsoWithVertex(recophoRef, 0.1, 0.02, 0.02, 0.0, 0.2, 0.1, reco::PFCandidate::h);

	phoCiCPF4chgpfIso005_.push_back(vtxIsolations005);
	phoCiCPF4chgpfIso01_ .push_back(vtxIsolations01);
      }
      std::vector<float> vtxIsolations02 = cicPhotonId_->pfTkIsoWithVertex(recophoRef, 0.2, 0.02, 0.02, 0.0, 0.2, 0.1, reco::PFCandidate::h);
      std::vector<float> vtxIsolations03 = cicPhotonId_->pfTkIsoWithVertex(recophoRef, 0.3, 0.02, 0.02, 0.0, 0.2, 0.1, reco::PFCandidate::h);
      std::vector<float> vtxIsolations04 = cicPhotonId_->pfTkIsoWithVertex(recophoRef, 0.4, 0.02, 0.02, 0.0, 0.2, 0.1, reco::PFCandidate::h);
      phoCiCPF4chgpfIso02_.push_back(vtxIsolations02);
      phoCiCPF4chgpfIso03_.push_back(vtxIsolations03);
      phoCiCPF4chgpfIso04_.push_back(vtxIsolations04);
      if (dumpESClusterInfo_) {
	std::vector<float> vtxIsolations05    = cicPhotonId_->pfTkIsoWithVertex(recophoRef, 0.5, 0.02, 0.02, 0.0, 0.2, 0.1, reco::PFCandidate::h);
	std::vector<float> vtxIsolations06    = cicPhotonId_->pfTkIsoWithVertex(recophoRef, 0.6, 0.02, 0.02, 0.0, 0.2, 0.1, reco::PFCandidate::h);
	std::vector<float> vtxIsolations07    = cicPhotonId_->pfTkIsoWithVertex(recophoRef, 0.7, 0.02, 0.02, 0.0, 0.2, 0.1, reco::PFCandidate::h);
	std::vector<float> vtxIsolations08    = cicPhotonId_->pfTkIsoWithVertex(recophoRef, 0.8, 0.02, 0.02, 0.0, 0.2, 0.1, reco::PFCandidate::h);
	
	std::vector<float> vtxIsolations005NoV   = cicPhotonId_->pfTkIsoWithVertex(recophoRef, 0.05, 0.0, 0.0, 0.0, 0.2, 0.1, reco::PFCandidate::h);
	std::vector<float> vtxIsolations01NoV    = cicPhotonId_->pfTkIsoWithVertex(recophoRef, 0.1, 0.0, 0.0, 0.0, 0.2, 0.1, reco::PFCandidate::h);
	std::vector<float> vtxIsolations02NoV    = cicPhotonId_->pfTkIsoWithVertex(recophoRef, 0.2, 0.0, 0.0, 0.0, 0.2, 0.1, reco::PFCandidate::h);
	std::vector<float> vtxIsolations03NoV    = cicPhotonId_->pfTkIsoWithVertex(recophoRef, 0.3, 0.0, 0.0, 0.0, 0.2, 0.1, reco::PFCandidate::h);
	std::vector<float> vtxIsolations04NoV    = cicPhotonId_->pfTkIsoWithVertex(recophoRef, 0.4, 0.0, 0.0, 0.0, 0.2, 0.1, reco::PFCandidate::h);
	std::vector<float> vtxIsolations05NoV    = cicPhotonId_->pfTkIsoWithVertex(recophoRef, 0.5, 0.0, 0.0, 0.0, 0.2, 0.1, reco::PFCandidate::h);
	std::vector<float> vtxIsolations06NoV    = cicPhotonId_->pfTkIsoWithVertex(recophoRef, 0.6, 0.0, 0.0, 0.0, 0.2, 0.1, reco::PFCandidate::h);
	std::vector<float> vtxIsolations07NoV    = cicPhotonId_->pfTkIsoWithVertex(recophoRef, 0.7, 0.0, 0.0, 0.0, 0.2, 0.1, reco::PFCandidate::h);
	std::vector<float> vtxIsolations08NoV    = cicPhotonId_->pfTkIsoWithVertex(recophoRef, 0.8, 0.0, 0.0, 0.0, 0.2, 0.1, reco::PFCandidate::h);

	phoCiCPF4chgpfIso05_.push_back(vtxIsolations05);
	phoCiCPF4chgpfIso06_.push_back(vtxIsolations06);
	phoCiCPF4chgpfIso07_.push_back(vtxIsolations07);
	phoCiCPF4chgpfIso08_.push_back(vtxIsolations08);
	
	phoCiCPF4chgpfIsoNoVETO005_.push_back(vtxIsolations005NoV);
	phoCiCPF4chgpfIsoNoVETO01_.push_back(vtxIsolations01NoV);
	phoCiCPF4chgpfIsoNoVETO02_.push_back(vtxIsolations02NoV);
	phoCiCPF4chgpfIsoNoVETO03_.push_back(vtxIsolations03NoV);
	phoCiCPF4chgpfIsoNoVETO04_.push_back(vtxIsolations04NoV);
	phoCiCPF4chgpfIsoNoVETO05_.push_back(vtxIsolations05NoV);
	phoCiCPF4chgpfIsoNoVETO06_.push_back(vtxIsolations06NoV);
	phoCiCPF4chgpfIsoNoVETO07_.push_back(vtxIsolations07NoV);
	phoCiCPF4chgpfIsoNoVETO08_.push_back(vtxIsolations08NoV);
      }
      
      vector<float> phoEtaVtx;
      vector<float> phoPhiVtx;
      vector<float> phoEtVtx;
      vector<float> phoCiCTrkIsoDR03;
      vector<float> phoCiCTrkIsoDR04;

      for (size_t iv=0; iv<recVtxsBS_->size(); ++iv) {
        TVector3 *v3 = new TVector3(iPho->caloPosition().x()-(*recVtxsBS_)[iv].x(), iPho->caloPosition().y()-(*recVtxsBS_)[iv].y(), iPho->caloPosition().z()-(*recVtxsBS_)[iv].z());
        phoEtaVtx.push_back(v3->Eta());
        phoPhiVtx.push_back(v3->Phi());
        phoEtVtx.push_back(iPho->energy() * TMath::Sin(2*TMath::ATan(TMath::Exp( - v3->Eta() ))));

	//phoCiCTrkIsoDR03.push_back(getPhotonTrkIso((*recVtxsBS_)[iv].z(), (*recVtxsBS_)[iv].x(), (*recVtxsBS_)[iv].y(), (*recVtxsBS_)[iv].z(), v3->Eta(), v3->Phi(), tracksHandle_, 0.3, 0.02, 0.0, 1., 0.1, 1));
	//phoCiCTrkIsoDR04.push_back(getPhotonTrkIso((*recVtxsBS_)[iv].z(), (*recVtxsBS_)[iv].x(), (*recVtxsBS_)[iv].y(), (*recVtxsBS_)[iv].z(), v3->Eta(), v3->Phi(), tracksHandle_, 0.4, 0.02, 0.0, 1., 0.1, 1));
	
	delete v3;
      }
      phoEtaVtx_.push_back(phoEtaVtx);
      phoPhiVtx_.push_back(phoPhiVtx);
      phoEtVtx_ .push_back(phoEtVtx);
      phoCiCTrkIsoDR03_.push_back(phoCiCTrkIsoDR03);
      phoCiCTrkIsoDR04_.push_back(phoCiCTrkIsoDR04);

      phoCiCdRtoTrk_.push_back(getPhotondRtoTrk(gsfElectronHandle_, *iPho, 2.5, 0));

      // where is photon ? (0: EB, 1: EE, 2: EBGap, 3: EEGap, 4: EBEEGap)
      int phoPos = -1;
      if (iPho->isEB() == true) phoPos = 0;
      if (iPho->isEE() == true) phoPos = 1;
      if (iPho->isEBGap() == true) phoPos = 2;
      if (iPho->isEEGap() == true) phoPos = 3;
      if (iPho->isEBEEGap() == true) phoPos = 4;
      phoPos_.push_back(phoPos);

      float phoSeedTime = -999.;
      int phoRecoFlag   = -999;
      int phoSeedDetId1 = -999;
      int phoSeedDetId2 = -999;
      const reco::CaloClusterPtr phoSeed = (*iPho).superCluster()->seed();
      DetId phoSeedDetId = lazyTool->getMaximum(*phoSeed).first;
   
      vector<float> phoCov;
      phoCov = lazyTool->localCovariances(*phoSeed);
      phoSigmaIEtaIEta_.push_back(iPho->sigmaIetaIeta());
      phoSigmaIEtaIPhi_.push_back(phoCov[1]);
      phoSigmaIPhiIPhi_.push_back(phoCov[2]);
      
      if ( iPho->isEB() && EBReducedRecHits_.isValid() ) {
        EcalRecHitCollection::const_iterator ebrhit = EBReducedRecHits_->find(phoSeedDetId);
        if ( ebrhit != EBReducedRecHits_->end() ) { 
	   phoSeedTime = ebrhit->time();
           phoRecoFlag = ebrhit->recoFlag();
	   EBDetId phoEBid = EBDetId(ebrhit->id());
	   phoSeedDetId1 = phoEBid.ieta();
	   phoSeedDetId2 = phoEBid.iphi();
        }
      } else if ( EEReducedRecHits_.isValid() ) {
        EcalRecHitCollection::const_iterator eerhit = EEReducedRecHits_->find(phoSeedDetId);
        if ( eerhit != EEReducedRecHits_->end() ) { 
	   phoSeedTime = eerhit->time(); 
           phoRecoFlag = eerhit->recoFlag();
	   EEDetId phoEEid = EEDetId(eerhit->id());
	   phoSeedDetId1 = phoEEid.ix();
	   phoSeedDetId2 = phoEEid.iy();
	}
      }
      phoSeedTime_.push_back(phoSeedTime);
      phoRecoFlag_.push_back(phoRecoFlag);
      phoSeedDetId1_.push_back(phoSeedDetId1);
      phoSeedDetId2_.push_back(phoSeedDetId2);

      // computing ECAL timing for spikes
      Float_t LICTD = 0.;
      vector< pair<DetId, float> > PhotonHit_DetIds  = iPho->superCluster()->hitsAndFractions();
      vector< std::pair<DetId, float> >::const_iterator detitr;
      for (detitr = PhotonHit_DetIds.begin(); detitr != PhotonHit_DetIds.end(); ++detitr) {
	if (((*detitr).first).det() == DetId::Ecal && ((*detitr).first).subdetId() == EcalBarrel) {
	  EcalRecHitCollection::const_iterator ebhit = EBReducedRecHits_->find(((*detitr).first));
	  EcalRecHitCollection::const_iterator thishit;
	  if (ebhit != EBReducedRecHits_->end()) thishit = ebhit;
	  if (ebhit == EBReducedRecHits_->end()) continue;
	  if (thishit->energy() > 1.) {
	    if (fabs(thishit->time() - phoSeedTime) > LICTD) 
	      LICTD = fabs(thishit->time() - phoSeedTime);   
	  }
	} else if (((*detitr).first).det() == DetId::Ecal && ((*detitr).first).subdetId() == EcalEndcap) {
	  if (((*detitr).first).det() == DetId::Ecal && ((*detitr).first).subdetId() == EcalEndcap) {
	    EcalRecHitCollection::const_iterator eehit = EEReducedRecHits_->find(((*detitr).first));
	    EcalRecHitCollection::const_iterator thishit;
	    if (eehit != EEReducedRecHits_->end()) thishit = eehit;
	    if (eehit == EEReducedRecHits_->end()) continue;
	    if (thishit->energy() > 1.) {
	      if (fabs(thishit->time() - phoSeedTime) > LICTD)
		LICTD = fabs(thishit->time() - phoSeedTime);
	    }
	  }
	}
      }
      phoLICTD_.push_back(LICTD);
      phoSeedEta_.push_back(phoSeed->eta());
      phoSeedPhi_.push_back(phoSeed->phi());
      phoSeedE_.push_back(phoSeed->energy());
      phoEmax_.push_back(iPho->maxEnergyXtal());
      phoE2ndMax_.push_back(lazyTool->e2nd(*phoSeed));
      phoETop_.push_back(lazyTool->eTop(*phoSeed));
      phoEBottom_.push_back(lazyTool->eBottom(*phoSeed));
      phoELeft_.push_back(lazyTool->eLeft(*phoSeed));
      phoERight_.push_back(lazyTool->eRight(*phoSeed));
      phoE3x3_.push_back(iPho->e3x3());
      phoE5x5_.push_back(iPho->e5x5());
      phoE1x5_.push_back(iPho->e1x5());
      phoE3x1_.push_back(lazyTool->e3x1(*phoSeed));
      phoE1x3_.push_back(lazyTool->e1x3(*phoSeed));
      phoE2x2_.push_back(lazyTool->e2x2(*phoSeed));
      phoE2x5Max_.push_back(iPho->e2x5());
      phoE2x5Top_.push_back(lazyTool->e2x5Top(*phoSeed));
      phoE2x5Bottom_.push_back(lazyTool->e2x5Bottom(*phoSeed));
      phoE2x5Top_.push_back(lazyTool->e2x5Top(*phoSeed));
      phoE2x5Right_.push_back(lazyTool->e2x5Right(*phoSeed));
      phoE2x5Left_.push_back(lazyTool->e2x5Left(*phoSeed));
      phoE2x5Max_.push_back(lazyTool->e2x5Max(*phoSeed));
       if(iPho->isEB()){
      float betacry, bphicry, bthetatilt, bphitilt;
      int bieta, biphi;
         _ecalLocal.localCoordsEB(*phoSeed,es,betacry,bphicry,bieta,biphi,bthetatilt,bphitilt);
         phoCrysEta_.push_back(betacry);
         phoCrysPhi_.push_back(bphicry);
         phoCrysIEta_.push_back(bieta);
         phoCrysIPhi_.push_back(biphi);
      }
       else{
         phoCrysEta_.push_back(-99);
         phoCrysPhi_.push_back(-99);
         phoCrysIEta_.push_back(-9999);
         phoCrysIPhi_.push_back(-9999);
      }
       
      // Regression Correction
      if (!egCorrPho_.IsInitialized()) {
	//std::string filenamePho = "http://homepages.spa.umn.edu/~rekovic/cms/regweights52xV3/gbrv3ph_52x.root";
	std::string filenamePho = "gbrv3ph_52x.root";
	egCorrPho_.Initialize(es, filenamePho);
      }

      std::pair<double,double> regrCorPho = egCorrPho_.CorrectedEnergyWithErrorV3(*phoRef, *(recVtxsBS_.product()), *rhoHandle_enReg, *lazyTool, es);
      phoRegrE_.push_back(regrCorPho.first);
      phoRegrEerr_.push_back(regrCorPho.second);
      
      // Gen Particle
      float phoGenMomPID  = -999;
      float phoGenMomPt   = -999;
      int   phoGenGMomPID = -999;
      int   phoGenIndex = -1;
      int   PhoGenIndex_count = 0;
      if ( !isData_ && genParticlesHandle_.isValid() ) {
        if ((*iPho).genPhoton()) {
          for (vector<GenParticle>::const_iterator iGen = genParticlesHandle_->begin(); 
	       iGen != genParticlesHandle_->end(); ++iGen) {

            if (iGen->p4() == (*iPho).genPhoton()->p4() && 
		iGen->pdgId() == (*iPho).genPhoton()->pdgId() && 
		iGen->status() == (*iPho).genPhoton()->status()) {

              phoGenIndex = PhoGenIndex_count;

              const Candidate *phop = (const Candidate*)&(*iGen);

              for (size_t j=0; j<phop->numberOfMothers(); ++j) {
                phomom = phop->mother(j);
                phoGenMomPID = phomom->pdgId();
                phoGenMomPt  = phomom->pt();
                if (phomom->mother()) phoGenGMomPID = phomom->mother()->pdgId();
              }
            }

            PhoGenIndex_count++;
          }
        }
      }
      phoGenMomPID_.push_back(phoGenMomPID);
      phoGenMomPt_ .push_back(phoGenMomPt);
      phoGenGMomPID_.push_back(phoGenGMomPID);
      phoGenIndex_.push_back(phoGenIndex);
      
      // Super Cluster
      phoSCE_   .push_back((*iPho).superCluster()->energy());
      phoSCRawE_.push_back((*iPho).superCluster()->rawEnergy());
      phoESEn_  .push_back((*iPho).superCluster()->preshowerEnergy());
      phoSCEta_ .push_back((*iPho).superCluster()->eta());
      phoSCPhi_ .push_back((*iPho).superCluster()->phi());
      phoSCEt_  .push_back((*iPho).superCluster()->energy()/cosh((*iPho).superCluster()->eta()));
      phoSCEtaWidth_.push_back((*iPho).superCluster()->etaWidth());
      phoSCPhiWidth_.push_back((*iPho).superCluster()->phiWidth());
      phoSCBrem_    .push_back((*iPho).superCluster()->phiWidth()/(*iPho).superCluster()->etaWidth());
      
      // supercluster removal PF isolations
      phoSCRChIso_ .push_back(scRemover03.PFIsolation("charged", (*iPho).superCluster(), firstGoodVtx));// not sure about vertex index
      phoSCRPhoIso_.push_back(scRemover03.PFIsolation("photon",  (*iPho).superCluster()));
      phoSCRNeuIso_.push_back(scRemover03.PFIsolation("neutral", (*iPho).superCluster()));
      
      phoSCRChIso04_ .push_back(scRemover04.PFIsolation("charged", (*iPho).superCluster(), firstGoodVtx));// not sure about vertex index
      phoSCRPhoIso04_.push_back(scRemover04.PFIsolation("photon",  (*iPho).superCluster()));
      phoSCRNeuIso04_.push_back(scRemover04.PFIsolation("neutral", (*iPho).superCluster()));

      phoRandConeChIso_.push_back(scRemover03.RandomConeIsolation("charged", (*iPho).superCluster(), firstGoodVtx));// not sure about vertex index
      phoRandConePhoIso_.push_back(scRemover03.RandomConeIsolation("photon",  (*iPho).superCluster()));
      phoRandConeNeuIso_.push_back(scRemover03.RandomConeIsolation("neutral", (*iPho).superCluster()));
      
      phoRandConeChIso04_.push_back(scRemover04.RandomConeIsolation("charged", (*iPho).superCluster(), firstGoodVtx));// not sure about vertex index
      phoRandConePhoIso04_.push_back(scRemover04.RandomConeIsolation("photon",  (*iPho).superCluster()));
      phoRandConeNeuIso04_.push_back(scRemover04.RandomConeIsolation("neutral", (*iPho).superCluster()));
      
      // skim on good photons
      float phoHoverE = iPho->hadronicOverEm();
      float phoSCEta = (*iPho).superCluster()->eta();
      float phoSigmaIEtaIEta = iPho->sigmaIetaIeta();
      if (phoHoverE < 0.15 && fabs(phoSCEta) < 2.5) {
	if (fabs(phoSCEta) < 1.4442) {
	  if (phoSigmaIEtaIEta < 0.017) nGoodPho++;
	} else {
	  if (phoSigmaIEtaIEta < 0.04)  nGoodPho++;
	}
      }

      // Conversion
      phoIsConv_.push_back(iPho->hasConversionTracks());
      phoNConv_.push_back(iPho->conversions().size());
 
      float phoConvInvMass       = 0;
      float phoConvCotTheta      = 0;
      float phoConvEoverP        = 0;
      float phoConvZofPVfromTrks = 0;
      float phoConvMinDist       = 0;
      float phoConvdPhiAtVtx     = 0;
      float phoConvdPhiAtCalo    = 0;
      float phoConvdEtaAtCalo    = 0;
      float phoConvTrkd0_x       = 0;
      float phoConvTrkd0_y       = 0;
      float phoConvTrkPin_x      = 0;
      float phoConvTrkPin_y      = 0;
      float phoConvTrkPout_x     = 0;
      float phoConvTrkPout_y     = 0;
      float phoConvTrkdz_x       = 0;
      float phoConvTrkdz_y       = 0;
      float phoConvTrkdzErr_x    = 0;
      float phoConvTrkdzErr_y    = 0;
      float phoConvChi2          = 0;
      float phoConvChi2Prob      = 0;
      float phoConvNTrks         = 0;
      float phoConvLikeLihood    = 0;
      float phoConvP4_0 = 0.;
      float phoConvP4_1 = 0.;
      float phoConvP4_2 = 0.;
      float phoConvP4_3 = 0.;
      float phoConvVtx_x = 0.;
      float phoConvVtx_y = 0.;
      float phoConvVtx_z = 0.;
      float phoConvVtxErr_x = 0.;
      float phoConvVtxErr_y = 0.;
      float phoConvVtxErr_z = 0.;
      float phoConvPairMomentum_x = 0.;
      float phoConvPairMomentum_y = 0.;
      float phoConvPairMomentum_z = 0.;
      float phoConvRefittedMomentum_x = 0.;
      float phoConvRefittedMomentum_y = 0.;
      float phoConvRefittedMomentum_z = 0.;
      bool  phoConvValidVtx = 0;

      float phoConvCharge1 = 0;
      float phoConvCharge2 = 0;

      if (iPho->hasConversionTracks()) {
 
        reco::ConversionRefVector conversions = iPho->conversions();
        reco::ConversionRef conv = conversions[0];
 
        phoConvValidVtx = conv->conversionVertex().isValid();
	
        for (unsigned int j=0; j<conversions.size(); ++j) {
          conv = conversions[j];
          reco::Vertex vtx=conv->conversionVertex();

          phoConvP4_0 = conv->refittedPair4Momentum().px();
          phoConvP4_1 = conv->refittedPair4Momentum().py();
          phoConvP4_2 = conv->refittedPair4Momentum().pz();
          phoConvP4_3 = conv->refittedPair4Momentum().energy();

          phoConvVtx_x = vtx.x();
          phoConvVtx_y = vtx.y();
          phoConvVtx_z = vtx.z();
          phoConvVtxErr_x = vtx.xError();
          phoConvVtxErr_y = vtx.yError();
          phoConvVtxErr_z = vtx.zError();
          phoConvPairMomentum_x = conv->pairMomentum().x();
          phoConvPairMomentum_y = conv->pairMomentum().y();
          phoConvPairMomentum_z = conv->pairMomentum().z();
          phoConvRefittedMomentum_x = conv->refittedPairMomentum().x();
          phoConvRefittedMomentum_y = conv->refittedPairMomentum().y();
          phoConvRefittedMomentum_z = conv->refittedPairMomentum().z();
 
          phoConvChi2       = vtx.chi2();
          phoConvChi2Prob   = ChiSquaredProbability(vtx.chi2(), vtx.ndof());
          phoConvNTrks      = conv->nTracks();
          phoConvLikeLihood = conv->MVAout();

          const std::vector<edm::RefToBase<reco::Track> > tracks = conv->tracks();
          for (unsigned int k=0; k<tracks.size(); ++k) {
            if (k==0) {
              phoConvTrkdz_x    = tracks[k]->dz();
              phoConvTrkdzErr_x = tracks[k]->dzError();
	      phoConvCharge1    = tracks[k]->charge();
	    } 
	    if (k==1) {
              phoConvTrkdz_y    = tracks[k]->dz();
              phoConvTrkdzErr_y = tracks[k]->dzError();
	      phoConvCharge2    = tracks[k]->charge();
	    }
          }
	}
 
        phoConvInvMass       = conv->pairInvariantMass();
        phoConvCotTheta      = conv->pairCotThetaSeparation();
        phoConvEoverP        = conv->EoverPrefittedTracks();
        phoConvZofPVfromTrks = conv->zOfPrimaryVertexFromTracks();
        phoConvMinDist       = conv->distOfMinimumApproach();
        phoConvdPhiAtVtx     = conv->dPhiTracksAtVtx();
	phoConvdPhiAtCalo    = conv->dPhiTracksAtEcal();
	phoConvdEtaAtCalo    = conv->dEtaTracksAtEcal();

        if (conv->tracks().size() > 0) {
          phoConvTrkd0_x   = conv->tracksSigned_d0()[0];
          phoConvTrkPin_x  = sqrt(conv->tracksPin()[0].Mag2());
          phoConvTrkPout_x = sqrt(conv->tracksPout()[0].Mag2());
        }
 
        if (conv->tracks().size() > 1) {
          phoConvTrkd0_y   = conv->tracksSigned_d0()[1];
          phoConvTrkPin_y  = sqrt(conv->tracksPin()[1].Mag2());
          phoConvTrkPout_y = sqrt(conv->tracksPout()[1].Mag2());
        }
 
      }
      // fill the vectors
      phoConvInvMass_      .push_back(phoConvInvMass);
      phoConvCotTheta_     .push_back(phoConvCotTheta);
      phoConvEoverP_       .push_back(phoConvEoverP);
      phoConvZofPVfromTrks_.push_back(phoConvZofPVfromTrks);
      phoConvMinDist_      .push_back(phoConvMinDist);
      phoConvdPhiAtVtx_    .push_back(phoConvdPhiAtVtx);
      phoConvdPhiAtCalo_   .push_back(phoConvdPhiAtCalo);
      phoConvdEtaAtCalo_   .push_back(phoConvdEtaAtCalo);
      phoConvTrkd0_x_      .push_back(phoConvTrkd0_x);
      phoConvTrkd0_y_      .push_back(phoConvTrkd0_y);
      phoConvTrkPin_x_     .push_back(phoConvTrkPin_x);
      phoConvTrkPin_y_     .push_back(phoConvTrkPin_y);
      phoConvTrkPout_x_    .push_back(phoConvTrkPout_x);
      phoConvTrkPout_y_    .push_back(phoConvTrkPout_y);
      phoConvTrkdz_x_      .push_back(phoConvTrkdz_x);
      phoConvTrkdz_y_      .push_back(phoConvTrkdz_y);
      phoConvTrkdzErr_x_   .push_back(phoConvTrkdzErr_x);
      phoConvTrkdzErr_y_   .push_back(phoConvTrkdzErr_y);
      phoConvChi2_         .push_back(phoConvChi2);
      phoConvChi2Prob_     .push_back(phoConvChi2Prob);
      phoConvNTrks_        .push_back(phoConvNTrks);
      phoConvCharge1_      .push_back(phoConvCharge1);
      phoConvCharge2_      .push_back(phoConvCharge2);
      phoConvValidVtx_     .push_back(phoConvValidVtx);
      phoConvLikeLihood_   .push_back(phoConvLikeLihood);
      phoConvP4_0_         .push_back(phoConvP4_0);
      phoConvP4_1_         .push_back(phoConvP4_1);
      phoConvP4_2_         .push_back(phoConvP4_2);
      phoConvP4_3_         .push_back(phoConvP4_3);
      phoConvVtx_x_        .push_back(phoConvVtx_x);
      phoConvVtx_y_        .push_back(phoConvVtx_y);
      phoConvVtx_z_        .push_back(phoConvVtx_z);
      phoConvVtxErr_x_     .push_back(phoConvVtxErr_x);
      phoConvVtxErr_y_     .push_back(phoConvVtxErr_y);
      phoConvVtxErr_z_     .push_back(phoConvVtxErr_z);
      phoConvPairMomentum_x_.push_back(phoConvPairMomentum_x);
      phoConvPairMomentum_y_.push_back(phoConvPairMomentum_y);
      phoConvPairMomentum_z_.push_back(phoConvPairMomentum_z);
      phoConvRefittedMomentum_x_.push_back(phoConvRefittedMomentum_x);
      phoConvRefittedMomentum_y_.push_back(phoConvRefittedMomentum_y);
      phoConvRefittedMomentum_z_.push_back(phoConvRefittedMomentum_z);
      //phoConvInvMass_      .push_back(phoConvInvMass);// EXTRA
      //phoConvCotTheta_     .push_back(phoConvCotTheta);
      //phoConvEoverP_       .push_back(phoConvEoverP);
      //phoConvZofPVfromTrks_.push_back(phoConvZofPVfromTrks);
      //phoConvMinDist_      .push_back(phoConvMinDist);
      //phoConvdPhiAtVtx_    .push_back(phoConvdPhiAtVtx);
      //phoConvdPhiAtCalo_   .push_back(phoConvdPhiAtCalo);
      //phoConvdEtaAtCalo_   .push_back(phoConvdEtaAtCalo);

      //do the same for Single Legs:
      int SingleLegConv = 0;
      vector<float> phoPFConvVtx_x;
      vector<float> phoPFConvVtx_y;
      vector<float> phoPFConvVtx_z;
      vector<float> phoPFConvMom_x;
      vector<float> phoPFConvMom_y;
      vector<float> phoPFConvMom_z;
      reco::PhotonCollection::const_iterator pfphot=pfPhoTranslator_->begin();
      for(;pfphot!=pfPhoTranslator_->end();++pfphot){
	
	if(recophoRef->superCluster()!=pfphot->superCluster())continue;
	//	cout<<"MATCHED PHOTON "<<endl;
	reco::ConversionRefVector SLconversions=pfphot->conversionsOneLeg();
	SingleLegConv = SLconversions.size();
	for(unsigned int SL=0; SL<SLconversions.size(); ++SL){
	  const std::vector<edm::RefToBase<reco::Track> > tracks = SLconversions[SL]->tracks();
	  phoPFConvVtx_x.push_back(SLconversions[SL]->conversionVertex().x());
	  phoPFConvVtx_y.push_back(SLconversions[SL]->conversionVertex().y());
	  phoPFConvVtx_z.push_back(SLconversions[SL]->conversionVertex().z());
	  phoPFConvMom_x.push_back(tracks[0]->px());
	  phoPFConvMom_y.push_back(tracks[0]->py());
	  phoPFConvMom_z.push_back(tracks[0]->pz());
	}
      }
      SingleLegConv_ .push_back(SingleLegConv);
      phoPFConvVtx_x_.push_back(phoPFConvVtx_x);
      phoPFConvVtx_y_.push_back(phoPFConvVtx_y);
      phoPFConvVtx_z_.push_back(phoPFConvVtx_z);
      phoPFConvMom_x_.push_back(phoPFConvMom_x);
      phoPFConvMom_y_.push_back(phoPFConvMom_y);
      phoPFConvMom_z_.push_back(phoPFConvMom_z);

      //init PhotonInfo and do Vertex Selection by Rishi:
      
      //reco::PhotonRef recophoRef(recoPhotonHandle_,nPho_);
      ggPFPhotons ggPFPhoton(*recophoRef, pfPhoTranslator_,
			     gsfElectronHandle_,
			     pfCandidates,
			     EBReducedRecHits_,
			     EEReducedRecHits_,
			     ESRecHits_,
			     geomBar,
			     geomEnd,
			     beamSpotHandle_
			     );	
      
      float MustacheEin       = -9999;
      float MustacheEOut      = -9999;
      float MustacheEtOut     = -9999;
      float PFLowestClustE    = -9999;
      float PFClustdEta       = -9999;
      float PFClustdPhi       = -9999;
      float PFClustRMSPhi     = -9999;
      float PFClustRMSPhiMust = -9999;
      float pho_pfconvVtxZ    = -9999;
      float pho_pfconvVtxZErr = -9999;
      float PFClustEneCorr    = -9999;
      int pho_hasConvPf   = 0;
      int pho_hasSLConvPf = 0;
      int PFEleMatch      = 0;
      int PFEleVeto       = 0;
      int pho_nSLConv     = 0;
      int PFRecoMatch     = 0;

      vector<float> pho_pfSLConvVtxZ;
      vector<float> pho_pfSLConvPos_x;
      vector<float> pho_pfSLConvPos_y;
      vector<float> pho_pfSLConvPos_z;

      if(ggPFPhoton.MatchPFReco()){
	PFRecoMatch=1;

	if(ggPFPhoton.isConv()){
	  pho_hasConvPf = 1;
	  
	  std::pair<float, float>VertexZ=ggPFPhoton.SLPoint();
	  pho_pfconvVtxZ = VertexZ.first;
	  pho_pfconvVtxZErr = VertexZ.second;
	}
	
	if(ggPFPhoton.hasSLConv()){
	  pho_hasSLConvPf = 1;
	  reco::ConversionRefVector SLconversions;
	  std::vector<float> Zint;

	  if(SLPoint(*recophoRef, SLconversions, Zint)){
	    for(unsigned int z=0; z<SLconversions.size(); ++z){
	      pho_nSLConv=z+1;
	      pho_pfSLConvVtxZ.push_back(Zint[z]);
	      pho_pfSLConvPos_x.push_back(SLconversions[z]->conversionVertex().x());
	      pho_pfSLConvPos_y.push_back(SLconversions[z]->conversionVertex().y());
	      pho_pfSLConvPos_z.push_back(SLconversions[z]->conversionVertex().z());
	    }
	  }
	}
	
	if(ggPFPhoton.isPFEle())PFEleMatch=1;
	
	if(ggPFPhoton.PFElectronVeto(convH_, gsfElectronHandle_))PFEleVeto = 1;
	MustacheEin    = ggPFPhoton.MustE();
	MustacheEOut   = ggPFPhoton.MustEOut();
	MustacheEtOut  = ggPFPhoton.MustEtOut();
	PFLowestClustE = ggPFPhoton.PFLowE();
	PFClustdEta    = ggPFPhoton.PFdEta();
	PFClustdPhi    = ggPFPhoton.PFdPhi();
	PFClustRMSPhi  = ggPFPhoton.PFClusRMSTot();
	PFClustRMSPhiMust = ggPFPhoton.PFClusRMSMust();	
	//std::vector<reco::CaloCluster>PFC=ggPFPhoton.PFClusters();
	//PFClustEneCorr_[nPho_]=ggPFPhoton.getPFPhoECorr(PFC, ReaderLCEB_, ReaderLCEE_);
      }
      else{	
	PFRecoMatch = 0;	
	ggPFPhoton.recoPhotonClusterLink(*iPho, pfCandidates);  
	MustacheEin = ggPFPhoton.MustE();
	MustacheEOut = ggPFPhoton.MustEOut();
	MustacheEtOut = ggPFPhoton.MustEtOut();
	PFLowestClustE = ggPFPhoton.PFLowE();
	PFClustdEta = ggPFPhoton.PFdEta();
	PFClustdPhi = ggPFPhoton.PFdPhi();
	PFClustRMSPhi = ggPFPhoton.PFClusRMSTot();
	PFClustRMSPhiMust = ggPFPhoton.PFClusRMSMust();	
      }
      // fill the variables
      MustacheEin_      .push_back(MustacheEin);
      MustacheEOut_     .push_back(MustacheEOut);
      MustacheEtOut_    .push_back(MustacheEtOut);
      PFLowestClustE_   .push_back(PFLowestClustE);
      PFClustdEta_      .push_back(PFClustdEta);
      PFClustdPhi_      .push_back(PFClustdPhi);
      PFClustRMSPhi_    .push_back(PFClustRMSPhi);
      PFClustRMSPhiMust_.push_back(PFClustRMSPhiMust);
      pho_pfconvVtxZ_   .push_back(pho_pfconvVtxZ);
      pho_pfconvVtxZErr_.push_back(pho_pfconvVtxZErr);
      PFClustEneCorr_   .push_back(PFClustEneCorr);
      pho_hasConvPf_    .push_back(pho_hasConvPf);
      pho_hasSLConvPf_  .push_back(pho_hasSLConvPf);
      PFEleMatch_       .push_back(PFEleMatch);
      PFEleVeto_        .push_back(PFEleVeto);
      pho_nSLConv_      .push_back(pho_nSLConv);
      PFRecoMatch_      .push_back(PFRecoMatch);

      pho_pfSLConvVtxZ_ .push_back(pho_pfSLConvVtxZ);
      pho_pfSLConvPos_x_.push_back(pho_pfSLConvPos_x);
      pho_pfSLConvPos_y_.push_back(pho_pfSLConvPos_y);
      pho_pfSLConvPos_z_.push_back(pho_pfSLConvPos_z);

      //------------ Applying corrections ----------      
      if (develop_) {
	float absEta = fabs((*iPho).superCluster()->eta());
	float phoSCRawE = (*iPho).superCluster()->rawEnergy();
	float phoSCBrem = (*iPho).superCluster()->phiWidth() / (*iPho).superCluster()->etaWidth();

	// C(eta) for EB only
	if (iPho->isEB()){
	  float energyWithEtaCorrection = phoSCRawE/fcorrs::f5x5((int)(absEta*(5/0.087)));
	  float pTWithEtaCorrection     = energyWithEtaCorrection / cosh(phoSCEta);
	  
	  phoCetaCorrE_.push_back(energyWithEtaCorrection);
	  phoCetaCorrEt_.push_back(pTWithEtaCorrection);
	  
	  // f(brem)-corrected energies
	  float phoBremCorrE = fcorrs::fBrem(phoSCBrem, energyWithEtaCorrection);
	  phoBremCorrE_.push_back(phoBremCorrE);
	  phoBremCorrEt_.push_back(phoBremCorrE / cosh(phoSCEta));
	  
	  // fully-corrected energies
	  float phoFullCorrEt = fcorrs::fullCorr(phoBremCorrE/cosh(phoSCEta), absEta);
	  phoFullCorrEt_.push_back(phoFullCorrEt);
	  phoFullCorrE_.push_back(phoFullCorrEt * cosh(phoSCEta)); 
	}
	else{
	  float energyWithEtaCorrection = phoSCRawE + (*iPho).superCluster()->preshowerEnergy();
	  float pTWithEtaCorrection     = energyWithEtaCorrection / cosh(phoSCEta);
	  
	  phoCetaCorrE_.push_back(energyWithEtaCorrection);
	  phoCetaCorrEt_.push_back(pTWithEtaCorrection);
	  
	  // f(brem)-corrected energies
	  float phoBremCorrE = fcorrs::fBrem_ee(phoSCBrem, energyWithEtaCorrection);
	  phoBremCorrE_.push_back(phoBremCorrE);
	  phoBremCorrEt_.push_back(phoBremCorrE / cosh(phoSCEta));
	  
	  // fully-corrected energies
	  float phoFullCorrEt = fcorrs::fullCorr_ee(phoBremCorrE/cosh(phoSCEta), absEta);
	  phoFullCorrEt_.push_back(phoFullCorrEt);
	  phoFullCorrE_.push_back(phoFullCorrEt * cosh(phoSCEta));
	}
      }
      //--------------------------------------------
      
      for (unsigned int i=0; i<2; ++i) phoESDetId_[nPho_][i] = 0;
      for (unsigned int i=0; i<3; ++i) 
	for (unsigned int j=0; j<62; ++j) 
	  phoESHits_[nPho_][i][j] = 0;
      if (dumpESClusterInfo_) {
	for (unsigned int i=0; i<2; ++i) {
	  phoESE1_[nPho_][i] = 0.;
	  phoESE3_[nPho_][i] = 0.;
	  phoESE5_[nPho_][i] = 0.;
	  phoESE7_[nPho_][i] = 0.;
	  phoESE11_[nPho_][i] = 0.;
	  phoESE21_[nPho_][i] = 0.;
	}
      }
    
      float phoESEffSigmaRR_x = 0.;
      float phoESEffSigmaRR_y = 0.;
      float phoESEffSigmaRR_z = 0.;
      
      if (ESRecHits_.isValid() && (fabs(phoSCEta_[nPho_]) > 1.6 && fabs(phoSCEta_[nPho_]) < 3)) {
        const GlobalPoint scPoint((*iPho).superCluster()->x(), (*iPho).superCluster()->y(), (*iPho).superCluster()->z());
        DetId extrapolatedESId1 = (dynamic_cast<const EcalPreshowerGeometry*>(geometry_p))->getClosestCellInPlane(scPoint, 1);
        DetId extrapolatedESId2 = (dynamic_cast<const EcalPreshowerGeometry*>(geometry_p))->getClosestCellInPlane(scPoint, 2);
        phoESDetId_[nPho_][0] = extrapolatedESId1.rawId();
        phoESDetId_[nPho_][1] = extrapolatedESId2.rawId();

        vector<float> phoESHits0 = getESHits((*iPho).superCluster()->x(), (*iPho).superCluster()->y(), (*iPho).superCluster()->z(), rechits_map_, geometry_p, topology_p, 0);
        vector<float> phoESHits1 = getESHits((*iPho).superCluster()->x(), (*iPho).superCluster()->y(), (*iPho).superCluster()->z(), rechits_map_, geometry_p, topology_p, 1);
        vector<float> phoESHits2 = getESHits((*iPho).superCluster()->x(), (*iPho).superCluster()->y(), (*iPho).superCluster()->z(), rechits_map_, geometry_p, topology_p, -1);
        for (unsigned int i=0; i<62; ++i) {
          phoESHits_[nPho_][0][i] = phoESHits0[i];
          phoESHits_[nPho_][1][i] = phoESHits1[i];
          phoESHits_[nPho_][2][i] = phoESHits2[i];
        }

        vector<float> phoESShape = getESEffSigmaRR(phoESHits0);
        phoESEffSigmaRR_x = phoESShape[0];
        phoESEffSigmaRR_y = phoESShape[1];
        phoESEffSigmaRR_z = phoESShape[2];
	if (dumpESClusterInfo_) {
	  vector<float> phoESEn = getESEn(phoESHits0);
	  phoESE1_[nPho_][0] =  phoESEn[0];   phoESE1_[nPho_][1] =  phoESEn[1];
	  phoESE3_[nPho_][0] =  phoESEn[2];   phoESE3_[nPho_][1] =  phoESEn[3];
	  phoESE5_[nPho_][0] =  phoESEn[4];   phoESE5_[nPho_][1] =  phoESEn[5];
	  phoESE7_[nPho_][0] =  phoESEn[6];   phoESE7_[nPho_][1] =  phoESEn[7];
	  phoESE11_[nPho_][0] = phoESEn[8];   phoESE11_[nPho_][1] = phoESEn[9];
	  phoESE21_[nPho_][0] = phoESEn[10];  phoESE21_[nPho_][1] = phoESEn[11];
	}
      }
      phoESEffSigmaRR_x_.push_back(phoESEffSigmaRR_x);
      phoESEffSigmaRR_y_.push_back(phoESEffSigmaRR_y);
      phoESEffSigmaRR_z_.push_back(phoESEffSigmaRR_z);

      nPho_++;
    }
  
  // muon trigger matching
  const TriggerObjectMatch * muTriggerMatch1( triggerEvent_->triggerObjectMatchResult( "muonTriggerMatchHLTIsoMu24eta2p1" ) );
  const TriggerObjectMatch * muTriggerMatch2( triggerEvent_->triggerObjectMatchResult( "muonTriggerMatchHLTIsoMu24" ) );
  const TriggerObjectMatch * muTriggerMatch3( triggerEvent_->triggerObjectMatchResult( "muonTriggerMatchHLTMu17Mu8" ) );
  const TriggerObjectMatch * muTriggerMatch4( triggerEvent_->triggerObjectMatchResult( "muonTriggerMatchHLTMu17TkMu8" ) );
  const TriggerObjectMatch * muTriggerMatch5( triggerEvent_->triggerObjectMatchResult( "muonTriggerMatchHLTMu22TkMu8" ) );
  const TriggerObjectMatch * muTriggerMatch6( triggerEvent_->triggerObjectMatchResult( "muonTriggerMatchHLTMu22Mu8" ) );
  const TriggerObjectMatch * muTriggerMatch7( triggerEvent_->triggerObjectMatchResult( "muonTriggerMatchHLTMu17forMu17Mu8" ) );
  const TriggerObjectMatch * muTriggerMatch8( triggerEvent_->triggerObjectMatchResult( "muonTriggerMatchHLTMu8forMu17Mu8" ) );
  const TriggerObjectMatch * muTriggerMatch9( triggerEvent_->triggerObjectMatchResult( "muonTriggerMatchHLTMu17forMu17TkMu8" ) );
  const TriggerObjectMatch * muTriggerMatch10( triggerEvent_->triggerObjectMatchResult( "muonTriggerMatchHLTMu8forMu17TkMu8" ) );
  
  // Muon
  nMu_ = 0;
  if ( muonHandle_.isValid() ) {
    for (View<pat::Muon>::const_iterator iMu = muonHandle_->begin(); iMu != muonHandle_->end(); ++iMu) {

      if (iMu->pt() < 5) continue;
      if (! (iMu->isPFMuon() || iMu->isGlobalMuon() || iMu->isTrackerMuon())) continue; 
      
      edm::RefToBase<pat::Muon> muRef = muonHandle_->refAt(nMu_);
      reco::CandidateBaseRef muBaseRef(muRef);
      const TriggerObjectRef muTrigRef1( matchHelper.triggerMatchObject( muBaseRef, muTriggerMatch1, e, *triggerEvent_ ) );
      const TriggerObjectRef muTrigRef2( matchHelper.triggerMatchObject( muBaseRef, muTriggerMatch2, e, *triggerEvent_ ) );
      const TriggerObjectRef muTrigRef3( matchHelper.triggerMatchObject( muBaseRef, muTriggerMatch3, e, *triggerEvent_ ) );
      const TriggerObjectRef muTrigRef4( matchHelper.triggerMatchObject( muBaseRef, muTriggerMatch4, e, *triggerEvent_ ) );
      const TriggerObjectRef muTrigRef5( matchHelper.triggerMatchObject( muBaseRef, muTriggerMatch5, e, *triggerEvent_ ) );
      const TriggerObjectRef muTrigRef6( matchHelper.triggerMatchObject( muBaseRef, muTriggerMatch6, e, *triggerEvent_ ) );
      const TriggerObjectRef muTrigRef7( matchHelper.triggerMatchObject( muBaseRef, muTriggerMatch7, e, *triggerEvent_ ) );
      const TriggerObjectRef muTrigRef8( matchHelper.triggerMatchObject( muBaseRef, muTriggerMatch8, e, *triggerEvent_ ) );
      const TriggerObjectRef muTrigRef9( matchHelper.triggerMatchObject( muBaseRef, muTriggerMatch9, e, *triggerEvent_ ) );
      const TriggerObjectRef muTrigRef10( matchHelper.triggerMatchObject( muBaseRef, muTriggerMatch10, e, *triggerEvent_ ) );
      ULong_t results = 0;
      results += muTrigRef1 .isAvailable() ? powerOfTwo[0] : 0;
      results += muTrigRef2 .isAvailable() ? powerOfTwo[1] : 0;
      results += muTrigRef3 .isAvailable() ? powerOfTwo[2] : 0;      
      results += muTrigRef4 .isAvailable() ? powerOfTwo[3] : 0;
      results += muTrigRef5 .isAvailable() ? powerOfTwo[4] : 0;
      results += muTrigRef6 .isAvailable() ? powerOfTwo[5] : 0;
      results += muTrigRef7 .isAvailable() ? powerOfTwo[6] : 0;
      results += muTrigRef8 .isAvailable() ? powerOfTwo[7] : 0;
      results += muTrigRef9 .isAvailable() ? powerOfTwo[8] : 0;
      results += muTrigRef10.isAvailable() ? powerOfTwo[9] : 0;
      muTrg_.push_back(results);

      Int_t goodMuonTrack = 1;
      const reco::TrackRef trkr   = iMu->globalTrack();
      const reco::TrackRef inntrk = iMu->innerTrack();
      vector<float> muD0Vtx;
      vector<float> muDzVtx;
      if (trkr.isNull()) {
	goodMuonTrack = 0;
        muD0_     .push_back(-99.);
	muDz_     .push_back(-99.);
        muD0GV_   .push_back(-99.);
	muDzGV_   .push_back(-99.);
	muChi2NDF_.push_back(-99.);
	muNumberOfValidTrkHits_ .push_back(-99);
	muNumberOfValidMuonHits_.push_back(-99);
	muVtxGlb_x_.push_back(-99);
	muVtxGlb_y_.push_back(-99);
	muVtxGlb_z_.push_back(-99);
      } else {
        muD0_     .push_back(trkr->dxy(pv));
        muDz_     .push_back(trkr->dz(pv));
        muD0GV_   .push_back(trkr->dxy(gv));
	muDzGV_   .push_back(trkr->dz(gv));
  	muChi2NDF_.push_back(trkr->normalizedChi2());
        muNumberOfValidTrkHits_.push_back(trkr->hitPattern().numberOfValidTrackerHits());
        muNumberOfValidMuonHits_.push_back(trkr->hitPattern().numberOfValidMuonHits());
	muVtxGlb_x_.push_back(trkr->vx());
	muVtxGlb_y_.push_back(trkr->vy());
	muVtxGlb_z_.push_back(trkr->vz());
      }

      if (inntrk.isNull()) {
	goodMuonTrack = 0;
        muInnerD0_   .push_back(-99.);
	muInnerDz_   .push_back(-99.);
        muInnerD0GV_ .push_back(-99.);
	muInnerDzGV_ .push_back(-99.);
	muInnerChi2NDF_ .push_back(-99.);
	muInnerPt_      .push_back(-99.);
	muInnerPtErr_   .push_back(-99.);
	muNumberOfValidTrkLayers_   .push_back(-99);
	muNumberOfValidPixelLayers_ .push_back(-99);
	muNumberOfValidPixelHits_   .push_back(-99);
      } else {
        muInnerD0_   .push_back(inntrk->dxy(pv));
        muInnerDz_   .push_back(inntrk->dz(pv));
        muInnerD0GV_ .push_back(inntrk->dxy(gv));
        muInnerDzGV_ .push_back(inntrk->dz(gv));
  	muInnerChi2NDF_ .push_back(inntrk->normalizedChi2());
	muInnerPt_      .push_back(inntrk->pt());
	muInnerPtErr_   .push_back(inntrk->ptError());
	muNumberOfValidTrkLayers_   .push_back(inntrk->hitPattern().trackerLayersWithMeasurement());
	muNumberOfValidPixelLayers_ .push_back(inntrk->hitPattern().pixelLayersWithMeasurement());
        muNumberOfValidPixelHits_   .push_back(inntrk->hitPattern().numberOfValidPixelHits());
	math::XYZPoint vtxPoint[100];
	for (size_t iv=0; iv<recVtxsBS_->size(); ++iv) {   
	  vtxPoint[iv] = math::XYZPoint(vtxbs_x_[iv], vtxbs_y_[iv], vtxbs_z_[iv]);	 
	  muD0Vtx.push_back(inntrk->dxy(vtxPoint[iv]));
	  muDzVtx.push_back(inntrk->dz(vtxPoint[iv]));	  
	}
      }
      muD0Vtx_.push_back(muD0Vtx);
      muDzVtx_.push_back(muDzVtx);

      muStations_.push_back(iMu->numberOfMatchedStations());
      muChambers_.push_back(iMu->numberOfMatches());

      // IP3D
      float muIP3D    = -999.;
      float muIP3DErr = -99.;
      float muVtx_x  = -99; 
      float muVtx_y  = -99; 
      float muVtx_z  = -99; 
      if (! iMu->track().isNull()) {
	const double musign = ((- iMu->track()->dxy(pv)) >= 0) ? 1. : -1.;
	const reco::TransientTrack tt_mu = thebuilder.build(iMu->track());
	const std::pair<bool,Measurement1D> ip3dpv_mu = IPTools::absoluteImpactParameter3D(tt_mu, *recVtxs_->begin());
	if (ip3dpv_mu.first) {
	  muIP3D    = musign*ip3dpv_mu.second.value();
	  muIP3DErr = ip3dpv_mu.second.error();
	}
	muVtx_x = iMu->track()->vx();
	muVtx_y = iMu->track()->vy();
	muVtx_z = iMu->track()->vz();
      }
      muIP3D_.push_back(muIP3D);
      muIP3DErr_.push_back(muIP3DErr);
      muVtx_x_.push_back(muVtx_x);
      muVtx_y_.push_back(muVtx_y);
      muVtx_z_.push_back(muVtx_z);

      muType_.push_back(iMu->type());
      muEta_ .push_back(iMu->eta());
      muPhi_ .push_back(iMu->phi());
      muCharge_.push_back(iMu->charge());
      muPt_  .push_back(iMu->pt());
      muPz_  .push_back(iMu->pz());

      // muon cocktail
      edm::Ptr<reco::Candidate> recoMuRef = iMu->originalObjectRef();
      const reco::Muon *recoMu = dynamic_cast<const reco::Muon *>(recoMuRef.get());

      if (iMu->tpfmsTrack().isNull()) goodMuonTrack = 0;
      if (iMu->pickyTrack().isNull()) goodMuonTrack = 0;
      reco::TrackRef cktTrack;
      if (goodMuonTrack == 1) cktTrack = (muon::tevOptimized(*recoMu, 200, 17., 40., 0.25)).first;
      mucktPt_   .push_back((! cktTrack.isNull()) ? cktTrack->pt() : -99);
      mucktPtErr_.push_back((! cktTrack.isNull()) ? cktTrack->ptError() : -99);
      mucktEta_  .push_back((! cktTrack.isNull()) ? cktTrack->eta() : -99);
      mucktPhi_  .push_back((! cktTrack.isNull()) ? cktTrack->phi() : -99);
      mucktdxy_  .push_back((! cktTrack.isNull()) ? cktTrack->dxy(pv) : -99);
      mucktdz_   .push_back((! cktTrack.isNull()) ? cktTrack->dz(pv) : -99);

      muIsoTrk_ .push_back(iMu->trackIso());
      muIsoCalo_.push_back(iMu->caloIso());
      muIsoEcal_.push_back(iMu->ecalIso());
      muIsoHcal_.push_back(iMu->hcalIso());

      muPFIsoR03_CH_   .push_back(iMu->pfIsolationR03().sumChargedHadronPt);
      muPFIsoR03_NH_   .push_back(iMu->pfIsolationR03().sumNeutralHadronEt);
      muPFIsoR03_Pho_  .push_back(iMu->pfIsolationR03().sumPhotonEt);
      
      muPFIsoR03_PU_   .push_back(iMu->pfIsolationR03().sumPUPt);
      muPFIsoR03_CPart_.push_back(iMu->pfIsolationR03().sumChargedParticlePt);
      muPFIsoR03_NHHT_ .push_back(iMu->pfIsolationR03().sumNeutralHadronEtHighThreshold);
      muPFIsoR03_PhoHT_.push_back(iMu->pfIsolationR03().sumPhotonEtHighThreshold); 

      muPFIsoR04_CH_   .push_back(iMu->pfIsolationR04().sumChargedHadronPt);
      muPFIsoR04_NH_   .push_back(iMu->pfIsolationR04().sumNeutralHadronEt);
      muPFIsoR04_Pho_  .push_back(iMu->pfIsolationR04().sumPhotonEt);
      
      muPFIsoR04_PU_   .push_back(iMu->pfIsolationR04().sumPUPt);
      muPFIsoR04_CPart_.push_back(iMu->pfIsolationR04().sumChargedParticlePt);
      muPFIsoR04_NHHT_ .push_back(iMu->pfIsolationR04().sumNeutralHadronEtHighThreshold);
      muPFIsoR04_PhoHT_.push_back(iMu->pfIsolationR04().sumPhotonEtHighThreshold); 

      int muGenIndex = -1;
      int MuGenIndex = 0;
      if (!isData_) {
        if ( (*iMu).genLepton() && genParticlesHandle_.isValid() ) {
          if (fabs((*iMu).genLepton()->pdgId())==13) {
            for (vector<GenParticle>::const_iterator iGen = genParticlesHandle_->begin(); 
		 iGen != genParticlesHandle_->end(); ++iGen) {

              if (iGen->p4() == (*iMu).genLepton()->p4() && 
		  iGen->pdgId() == (*iMu).genLepton()->pdgId() && 
		  iGen->status() == (*iMu).genLepton()->status()) 
	        muGenIndex = MuGenIndex;

              MuGenIndex++;
            }
          }
        }
      }
      muGenIndex_.push_back(muGenIndex);
      nMu_++;
    }  // iMu loop
  } // muonHandle is valid
 

  // Taus
  nTau_ = 0;
  if ( tauHandle_.isValid() ) {
    for(vector<pat::Tau>::const_iterator itau = tauHandle_->begin(); itau != tauHandle_->end(); ++itau) {
      // ID related variables
      tauDecayModeFinding_.push_back(itau->tauID("decayModeFinding"));
      tauAgainstElectronLooseMVA3_.push_back(itau->tauID("againstElectronLooseMVA3"));
      tauAgainstElectronMediumMVA3_.push_back(itau->tauID("againstElectronMediumMVA3"));
      tauAgainstElectronTightMVA3_.push_back(itau->tauID("againstElectronTightMVA3"));
      tauAgainstElectronVTightMVA3_.push_back(itau->tauID("againstElectronVTightMVA3"));
      tauAgainstElectronDeadECAL_.push_back(itau->tauID("againstElectronDeadECAL"));
      tauAgainstMuonLoose2_.push_back(itau->tauID("againstMuonLoose2"));
      tauAgainstMuonMedium2_.push_back(itau->tauID("againstMuonMedium2"));
      tauAgainstMuonTight2_.push_back(itau->tauID("againstMuonTight2"));
      tauCombinedIsolationDeltaBetaCorrRaw3Hits_.push_back(itau->tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits"));
      tauLooseCombinedIsolationDeltaBetaCorr3Hits_.push_back(itau->tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits"));
      tauMediumCombinedIsolationDeltaBetaCorr3Hits_.push_back(itau->tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits"));
      tauTightCombinedIsolationDeltaBetaCorr3Hits_.push_back(itau->tauID("byTightCombinedIsolationDeltaBetaCorr3Hits"));

      // kinematics + isolation
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
      
      reco::PFCandidateRef leadPFChargedHadrCand_Ref = itau->leadPFChargedHadrCand();
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

    // PFPhotons
  nPFPho_  = 0;
  nPFEle_  = 0;
  nPFchad_ = 0;
  nPFnhad_ = 0;    
  
  float PFPhoE = 0;
  float PFPhoEt = 0;
  float PFPhoEta = 0;
  float PFPhoPhi = 0;
  float PFPhoIso = 0;
  int PFPhoType = 999;

  if (develop_) { // loop all PF candidates

    for (unsigned iCand=0; iCand< pfAllCandidates->size(); ++iCand) {
      
      const reco::PFCandidate  & pfParticleT((*pfAllCandidates)[iCand]);
      const reco::PFCandidate  *pfParticle = (&pfParticleT);
      
      if (pfParticle->pdgId() == 22 && pfParticle->pt() > 2 && fabs(pfParticle->eta()) < 2.5) {
	
	PFPhoE    = pfParticle->energy();
	PFPhoEt   = pfParticle->pt();
	PFPhoEta  = pfParticle->eta();
	PFPhoPhi  = pfParticle->phi();
	PFPhoIso  = fsrPhotonIso03(pfAllCandidates.product(), pfCandidates.product(), *pfParticle);
	PFPhoType = 1;

	// non-develop arrays
	PFPhoE_.push_back(PFPhoE);
	PFPhoEt_.push_back(PFPhoEt);
	PFPhoEta_.push_back(PFPhoEta);
	PFPhoPhi_.push_back(PFPhoPhi);
	PFPhoType_.push_back(PFPhoType);
	PFPhoIso_.push_back(PFPhoIso);

	nPFPho_++;
      }
   
      if (abs(pfParticle->pdgId()) == 13 && pfParticle->muonRef().isNonnull() && pfParticle->ecalEnergy() > 2) {
	if ((pfParticle->ecalEnergy() * pfParticle->pt() / pfParticle-> energy()) > 2) {
	  
	  PFPhoE    = pfParticle->ecalEnergy();
	  PFPhoEt   = pfParticle->ecalEnergy() * pfParticle->pt() / pfParticle->energy();
	  PFPhoEta  = pfParticle->eta();
	  PFPhoPhi  = pfParticle->phi();
	  PFPhoIso  = fsrPhotonIso03(pfAllCandidates.product(), pfCandidates.product(), *pfParticle);
	  PFPhoType = 2;
	  
	  // non-develop arrays
	  PFPhoE_.push_back(PFPhoE);
	  PFPhoEt_.push_back(PFPhoEt);
	  PFPhoEta_.push_back(PFPhoEta);
	  PFPhoPhi_.push_back(PFPhoPhi);
	  PFPhoType_.push_back(PFPhoType);
	  PFPhoIso_.push_back(PFPhoIso);

	  nPFPho_++;
	}
      }
    }
  }
  
  if (develop_ ) {  
    //PFElectron-PFPhoton disambiguation:
    std::vector<reco::PFCandidate>matchPFEle;
    //std::vector<reco::PFCandidate>ConvPFEle;
    std::vector<reco::PFCandidate>PFEle;  
    int nPFE=0;
    unsigned ncandidates= pfCandidatePhotons->size();
    for(unsigned iCand=0;iCand<ncandidates;++iCand) {
      const reco::PFCandidate  & pfParticleT((*pfCandidatePhotons)[iCand]);
      const reco::PFCandidate  *pfParticle = (&pfParticleT);
      
      if(abs(pfParticle->pdgId())==11){
	//match to GSFElectron:
	nPFE++;    
	
	reco::GsfElectronCollection::const_iterator pfele=gsfElectronHandle_->begin();
	bool EleVeto=false;
	bool Conv=false;
	for(;pfele!=gsfElectronHandle_->end();++pfele){
	  if(pfele->superCluster().isNull())continue;
	  if(pfele->gsfTrack()==pfParticle->gsfTrackRef()){
	    if(pfele->gsfTrack()==pfParticle->gsfTrackRef()){
	      
	      //Conversion Veto Based on Matched GsfElectron:
	      EleVeto= ConversionTools::hasMatchedPromptElectron(pfele->superCluster(), gsfElectronHandle_, convH_, beamSpotHandle_->position());
	      if(!EleVeto)Conv=true;
	      break;
	    }
	  }
	}
	if(Conv)matchPFEle.push_back(*pfParticle); 
	if(!Conv)PFEle.push_back(*pfParticle);
      }
    }
    
    ncandidates = pfCandidatePhotons->size();
    for (unsigned iCand=0; iCand<ncandidates; ++iCand) {
      
      const reco::PFCandidate  & pfParticleT((*pfCandidatePhotons)[iCand]);
      const reco::PFCandidate  *pfParticle = (&pfParticleT);
      
      if(pfParticle->pdgId()==22 && pfParticle->mva_nothing_gamma()>0.1) {
	reco::SuperClusterRef PhotscRef= pfParticle->superClusterRef();
	PFPhoE    = pfParticle->energy();
	PFPhoEt   = pfParticle->pt();
	PFPhoEta  = pfParticle->eta();
	PFPhoPhi  = pfParticle->phi();
	PFPhoType = 1;
      
	PFCandidate::ElementsInBlocks eleInBlocks = pfParticle->elementsInBlocks();
	nPFPho_tks=0;
	nPFPhoTks_[nPFPho_]=nPFPho_tks;
	nPFPho_clust=0;
	nPFPho_ES1clust=0;
	nPFPho_ES2clust=0;
	fillBlockInfo(*pfParticle, true, false,geomBar, geomEnd, geometry_p); 

	// non-develop arrays
	PFPhoE_.push_back(PFPhoE);
	PFPhoEt_.push_back(PFPhoEt);
	PFPhoEta_.push_back(PFPhoEta);
	PFPhoPhi_.push_back(PFPhoPhi);
	PFPhoType_.push_back(PFPhoType);
	PFPhoIso_.push_back(PFPhoIso);
	
	nPFPho_++;
      }
    }
  
    for (unsigned int pfe=0; pfe<matchPFEle.size(); ++pfe) {
      //for AOD match to GSF Electron:
      if ( electronHandle_.isValid() )
	for (View<pat::Electron>::const_iterator iEle = electronHandle_->begin(); iEle != electronHandle_->end(); ++iEle){
	  float deta=fabs(iEle->eta()- matchPFEle[pfe].eta());
	  float dphi=acos(cos(iEle->phi()- matchPFEle[pfe].phi()));
	  if(deta<0.1 && dphi<0.4){
	    PFPhoGsfeta_In_[nPFPho_]=iEle->trackMomentumAtVtx().eta();
	    PFPhoGsfeta_Out_[nPFPho_]=iEle->trackMomentumOut().eta();	      
	    PFPhoGsfphi_In_[nPFPho_]=iEle->trackMomentumAtVtx().phi();
	    PFPhoGsfphi_Out_[nPFPho_]=iEle->trackMomentumOut().phi();;
	    PFPhoGSFPin_[nPFPho_][0]=iEle->trackMomentumAtVtx().X();
	    PFPhoGSFPin_[nPFPho_][1]=iEle->trackMomentumAtVtx().Y();;
	    PFPhoGSFPin_[nPFPho_][2]=iEle->trackMomentumAtVtx().Z();;
	    PFPhoGsf_In_[nPFPho_][0]=iEle->trackPositionAtVtx().x();
	    PFPhoGsf_In_[nPFPho_][1]=iEle->trackPositionAtVtx().y();;
	    PFPhoGsf_In_[nPFPho_][2]=iEle->trackPositionAtVtx().z();;
	    PFGsf_Theta_[nPFPho_]=iEle->trackMomentumAtVtx().Theta();;
	    PFPhoGsf_ThetaErr_[nPFPho_]=iEle->gsfTrack()->thetaError();
	    PFPhoGSFPout_[nPFPho_]=iEle->trackMomentumOut().R();
	    break;
	  }
	}
    
      PFPhoE=matchPFEle[pfe].energy();
      PFPhoEt=matchPFEle[pfe].pt();
      PFPhoEta=matchPFEle[pfe].eta();
      PFPhoPhi=matchPFEle[pfe].phi();
      
      PFCandidate::ElementsInBlocks eleInBlocks = matchPFEle[pfe].elementsInBlocks();
      nPFPho_tks=0;
      nPFPhoTks_[nPFPho_]=nPFPho_tks;
      nPFPho_clust=0;
      nPFPho_ES1clust=0;
      nPFPho_ES2clust=0;
      //cout << endl << "Called by photons (fill)";
      if (develop_) fillBlockInfo(matchPFEle[pfe], true,true, geomBar, geomEnd, geometry_p);

      // non-develop arrays
      PFPhoE_.push_back(PFPhoE);
      PFPhoEt_.push_back(PFPhoEt);
      PFPhoEta_.push_back(PFPhoEta);
      PFPhoPhi_.push_back(PFPhoPhi);
      PFPhoType_.push_back(PFPhoType);
      PFPhoIso_.push_back(PFPhoIso);
      
      nPFPho_++;
    }
  
    matchPFEle.clear();
    for(unsigned int pfe=0; pfe<PFEle.size(); ++pfe){    
      PFElePt_[nPFEle_]=PFEle[pfe].pt();
      PFEleEta_[nPFEle_]=PFEle[pfe].eta();
      PFElePhi_[nPFEle_]=PFEle[pfe].phi();
      PFEleEn_[nPFEle_]=PFEle[pfe].energy();
      //for AOD match to GSF Electron:
      if ( electronHandle_.isValid() )
	for (View<pat::Electron>::const_iterator iEle = electronHandle_->begin(); iEle != electronHandle_->end(); ++iEle) {
	  if(iEle->gsfTrack()==PFEle[pfe].gsfTrackRef()){
	    PFGsfeta_In_[nPFEle_]=iEle->trackMomentumAtVtx().eta();
	    PFGsfeta_Out_[nPFEle_]=iEle->trackMomentumOut().eta();	      
	    PFGsfphi_In_[nPFEle_]=iEle->trackMomentumAtVtx().phi();
	    PFGsfphi_Out_[nPFEle_]=iEle->trackMomentumOut().phi();;
	    PFElePin_[nPFEle_][0]=iEle->trackMomentumAtVtx().X();
	    PFElePin_[nPFEle_][1]=iEle->trackMomentumAtVtx().Y();;
	    PFElePin_[nPFEle_][2]=iEle->trackMomentumAtVtx().Z();;
	    PFGsf_In_[nPFEle_][0]=iEle->trackPositionAtVtx().x();
	    PFGsf_In_[nPFEle_][1]=iEle->trackPositionAtVtx().y();;
	    PFGsf_In_[nPFEle_][2]=iEle->trackPositionAtVtx().z();;
	    PFGsf_Theta_[nPFEle_]=iEle->trackMomentumAtVtx().Theta();;
	    PFGsf_ThetaErr_[nPFEle_]=iEle->gsfTrack()->thetaError();
	    PFElePout_[nPFEle_]=iEle->trackMomentumOut().R();
	    PFEleChi2NDF_[nPFEle_]=iEle->gsfTrack()->normalizedChi2();
	    break;
	  }
	}
      
      PFCandidate::ElementsInBlocks eleInBlocks = PFEle[pfe].elementsInBlocks();
      
      nPFEle_clust=0;
      //cout << endl << "Called by electrons (fill)";
      if (develop_) fillBlockInfo(PFEle[pfe], false, false, geomBar, geomEnd, geometry_p);
      nPFEle_++;
    }
    PFEle.clear();
      
    for (PFCandidateCollection::const_iterator pfParticle =pfCandidates->begin(); pfParticle!=pfCandidates->end(); ++pfParticle) {
      if (useAllPF_) {	  
	if (abs(pfParticle->pdgId()) == 111) {
	  nPFchad_++; 
	  PFhad_charge_[nPFchad_]=pfParticle->charge();
	  PFchad_Eta_[nPFchad_]=pfParticle->eta();
	  PFchad_Phi_[nPFchad_]=pfParticle->phi();
	  PFchad_Pt_[nPFchad_]=pfParticle->pt();
	  PFchad_E_[nPFchad_]=pfParticle->energy();
	  PFchad_P_[nPFchad_]=pfParticle->p();
	  PFchad_Ecal_[nPFchad_]=pfParticle->ecalEnergy();
	  PFchad_Hcal_[nPFchad_]=pfParticle->hcalEnergy();
	}
	if (abs(pfParticle->pdgId()) == 211) {
	  nPFnhad_++; 
	  PFnhad_Eta_[nPFchad_]=pfParticle->eta();
	  PFnhad_Phi_[nPFchad_]=pfParticle->phi();
	  PFnhad_Pt_[nPFchad_]=pfParticle->pt();
	  PFnhad_E_[nPFchad_]=pfParticle->energy();
	  PFnhad_P_[nPFchad_]=pfParticle->p();
	  PFnhad_Hcal_[nPFnhad_]=pfParticle->hcalEnergy();
	}
      }
    }
  }

  // Jet                                                                            
  const TriggerObjectMatch *jetTriggerMatch1(triggerEvent_->triggerObjectMatchResult("jetTriggerMatchHLTJet30"));
  const TriggerObjectMatch *jetTriggerMatch2(triggerEvent_->triggerObjectMatchResult("jetTriggerMatchHLTJet60"));
  const TriggerObjectMatch *jetTriggerMatch3(triggerEvent_->triggerObjectMatchResult("jetTriggerMatchHLTJet80"));
  const TriggerObjectMatch *jetTriggerMatch4(triggerEvent_->triggerObjectMatchResult("jetTriggerMatchHLTJet110"));
  const TriggerObjectMatch *jetTriggerMatch5(triggerEvent_->triggerObjectMatchResult("jetTriggerMatchHLTJet150"));
  const TriggerObjectMatch *jetTriggerMatch6(triggerEvent_->triggerObjectMatchResult("jetTriggerMatchHLTJet190"));
  const TriggerObjectMatch *jetTriggerMatch7(triggerEvent_->triggerObjectMatchResult("jetTriggerMatchHLTJet240"));
  const TriggerObjectMatch *jetTriggerMatch8(triggerEvent_->triggerObjectMatchResult("jetTriggerMatchHLTJet300"));
  const TriggerObjectMatch *jetTriggerMatch9(triggerEvent_->triggerObjectMatchResult("jetTriggerMatchHLTJet370"));

  edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
  es.get<JetCorrectionsRecord>().get("AK5PFchs",JetCorParColl);
  JetCorrectionUncertainty *jecUnc=0;
  JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
  jecUnc = new JetCorrectionUncertainty(JetCorPar);

  PileupJetIdAlgo* jetMVACalculator = 0;
  if (pujetIDalgos_.size() > 0) jetMVACalculator = pujetIDalgos_[0];
  nLowPtJet_ = 0;
  nJet_ = 0;
  if (jetHandle_.isValid()) 
    for (View<pat::Jet>::const_iterator iJet = jetHandle_->begin(); iJet != jetHandle_->end(); ++iJet) {
     
      if ((*iJet).correctedJet("Uncorrected").pt() >= 5. && (*iJet).correctedJet("Uncorrected").pt() < 12. && nLowPtJet_<600 && dumpTrks_) {
        jetLowPtEn_       .push_back(iJet->energy());
        jetLowPtPt_       .push_back(iJet->pt());
        jetLowPtEta_      .push_back(iJet->eta());
        jetLowPtPhi_      .push_back(iJet->phi());
        jetLowPtCharge_   .push_back(iJet->jetCharge());
        jetLowPtEt_       .push_back(iJet->et());
        jetLowPtRawPt_    .push_back((*iJet).correctedJet("Uncorrected").pt());
        jetLowPtRawEn_    .push_back((*iJet).correctedJet("Uncorrected").energy());
        jetLowPtArea_     .push_back(iJet->jetArea());
        jetLowPtPartonID_ .push_back(iJet->partonFlavour());
	
        int jetLowPtGenPartonID = -99;
	float jetLowPtGenEn  = 0;
	float jetLowPtGenPt  = 0;
	float jetLowPtGenEta = 0;
	float jetLowPtGenPhi = 0;
        if (!isData_ && genParticlesHandle_.isValid() ) {
          if ((*iJet).genParton()) {
            jetLowPtGenPartonID = (*iJet).genParton()->pdgId();
            jetLowPtGenEn       = (*iJet).genParton()->energy();
            jetLowPtGenPt       = (*iJet).genParton()->pt();
            jetLowPtGenEta      = (*iJet).genParton()->eta();
            jetLowPtGenPhi      = (*iJet).genParton()->phi();
          }
        }
	jetLowPtGenPartonID_.push_back(jetLowPtGenPartonID);
	jetLowPtGenEn_.push_back(jetLowPtGenEn);
	jetLowPtGenPt_.push_back(jetLowPtGenPt);
	jetLowPtGenEta_.push_back(jetLowPtGenEta);
	jetLowPtGenPhi_.push_back(jetLowPtGenPhi);

        float jetLowPtGenJetEn   = -1; 
        float jetLowPtGenJetPt   = -999;
        float jetLowPtGenJetEta  = -999;
        float jetLowPtGenJetPhi  = -999;
        if (!isData_ && genParticlesHandle_.isValid() ) {
          if ((*iJet).genJet()) {
            jetLowPtGenJetEn  = (*iJet).genJet()->energy();
            jetLowPtGenJetPt  = (*iJet).genJet()->pt();
            jetLowPtGenJetEta = (*iJet).genJet()->eta();
            jetLowPtGenJetPhi = (*iJet).genJet()->phi();
          }
        }
	jetLowPtGenJetEn_ .push_back(jetLowPtGenJetEn);
	jetLowPtGenJetPt_ .push_back(jetLowPtGenJetPt);
	jetLowPtGenJetEta_.push_back(jetLowPtGenJetEta);
	jetLowPtGenJetPhi_.push_back(jetLowPtGenJetPhi);
        nLowPtJet_++;
      }

      if ((*iJet).correctedJet("Uncorrected").pt() < 12. || nJet_>=600) continue;

      edm::RefToBase<pat::Jet> jetRef = jetHandle_->refAt(nJet_);
      reco::CandidateBaseRef jetBaseRef(jetRef);
      const TriggerObjectRef jetTrigRef1( matchHelper.triggerMatchObject( jetBaseRef, jetTriggerMatch1, e, *triggerEvent_ ) );
      const TriggerObjectRef jetTrigRef2( matchHelper.triggerMatchObject( jetBaseRef, jetTriggerMatch2, e, *triggerEvent_ ) );
      const TriggerObjectRef jetTrigRef3( matchHelper.triggerMatchObject( jetBaseRef, jetTriggerMatch3, e, *triggerEvent_ ) );
      const TriggerObjectRef jetTrigRef4( matchHelper.triggerMatchObject( jetBaseRef, jetTriggerMatch4, e, *triggerEvent_ ) );
      const TriggerObjectRef jetTrigRef5( matchHelper.triggerMatchObject( jetBaseRef, jetTriggerMatch5, e, *triggerEvent_ ) );
      const TriggerObjectRef jetTrigRef6( matchHelper.triggerMatchObject( jetBaseRef, jetTriggerMatch6, e, *triggerEvent_ ) );
      const TriggerObjectRef jetTrigRef7( matchHelper.triggerMatchObject( jetBaseRef, jetTriggerMatch7, e, *triggerEvent_ ) );
      const TriggerObjectRef jetTrigRef8( matchHelper.triggerMatchObject( jetBaseRef, jetTriggerMatch8, e, *triggerEvent_ ) );
      const TriggerObjectRef jetTrigRef9( matchHelper.triggerMatchObject( jetBaseRef, jetTriggerMatch9, e, *triggerEvent_ ) );
      ULong_t results = 0;
      results += jetTrigRef1.isAvailable() ? powerOfTwo[0] : 0;
      results += jetTrigRef2.isAvailable() ? powerOfTwo[1] : 0;
      results += jetTrigRef3.isAvailable() ? powerOfTwo[2] : 0;
      results += jetTrigRef4.isAvailable() ? powerOfTwo[3] : 0;
      results += jetTrigRef5.isAvailable() ? powerOfTwo[4] : 0;
      results += jetTrigRef6.isAvailable() ? powerOfTwo[5] : 0;
      results += jetTrigRef7.isAvailable() ? powerOfTwo[6] : 0;
      results += jetTrigRef8.isAvailable() ? powerOfTwo[7] : 0;
      results += jetTrigRef9.isAvailable() ? powerOfTwo[8] : 0;
      jetTrg_.push_back(results);

      jetEn_.push_back(    iJet->energy());
      jetPt_.push_back(    iJet->pt());
      jetEta_.push_back(   iJet->eta());
      jetPhi_.push_back(   iJet->phi());
      jetCharge_.push_back(iJet->jetCharge());
      jetEt_.push_back(    iJet->et());
      jetMt_.push_back(    iJet->mt());
      jetRawPt_.push_back( (*iJet).correctedJet("Uncorrected").pt());
      jetRawEn_.push_back( (*iJet).correctedJet("Uncorrected").energy());

      jetArea_.push_back(  iJet->jetArea());
      jetCEF_.push_back(   iJet->chargedEmEnergyFraction());
      jetNEF_.push_back(   iJet->neutralEmEnergyFraction());
      jetCHF_.push_back(   iJet->chargedHadronEnergyFraction());
      jetNHF_.push_back(   iJet->neutralHadronEnergyFraction());
      jetHFHAE_.push_back( iJet->HFHadronEnergy());
      jetHFEME_.push_back( iJet->HFEMEnergy());
      jetNCH_.push_back(   iJet->chargedMultiplicity());
      jetNConstituents_.push_back(iJet->getPFConstituents().size());

      if (fabs(iJet->eta()) < 5.2) {
	jecUnc->setJetEta(iJet->eta());
	jecUnc->setJetPt(iJet->pt()); // here you must use the CORRECTED jet pt
	jetJECUnc_.push_back(jecUnc->getUncertainty(true));
      } else {
	jetJECUnc_.push_back(-1.); 
      }

      float PI = 3.1415927;
      float dphiTemp=fabs(pfMETPhi_-iJet->phi());
      if(dphiTemp>PI) dphiTemp = 2*PI-dphiTemp;
      jetDPhiMETJet_.push_back(dphiTemp);

      float jetLeadTrackPt = -9999.;
      if (iJet->isPFJet() == true) {
	std::vector <reco::PFCandidatePtr> constituents = iJet->getPFConstituents ();
	for (unsigned ic = 0; ic < constituents.size (); ++ic) {
	  if ( constituents[ic]->particleId() > 3 ) continue;
	  reco::TrackRef trackRef = constituents[ic]->trackRef();
	  if ( trackRef.isNonnull() ) { if(trackRef->pt() > jetLeadTrackPt) jetLeadTrackPt = trackRef->pt(); }
	}
      }
      jetLeadTrackPt_.push_back(jetLeadTrackPt);

      float jetVtxPt   = -99;
      float jetVtxMass = -99;
      float jetVtx3dL  = -99;
      float jetVtx3deL = -99;
      const reco::SecondaryVertexTagInfo * tf = iJet->tagInfoSecondaryVertex("secondaryVertex");
      if (tf){
      	math::XYZTLorentzVectorD vertexSum;
      	for(size_t vi=0;vi< tf->nVertices(); ++vi)
      	  {
      	    vertexSum+=tf->secondaryVertex(vi).p4();
      	  }
      	jetVtxPt = vertexSum.Pt();
      	if (tf->nVertices() >0){
      	  jetVtxMass =  tf->secondaryVertex(0).p4().mass();
      	  Measurement1D m = tf->flightDistance(0);
      	  jetVtx3dL = m.value();
      	  jetVtx3deL = m.error();
      	}
      }
      jetVtxPt_  .push_back(jetVtxPt);
      jetVtxMass_.push_back(jetVtxMass);
      jetVtx3dL_ .push_back(jetVtx3dL);
      jetVtx3deL_.push_back(jetVtx3deL);

      int isSemiLept=0;
      float jetSoftLeptdR        = -99.;
      float jetSoftLeptPt        = -99.;
      float jetSoftLeptPtRel     = -99.;
      float jetSoftLeptIdlooseMu = -99.;
      float jetSoftLeptIdEle95   = -99.;
      edm::Handle<edm::View<reco::Candidate> > muonNoCutsHandle;
      e.getByLabel(muonNoCutsLabel_, muonNoCutsHandle);
      edm::View<reco::Candidate> muonsNoCuts = *muonNoCutsHandle;
      for(edm::View<reco::Candidate>::const_iterator mu = muonsNoCuts.begin(); mu!=muonsNoCuts.end() && isSemiLept!=1; ++mu){
	const pat::Muon& m = static_cast <const pat::Muon&> (*mu);
	float Smpt = m.pt();
	float Smeta = m.eta();
	float Smphi = m.phi();
	float SmJdR = deltaR(Smeta, Smphi, iJet->eta(), iJet->phi());
	if   ( Smpt > 5 && SmJdR < 0.5) { //lep_ptCutForBjets_ = 5 GeV
	  isSemiLept=1;
	  jetSoftLeptdR = SmJdR;
	  jetSoftLeptPt =Smpt;
	  TLorentzVector jvec (iJet->p4().X(), iJet->p4().Y(), iJet->p4().Z(), iJet->p4().T());
	  TVector3 mvec ( m.p4().Vect().X(), m.p4().Vect().Y(), m.p4().Vect().Z()  );
	  jetSoftLeptPtRel = jvec.Perp(  mvec );
	  jetSoftLeptIdlooseMu = m.muonID("TMLastStationLoose");
	}
      }
      jetSoftLeptIdlooseMu_.push_back(jetSoftLeptIdlooseMu);

      edm::Handle<edm::View<reco::Candidate> > eleNoCutsHandle;
      e.getByLabel(eleNoCutsLabel_, eleNoCutsHandle);
      edm::View<reco::Candidate> elesNoCuts = *eleNoCutsHandle;
      for(edm::View<reco::Candidate>::const_iterator ele = elesNoCuts.begin(); 
	  ele!=elesNoCuts.end() && isSemiLept!=1; ++ele){
        const pat::Electron& e = static_cast <const pat::Electron&> (*ele);
        float Smpt = e.pt();
        float Smeta = e.eta();
        float Smphi = e.phi();
        float SmJdR = deltaR(Smeta, Smphi, iJet->eta(), iJet->phi());
        if   ( Smpt> 5 && SmJdR <0.5) { //lep_ptCutForBjets_ = 5 GeV
	  isSemiLept=1;
	  jetSoftLeptdR = SmJdR;
	  jetSoftLeptPt = Smpt;
	  TLorentzVector jvec (iJet->p4().X(), iJet->p4().Y(), iJet->p4().Z(), iJet->p4().T());
	  TVector3 mvec ( e.p4().Vect().X(), e.p4().Vect().Y(), e.p4().Vect().Z()  );
	  jetSoftLeptPtRel = jvec.Perp(  mvec );
	  if (
	      ( fabs(Smeta)<2.5 && !( fabs(Smeta)>1.4442 && fabs(Smeta)<1.566))  &&
	      (( fabs(Smeta)>1.566  && e.sigmaIetaIeta()<0.01 && fabs(e.deltaPhiSuperClusterTrackAtVtx())<0.8 && fabs(e.deltaEtaSuperClusterTrackAtVtx())<0.007 ) ||
	       ( fabs(Smeta)<1.4442 && e.sigmaIetaIeta()<0.03 && fabs(e.deltaPhiSuperClusterTrackAtVtx())<0.7 && fabs(e.deltaEtaSuperClusterTrackAtVtx())<0.01 ))     )
	    jetSoftLeptIdEle95 = 1;
	}
      }
      jetSoftLeptdR_       .push_back(jetSoftLeptdR);
      jetSoftLeptPt_       .push_back(jetSoftLeptPt);
      jetSoftLeptPtRel_    .push_back(jetSoftLeptPtRel);
      jetSoftLeptIdEle95_  .push_back(jetSoftLeptIdEle95);

      // b-tagging
      jetCombinedSecondaryVtxBJetTags_.push_back(iJet->bDiscriminator("combinedSecondaryVertexBJetTags"));
      jetCombinedSecondaryVtxMVABJetTags_.push_back(iJet->bDiscriminator("combinedSecondaryVertexMVABJetTags"));
      jetJetProbabilityBJetTags_.push_back(iJet->bDiscriminator("jetProbabilityBJetTags"));  
      jetJetBProbabilityBJetTags_.push_back(iJet->bDiscriminator("jetBProbabilityBJetTags"));

      // betastar
      vector<float> jetBetaStar;
      for (int iVTX = 0; iVTX < nVtx_; ++iVTX) {
        double tracks_x = 0.;
        double tracks_y = 0.;
        double tracks_x_tot = 0.;
        double tracks_y_tot = 0.;
        for (unsigned i = 0;  i <  iJet->numberOfDaughters (); ++i) {	     
          const reco::PFCandidatePtr pfcand = iJet->getPFConstituent(i);
          reco::TrackRef trackref = pfcand->trackRef();
          if( trackref.isNonnull()) {
            // track_all
            tracks_x_tot += (trackref)->px();
            tracks_y_tot += (trackref)->py();
            for (int jVTX = 0; jVTX < nVtx_; ++jVTX) {
              if (jVTX == iVTX) continue;
 
              if (fabs((trackref)->vz()-vtx_z_[jVTX]) < 0.1) {        	  
                tracks_x += (trackref)->px();
                tracks_y += (trackref)->py();	
                break;	
              } // track_PU
            } // non-PV loop (assume that iVTX is PV)
          }
        } // jet associated tracks

        if (tracks_x_tot!=0. || tracks_y_tot!=0.) 
	  jetBetaStar.push_back(sqrt(tracks_x*tracks_x+tracks_y*tracks_y)/sqrt(tracks_x_tot*tracks_x_tot+tracks_y_tot*tracks_y_tot));
      } // vtx loop
      jetBetaStar_.push_back(jetBetaStar);

      //jet PF Loose ID
      pat::strbitset retjet = pfLooseId_.getBitTemplate();
      jetPFLooseId_.push_back(pfLooseId_(*iJet, retjet));

      //PileupJet ID variables
      float jetDRMean = 0; 
      float jetDR2Mean= 0; 
      float jetDZ     = 0; 
      float jetFrac01 = 0; 
      float jetFrac02 = 0; 
      float jetFrac03 = 0; 
      float jetFrac04 = 0;  
      float jetFrac05 = 0; 
      float jetFrac06 = 0;  
      float jetFrac07 = 0;
      float jetBeta = 0;    
      float jetBetaStarCMG = 0;
      float jetBetaStarClassic = 0;
      float jetNNeutrals = 0;
      float jetNCharged  = 0;
      
      float jetPuJetIdL=0;
      float jetPuJetIdM=0;
      float jetPuJetIdT=0;
    
      // for [nJet][imva]
      vector<float> jetMVAs;
      vector<int> jetWPLevels;

      vector<float> jetBetaExt;
      vector<float> jetBetaStarCMGExt;
      vector<float> jetBetaStarClassicExt;

      jetMVAsExt_simple_.clear();
      jetMVAsExt_full_.clear();
      jetMVAsExt_cutBased_.clear();
      jetMVAsExt_philv1_.clear();
      jetWPLevelsExt_simple_.clear();
      jetWPLevelsExt_full_.clear();
      jetWPLevelsExt_cutBased_.clear();
      jetWPLevelsExt_philv1_.clear();

//       for(unsigned int iVtx = 0; iVtx < 100; ++iVtx) {
//         for(unsigned int imva = 0; imva < jetMVAAlgos_.size(); ++imva) {
//           jetMVAsExt_[nJet_][imva][iVtx]  = -3.;
//           jetWPLevelsExt_[nJet_][imva][iVtx] = -1;
//         }
//       }

      if (pujetIDalgos_.size()>0) {
        float jecSetup = 0.0;
        const reco::VertexCollection vertexCollection = *(recVtxsBS_.product());
        const reco::Vertex* selectedVtx  = &(*vertexCollection.begin());;
        const pat::Jet* thisjet = &(*iJet);
  
        PileupJetIdentifier jetMVAinputs = jetMVACalculator->computeIdVariables( thisjet, jecSetup, selectedVtx, vertexCollection);
  
        jetDRMean  = jetMVAinputs.dRMean();
        jetDR2Mean = jetMVAinputs.dR2Mean();
        jetDZ      = jetMVAinputs.dZ();
        jetFrac01  = jetMVAinputs.frac01();
        jetFrac02  = jetMVAinputs.frac02();
        jetFrac03  = jetMVAinputs.frac03();
        jetFrac04  = jetMVAinputs.frac04();
        jetFrac05  = jetMVAinputs.frac05();
        jetFrac06  = jetMVAinputs.frac06();
        jetFrac07  = jetMVAinputs.frac07();
        jetBeta    = jetMVAinputs.beta();
        jetBetaStarCMG = jetMVAinputs.betaStar();
        jetBetaStarClassic = jetMVAinputs.betaStarClassic();
        jetNNeutrals = jetMVAinputs.nNeutrals();
        jetNCharged  = jetMVAinputs.nCharged();

	for(unsigned int imva=0; imva<jetMVAAlgos_.size(); ++imva){
          PileupJetIdAlgo* ialgo = (pujetIDalgos_[imva]);
          ialgo->set(jetMVAinputs);
          PileupJetIdentifier id = ialgo->computeMva();
          jetMVAs.push_back(id.mva());
          jetWPLevels.push_back(id.idFlag());

	  if(imva==1){
	    int    idflag = id.idFlag();
	    if( PileupJetIdentifier::passJetId( idflag, PileupJetIdentifier::kLoose )) {
	      jetPuJetIdL=1;
	    }
	    if( PileupJetIdentifier::passJetId( idflag, PileupJetIdentifier::kMedium )) {
	      jetPuJetIdM=1;
	    }
	    if( PileupJetIdentifier::passJetId( idflag, PileupJetIdentifier::kTight )) {
	      jetPuJetIdT=1;
	    }
	  }
        } // loop over algos

        vector<float> tmp_MVAsExt[4];
        vector<int> tmp_WPsExt[4];

        for(size_t iVtx=0; 
	    iVtx<vertexCollection.size() && iVtx < 100; ++iVtx) {
          PileupJetIdentifier jetMVAinputsExt = jetMVACalculator->computeIdVariables( thisjet, jecSetup, &vertexCollection[iVtx], vertexCollection);
          jetBetaExt           .push_back(jetMVAinputsExt.beta());
          jetBetaStarCMGExt    .push_back(jetMVAinputsExt.betaStar());
          jetBetaStarClassicExt.push_back(jetMVAinputsExt.betaStarClassic());
	  
	  vector<float> tmp_MVAs;
	  vector<int> tmp_WPs;
	  
          for(unsigned int imva=0; imva<jetMVAAlgos_.size(); ++imva){
            PileupJetIdAlgo* ialgo = (pujetIDalgos_[imva]);
            ialgo->set(jetMVAinputsExt);
            PileupJetIdentifier id = ialgo->computeMva();
            tmp_MVAs.push_back(id.mva());
            tmp_WPs .push_back(id.idFlag());

            tmp_MVAsExt[imva].push_back(id.mva());
            tmp_WPsExt[imva].push_back(id.idFlag());
          }
        } // loop over vtxs

        jetMVAsExt_simple_.push_back(tmp_MVAsExt[0]);
        jetMVAsExt_full_.push_back(tmp_MVAsExt[1]);
        jetMVAsExt_cutBased_.push_back(tmp_MVAsExt[2]);
        jetMVAsExt_philv1_.push_back(tmp_MVAsExt[3]);

        jetWPLevelsExt_simple_.push_back(tmp_WPsExt[0]);
        jetWPLevelsExt_full_.push_back(tmp_WPsExt[1]);
        jetWPLevelsExt_cutBased_.push_back(tmp_WPsExt[2]);
        jetWPLevelsExt_philv1_.push_back(tmp_WPsExt[3]);

      }
      // Fill the algo jet containers
      jetDRMean_ .push_back(jetDRMean);
      jetDR2Mean_.push_back(jetDR2Mean);
      jetDZ_     .push_back(jetDZ);      
      jetFrac01_ .push_back(jetFrac01);
      jetFrac02_ .push_back(jetFrac02);
      jetFrac03_ .push_back(jetFrac03);
      jetFrac04_ .push_back(jetFrac04);
      jetFrac05_ .push_back(jetFrac05);
      jetFrac06_ .push_back(jetFrac06);
      jetFrac07_ .push_back(jetFrac07);
      jetBeta_   .push_back(jetBeta);
      jetBetaStarCMG_.push_back(jetBetaStarCMG);
      jetBetaStarClassic_.push_back(jetBetaStarClassic);
      jetNNeutrals_.push_back(jetNNeutrals);
      jetNCharged_ .push_back(jetNCharged);

      jetPuJetIdL_.push_back(jetPuJetIdL);
      jetPuJetIdM_.push_back(jetPuJetIdM);
      jetPuJetIdT_.push_back(jetPuJetIdT);

      // vectors
      jetMVAs_              .push_back(jetMVAs);
      jetWPLevels_          .push_back(jetWPLevels);
      jetBetaExt_           .push_back(jetBetaExt);
      jetBetaStarCMGExt_    .push_back(jetBetaStarCMGExt);
      jetBetaStarClassicExt_.push_back(jetBetaStarClassicExt);

      // gen jet and parton
      jetPartonID_.push_back(iJet->partonFlavour());

      int jetGenPartonID = -99;
      int jetGenPartonMomID = -99;
      float jetGenEn  = 0;
      float jetGenPt  = 0;
      float jetGenEta = 0;
      float jetGenPhi = 0;
      if (!isData_ && genParticlesHandle_.isValid() ) {
        if ((*iJet).genParton()) {
          jetGenPartonID = (*iJet).genParton()->pdgId();
          jetGenEn       = (*iJet).genParton()->energy();
          jetGenPt       = (*iJet).genParton()->pt();
          jetGenEta      = (*iJet).genParton()->eta();
          jetGenPhi      = (*iJet).genParton()->phi();
	  
	  if ((*iJet).genParton()->mother()) {
	    jetGenPartonMomID = (*iJet).genParton()->mother()->pdgId();
	  }
	}
      }
      jetGenPartonID_.push_back(jetGenPartonID);
      jetGenPartonMomID_.push_back(jetGenPartonMomID);
      jetGenEn_      .push_back(jetGenEn);
      jetGenPt_      .push_back(jetGenPt);
      jetGenEta_     .push_back(jetGenEta);
      jetGenPhi_     .push_back(jetGenPhi);

      int jetGenJetIndex = -1;
      float jetGenJetEn  = -1;
      float jetGenJetPt  = -999;
      float jetGenJetEta = -999;
      float jetGenJetPhi = -999;      
      if (!isData_ && genParticlesHandle_.isValid() ) {
	if ((*iJet).genJet()) {
	  jetGenJetIndex = 1;
	  jetGenJetEn = (*iJet).genJet()->energy();
	  jetGenJetPt = (*iJet).genJet()->pt();
	  jetGenJetEta = (*iJet).genJet()->eta();
	  jetGenJetPhi = (*iJet).genJet()->phi();
	}
      }
      jetGenJetIndex_.push_back(jetGenJetIndex);
      jetGenJetEn_.push_back(jetGenJetEn);
      jetGenJetPt_.push_back(jetGenJetPt);
      jetGenJetEta_.push_back(jetGenJetEta);
      jetGenJetPhi_.push_back(jetGenJetPhi);

      nJet_++;
    }

  //#Lvdp
  // gluon tagger
  edm::Handle<edm::ValueMap<float> >  QGTagsHandleMLP;
  edm::Handle<edm::ValueMap<float> >  QGTagsHandleLikelihood;
  e.getByLabel(QGTagsHandleMLPLabel_, QGTagsHandleMLP);
  e.getByLabel(QGTagsHandleLikelihoodLabel_, QGTagsHandleLikelihood);
  
  QGTag_MLP_ = -999.;
  QGTag_likelihood_ = -999.;
  edm::Handle<edm::View<pat::Jet> > QGtagjetHandle;
  //  Handle<PFJetCollection> jets;
  e.getByLabel(QGtagjetLabel_, QGtagjetHandle);
  for (View<pat::Jet>::const_iterator iQGJet = QGtagjetHandle->begin(); iQGJet != QGtagjetHandle->end(); ++iQGJet){
    int ijet = iQGJet - QGtagjetHandle->begin();
    edm::RefToBase<pat::Jet> jetRef = QGtagjetHandle->refAt(ijet);//(edm::Ref<pat::Jet>(QGtagjetHandle,iQGJet));
    //cout<<ijet<<'\t'<<iQGJet->pt()<<'\t'<<iQGJet->eta()<<endl;
    if(iQGJet->pt()>20.){
      //  if (QGTagsHandleMLP.isValid()) cout << "MLP: " << (*QGTagsHandleMLP)[jetRef] << endl;
      //if (QGTagsHandleLikelihood.isValid()) cout << "Likelihood: " << (*QGTagsHandleLikelihood)[jetRef] << endl;
      if (QGTagsHandleMLP.isValid()) QGTag_MLP_ = (*QGTagsHandleMLP)[jetRef] ;
      if (QGTagsHandleLikelihood.isValid()) QGTag_likelihood_ = (*QGTagsHandleLikelihood)[jetRef] ;
    }
  }
  if(dumpSubJets_){
  nCA8Jet_ = 0;
  //jet substructure
  edm::Handle<edm::View<pat::Jet> > jetsCHSpruned;
  e.getByLabel(jetsCHSprunedLabel_, jetsCHSpruned);

  edm::Handle<edm::View<pat::Jet> > jetsCHS;
  e.getByLabel(jetsCHSLabel_, jetsCHS);

  View<pat::Jet>::const_iterator beginCHS = jetsCHS->begin();
  View<pat::Jet>::const_iterator endCHS = jetsCHS->end();
  View<pat::Jet>::const_iterator ijetCHS = beginCHS;
    
  View<pat::Jet>::const_iterator beginCHSpruned = jetsCHSpruned->begin();
  View<pat::Jet>::const_iterator endCHSpruned = jetsCHSpruned->end();
  View<pat::Jet>::const_iterator ijetCHSpruned = beginCHSpruned;

  // Loop over the "hard" jets
  for(ijetCHS = beginCHS; ijetCHS != endCHS; ++ijetCHS ){
    if( ijetCHS->pt() < 50.0 ) continue;
    nCA8Jet_++;
    CA8JetPt_.push_back( ijetCHS->pt() );
    CA8JetEta_.push_back( ijetCHS->eta() );
    CA8JetPhi_.push_back( ijetCHS->phi() );
    CA8JetMass_.push_back( ijetCHS->mass() );
    CA8JetArea_.push_back( ijetCHS->jetArea() );
    CA8Jet_tau1_.push_back( ijetCHS->userFloat("tau1") );
    CA8Jet_tau2_.push_back( ijetCHS->userFloat("tau2") );
    CA8Jet_tau3_.push_back( ijetCHS->userFloat("tau3") );
    
    View<pat::Jet>::const_iterator ijetCHSmatch;
    bool CA8matched = false;
    // iterate over pruned jets, choose the mathed one
    for( ijetCHSpruned = beginCHSpruned; ijetCHSpruned != endCHSpruned; ++ijetCHSpruned ){
      if(reco::deltaR(ijetCHS->momentum(), ijetCHSpruned->momentum()) < 0.25) {
	ijetCHSmatch = ijetCHSpruned;
	CA8matched = true;
	break;
      }
    }
    
    CA8prunedJetMass_.push_back( !CA8matched ? -999. : ijetCHSmatch->mass() );
    CA8prunedJet_nSubJets_.push_back( !CA8matched ? 0 : ijetCHSmatch->numberOfDaughters() );
    std::vector<TLorentzVector> subJets;
    if( CA8matched ){
      for(unsigned int idaughter=0; idaughter!=ijetCHSmatch->numberOfDaughters(); idaughter++){
	if( ijetCHSmatch->daughter(idaughter) != 0 && ijetCHSmatch->daughter(idaughter)->pt()>0.)
	  subJets.push_back( TLorentzVector(ijetCHSmatch->daughter(idaughter)->px(),
					    ijetCHSmatch->daughter(idaughter)->py(),
					    ijetCHSmatch->daughter(idaughter)->pz(),
					    ijetCHSmatch->daughter(idaughter)->energy())
			     );
      }
      
      //Now sort the vector.
      std::sort(subJets.begin(),subJets.end(),compareMass);
    }
    std::vector<float> SubjetPt;
    std::vector<float> SubjetEta;
    std::vector<float> SubjetPhi;
    std::vector<float> SubjetMass;

    for(std::vector<TLorentzVector>::iterator isubJet = subJets.begin(); isubJet!=subJets.end(); isubJet++){
      SubjetPt.push_back(isubJet->Perp());
      SubjetEta.push_back(isubJet->Eta());
      SubjetPhi.push_back(isubJet->Phi());
      SubjetMass.push_back(isubJet->M());
    }
    CA8prunedJet_SubjetPt_.push_back(SubjetPt);
    CA8prunedJet_SubjetEta_.push_back(SubjetEta);
    CA8prunedJet_SubjetPhi_.push_back(SubjetPhi);
    CA8prunedJet_SubjetMass_.push_back(SubjetMass);
  } // end loop over "hard" jets
  
  // end Lvdp
  }
  
  nConv_=0;
  
  if (convH_.isValid()) {
    for( reco::ConversionCollection::const_iterator  iConv = convH_->begin(); iConv != convH_->end(); ++iConv) {
      if (nConv_ >= 500) {
        cout << "WARNING :TOO MANY CONVERTED CANDIDATES: " << convH_->size() << endl;
        break;
      }
      reco::Conversion localConv = reco::Conversion(*iConv);

      convP4_x_.push_back(localConv.refittedPair4Momentum().px());
      convP4_y_.push_back(localConv.refittedPair4Momentum().py());
      convP4_z_.push_back(localConv.refittedPair4Momentum().pz());
      convP4_E_.push_back(localConv.refittedPair4Momentum().energy());

      convValidVtx_.push_back(localConv.conversionVertex().isValid());
      if ( !localConv.conversionVertex().isValid() ) {
	// put dummy values
	convVtx_x_.push_back(0);
	convVtx_y_.push_back(0);
	convVtx_z_.push_back(0);
	convVtxErr_x_.push_back(0);
	convVtxErr_y_.push_back(0);
	convVtxErr_z_.push_back(0);
	convPairMomentum_x_.push_back(0);
	convPairMomentum_y_.push_back(0);
	convPairMomentum_z_.push_back(0);
	convRefittedMomentum_x_.push_back(0);
	convRefittedMomentum_y_.push_back(0);
	convRefittedMomentum_z_.push_back(0);
	
	convChi2_.push_back(0);
	convChi2Probability_.push_back(0);
	convNTracks_.push_back(0);
	convMVALikelihood_.push_back(0);
	convTk1Dz_.push_back(0);
	convTk1DzErr_.push_back(0);
	convCh1Ch2_.push_back(0);
	//convTk2Dz_.push_back(0);// EXTRA
	//convTk2DzErr_.push_back(0);
	//convCh1Ch2_.push_back(0);

	convPairInvMass_.push_back(0);
	convPairCotThetaSep_.push_back(0);
	convEoverP_.push_back(0);
	convZofPrimVtxFromTrks_.push_back(0);
	convDistOfMinApproach_.push_back(0);
	convDPhiTrksAtVtx_.push_back(0);
	convDxy_.push_back(0);
	convDz_.push_back(0);
	convLxy_.push_back(0);
	convLz_.push_back(0);
	
	convNHitsBeforeVtx_0_.push_back(0);
	convNHitsBeforeVtx_1_.push_back(0);

	convNSharedHits_.push_back(0);
	
	convTk1D0_.push_back(0);
	convTk1Pout_.push_back(0);
	convTk1Pin_.push_back(0);
	
	convTk2D0_.push_back(0);
	convTk2Pout_.push_back(0);
	convTk2Pin_.push_back(0);
	
	continue;

      } else {
	reco::Vertex vtx=localConv.conversionVertex();

	convVtx_x_.push_back(vtx.x());
	convVtx_y_.push_back(vtx.y());
	convVtx_z_.push_back(vtx.z());
	convVtxErr_x_.push_back(vtx.xError());
	convVtxErr_y_.push_back(vtx.yError());
	convVtxErr_z_.push_back(vtx.zError());
	convPairMomentum_x_.push_back(localConv.pairMomentum().x());
	convPairMomentum_y_.push_back(localConv.pairMomentum().y());
	convPairMomentum_z_.push_back(localConv.pairMomentum().z());
	convRefittedMomentum_x_.push_back(localConv.refittedPairMomentum().x());
	convRefittedMomentum_y_.push_back(localConv.refittedPairMomentum().y());
	convRefittedMomentum_z_.push_back(localConv.refittedPairMomentum().z());
	
	convChi2_.push_back(vtx.chi2());
	convChi2Probability_.push_back(ChiSquaredProbability(vtx.chi2(), vtx.ndof()));
	convNTracks_.push_back(localConv.nTracks());
	convMVALikelihood_.push_back(localConv.MVAout());
	
	if( localConv.nTracks()) {
	  const std::vector<edm::RefToBase<reco::Track> > tracks = localConv.tracks();
	  for (unsigned int i=0; i<tracks.size(); ++i) {
	    if(i==0) {
	      convTk1Dz_.push_back(tracks[i]->dz());
	      convTk1DzErr_.push_back(tracks[i]->dzError());
	      convCh1Ch2_.push_back(tracks[i]->charge());
	    }
	    else if(i==1) {
	      convTk2Dz_.push_back(tracks[i]->dz());
	      convTk2DzErr_.push_back(tracks[i]->dzError());
	      convCh1Ch2_.back()*=(tracks[i]->charge());//EXTRA
	    }
	  }
	}
	
	convPairInvMass_.push_back(localConv.pairInvariantMass());
	convPairCotThetaSep_.push_back(localConv.pairCotThetaSeparation());
	convEoverP_.push_back(localConv.EoverPrefittedTracks());
	convZofPrimVtxFromTrks_.push_back(localConv.zOfPrimaryVertexFromTracks());
	convDistOfMinApproach_.push_back(localConv.distOfMinimumApproach());
	convDPhiTrksAtVtx_.push_back(localConv.dPhiTracksAtVtx());
	convDxy_.push_back(localConv.dxy());
	convDz_.push_back(localConv.dz());
	convLxy_.push_back(localConv.lxy());
	convLz_.push_back(localConv.lz());
	
	for (unsigned int i=0; i<localConv.nHitsBeforeVtx().size(); ++i) {
	  if (i==0) convNHitsBeforeVtx_0_.push_back(localConv.nHitsBeforeVtx()[i]);
	  if (i==1) convNHitsBeforeVtx_1_.push_back(localConv.nHitsBeforeVtx()[i]);
	}
	convNSharedHits_.push_back(localConv.nSharedHits());
	
	if(localConv.tracks().size() > 0) {
	  convTk1D0_.push_back(localConv.tracksSigned_d0()[0]);
	  convTk1Pout_.push_back(sqrt(localConv.tracksPout()[0].Mag2()));
	  convTk1Pin_.push_back(sqrt(localConv.tracksPin()[0].Mag2()));
	}
	
	if(localConv.tracks().size() > 1) {
          convTk2D0_.push_back(localConv.tracksSigned_d0()[1]);
          convTk2Pout_.push_back(sqrt(localConv.tracksPout()[1].Mag2()));
          convTk2Pin_.push_back(sqrt(localConv.tracksPin()[1].Mag2()));
	}
      } // conv vertex is valid
      
      nConv_++;
    }
    
  } // End of Converted Photon Collection
  
  hEvents_->Fill(1.5);
  if (doSkim_) {
    if (nEle_ >= 2 || nMu_ >= 2 || nGoodPho >= 2) tree_->Fill();
  } else {
    tree_->Fill();
  }

  delete lazyTool;
  delete topology_p;
  delete jecUnc;

}

void ggNtuplizer::beginJob() {
}

void ggNtuplizer::endJob() {
}

double ggNtuplizer::eT(double pt1, double pt2) const {
  double et = pt1 + pt2;
  return et;
}

double ggNtuplizer::massT(double pt1, double pt2, double wpx, double wpy) const {
  double mt = eT(pt1, pt2)*eT(pt1, pt2) - wpx*wpx - wpy*wpy;
  mt = (mt>0) ? sqrt(mt) : 0;
  return mt;
}

double ggNtuplizer::acop(double phi1, double phi2) const {
  Geom::Phi<double> deltaphi(phi1-phi2);
  double acop = deltaphi.value();
  if (acop<0) acop = - acop;
  acop = M_PI - acop;
  return acop;
}

float ggNtuplizer::getPhotonTrkIso(double photonVz, double vbsx, double vbsy, double vbsz, double photonEta, double photonPhi, edm::Handle<reco::TrackCollection> tracksHandle, double outercone, double innercone, double etastrip, double dZCut, double dxy, int Option) {
  
  float myiso = 0;
  math::XYZPoint vtx(vbsx, vbsy, vbsz);

  for (reco::TrackCollection::const_iterator trItr = tracksHandle.product()->begin(); trItr != tracksHandle_.product()->end(); ++trItr) {

    double dz = 0;
    if (Option == 0) dz = fabs( (*trItr).vz() - photonVz );
    else dz = (*trItr).dz(vtx);
    if (fabs(dz) > dZCut) continue;

    double this_dxy = (*trItr).dxy(vtx);
    if (fabs(this_dxy) > dxy) continue;

    double dr = deltaR(photonEta, photonPhi, (*trItr).eta(), (*trItr).phi());
    double deta = (*trItr).eta() - photonEta;

    if (dr < outercone && dr >= innercone && fabs(deta) >= etastrip) myiso += (*trItr).pt();
  }

  return myiso;
}

float ggNtuplizer::getPhotondRtoTrk(edm::Handle<reco::GsfElectronCollection> gsfHandle, const pat::Photon & pho, float minPt, int maxMissingHits) {

  float dr = 99;
  for (reco::GsfElectronCollection::const_iterator igsf = gsfHandle->begin(); igsf != gsfHandle->end(); ++igsf) {

    reco::GsfElectron ele = reco::GsfElectron(*igsf);

    if (ele.gsfTrack()->trackerExpectedHitsInner().numberOfHits() > maxMissingHits) continue;
    //if (ele.pt() < minPt) continue;                                                                                                                                              

    if (fabs(ele.superCluster()->eta() - pho.superCluster()->eta()) == 0.)
      if (fabs(ele.superCluster()->phi() - pho.superCluster()->phi()) == 0.) {
        dr = sqrt(pow(ele.deltaEtaSuperClusterTrackAtVtx(),2)+pow(ele.deltaPhiSuperClusterTrackAtVtx(),2));
      }
  }

  return dr;
}

float ggNtuplizer::getGenCalIso(edm::Handle<reco::GenParticleCollection> handle, reco::GenParticleCollection::const_iterator thisPho,
				const Float_t dRMax, bool removeMu, bool removeNu) {
  
  const Float_t etMin = 0.0;
  Float_t genCalIsoSum = 0.0;
  if (!doGenParticles_)  return genCalIsoSum;
  if (!handle.isValid()) return genCalIsoSum;

  for (reco::GenParticleCollection::const_iterator it_gen=handle->begin(); it_gen!=handle->end(); ++it_gen) {

    if (it_gen == thisPho) continue;        // can't be the original photon
    if (it_gen->status() != 1) continue;    // need to be a stable particle
    if (thisPho->collisionId() != it_gen->collisionId()) continue; // has to come from the same collision
   
    Int_t pdgCode = abs(it_gen->pdgId());
    // we should not count neutrinos, muons
    if (removeMu && pdgCode == 13 ) continue;
    if (removeNu && (pdgCode == 12 || pdgCode == 14 || pdgCode == 16)) continue;

    Float_t et = it_gen->et();
    if (et < etMin) continue; // pass a minimum et threshold, default 0

    Float_t dR = reco::deltaR(thisPho->momentum(), it_gen->momentum());
    if (dR > dRMax) continue; // within deltaR cone

    genCalIsoSum += et;
    
  } 

  return genCalIsoSum;
}

float ggNtuplizer::getGenTrkIso(edm::Handle<reco::GenParticleCollection> handle, reco::GenParticleCollection::const_iterator thisPho, const Float_t dRMax) {

  const Float_t ptMin = 0.0;
  Float_t genTrkIsoSum = 0.0;
  if (!doGenParticles_)  return genTrkIsoSum;
  if (!handle.isValid()) return genTrkIsoSum;

  for (reco::GenParticleCollection::const_iterator it_gen=handle->begin(); it_gen!=handle->end(); ++it_gen){

    if (it_gen == thisPho) continue;        // can't be the original photon
    if (it_gen->status() != 1) continue;    // need to be a stable particle
    if (thisPho->collisionId() != it_gen->collisionId()) continue; // has to come from the same collision
   
    if (it_gen->charge() == 0) continue;    // we should not count neutral particles
   
    Float_t pt = it_gen->pt();
    if (pt < ptMin) continue; // pass a minimum pt threshold, default 0

    Float_t dR = reco::deltaR(thisPho->momentum(), it_gen->momentum());
    if (dR > dRMax) continue; // within deltaR cone
    genTrkIsoSum += pt;
    
  }// end of loop over gen particles

  return genTrkIsoSum;
}



vector<float> ggNtuplizer::getESHits(double X, double Y, double Z, map<DetId, EcalRecHit> rechits_map, const CaloSubdetectorGeometry*& geometry_p, CaloSubdetectorTopology *topology_p, int row) {

  //cout<<row<<endl;

  vector<float> esHits;

  //double X = bcPtr->x();
  //double Y = bcPtr->y();
  //double Z = bcPtr->z();
  const GlobalPoint point(X,Y,Z);

  DetId esId1 = (dynamic_cast<const EcalPreshowerGeometry*>(geometry_p))->getClosestCellInPlane(point, 1);
  DetId esId2 = (dynamic_cast<const EcalPreshowerGeometry*>(geometry_p))->getClosestCellInPlane(point, 2);
  ESDetId esDetId1 = (esId1 == DetId(0)) ? ESDetId(0) : ESDetId(esId1);
  ESDetId esDetId2 = (esId2 == DetId(0)) ? ESDetId(0) : ESDetId(esId2);  

  map<DetId, EcalRecHit>::iterator it;
  ESDetId next;
  ESDetId strip1;
  ESDetId strip2;

  strip1 = esDetId1;
  strip2 = esDetId2;

  EcalPreshowerNavigator theESNav1(strip1, topology_p);
  theESNav1.setHome(strip1);
  
  EcalPreshowerNavigator theESNav2(strip2, topology_p);
  theESNav2.setHome(strip2);

  if (row == 1) {
    if (strip1 != ESDetId(0)) strip1 = theESNav1.north();
    if (strip2 != ESDetId(0)) strip2 = theESNav2.east();
  } else if (row == -1) {
    if (strip1 != ESDetId(0)) strip1 = theESNav1.south();
    if (strip2 != ESDetId(0)) strip2 = theESNav2.west();
  }

  // Plane 1 
  if (strip1 == ESDetId(0)) {
    for (unsigned int i=0; i<31; ++i) esHits.push_back(0);
  } else {
    
    it = rechits_map.find(strip1);
    if (it->second.energy() > 1.0e-10 && it != rechits_map.end()) esHits.push_back(it->second.energy());  
    else esHits.push_back(0);
    //cout<<"center : "<<strip1<<" "<<it->second.energy()<<endl;      

    // east road 
    for (unsigned int i=0; i<15; ++i) {
      next = theESNav1.east();
      if (next != ESDetId(0)) {
	it = rechits_map.find(next);
	if (it->second.energy() > 1.0e-10 && it != rechits_map.end()) esHits.push_back(it->second.energy());  
	else esHits.push_back(0);
	//cout<<"east "<<i<<" : "<<next<<" "<<it->second.energy()<<endl;
      } else {
        for (unsigned int j=i; j<15; ++j) esHits.push_back(0);
        break;
	//cout<<"east "<<i<<" : "<<next<<" "<<0<<endl;
      }
    }

    // west road 
    theESNav1.setHome(strip1);
    theESNav1.home();
    for (unsigned int i=0; i<15; ++i) {
      next = theESNav1.west();
      if (next != ESDetId(0)) {
	it = rechits_map.find(next);
	if (it->second.energy() > 1.0e-10 && it != rechits_map.end()) esHits.push_back(it->second.energy());  
	else esHits.push_back(0);
	//cout<<"west "<<i<<" : "<<next<<" "<<it->second.energy()<<endl;
      } else {
        for (unsigned int j=i; j<15; ++j) esHits.push_back(0);
        break;
	//cout<<"west "<<i<<" : "<<next<<" "<<0<<endl;
      }
    }
  }

  if (strip2 == ESDetId(0)) {
    for (unsigned int i=0; i<31; ++i) esHits.push_back(0);
  } else {

    it = rechits_map.find(strip2);
    if (it->second.energy() > 1.0e-10 && it != rechits_map.end()) esHits.push_back(it->second.energy());  
    else esHits.push_back(0);
    //cout<<"center : "<<strip2<<" "<<it->second.energy()<<endl;      

    // north road 
    for (unsigned int i=0; i<15; ++i) {
      next = theESNav2.north();
      if (next != ESDetId(0)) {
	it = rechits_map.find(next);
	if (it->second.energy() > 1.0e-10 && it != rechits_map.end()) esHits.push_back(it->second.energy());
	else esHits.push_back(0);
	//cout<<"north "<<i<<" : "<<next<<" "<<it->second.energy()<<endl;  
      } else {
        for (unsigned int j=i; j<15; ++j) esHits.push_back(0);
        break;
	//cout<<"north "<<i<<" : "<<next<<" "<<0<<endl;
      }
    }

    // south road 
    theESNav2.setHome(strip2);
    theESNav2.home();
    for (unsigned int i=0; i<15; ++i) {
      next = theESNav2.south();
      if (next != ESDetId(0)) {
	it = rechits_map.find(next);
	if (it->second.energy() > 1.0e-10 && it != rechits_map.end()) esHits.push_back(it->second.energy());  
	else esHits.push_back(0);
	//cout<<"south "<<i<<" : "<<next<<" "<<it->second.energy()<<endl;
      } else {
        for (unsigned int j=i; j<15; ++j) esHits.push_back(0);
        break;
	//cout<<"south "<<i<<" : "<<next<<" "<<0<<endl;
      }
    }
  }

  return esHits;
}

vector<float> ggNtuplizer::getESEffSigmaRR(vector<float> ESHits0)
{
  const int nBIN = 21;
  vector<float> esShape;

  TH1F *htmpF = new TH1F("htmpF","",nBIN,0,nBIN);
  TH1F *htmpR = new TH1F("htmpR","",nBIN,0,nBIN);
  htmpF->Reset(); htmpR->Reset();

  Float_t effsigmaRR=0.;

  for(int ibin=0; ibin<((nBIN+1)/2); ++ibin) {
    if (ibin==0) {
      htmpF->SetBinContent((nBIN+1)/2,ESHits0[ibin]);
      htmpR->SetBinContent((nBIN+1)/2,ESHits0[ibin+31]);
    } else { // hits sourd the seed
      htmpF->SetBinContent((nBIN+1)/2+ibin,ESHits0[ibin]);
      htmpF->SetBinContent((nBIN+1)/2-ibin,ESHits0[ibin+15]);
      htmpR->SetBinContent((nBIN+1)/2+ibin,ESHits0[ibin+31]);
      htmpR->SetBinContent((nBIN+1)/2-ibin,ESHits0[ibin+31+15]);
    }
  }

  // ---- Effective Energy Deposit Width ---- //
  double EffWidthSigmaXX = 0.;
  double EffWidthSigmaYY = 0.;
  double totalEnergyXX   = 0.;
  double totalEnergyYY   = 0.;
  double EffStatsXX      = 0.;
  double EffStatsYY      = 0.;
  for (int id_X=1; id_X<=21; ++id_X) {
    totalEnergyXX  += htmpF->GetBinContent(id_X);
    EffStatsXX     += htmpF->GetBinContent(id_X)*(id_X-11)*(id_X-11);
    totalEnergyYY  += htmpR->GetBinContent(id_X);
    EffStatsYY     += htmpR->GetBinContent(id_X)*(id_X-11)*(id_X-11);
  }
  // If denominator == 0, effsigmaRR = 0;
  EffWidthSigmaXX  = (totalEnergyXX>0.)  ? sqrt(fabs(EffStatsXX  / totalEnergyXX))   : 0.;
  EffWidthSigmaYY  = (totalEnergyYY>0.)  ? sqrt(fabs(EffStatsYY  / totalEnergyYY))   : 0.;
  effsigmaRR =  ((totalEnergyXX  + totalEnergyYY) >0.) ? sqrt(EffWidthSigmaXX  * EffWidthSigmaXX  + EffWidthSigmaYY  * EffWidthSigmaYY)  : 0.;
  esShape.push_back(effsigmaRR);
  esShape.push_back(EffWidthSigmaXX);
  esShape.push_back(EffWidthSigmaYY);

  delete htmpF;
  delete htmpR;

  return esShape;
}

vector<float> ggNtuplizer::getESEn(vector<float> ESHits0)
{
  const int nBIN = 21;
  vector<float> esEn;

  float esf_e1 = 0.;    float esr_e1 = 0.;    // Sum of energy deposited in p/m  0 strip
  float esf_e3 = 0.;    float esr_e3 = 0.;    // Sum of energy deposited in p/m  1 strip
  float esf_e5 = 0.;    float esr_e5 = 0.;    // Sum of energy deposited in p/m  2 strips
  float esf_e7 = 0.;    float esr_e7 = 0.;    // Sum of energy deposited in p/m  3 strips
  float esf_e11= 0.;    float esr_e11= 0.;    // Sum of energy deposited in p/m  5 strips
  float esf_e21= 0.;    float esr_e21= 0.;    // Sum of energy deposited in p/m 10 strips

  // Insert energy from central strip
  esf_e1  += ESHits0[0];     esr_e1  += ESHits0[31];
  esf_e3  += ESHits0[0];     esr_e3  += ESHits0[31];
  esf_e5  += ESHits0[0];     esr_e5  += ESHits0[31];
  esf_e7  += ESHits0[0];     esr_e7  += ESHits0[31];
  esf_e11 += ESHits0[0];     esr_e11 += ESHits0[31];
  esf_e21 += ESHits0[0];     esr_e21 += ESHits0[31];
  // hits sround the central strip
  for(int ibin=1; ibin<((nBIN+1)/2); ++ibin) {
    if (ibin<=1) {
      esf_e3 += ESHits0[ibin];
      esf_e3 += ESHits0[ibin+15];
      esr_e3 += ESHits0[ibin+31];
      esr_e3 += ESHits0[ibin+31+15];
    }
    if (ibin<=2) {
      esf_e5 += ESHits0[ibin];
      esf_e5 += ESHits0[ibin+15];
      esr_e5 += ESHits0[ibin+31];
      esr_e5 += ESHits0[ibin+31+15];
    }
    if (ibin<=3) {
      esf_e7 += ESHits0[ibin];
      esf_e7 += ESHits0[ibin+15];
      esr_e7 += ESHits0[ibin+31];
      esr_e7 += ESHits0[ibin+31+15];
    }
    if (ibin<=5) {
      esf_e11 += ESHits0[ibin];
      esf_e11 += ESHits0[ibin+15];
      esr_e11 += ESHits0[ibin+31];
      esr_e11 += ESHits0[ibin+31+15];
    }
    esf_e21 += ESHits0[ibin];
    esf_e21 += ESHits0[ibin+15];
    esr_e21 += ESHits0[ibin+31];
    esr_e21 += ESHits0[ibin+31+15];
  }

  esEn.push_back(esf_e1);  esEn.push_back(esr_e1);
  esEn.push_back(esf_e3);  esEn.push_back(esr_e3);
  esEn.push_back(esf_e5);  esEn.push_back(esr_e5);
  esEn.push_back(esf_e7);  esEn.push_back(esr_e7);
  esEn.push_back(esf_e11); esEn.push_back(esr_e11);
  esEn.push_back(esf_e21); esEn.push_back(esr_e21);

  return esEn;
}

double ggNtuplizer::fsrPhotonIso03(const reco::PFCandidateCollection* pfAllCandidates, const reco::PFCandidateCollection* pfCandidates, const reco::PFCandidate pfPhoton) {

  double chIso = 0;
  double phIso = 0;
  double nhIso = 0;
  double puIso = 0;
  double dr    = 0;

  for (unsigned iCand=0; iCand< pfCandidates->size(); ++iCand) {

    const reco::PFCandidate  & pfParticleT((*pfCandidates)[iCand]);
    const reco::PFCandidate  *pfParticle = (&pfParticleT);
   
    dr = deltaR(pfParticle->eta(), pfParticle->phi(), pfPhoton.eta(), pfPhoton.phi()); 

    if (dr > 0.3 || dr < 0.01 || pfParticle->pt() < 0.2) continue;

    if (PFCandidate::ParticleType(pfParticle->particleId()) == PFCandidate::h)
      if (dr > 0.01 && pfParticle->pt() > 0.2) chIso += pfParticle->pt();

    if (dr > 0.01 && pfParticle->pt() > 0.5 && PFCandidate::ParticleType(pfParticle->particleId()) == PFCandidate::gamma) phIso += pfParticle->pt(); 

    if (PFCandidate::ParticleType(pfParticle->particleId()) == PFCandidate::h0)
      if (dr > 0.01 && pfParticle->pt() > 0.5) nhIso += pfParticle->pt(); 
    
  }
  
  for (unsigned jCand=0; jCand< pfAllCandidates->size(); ++jCand) {
    const reco::PFCandidate  & pfAllParticleT((*pfAllCandidates)[jCand]);
    const reco::PFCandidate  *pfAllParticle = (&pfAllParticleT);
    
    if (pfAllParticle->charge() == 0) continue;
    int puFlag = 1;

    for (unsigned iCand=0; iCand< pfCandidates->size(); ++iCand) {
      
      const reco::PFCandidate  & pfParticleB((*pfCandidates)[iCand]);
      const reco::PFCandidate  *pfParticleBB = (&pfParticleB);

      if (pfParticleBB->charge() == 0) continue;
      dr = deltaR(pfAllParticle->eta(), pfAllParticle->phi(), pfParticleBB->eta(), pfParticleBB->phi());
      if (dr == 0) puFlag = 0;
    }    

    if (puFlag == 0) continue;
    dr = deltaR(pfAllParticle->eta(), pfAllParticle->phi(), pfPhoton.eta(), pfPhoton.phi());
    if (dr > 0.3) continue;
    if (dr > 0.01 && pfAllParticle->pt() > 0.2) puIso += pfAllParticle->pt();
  }

  return (chIso + phIso + nhIso + puIso) / (pfPhoton.ecalEnergy() * pfPhoton.pt() / pfPhoton.energy());
}

void ggNtuplizer::fillBlockInfo(reco::PFCandidate PFCand, bool doPhotons, bool reclassify, const CaloSubdetectorGeometry* geomBar, const CaloSubdetectorGeometry* geomEnd, const CaloSubdetectorGeometry* geometry_p
				){
  if(develop_){
    //first match to Photon in the Translator: 
    
    ggPFClusters PFClusterCollection(EBReducedRecHits_, EEReducedRecHits_, geomBar, geomEnd);
    
    if(doPhotons){
      reco::SuperClusterRef pfSuperCluster;
      vector<reco::CaloCluster>PF;
      reco::ConversionRefVector  conversions;
      reco::ConversionRefVector  SLconversions;
      if(!reclassify){
	reco::PhotonCollection::const_iterator pfphot=pfPhoTranslator_->begin();
	for(;pfphot!=pfPhoTranslator_->end();++pfphot){
	  
	  if(PFCand.superClusterRef()==pfphot->superCluster()){
	    pfSuperCluster=pfphot->pfSuperCluster();
	    PF=PFClusterCollection.getPFClusters(*pfSuperCluster);
	    conversions=pfphot->conversions();
	    SLconversions=pfphot->conversionsOneLeg();
	    break;
	  }
	}        
      }
      else{
	reco::GsfElectronCollection::const_iterator pfele=gsfElectronHandle_->begin();
	for(;pfele!=gsfElectronHandle_->end();++pfele){
	  if(pfele->superCluster().isNull())continue;	
	  if(pfele->pflowSuperCluster().isNull())continue;
	  if(pfele->gsfTrack()==PFCand.gsfTrackRef()){
	    pfSuperCluster=pfele->pflowSuperCluster();
	    PF=PFClusterCollection.getPFClusters(*pfSuperCluster);
	    
	    break;
	  }
	}
      	
      }
     
      std::map<unsigned int,unsigned int> ESEELinks;
      ggPFESClusters PFPreShowerColl(ESRecHits_, geomEnd);
      vector<reco::PreshowerCluster> PFPS;
      if (dumpESClusterInfo_) {
	if (pfSuperCluster->seed()->seed().subdetId()!=1){
	  PFPS = PFPreShowerColl.getPFESClusters(*pfSuperCluster);
	  // If there are no ES clusters, ESEELinks is empty; if for a ES cluster there are no linked EE clusters, there is no entry.
	  //ESEELinks = PFPreShowerColl.getClosestEECls(PFPS, PF);
	}
      }


      nPFPho_ES1clust=0;
      nPFPho_ES2clust=0; 
      for(unsigned int clust=0; clust<PF.size();++clust){
	float tmp_ES1dist = 0., tmp_ES2dist = 0.;
	float tmp_distEE_ES1 = 0., tmp_distEE_ES2 = 0.;
	int tmp_es1s = 0, tmp_es2s = 0;
        float tmp_ESsieta = 0., tmp_ESetaln = 0., tmp_ESsiphi = 0., tmp_ESphiln = 0.;
	
	PFPho_clustE_[nPFPho_][clust]=PF[clust].energy();
	PFPho_clustEt_[nPFPho_][clust]=PF[clust].energy()*sin(PF[clust].position().theta());
	PFPho_clusteta_[nPFPho_][clust]=PF[clust].eta();
	PFPho_clustphi_[nPFPho_][clust]=PF[clust].phi();
	
	//find PFSeed and get indices
	
      std::vector< std::pair<DetId, float> >bcCells=PF[clust].hitsAndFractions();
      bool isEb;
      DetId seedXtalId = bcCells[0].first ;
      int detector = seedXtalId.subdetId();
      if(detector==1)isEb=true;
      else isEb=false;
            
      PFPho_clustEseed_[nPFPho_][clust]=PFClusterCollection.get5x5Element(0, 0, bcCells, isEb);
      
      PFPho_clustEtop_[nPFPho_][clust]=PFClusterCollection.get5x5Element(0, 1, bcCells, isEb);
      PFPho_clustEbottom_[nPFPho_][clust]=PFClusterCollection.get5x5Element(0, -1, bcCells, isEb);
      PFPho_clustEleft_[nPFPho_][clust]=PFClusterCollection.get5x5Element(-1, 0, bcCells, isEb);
      PFPho_clustEright_[nPFPho_][clust]=PFClusterCollection.get5x5Element(1, 0, bcCells, isEb);
      //float e3x1=0; float e1x3=0; float e3x3=0;
      
      float e5x5Map[5][5];
      float e5x5=0;
      for(int i=-2; abs(i)<3; ++i)
	for(int j=-2; abs(j)<3; ++j){
	int ind1=i+2;
	int ind2=j+2;
	  e5x5Map[ind1][ind2]=PFClusterCollection.get5x5Element(i, j, bcCells, isEb);
	  e5x5=e5x5+PFClusterCollection.get5x5Element(i, j, bcCells, isEb);
	  //cout<<"e "<<i<<", "<<j<<" "<<e5x5Map[i+2][j+2]<<endl;
	}
      
      PFPho_clustE3x3_[nPFPho_][clust]=e5x5Map[2][2]+e5x5Map[2][1]+e5x5Map[2][3]+e5x5Map[3][2]+e5x5Map[1][2]+e5x5Map[3][3]+e5x5Map[1][1]+e5x5Map[3][1]+e5x5Map[1][3];
      
      PFPho_clustE3x1_[nPFPho_][clust]=e5x5Map[2][2]+e5x5Map[1][2]+e5x5Map[3][2];
      PFPho_clustE1x3_[nPFPho_][clust]=e5x5Map[2][2]+e5x5Map[2][1]+e5x5Map[2][3];
      PFPho_clustE5x1_[nPFPho_][clust]=e5x5Map[2][2]+e5x5Map[0][2]+e5x5Map[1][2]+e5x5Map[3][2]+e5x5Map[4][2];
      PFPho_clustE1x5_[nPFPho_][clust]=e5x5Map[2][2]+e5x5Map[2][0]+e5x5Map[2][1]+e5x5Map[2][3]+e5x5Map[2][4];
      
      PFPho_clustE5x5_[nPFPho_][clust]=e5x5;
      PFPho_clustE2x5Top_[nPFPho_][clust]=e5x5Map[0][4]+e5x5Map[1][4]+e5x5Map[2][4]+e5x5Map[3][4]+e5x5Map[4][4]
	+e5x5Map[0][3]+e5x5Map[1][3]+e5x5Map[2][3]+e5x5Map[3][3]+e5x5Map[4][3];
      //add up bottom edge of 5x5 2 rows
      PFPho_clustE2x5Bottom_[nPFPho_][clust]=e5x5Map[0][0]+e5x5Map[1][0]+e5x5Map[2][0]+e5x5Map[3][0]+e5x5Map[4][0]
	+e5x5Map[0][1]+e5x5Map[1][1]+e5x5Map[2][1]+e5x5Map[3][1]+e5x5Map[4][1];
      //add up left edge of 5x5 2 rows
      PFPho_clustE2x5Left_[nPFPho_][clust]=e5x5Map[0][0]+e5x5Map[0][1]+e5x5Map[0][2]+e5x5Map[0][3]+e5x5Map[0][4]
	+e5x5Map[1][0]+e5x5Map[1][1]+e5x5Map[1][2]+e5x5Map[1][3]+e5x5Map[1][4];
      //add up right edge of 5x5 2 rows
      PFPho_clustE2x5Right_[nPFPho_][clust]=e5x5Map[4][0]+e5x5Map[4][1]+e5x5Map[4][2]+e5x5Map[4][3]+e5x5Map[4][4]
	+e5x5Map[3][0]+e5x5Map[3][1]+e5x5Map[3][2]+e5x5Map[3][3]+e5x5Map[3][4];
      //find max 2x5 from the center
      float centerstrip=e5x5Map[2][2]+e5x5Map[2][0]+e5x5Map[2][1]+e5x5Map[2][3]+e5x5Map[2][4];
      float rightstrip=e5x5Map[3][2]+e5x5Map[3][0]+e5x5Map[3][1]+e5x5Map[3][3]+e5x5Map[3][4];
      float leftstrip=e5x5Map[1][2]+e5x5Map[1][0]+e5x5Map[1][1]+e5x5Map[1][3]+e5x5Map[1][4];
      if(rightstrip>leftstrip)PFPho_clustE2x5Max_[nPFPho_][clust]=rightstrip+centerstrip;
      else PFPho_clustE2x5Max_[nPFPho_][clust]=leftstrip+centerstrip;
      if(isEb){
	float etacry, phicry,thetatilt,  phitilt; 
	int ieta,iphi;
	PFClusterCollection.localCoordsEB(PF[clust], etacry, phicry, ieta, iphi, thetatilt, phitilt);
	PFPho_crysphifix_[nPFPho_][clust]=phicry;
	PFPho_crysetafix_[nPFPho_][clust]=etacry;
	PFPho_crysIeta_[nPFPho_][clust]=ieta; 
	PFPho_crysIphi_[nPFPho_][clust]=iphi;
	//PFPho_crysphiTilt_[nPFPho_][clust]=phitilt; 	 
	//PFPho_crysthetaTilt_[nPFPho_][clust]=thetatilt;
	//cout<<"From Photon Fix "<< PFClustFix.etaC()<<endl;
      }
      else{		    
	float xcry, ycry, thetatilt,  phitilt;
	int ix, iy;
	PFClusterCollection.localCoordsEE(PF[clust], xcry, ycry, ix, iy, thetatilt, phitilt);
	PFPho_crysxfix_[nPFPho_][clust]=xcry;
	PFPho_crysyfix_[nPFPho_][clust]=ycry;
	PFPho_crysIX_[nPFPho_][clust]=ix; 
	PFPho_crysIY_[nPFPho_][clust]=iy;
	PFPho_crysphiTilt_[nPFPho_][clust]=phitilt; 	 
	PFPho_crysthetaTilt_[nPFPho_][clust]=thetatilt;
      }
      //cout<<"nC "<<clust<<endl;

      // PreShower 
      float tmp_PS1TotE=0, tmp_PS2TotE=0;
      float max1=0., max2=0., totE1=0., totE2=0.;
      int index1=-1, index2=-1;
      if (dumpESClusterInfo_) {
	if(!isEb) {
	  if(PFPS.size()!=0){
	    if(ESEELinks.size()==0){
	      cout << endl << ">>>>>> ES clusters without links map! <<<<<<";
	      continue;
	    }
	    for(unsigned int p=0; p<PFPS.size(); ++p){
	      std::map<unsigned int,unsigned int>::iterator it = ESEELinks.find(p);
	      if(it==ESEELinks.end()) continue;
	      unsigned int closestEECl = (*it).second;
	      if(closestEECl!=clust) continue;

	      double linkDist = PFPreShowerColl.getLinkDist(PFPS[p], PF[clust]);
	      if(linkDist==-1.) {
		continue;
	      }

  	      // Quick loops for catching the clusters with max energy 
	      if(PFPS[p].plane()==1){
	        if(PFPS[p].energy()>max1){
	          max1 = PFPS[p].energy();
	          index1 = p;
	        }
	        totE1 += PFPS[p].energy();
	      }
	      if(PFPS[p].plane()==2){
	        if(PFPS[p].energy()>max2){
		  max2 = PFPS[p].energy();
		  index2 = p;
	        }
	        totE2 += PFPS[p].energy();
	      }
 
	      if(PFPS[p].plane()==1){
		float tmp_ES1stripdist = 0.;
		tmp_PS1TotE += PFPS[p].energy();
		PFPho_ES1clusteta_[nPFPho_][nPFPho_ES1clust] = PFPS[p].eta();
		PFPho_ES1clustphi_[nPFPho_][nPFPho_ES1clust] = PFPS[p].phi();
		PFPho_ES1clustx_[nPFPho_][nPFPho_ES1clust] = PFPS[p].x();
		PFPho_ES1clusty_[nPFPho_][nPFPho_ES1clust] = PFPS[p].y();
		PFPho_ES1clustz_[nPFPho_][nPFPho_ES1clust] = PFPS[p].z();
		PFPho_ES1clustE_[nPFPho_][nPFPho_ES1clust] = PFPS[p].energy();
		PFPho_ES1clustLinkDist_[nPFPho_][nPFPho_ES1clust] = linkDist;
		tmp_ES1dist += PFPS[p].x()*PFPS[p].energy();
  	        if ( (PFPS[p].phi()>3.14159/4. && PFPS[p].phi()<3.14159*0.75) || (PFPS[p].phi()>-3.14159*0.75 && PFPS[p].phi()<-3.14159/4.) ){
		  tmp_ESsiphi += TMath::Log(PFPS[p].energy()/totE1) * (PFPS[p].phi()-PFPS[index1].phi()) * (PFPS[p].phi()-PFPS[index1].phi());
                  tmp_ESphiln += TMath::Log(PFPS[p].energy()/totE1);
		}
		else {
                  tmp_ESsieta += TMath::Log(PFPS[p].energy()/totE1) * (PFPS[p].eta()-PFPS[index1].eta()) * (PFPS[p].eta()-PFPS[index1].eta());
                  tmp_ESetaln += TMath::Log(PFPS[p].energy()/totE1);
                }
		tmp_distEE_ES1 += linkDist;

		// ES rechits 
		std::vector< std::pair<DetId, float> > psCells = PFPS[p].hitsAndFractions();
		unsigned int s1 = 0;
		for(s1=0; s1<psCells.size(); ++s1){
		  if(s1>29) continue;
		  DetId hitid = DetId(psCells[s1].first.rawId());
		  PFPho_ES1_stripsDetId_[nPFPho_][nPFPho_ES1clust][s1] = hitid;
		  const CaloCellGeometry *strip = geometry_p->getGeometry(hitid);
		  if(!strip) continue;
		  const GlobalPoint strippos = strip->getPosition();
		  PFPho_ES1_stripsX_[nPFPho_][nPFPho_ES1clust][s1] = strippos.x();
		  PFPho_ES1_stripsY_[nPFPho_][nPFPho_ES1clust][s1] = strippos.y();
		  PFPho_ES1_stripsZ_[nPFPho_][nPFPho_ES1clust][s1] = strippos.z();
		  PFPho_ES1_stripsEta_[nPFPho_][nPFPho_ES1clust][s1] = strippos.eta();
		  PFPho_ES1_stripsPhi_[nPFPho_][nPFPho_ES1clust][s1] = strippos.phi();
		  PFPho_ES1_stripsFrac_[nPFPho_][nPFPho_ES1clust][s1] = psCells[s1].second;
		  for(EcalRecHitCollection::const_iterator es=ESRecHits_->begin(); es!= ESRecHits_->end(); ++es){
		    ESDetId sid = ESDetId(es->id().rawId());
		    if(sid.plane()!=1) continue;
		    if(es->id().rawId()!=hitid.rawId()) continue;
		    float E = es->energy() * psCells[s1].second;
		    PFPho_ES1_stripsE_[nPFPho_][nPFPho_ES1clust][s1] = E;
                   tmp_ES1stripdist += E*(PFPS[p].x()-strippos.x());
		    break;
		  }
		}
		PFPho_ES1size_[nPFPho_][nPFPho_ES1clust] = s1;
		tmp_es1s += s1;
               PFPho_weightedXst_perES1clust_[nPFPho_][nPFPho_ES1clust] = tmp_ES1stripdist/PFPS[p].energy();
		nPFPho_ES1clust++;
	      } // ES clusters plane 1 
	      
	      if(PFPS[p].plane()==2){
		float tmp_ES2stripdist = 0.;
		tmp_PS2TotE += PFPS[p].energy();
		PFPho_ES2clusteta_[nPFPho_][nPFPho_ES2clust] = PFPS[p].eta();
		PFPho_ES2clustphi_[nPFPho_][nPFPho_ES2clust] = PFPS[p].phi();
		PFPho_ES2clustx_[nPFPho_][nPFPho_ES2clust] = PFPS[p].x();
		PFPho_ES2clusty_[nPFPho_][nPFPho_ES2clust] = PFPS[p].y();
		PFPho_ES2clustz_[nPFPho_][nPFPho_ES2clust] = PFPS[p].z();
		PFPho_ES2clustE_[nPFPho_][nPFPho_ES2clust] = PFPS[p].energy();
		PFPho_ES2clustLinkDist_[nPFPho_][nPFPho_ES2clust] = linkDist;
		tmp_ES2dist += PFPS[p].y()*PFPS[p].energy();
                if ( (PFPS[p].phi()>3.14159/4. && PFPS[p].phi()<3.14159*0.75) || (PFPS[p].phi()>-3.14159*0.75 && PFPS[p].phi()<-3.14159/4.) ){
                  tmp_ESsieta += TMath::Log(PFPS[p].energy()/totE2) * (PFPS[p].eta()-PFPS[index2].eta()) * (PFPS[p].eta()-PFPS[index2].eta());
                  tmp_ESetaln += TMath::Log(PFPS[p].energy()/totE2);
		}
		else {
		  tmp_ESsiphi += TMath::Log(PFPS[p].energy()/totE2) * (PFPS[p].phi()-PFPS[index2].phi()) * (PFPS[p].phi()-PFPS[index2].phi());
		  tmp_ESphiln += TMath::Log(PFPS[p].energy()/totE2);
		}
		tmp_distEE_ES2 += linkDist;
		
		// ES rechits
		std::vector< std::pair<DetId, float> > psCells = PFPS[p].hitsAndFractions();
		unsigned int s2 = 0;
		for(s2=0; s2<psCells.size(); ++s2){
		  if(s2>29) continue;
		  DetId hitid = psCells[s2].first.rawId();
		  PFPho_ES2_stripsDetId_[nPFPho_][nPFPho_ES2clust][s2] = hitid;
		  const CaloCellGeometry *strip = geometry_p->getGeometry(hitid);
		  if(!strip) continue;
		  const GlobalPoint& strippos = strip->getPosition();
		  PFPho_ES2_stripsX_[nPFPho_][nPFPho_ES2clust][s2] = strippos.x();
		  PFPho_ES2_stripsY_[nPFPho_][nPFPho_ES2clust][s2] = strippos.y();
		  PFPho_ES2_stripsZ_[nPFPho_][nPFPho_ES2clust][s2] = strippos.z();
		  PFPho_ES2_stripsEta_[nPFPho_][nPFPho_ES2clust][s2] = strippos.eta();
		  PFPho_ES2_stripsPhi_[nPFPho_][nPFPho_ES2clust][s2] = strippos.phi();
		  PFPho_ES2_stripsFrac_[nPFPho_][nPFPho_ES2clust][s2] = psCells[s2].second;
		  for(EcalRecHitCollection::const_iterator es=ESRecHits_->begin(); es!= ESRecHits_->end(); ++es){
		    ESDetId sid = ESDetId(es->id().rawId());
		    if(sid.plane()!=2) continue;
		    if(es->id().rawId()!=hitid.rawId()) continue;
		    float E = es->energy() * psCells[s2].second;
		    PFPho_ES2_stripsE_[nPFPho_][nPFPho_ES2clust][s2] = E;
                    tmp_ES2stripdist += E*(PFPS[p].y()-strippos.y());
		    break;
		  }
		}
		
		PFPho_ES2size_[nPFPho_][nPFPho_ES2clust] = s2;
		tmp_es2s += s2;
                PFPho_weightedYst_perES2clust_[nPFPho_][nPFPho_ES2clust] = tmp_ES2stripdist/PFPS[p].energy();
		nPFPho_ES2clust++;
	      } // ES clusters plane 2 
	      
	    } // loop on ES planes 
	  } // check PFPS size 
	} // if !isEb
	
        PFPho_ES1Energy_[nPFPho_][clust] = tmp_PS1TotE;
        PFPho_ES2Energy_[nPFPho_][clust] = tmp_PS2TotE;
        nPFPhoClust_[nPFPho_] = clust;
        nPFPhoCrys_[nPFPho_][clust] = PF[clust].size();
        nPFPhoES1Clust_[nPFPho_] = nPFPho_ES1clust;
        nPFPhoES2Clust_[nPFPho_] = nPFPho_ES2clust;
      
        nPFPho_ES1Clust_perEEclust_[nPFPho_][clust] = nPFPho_ES1clust;
        nPFPho_ES2Clust_perEEclust_[nPFPho_][clust] = nPFPho_ES2clust;
        PFPho_ES1weightedX_perEEclust_[nPFPho_][clust] = PF[clust].x()-tmp_ES1dist/tmp_PS1TotE;
        PFPho_ES2weightedY_perEEclust_[nPFPho_][clust] = PF[clust].y()-tmp_ES2dist/tmp_PS2TotE;
        PFPho_ES_siphiiphi_[nPFPho_][clust] = tmp_ESsiphi / tmp_ESphiln;
        PFPho_ES_sietaieta_[nPFPho_][clust] = tmp_ESsieta / tmp_ESetaln;
        PFPho_ES1linkD_perEEclust_[nPFPho_][clust] = tmp_distEE_ES1;
        PFPho_ES2linkD_perEEclust_[nPFPho_][clust] = tmp_distEE_ES2;
        PFPho_ES1size_perEEclust_[nPFPho_][clust] = tmp_es1s;
        PFPho_ES2size_perEEclust_[nPFPho_][clust] = tmp_es2s;

       } // if dumpESClusterInfo_ 
      } // clust loop 
      
      std::pair<float, float>widths=PFClusterCollection.ClusterWidth(PF);
      
      PFPhoMaxEtaWidth_[nPFPho_]=widths.first;
      PFPhoMaxPhiWidth_[nPFPho_]=widths.second;
      nPFPhoClust_[nPFPho_]=PF.size();
      Mustache Must;
      Must.FillMustacheVar(PF);
      PFPhoMustEtOut_[nPFPho_]=Must.MustacheEtOut();
      PFPhoMustExcl_[nPFPho_]=Must.OutsideMust();     
      PFPhoMustEin_[nPFPho_]=Must.MustacheE();
      PFPhoMustEout_[nPFPho_]=Must.MustacheEOut();
      nPFPho_tks=0;
      nPFPhoTks_[nPFPho_]=nPFPho_tks;
      if(!reclassify){	
	for(unsigned int c=0; c<conversions.size(); ++c){
	  std::vector<math::XYZPointF> innerPos=conversions[c]->tracksInnerPosition();
	  std::vector<math::XYZVectorF>innerMom=conversions[c]->tracksPin();
	  const std::vector<edm::RefToBase<reco::Track> > tracks = conversions[c]->tracks();
	  for (unsigned int t=0; t<tracks.size(); ++t){
	    PFPho_tkerrthet_[nPFPho_][nPFPho_tks]=tracks[t]->thetaError();
		
	    PFPho_tkR_[nPFPho_][nPFPho_tks]=sqrt( innerPos[t].X() * innerPos[t].X() + innerPos[t].Y() * innerPos[t].Y());
	    PFPho_tkPos_[nPFPho_][nPFPho_tks][0]=  innerPos[t].X();
	    PFPho_tkPos_[nPFPho_][nPFPho_tks][1]=  innerPos[t].Y();
	    PFPho_tkPos_[nPFPho_][nPFPho_tks][2]=  innerPos[t].Z();
	    PFPho_tkTheta_[nPFPho_][nPFPho_tks]=   innerPos[t].Theta();
	    PFPho_tkpt_[nPFPho_][nPFPho_tks]=sqrt(innerMom[t].X()*innerMom[t].X() + innerMom[t].Y()* innerMom[t].Y() );
	    //cout<<"nPFPho_tks "<<nPFPho_tks<<endl;
	    nPFPho_tks++;
	    nPFPhoTks_[nPFPho_]=nPFPho_tks;
	    
	  }
	} 
	for(unsigned int SLc=0; SLc<conversions.size(); ++SLc){
	  std::vector<math::XYZPointF> innerPos=conversions[SLc]->tracksInnerPosition();
	  //std::vector<math::XYZVectorF>innerMom=conversions[SLc]->tracksPin();
	  const std::vector<edm::RefToBase<reco::Track> > tracks = conversions[SLc]->tracks();
	  for (unsigned int t=0; t<tracks.size(); ++t) {
	    //cout<<"tracks "<<tracks.size()<<endl;	   
	    PFPho_tkerreta_[nPFPho_][nPFPho_tks]=tracks[t]->etaError();
	    PFPho_tkerrphi_[nPFPho_][nPFPho_tks]=tracks[t]->phiError();
	    PFPho_tkerrthet_[nPFPho_][nPFPho_tks]=tracks[t]->thetaError();
	    
	    PFPho_tkR_[nPFPho_][nPFPho_tks]=sqrt(innerPos[t].X()* innerPos[t].X()+innerPos[t].Y() * innerPos[t].Y());
	    PFPho_tkPos_[nPFPho_][nPFPho_tks][0]=  tracks[t]->referencePoint().X();
	    PFPho_tkPos_[nPFPho_][nPFPho_tks][1]=  tracks[t]->referencePoint().Y();
	    PFPho_tkPos_[nPFPho_][nPFPho_tks][2]=  tracks[t]->referencePoint().Z();
	    PFPho_tkTheta_[nPFPho_][nPFPho_tks]=   innerPos[t].Theta();
	    PFPho_tkpt_[nPFPho_][nPFPho_tks]= tracks[t]->pt();//sqrt(innerMom[t].X()*innerMom[t].X() + innerMom[t].Y()* innerMom[t].Y() );
	    nPFPho_tks++;
	    nPFPhoTks_[nPFPho_]=nPFPho_tks;
	    
	  }
	}			
      }
    }  
    else{ //for Electrons:
      
      reco::SuperClusterRef pfSuperCluster;
      vector<reco::CaloCluster> PF;
      
      reco::GsfElectronCollection::const_iterator pfele=gsfElectronHandle_->begin();
      for(;pfele!=gsfElectronHandle_->end();++pfele){
	if(pfele->superCluster().isNull())continue;	
	if(pfele->pflowSuperCluster().isNull())continue;
	if(pfele->gsfTrack()==PFCand.gsfTrackRef()){
	  pfSuperCluster=pfele->pflowSuperCluster();
	  break;
	}
      }
      
      PF=PFClusterCollection.getPFClusters(*pfSuperCluster); 
      
      for(unsigned int clust=0; clust<PF.size();++clust){
	PFEle_clustE_[nPFEle_][clust]=PF[clust].energy();
	PFEle_clustEt_[nPFEle_][clust]=PF[clust].energy()*sin(PF[clust].position().theta());
	PFEle_clusteta_[nPFEle_][clust]=PF[clust].eta();
	PFEle_clustphi_[nPFEle_][clust]=PF[clust].phi();
	std::vector< std::pair<DetId, float> >bcCells=PF[clust].hitsAndFractions();
	bool isEb;
	DetId seedXtalId = bcCells[0].first ;
	int detector = seedXtalId.subdetId();
	if(detector==1)isEb=true;
	else isEb=false;
	PFEle_clustEseed_[nPFEle_][clust]=PFClusterCollection.get5x5Element(0, 0, bcCells, isEb);
	PFEle_clustEtop_[nPFEle_][clust]=PFClusterCollection.get5x5Element(0, 1, bcCells, isEb);
	PFEle_clustEbottom_[nPFEle_][clust]=PFClusterCollection.get5x5Element(0, -1, bcCells, isEb);
	PFEle_clustEleft_[nPFEle_][clust]=PFClusterCollection.get5x5Element(-1, 0, bcCells, isEb);
	PFEle_clustEright_[nPFEle_][clust]=PFClusterCollection.get5x5Element(1, 0, bcCells, isEb);
	//float e3x1=0; float e1x3=0; float e3x3=0;
	float e[5][5];
	float e5x5=0;
	for(int i=-2; abs(i)<3; ++i)
	  for(int j=-2; abs(j)<3; ++j){
	    e[i+2][j+2]=PFClusterCollection.get5x5Element(i, j, bcCells, isEb);
	    e5x5=e5x5+PFClusterCollection.get5x5Element(i, j, bcCells, isEb);
	    //cout<<"e "<<i<<", "<<j<<" "<<e[i+2][j+2]<<endl;
	}
	PFEle_clustE3x3_[nPFEle_][clust]=e[2][2]+e[2][1]+e[2][3]+e[3][2]+e[1][2]+e[3][3]+e[1][1]+e[3][1]+e[1][3];
	
	PFEle_clustE3x1_[nPFEle_][clust]=e[2][2]+e[1][2]+e[3][2];
	PFEle_clustE1x3_[nPFEle_][clust]=e[2][2]+e[2][1]+e[2][3];
	PFEle_clustE5x1_[nPFEle_][clust]=e[2][2]+e[0][2]+e[1][2]+e[3][2]+e[4][2];
	PFEle_clustE1x5_[nPFEle_][clust]=e[2][2]+e[2][0]+e[2][1]+e[2][3]+e[2][4];
	
	PFEle_clustE5x5_[nPFEle_][clust]=e5x5;
	PFEle_clustE2x5Top_[nPFEle_][clust]=e[0][4]+e[1][4]+e[2][4]+e[3][4]+e[4][4]
	  +e[0][3]+e[1][3]+e[2][3]+e[3][3]+e[4][3];
	//add up bottom edge of 5x5 2 rows
	PFEle_clustE2x5Bottom_[nPFEle_][clust]=e[0][0]+e[1][0]+e[2][0]+e[3][0]+e[4][0]
	  +e[0][1]+e[1][1]+e[2][1]+e[3][1]+e[4][1];
	//add up left edge of 5x5 2 rows
	PFEle_clustE2x5Left_[nPFEle_][clust]=e[0][0]+e[0][1]+e[0][2]+e[0][3]+e[0][4]
	  +e[1][0]+e[1][1]+e[1][2]+e[1][3]+e[1][4];
	//add up right edge of 5x5 2 rows
	PFEle_clustE2x5Right_[nPFEle_][clust]=e[4][0]+e[4][1]+e[4][2]+e[4][3]+e[4][4]
	  +e[3][0]+e[3][1]+e[3][2]+e[3][3]+e[3][4];
	//find max 2x5 from the center
	float centerstrip=e[2][2]+e[2][0]+e[2][1]+e[2][3]+e[2][4];
	float rightstrip=e[3][2]+e[3][0]+e[3][1]+e[3][3]+e[3][4];
	float leftstrip=e[1][2]+e[1][0]+e[1][1]+e[1][3]+e[1][4];
	if(rightstrip>leftstrip)PFEle_clustE2x5Max_[nPFEle_][clust]=rightstrip+centerstrip;
	else PFEle_clustE2x5Max_[nPFEle_][clust]=leftstrip+centerstrip;
	
	DetId idseed=PFClusterCollection.FindSeed(bcCells, isEb);
	if(idseed.subdetId()==EcalBarrel){
	  float etacry, phicry,thetatilt,  phitilt; 
	  int ieta,iphi;
	  PFClusterCollection.localCoordsEB(PF[clust], etacry, phicry, ieta, iphi, thetatilt, phitilt);
	  PFEle_crysphifix_[nPFEle_][clust]=phicry;
	  PFEle_crysetafix_[nPFEle_][clust]=etacry;
	  PFEle_crysIeta_[nPFEle_][clust]=ieta; 
	  PFEle_crysIphi_[nPFEle_][clust]=iphi;
	  //PFEle_crysphiTilt_[nPFEle_][clust]=phitilt; 	 
	  //PFEle_crysthetaTilt_[nPFEle_][clust]=thetatilt;
	  //cout<<"From Eleton Fix "<< PFClustFix.etaC()<<endl;
	}
	else{		    
	  
	  float xcry, ycry, thetatilt,  phitilt;
	  int ix, iy;
	  PFClusterCollection.localCoordsEE(PF[clust], xcry, ycry, ix, iy, thetatilt, phitilt);
	  PFEle_crysxfix_[nPFEle_][clust]=xcry;
	  PFEle_crysyfix_[nPFEle_][clust]=ycry;
	  PFEle_crysIX_[nPFEle_][clust]=ix; 
	  PFEle_crysIY_[nPFEle_][clust]=iy;
	  //PFEle_crysphiTilt_[nPFEle_][clust]=phitilt; 	 
	  //PFEle_crysthetaTilt_[nPFEle_][clust]=thetatilt;
	}
	
	// Preshower 
	float PS1TotE = 0;
	float PS2TotE = 0;
      	
	if(!isEb) {
	  ggPFESClusters PFPreShowerColl(ESRecHits_, geomEnd);
	  vector<reco::PreshowerCluster>PFPS=PFPreShowerColl.getPFESClusters(*pfSuperCluster);
	  if(PFPS.size()!=0){
	    for(unsigned int p=0; p<PFPS.size(); ++p){
	      double linkDist = PFPreShowerColl.getLinkDist(PFPS[p], PF[clust]);
	      if(linkDist==-1) continue;

	      if(PFPS[p].plane()==1){
		PS1TotE += PFPS[p].energy();
		PFEle_ES1clusteta_[nPFEle_][nPFEle_ES1clust] = PFPS[p].eta();
		PFEle_ES1clustphi_[nPFEle_][nPFEle_ES1clust] = PFPS[p].phi();
		PFEle_ES1clustx_[nPFEle_][nPFEle_ES1clust] = PFPS[p].x();
		PFEle_ES1clusty_[nPFEle_][nPFEle_ES1clust] = PFPS[p].y();
		PFEle_ES1clustz_[nPFEle_][nPFEle_ES1clust] = PFPS[p].z();
		PFEle_ES1clustE_[nPFEle_][nPFEle_ES1clust] = PFPS[p].energy();
		nPFEle_ES1clust++;
	      }
	      if(PFPS[p].plane()==2){
		PS2TotE += PFPS[p].energy();
		PFEle_ES2clusteta_[nPFEle_][nPFEle_ES2clust] = PFPS[p].eta();
		PFEle_ES2clustphi_[nPFEle_][nPFEle_ES2clust] = PFPS[p].phi();
		PFEle_ES2clustx_[nPFEle_][nPFEle_ES2clust] = PFPS[p].x();
		PFEle_ES2clusty_[nPFEle_][nPFEle_ES2clust] = PFPS[p].y();
		PFEle_ES2clustz_[nPFEle_][nPFEle_ES2clust] = PFPS[p].z();
		PFEle_ES2clustE_[nPFEle_][nPFEle_ES2clust] = PFPS[p].energy();
		nPFEle_ES2clust++;
	      }
	    }
	    
	  }
	}
	PFEle_ES1Energy_[nPFEle_][clust] = PS1TotE;
	PFEle_ES2Energy_[nPFEle_][clust] = PS2TotE;
	nPFEleES1Clust_[nPFEle_] = nPFEle_ES1clust;
	nPFEleES2Clust_[nPFEle_] = nPFEle_ES2clust;
	
      } // clust loop 
      
      
      std::pair<float, float>widths=PFClusterCollection.ClusterWidth(PF);
      PFEleMaxEtaWidth_[nPFEle_]=widths.first;
      PFEleMaxPhiWidth_[nPFEle_]=widths.second; 
      Mustache Must;
      Must.FillMustacheVar(PF);
      PFEleMustEtOut_[nPFEle_]=Must.MustacheEtOut();
      PFEleMustExcl_[nPFEle_]=Must.OutsideMust();     
      PFEleMustEin_[nPFEle_]=Must.MustacheE();
      PFEleMustEout_[nPFEle_]=Must.MustacheEOut();
    }
      
  }

}
bool ggNtuplizer::SLPoint(const reco::Photon phot, 
			  reco::ConversionRefVector  &SLconversions,
			  std::vector<float> &Zint 
			  ){
  bool hasSL=false;
  reco::PhotonCollection::const_iterator pfphot=pfPhoTranslator_->begin();
  for(;pfphot!=pfPhoTranslator_->end();++pfphot){
    
    if(phot.superCluster()==pfphot->superCluster()){
      SLconversions=pfphot->conversionsOneLeg();
      break;
    }
  }
  if(SLconversions.size()>0)hasSL=true;
  TVector3 beamspot(beamSpotHandle_->position().x(),beamSpotHandle_->position().y(),
		    beamSpotHandle_->position().z());
  TVector3 SCPos(phot.superCluster()->position().x()-beamspot[0], phot.superCluster()->position().y()-beamspot[1], phot.superCluster()->position().z()-beamspot[2]);
  for(unsigned int SL=0; SL<SLconversions.size(); ++SL){
    TVector3 TkPos(
		   SLconversions[SL]->conversionVertex().x()-beamspot.X(), 
		   SLconversions[SL]->conversionVertex().y()-beamspot.Y(),
		   SLconversions[SL]->conversionVertex().z()-beamspot.Z());
    //Intersection fromt the two points:
    float R1=sqrt(SCPos.X()* SCPos.X() + SCPos.Y()*SCPos.Y()); 
    float R2=sqrt(TkPos.X()* TkPos.X() + TkPos.Y()*TkPos.Y());
    float Z1=SCPos.Z();
    float Z2=TkPos.Z();
    float slope=(Z1-Z2)/(R1-R2);
    float Z=Z2 - R2*slope;
    Zint.push_back(Z); 
  }  
  return hasSL;  
}

DEFINE_FWK_MODULE(ggNtuplizer);
