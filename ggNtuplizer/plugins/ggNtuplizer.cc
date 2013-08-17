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
  runOnParticleGun_  = ps.getParameter<bool>("runOnParticleGun");
  develop_           = ps.getParameter<bool>("development");
  doCentrality_      = ps.getParameter<bool>("doCentrality");
  
  vtxlabel_          = ps.getParameter<InputTag>("VtxLabel");
  tracklabel_        = ps.getParameter<InputTag>("TrackLabel");
  gsfElectronlabel_  = ps.getParameter<InputTag>("gsfElectronLabel");
  pfMETlabel_        = ps.getParameter<InputTag>("pfMETLabel");
  recoPfMETlabel_    = ps.getParameter<InputTag>("recoPfMETLabel");

  genParticlesCollection_    = ps.getParameter<InputTag>("genParticleSrc");
  METCollection_             = ps.getParameter<InputTag>("METSrc");
  electronCollection_        = ps.getParameter<InputTag>("electronSrc");
  photonCollection_          = ps.getParameter<InputTag>("photonSrc");
  recophotonCollection_      = ps.getParameter<InputTag>("Photons");
  muonCollection_            = ps.getParameter<InputTag>("muonSrc");
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
  rhoCollection44_           = ps.getParameter<InputTag>("rhoCollection44");
  rho2011Label_		     = ps.getParameter<InputTag>("rho2011Label");
  pfAllParticles_            = ps.getParameter<InputTag>("PFAllCandidates");
  pfParticles_               = ps.getParameter<InputTag>("PFCandidates");
  pfPhotonParticles_         = ps.getParameter<InputTag>("PFCandidatePhotons");
  allConversionsColl_        = ps.getParameter<InputTag>("ConvertedPhotonColl");
  pfPhotonCollection_        = ps.getParameter<InputTag>("PFPhotons");
  jetMVAAlgos_               = ps.getUntrackedParameter<std::vector<edm::ParameterSet> >("puJetIDAlgos");
  pfLooseId_                 = ps.getParameter<edm::ParameterSet>("pfLooseId");

  // jetMVAs_.resize(jetMVAAlgos_.size());
  // jetWPLevels_.resize(jetMVAAlgos_.size());
  // jetMVAsExt_.resize(jetMVAAlgos_.size());
  // jetWPLevelsExt_.resize(jetMVAAlgos_.size());
  pujetIDalgos_.resize(jetMVAAlgos_.size());
  for(unsigned int imva=0; imva<jetMVAAlgos_.size(); imva++){
    pujetIDalgos_[imva] = new PileupJetIdAlgo((jetMVAAlgos_.at(imva)));
  }

  inputTagPhotonIsoDeposits_ = ps.getParameter<std::vector<edm::InputTag> >("PhotonIsoDeposits");

  // PF isolation
  inputTagIsoDepElectrons_     = ps.getParameter< std::vector<edm::InputTag> >("IsoDepElectron");
  inputTagIsoDepPhotons_       = ps.getParameter< std::vector<edm::InputTag> >("IsoDepPhoton");
  inputTagIsoValElectronsPFId_ = ps.getParameter< std::vector<edm::InputTag> >("IsoValElectronPF");
  inputTagIsoValPhotonsPFId_   = ps.getParameter< std::vector<edm::InputTag> >("IsoValPhoton");

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
  tree_->Branch("vtx", vtx_, "vtx[nVtx][3]/F");
  tree_->Branch("IsVtxGood", &IsVtxGood_, "IsVtxGood/I");
  tree_->Branch("nGoodVtx", &nGoodVtx_, "nGoodVtx/I");
  if (doCentrality_) tree_->Branch("centrality", centrality_, "centrality[5]/F");
  tree_->Branch("nVtxBS", &nVtxBS_, "nVtxBS/I");
  tree_->Branch("vtxbs", vtxbs_, "vtxbs[nVtxBS][3]/F");
  if (dumpTrks_) {
    tree_->Branch("vtxbsPtMod", vtxbsPtMod_, "vtxbsPtMod[nVtxBS]/F");
    tree_->Branch("vtxbsSumPt2", vtxbsSumPt2_, "vtxbsSumPt2[nVtxBS]/F");
    tree_->Branch("vtxbsTkIndex", "std::vector<std::vector<Int_t> >",   &vtxbsTkIndex_);
    tree_->Branch("vtxbsTkWeight","std::vector<std::vector<Float_t> >", &vtxbsTkWeight_);
    tree_->Branch("nTrk", &nTrk_, "nTrk/I");
    tree_->Branch("trkP", trkP_, "trkP[nTrk][3]/F");
    tree_->Branch("trkVtx", trkVtx_, "trkVtx[nTrk][3]/F");
    tree_->Branch("trkd0", trkd0_, "trkd0[nTrk]/F");
    tree_->Branch("trkd0Err", trkd0Err_, "trkd0Err[nTrk]/F");
    tree_->Branch("trkdz", trkdz_, "trkdz[nTrk]/F");
    tree_->Branch("trkdzErr", trkdzErr_, "trkdzErr[nTrk]/F");
    tree_->Branch("trkPtErr", trkPtErr_, "trkPtErr[nTrk]/F");
    tree_->Branch("trkQuality", trkQuality_, "trkQuality[nTrk]/I");
    tree_->Branch("nGoodTrk", &nGoodTrk_, "nGoodTrk/I");
    tree_->Branch("IsTracksGood", &IsTracksGood_, "IsTracksGood/I");
  }
  if (doGenParticles_) {
    tree_->Branch("pdf", pdf_, "pdf[7]/F");
    tree_->Branch("pthat", &pthat_, "pthat/F");
    tree_->Branch("processID", &processID_, "processID/F");
    // genParticle
    tree_->Branch("nMC", &nMC_, "nMC/I");
    tree_->Branch("mcPID", mcPID, "mcPID[nMC]/I");
    tree_->Branch("mcVtx", mcVtx, "mcVtx[nMC][3]/F");
    tree_->Branch("mcPt", mcPt, "mcPt[nMC]/F");
    tree_->Branch("mcMass", mcMass, "mcMass[nMC]/F");
    tree_->Branch("mcEta", mcEta, "mcEta[nMC]/F");
    tree_->Branch("mcPhi", mcPhi, "mcPhi[nMC]/F");
    tree_->Branch("mcE", mcE, "mcE[nMC]/F");
    tree_->Branch("mcEt", mcEt, "mcEt[nMC]/F");
    tree_->Branch("mcGMomPID", mcGMomPID, "mcGMomPID[nMC]/I");
    tree_->Branch("mcMomPID", mcMomPID, "mcMomPID[nMC]/I");
    tree_->Branch("mcMomPt", mcMomPt, "mcMomPt[nMC]/F");
    tree_->Branch("mcMomMass", mcMomMass, "mcMomMass[nMC]/F");
    tree_->Branch("mcMomEta", mcMomEta, "mcMomEta[nMC]/F");
    tree_->Branch("mcMomPhi", mcMomPhi, "mcMomPhi[nMC]/F");
    tree_->Branch("mcIndex", mcIndex, "mcIndex[nMC]/I");
    tree_->Branch("mcDecayType", mcDecayType, "mcDecayType[nMC]/I"); //-999:non W or Z, 1:hardronic, 2:e, 3:mu, 4:tau
    tree_->Branch("mcParentage", mcParentage, "mcParentage[nMC]/I"); // 16*lepton + 8*boson + 4*non-prompt + 2*qcd + exotics
    tree_->Branch("mcStatus", mcStatus, "mcStatus[nMC]/I"); // status of the particle
    // Gen MET
    tree_->Branch("genMET", &genMET_, "genMET/F");
    tree_->Branch("genMETPhi", &genMETPhi_, "genMETPhi/F");
    // PU Info
    tree_->Branch("nPUInfo", &nPUInfo_, "nPUInfo/I");
    tree_->Branch("nPU", nPU_, "nPU[nPUInfo]/I");
    tree_->Branch("puBX", puBX_, "puBX[nPUInfo]/I");
    tree_->Branch("puTrue", puTrue_, "puTrue[nPUInfo]/F");
  }
  // pfMET
  tree_->Branch("pfMET", &pfMET_, "pfMET/F");
  tree_->Branch("pfMETPhi", &pfMETPhi_, "pfMETPhi/F");
  tree_->Branch("pfMETsumEt", &pfMETsumEt_, "pfMETsumEt/F");
  tree_->Branch("pfMETmEtSig", &pfMETmEtSig_, "pfMETmEtSig/F");
  tree_->Branch("pfMETSig", &pfMETSig_, "pfMETSig/F");
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
  tree_->Branch("trkMETx",      trkMETx_,     "trkMETx[nVtxBS]/F");
  tree_->Branch("trkMETy",      trkMETy_,     "trkMETy[nVtxBS]/F");
  tree_->Branch("trkMETPhi",    trkMETPhi_,   "trkMETPhi[nVtxBS]/F");
  tree_->Branch("trkMET",       trkMET_,      "trkMET[nVtxBS]/F");
  // MET filters
  tree_->Branch("metFilters", metFilters_, "metFilters[10]/I");
  // Electron
  tree_->Branch("nEle", &nEle_, "nEle/I");
  tree_->Branch("eleTrg", eleTrg_, "eleTrg[nEle][16]/I");
  tree_->Branch("eleClass", eleClass_, "eleClass[nEle]/I");
  tree_->Branch("eleIsEcalDriven", eleIsEcalDriven_, "eleIsEcalDriven[nEle]/I");
  tree_->Branch("eleCharge", eleCharge_, "eleCharge[nEle]/I");
  tree_->Branch("eleChargeConsistent", eleChargeConsistent_, "eleChargeConsistent[nEle]/I");
  tree_->Branch("eleEn", eleEn_, "eleEn[nEle]/F");
  tree_->Branch("eleEcalEn", eleEcalEn_, "eleEcalEn[nEle]/F");
  tree_->Branch("eleSCRawEn", eleSCRawEn_, "eleSCRawEn[nEle]/F");
  tree_->Branch("eleSCEn", eleSCEn_, "eleSCEn[nEle]/F");
  tree_->Branch("eleESEn", eleESEn_, "eleESEn[nEle]/F");
  tree_->Branch("elePt", elePt_, "elePt[nEle]/F");
  tree_->Branch("eleEta", eleEta_, "eleEta[nEle]/F");
  tree_->Branch("elePhi", elePhi_, "elePhi[nEle]/F");
  tree_->Branch("eleEtaVtx", eleEtaVtx_, "eleEtaVtx[nEle][100]/F");
  tree_->Branch("elePhiVtx", elePhiVtx_, "elePhiVtx[nEle][100]/F");
  tree_->Branch("eleEtVtx", eleEtVtx_, "eleEtVtx[nEle][100]/F");
  tree_->Branch("eleSCEta", eleSCEta_, "eleSCEta[nEle]/F");
  tree_->Branch("eleSCPhi", eleSCPhi_, "eleSCPhi[nEle]/F");
  tree_->Branch("eleSCEtaWidth", eleSCEtaWidth_, "eleSCEtaWidth[nEle]/F");
  tree_->Branch("eleSCPhiWidth", eleSCPhiWidth_, "eleSCPhiWidth[nEle]/F");
  tree_->Branch("eleVtx", eleVtx_, "eleVtx[nEle][3]/F");
  tree_->Branch("eleD0", eleD0_, "eleD0[nEle]/F");
  tree_->Branch("eleDz", eleDz_, "eleDz[nEle]/F");
  tree_->Branch("eleD0GV", eleD0GV_, "eleD0GV[nEle]/F");
  tree_->Branch("eleDzGV", eleDzGV_, "eleDzGV[nEle]/F");
  tree_->Branch("eleD0Vtx", eleD0Vtx_, "eleD0Vtx[nEle][100]/F");
  tree_->Branch("eleDzVtx", eleDzVtx_, "eleDzVtx[nEle][100]/F");
  tree_->Branch("eleHoverE", eleHoverE_, "eleHoverE[nEle]/F");
  tree_->Branch("eleHoverE12", eleHoverE12_, "eleHoverE12[nEle]/F");
  tree_->Branch("eleEoverP", eleEoverP_, "eleEoverP[nEle]/F");
  tree_->Branch("elePin", elePin_, "elePin[nEle]/F");
  tree_->Branch("elePout", elePout_, "elePout[nEle]/F");
  tree_->Branch("eleTrkMomErr", eleTrkMomErr_, "eleTrkMomErr[nEle]/F");
  tree_->Branch("eleBrem", eleBrem_, "eleBrem[nEle]/F");
  tree_->Branch("eledEtaAtVtx", eledEtaAtVtx_, "eledEtaAtVtx[nEle]/F");
  tree_->Branch("eledPhiAtVtx", eledPhiAtVtx_, "eledPhiAtVtx[nEle]/F");
  tree_->Branch("eleSigmaIEtaIEta", eleSigmaIEtaIEta_, "eleSigmaIEtaIEta[nEle]/F");
  tree_->Branch("eleSigmaIEtaIPhi", eleSigmaIEtaIPhi_, "eleSigmaIEtaIPhi[nEle]/F");
  tree_->Branch("eleSigmaIPhiIPhi", eleSigmaIPhiIPhi_, "eleSigmaIPhiIPhi[nEle]/F");
  tree_->Branch("eleEmax", eleEmax_, "eleEmax[nEle]/F");
  tree_->Branch("eleE1x5", eleE1x5_, "eleE1x5[nEle]/F");
  tree_->Branch("eleE3x3", eleE3x3_, "eleE3x3[nEle]/F");
  tree_->Branch("eleE5x5", eleE5x5_, "eleE5x5[nEle]/F");
  tree_->Branch("eleE2x5Max", eleE2x5Max_, "eleE2x5Max[nEle]/F");
  tree_->Branch("eleRegrE", eleRegrE_, "eleRegrE[nEle]/F");
  tree_->Branch("eleRegrEerr", eleRegrEerr_, "eleRegrEerr[nEle]/F");
  tree_->Branch("elePhoRegrE", elePhoRegrE_, "elePhoRegrE[nEle]/F");
  tree_->Branch("elePhoRegrEerr", elePhoRegrEerr_, "elePhoRegrEerr[nEle]/F");
  tree_->Branch("eleSeedTime", eleSeedTime_, "eleSeedTime[nEle]/F");
  // If Flag == 2, it means that rechit is out of time
  // https://twiki.cern.ch/twiki/bin/viewauth/CMS/EcalFirstBeam09Anomalous#Spike_identification_in_collision
  tree_->Branch("eleRecoFlag", eleRecoFlag_, "eleRecoFlag[nEle]/I");
  tree_->Branch("elePos", elePos_, "elePos[nEle]/I");
  if (doGenParticles_) {
    tree_->Branch("eleGenIndex", eleGenIndex_, "eleGenIndex[nEle]/I");
    tree_->Branch("eleGenGMomPID", eleGenGMomPID_, "eleGenGMomPID[nEle]/I");
    tree_->Branch("eleGenMomPID", eleGenMomPID_, "eleGenMomPID[nEle]/I");
    tree_->Branch("eleGenMomPt", eleGenMomPt_, "eleGenMomPt[nEle]/F");
  }
  tree_->Branch("eleIsoTrkDR03", eleIsoTrkDR03_, "eleIsoTrkDR03[nEle]/F");
  tree_->Branch("eleIsoEcalDR03", eleIsoEcalDR03_, "eleIsoEcalDR03[nEle]/F");
  tree_->Branch("eleIsoHcalDR03", eleIsoHcalDR03_, "eleIsoHcalDR03[nEle]/F");
  tree_->Branch("eleIsoHcalDR0312", eleIsoHcalDR0312_, "eleIsoHcalDR0312[nEle]/F");
  tree_->Branch("eleIsoTrkDR04", eleIsoTrkDR04_, "eleIsoTrkDR04[nEle]/F");
  tree_->Branch("eleIsoEcalDR04", eleIsoEcalDR04_, "eleIsoEcalDR04[nEle]/F");
  tree_->Branch("eleIsoHcalDR04", eleIsoHcalDR04_, "eleIsoHcalDR04[nEle]/F");
  tree_->Branch("eleIsoHcalDR0412", eleIsoHcalDR0412_, "eleIsoHcalDR0412[nEle]/F");
  tree_->Branch("eleModIsoTrk", eleModIsoTrk_, "eleModIsoTrk[nEle]/F");
  tree_->Branch("eleModIsoEcal", eleModIsoEcal_, "eleModIsoEcal[nEle]/F");
  tree_->Branch("eleModIsoHcal", eleModIsoHcal_, "eleModIsoHcal[nEle]/F");
  tree_->Branch("eleMissHits", eleMissHits_, "eleMissHits[nEle]/I");
  tree_->Branch("eleConvDist", eleConvDist_, "eleConvDist[nEle]/F");
  tree_->Branch("eleConvDcot", eleConvDcot_, "eleConvDcot[nEle]/F");
  tree_->Branch("eleConvVtxFit", eleConvVtxFit_, "eleConvVtxFit[nEle]/I");
  tree_->Branch("eleIP3D", eleIP3D_, "eleIP3D[nEle]/F");
  tree_->Branch("eleIP3DErr", eleIP3DErr_, "eleIP3DErr[nEle]/F");
  tree_->Branch("eleIDMVANonTrig", eleIDMVANonTrig_, "eleIDMVANonTrig[nEle]/F");
  tree_->Branch("eleIDMVATrig", eleIDMVATrig_, "eleIDMVATrig[nEle]/F");
  tree_->Branch("elePFChIso03", elePFChIso03_, "elePFChIso03[nEle]/F");
  tree_->Branch("elePFPhoIso03", elePFPhoIso03_, "elePFPhoIso03[nEle]/F");
  tree_->Branch("elePFNeuIso03", elePFNeuIso03_, "elePFNeuIso03[nEle]/F");
  tree_->Branch("elePFChIso04", elePFChIso04_, "elePFChIso04[nEle]/F");
  tree_->Branch("elePFPhoIso04", elePFPhoIso04_, "elePFPhoIso04[nEle]/F");
  tree_->Branch("elePFNeuIso04", elePFNeuIso04_, "elePFNeuIso04[nEle]/F");
  tree_->Branch("eleESEffSigmaRR", eleESEffSigmaRR_, "eleESEffSigmaRR_[nEle][3]/F");
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
    tree_->Branch("eleNBC", eleNBC_, "eleNBC[nEle]/I");
    tree_->Branch("eleBrLinear"  , eleBrLinear_  , "eleBrLinear[nEle]/F");
    tree_->Branch("eleCetaCorrE" , eleCetaCorrE_ , "eleCetaCorrE[nEle]/F");
    tree_->Branch("eleCetaCorrEt", eleCetaCorrEt_, "eleCetaCorrEt[nEle]/F");
    tree_->Branch("eleBremCorrE" , eleBremCorrE_ , "eleBremCorrE[nEle]/F");
    tree_->Branch("eleBremCorrEt", eleBremCorrEt_, "eleBremCorrEt[nEle]/F");
    tree_->Branch("eleFullCorrE" , eleFullCorrE_ , "eleFullCorrE[nEle]/F");
    tree_->Branch("eleFullCorrEt", eleFullCorrEt_, "eleFullCorrEt[nEle]/F");  
  }
  // Photon
  tree_->Branch("nPho", &nPho_, "nPho/I");
  tree_->Branch("phoTrg", phoTrg_, "phoTrg[nPho][8]/I");
  tree_->Branch("phoTrgFilter", phoTrgFilter_, "phoTrgFilter[nPho][50]/I");
  tree_->Branch("phoIsPhoton", phoIsPhoton_, "phoIsPhoton[nPho]/O");
  tree_->Branch("phoSCPos", phoSCPos_, "phoSCPos[nPho][3]/F");
  tree_->Branch("phoCaloPos", phoCaloPos_, "phoCaloPos[nPho][3]/F");
  tree_->Branch("phoE", phoE_, "phoE[nPho]/F");
  tree_->Branch("phoEt", phoEt_, "phoEt[nPho]/F");
  tree_->Branch("phoEta", phoEta_, "phoEta[nPho]/F");
  tree_->Branch("phoVtx", phoVtx_, "phoVtx[nPho][3]/F");
  tree_->Branch("phoPhi", phoPhi_, "phoPhi[nPho]/F");
  tree_->Branch("phoEtVtx", phoEtVtx_, "phoEtVtx[nPho][100]/F");
  tree_->Branch("phoEtaVtx", phoEtaVtx_, "phoEtaVtx[nPho][100]/F");
  tree_->Branch("phoPhiVtx", phoPhiVtx_, "phoPhiVtx[nPho][100]/F");
  tree_->Branch("phoR9", phoR9_, "phoR9[nPho]/F");
  if (develop_) {
    tree_->Branch("phoCetaCorrE" , phoCetaCorrE_ , "phoCetaCorrE[nPho]/F");
    tree_->Branch("phoCetaCorrEt", phoCetaCorrEt_, "phoCetaCorrEt[nPho]/F");
    tree_->Branch("phoBremCorrE" , phoBremCorrE_ , "phoBremCorrE[nPho]/F");
    tree_->Branch("phoBremCorrEt", phoBremCorrEt_, "phoBremCorrEt[nPho]/F");
    tree_->Branch("phoFullCorrE" , phoFullCorrE_ , "phoFullCorrE[nPho]/F");
    tree_->Branch("phoFullCorrEt", phoFullCorrEt_, "phoFullCorrEt[nPho]/F");
  }
  //tree_->Branch("phoTrkIsoSolidDR03", phoTrkIsoSolidDR03_, "phoTrkIsoSolidDR03[nPho]/F");
  tree_->Branch("phoTrkIsoHollowDR03", phoTrkIsoHollowDR03_, "phoTrkIsoHollowDR03[nPho]/F");
  tree_->Branch("phoEcalIsoDR03", phoEcalIsoDR03_, "phoEcalIsoDR03[nPho]/F");
  tree_->Branch("phoHcalIsoDR03", phoHcalIsoDR03_, "phoHcalIsoDR03[nPho]/F");
  tree_->Branch("phoHcalIsoDR0312", phoHcalIsoDR0312_, "phoHcalIsoDR0312[nPho]/F");
  //tree_->Branch("phoHcalIsoSolidDR03", phoHcalIsoSolidDR03_, "phoHcalIsoSolidDR03[nPho]/F");
  //tree_->Branch("phoTrkIsoSolidDR04", phoTrkIsoSolidDR04_, "phoTrkIsoSolidDR04[nPho]/F");
  tree_->Branch("phoTrkIsoHollowDR04", phoTrkIsoHollowDR04_, "phoTrkIsoHollowDR04[nPho]/F");
  tree_->Branch("phoCiCTrkIsoDR03", phoCiCTrkIsoDR03_, "phoCiCTrkIsoDR03[nPho][100]/F");
  tree_->Branch("phoCiCTrkIsoDR04", phoCiCTrkIsoDR04_, "phoCiCTrkIsoDR04[nPho][100]/F");
  tree_->Branch("phoCiCdRtoTrk", phoCiCdRtoTrk_, "phoCiCdRtoTrk[nPho]/F");
  tree_->Branch("phoEcalIsoDR04", phoEcalIsoDR04_, "phoEcalIsoDR04[nPho]/F");
  tree_->Branch("phoHcalIsoDR04", phoHcalIsoDR04_, "phoHcalIsoDR04[nPho]/F");
  tree_->Branch("phoHcalIsoDR0412", phoHcalIsoDR0412_, "phoHcalIsoDR0412[nPho]/F");
  //tree_->Branch("phoHcalIsoSolidDR04", phoHcalIsoSolidDR04_, "phoHcalIsoSolidDR04[nPho]/F");
  tree_->Branch("phoHoverE", phoHoverE_, "phoHoverE[nPho]/F");
  tree_->Branch("phoHoverE12", phoHoverE12_, "phoHoverE12[nPho]/F");
  tree_->Branch("phoEleVeto", phoEleVeto_, "phoEleVeto[nPho]/I");
  tree_->Branch("phoSigmaIEtaIEta", phoSigmaIEtaIEta_, "phoSigmaIEtaIEta[nPho]/F");
  tree_->Branch("phoSigmaIEtaIPhi", phoSigmaIEtaIPhi_, "phoSigmaIEtaIPhi[nPho]/F");
  tree_->Branch("phoSigmaIPhiIPhi", phoSigmaIPhiIPhi_, "phoSigmaIPhiIPhi[nPho]/F");
  if (dumpESClusterInfo_) {
    tree_->Branch("phoCiCPF4phopfIso005", phoCiCPF4phopfIso005_, "phoCiCPF4phopfIso005[nPho]/F");
    tree_->Branch("phoCiCPF4phopfIso01", phoCiCPF4phopfIso01_, "phoCiCPF4phopfIso01[nPho]/F");
    tree_->Branch("phoCiCPF4phopfIso02", phoCiCPF4phopfIso02_, "phoCiCPF4phopfIso02[nPho]/F");
  }
  tree_->Branch("phoCiCPF4phopfIso03", phoCiCPF4phopfIso03_, "phoCiCPF4phopfIso03[nPho]/F");
  tree_->Branch("phoCiCPF4phopfIso04", phoCiCPF4phopfIso04_, "phoCiCPF4phopfIso04[nPho]/F");
  if (dumpESClusterInfo_) {
    tree_->Branch("phoCiCPF4phopfIso05", phoCiCPF4phopfIso05_, "phoCiCPF4phopfIso05[nPho]/F");
    tree_->Branch("phoCiCPF4phopfIso06", phoCiCPF4phopfIso06_, "phoCiCPF4phopfIso06[nPho]/F");
    tree_->Branch("phoCiCPF4phopfIso07", phoCiCPF4phopfIso07_, "phoCiCPF4phopfIso07[nPho]/F");
    tree_->Branch("phoCiCPF4phopfIso08", phoCiCPF4phopfIso08_, "phoCiCPF4phopfIso08[nPho]/F");
    tree_->Branch("phoCiCPF4chgpfIso005", phoCiCPF4chgpfIso005_, "phoCiCPF4chgpfIso005[nPho][100]/F");
    tree_->Branch("phoCiCPF4chgpfIso01", phoCiCPF4chgpfIso01_, "phoCiCPF4chgpfIso01[nPho][100]/F");
  }
  tree_->Branch("phoCiCPF4chgpfIso02", phoCiCPF4chgpfIso02_, "phoCiCPF4chgpfIso02[nPho][100]/F");
  tree_->Branch("phoCiCPF4chgpfIso03", phoCiCPF4chgpfIso03_, "phoCiCPF4chgpfIso03[nPho][100]/F");
  tree_->Branch("phoCiCPF4chgpfIso04", phoCiCPF4chgpfIso04_, "phoCiCPF4chgpfIso04[nPho][100]/F");
  if (dumpESClusterInfo_) {
    tree_->Branch("phoCiCPF4chgpfIso05", phoCiCPF4chgpfIso05_, "phoCiCPF4chgpfIso05[nPho][100]/F");
    tree_->Branch("phoCiCPF4chgpfIso06", phoCiCPF4chgpfIso06_, "phoCiCPF4chgpfIso06[nPho][100]/F");
    tree_->Branch("phoCiCPF4chgpfIso07", phoCiCPF4chgpfIso07_, "phoCiCPF4chgpfIso07[nPho][100]/F");
    tree_->Branch("phoCiCPF4chgpfIso08", phoCiCPF4chgpfIso08_, "phoCiCPF4chgpfIso08[nPho][100]/F");
    tree_->Branch("phoCiCPF4phopfIsoNoVETO005", phoCiCPF4phopfIsoNoVETO005_, "phoCiCPF4phopfIsoNoVETO005[nPho]/F");
    tree_->Branch("phoCiCPF4phopfIsoNoVETO01", phoCiCPF4phopfIsoNoVETO01_, "phoCiCPF4phopfIsoNoVETO01[nPho]/F");
    tree_->Branch("phoCiCPF4phopfIsoNoVETO02", phoCiCPF4phopfIsoNoVETO02_, "phoCiCPF4phopfIsoNoVETO02[nPho]/F");
    tree_->Branch("phoCiCPF4phopfIsoNoVETO03", phoCiCPF4phopfIsoNoVETO03_, "phoCiCPF4phopfIsoNoVETO03[nPho]/F");
    tree_->Branch("phoCiCPF4phopfIsoNoVETO04", phoCiCPF4phopfIsoNoVETO04_, "phoCiCPF4phopfIsoNoVETO04[nPho]/F");
    tree_->Branch("phoCiCPF4phopfIsoNoVETO05", phoCiCPF4phopfIsoNoVETO05_, "phoCiCPF4phopfIsoNoVETO05[nPho]/F");
    tree_->Branch("phoCiCPF4phopfIsoNoVETO06", phoCiCPF4phopfIsoNoVETO06_, "phoCiCPF4phopfIsoNoVETO06[nPho]/F");
    tree_->Branch("phoCiCPF4phopfIsoNoVETO07", phoCiCPF4phopfIsoNoVETO07_, "phoCiCPF4phopfIsoNoVETO07[nPho]/F");
    tree_->Branch("phoCiCPF4phopfIsoNoVETO08", phoCiCPF4phopfIsoNoVETO08_, "phoCiCPF4phopfIsoNoVETO08[nPho]/F");
    tree_->Branch("phoCiCPF4chgpfIsoNoVETO005", phoCiCPF4chgpfIsoNoVETO005_, "phoCiCPF4chgpfIsoNoVETO005[nPho][100]/F");
    tree_->Branch("phoCiCPF4chgpfIsoNoVETO01", phoCiCPF4chgpfIsoNoVETO01_, "phoCiCPF4chgpfIsoNoVETO01[nPho][100]/F");
    tree_->Branch("phoCiCPF4chgpfIsoNoVETO02", phoCiCPF4chgpfIsoNoVETO02_, "phoCiCPF4chgpfIsoNoVETO02[nPho][100]/F");
    tree_->Branch("phoCiCPF4chgpfIsoNoVETO03", phoCiCPF4chgpfIsoNoVETO03_, "phoCiCPF4chgpfIsoNoVETO03[nPho][100]/F");
    tree_->Branch("phoCiCPF4chgpfIsoNoVETO04", phoCiCPF4chgpfIsoNoVETO04_, "phoCiCPF4chgpfIsoNoVETO04[nPho][100]/F");
    tree_->Branch("phoCiCPF4chgpfIsoNoVETO05", phoCiCPF4chgpfIsoNoVETO05_, "phoCiCPF4chgpfIsoNoVETO05[nPho][100]/F");
    tree_->Branch("phoCiCPF4chgpfIsoNoVETO06", phoCiCPF4chgpfIsoNoVETO06_, "phoCiCPF4chgpfIsoNoVETO06[nPho][100]/F");
    tree_->Branch("phoCiCPF4chgpfIsoNoVETO07", phoCiCPF4chgpfIsoNoVETO07_, "phoCiCPF4chgpfIsoNoVETO07[nPho][100]/F");
    tree_->Branch("phoCiCPF4chgpfIsoNoVETO08", phoCiCPF4chgpfIsoNoVETO08_, "phoCiCPF4chgpfIsoNoVETO08[nPho][100]/F");
  }
  tree_->Branch("phoEmax", phoEmax_, "phoEmax[nPho]/F");
  tree_->Branch("phoE3x3", phoE3x3_, "phoE3x3[nPho]/F");
  tree_->Branch("phoE3x1", phoE3x1_, "phoE3x1[nPho]/F");
  tree_->Branch("phoE1x3", phoE1x3_, "phoE1x3[nPho]/F");
  tree_->Branch("phoE5x5", phoE5x5_, "phoE5x5[nPho]/F");
  tree_->Branch("phoE1x5", phoE1x5_, "phoE1x5[nPho]/F");
  tree_->Branch("phoE2x2", phoE2x2_, "phoE2x2[nPho]/F");
  tree_->Branch("phoE2x5Max", phoE2x5Max_, "phoE2x5Max[nPho]/F");
  tree_->Branch("phoPFChIso", phoPFChIso_, "phoPFChIso[nPho]/F");
  tree_->Branch("phoPFPhoIso", phoPFPhoIso_, "phoPFPhoIso[nPho]/F");
  tree_->Branch("phoPFNeuIso", phoPFNeuIso_, "phoPFNeuIso[nPho]/F");
  tree_->Branch("phoSCRChIso", phoSCRChIso_, "phoSCRChIso[nPho]/F");
  tree_->Branch("phoSCRPhoIso", phoSCRPhoIso_, "phoSCRPhoIso[nPho]/F");
  tree_->Branch("phoSCRNeuIso", phoSCRNeuIso_, "phoSCRNeuIso[nPho]/F");
  tree_->Branch("phoRegrE", phoRegrE_, "phoRegrE[nPho]/F");
  tree_->Branch("phoRegrEerr", phoRegrEerr_, "phoRegrEerr[nPho]/F");
  tree_->Branch("phoSeedTime", phoSeedTime_, "phoSeedTime[nPho]/F");
  tree_->Branch("phoSeedDetId1", phoSeedDetId1_, "phoSeedDetId1[nPho]/I");
  tree_->Branch("phoSeedDetId2", phoSeedDetId2_, "phoSeedDetId2[nPho]/I");
  tree_->Branch("phoLICTD", phoLICTD_, "phoLICTD[nPho]/F");
  // If Flag == 2, it means that rechit is out of time
  tree_->Branch("phoRecoFlag", phoRecoFlag_, "phoRecoFlag[nPho]/I");
  tree_->Branch("phoPos", phoPos_, "phoPos[nPho]/I");
  if (doGenParticles_) {
    tree_->Branch("phoGenIndex", phoGenIndex_, "phoGenIndex[nPho]/I");
    tree_->Branch("phoGenGMomPID", phoGenGMomPID, "phoGenGMomPID[nPho]/I");
    tree_->Branch("phoGenMomPID", phoGenMomPID, "phoGenMomPID[nPho]/I");
    tree_->Branch("phoGenMomPt", phoGenMomPt, "phoGenMomPt[nPho]/F");
  }
  tree_->Branch("phoSCE", phoSCE_, "phoSCE[nPho]/F");
  tree_->Branch("phoSCRawE", phoSCRawE_, "phoSCRawE[nPho]/F");
  tree_->Branch("phoESEn", phoESEn_, "phoESEn[nPho]/F");
  tree_->Branch("phoSCEt", phoSCEt_, "phoSCEt[nPho]/F");
  tree_->Branch("phoSCEta", phoSCEta_, "phoSCEta[nPho]/F");
  tree_->Branch("phoSCPhi", phoSCPhi_, "phoSCPhi[nPho]/F");
  tree_->Branch("phoSCEtaWidth", phoSCEtaWidth_, "phoSCEtaWidth[nPho]/F");
  tree_->Branch("phoSCPhiWidth", phoSCPhiWidth_, "phoSCPhiWidth[nPho]/F");
  tree_->Branch("phoSCBrem", phoSCBrem_, "phoSCBrem[nPho]/F");
  tree_->Branch("phoOverlap", phoOverlap_, "phoOverlap[nPho]/I");
  tree_->Branch("phohasPixelSeed", phohasPixelSeed_, "phohasPixelSeed[nPho]/I");
  tree_->Branch("pho_hasConvPf", pho_hasConvPf_, "pho_hasConvPf[nPho]/I");
  tree_->Branch("pho_hasSLConvPf", pho_hasSLConvPf_, "pho_hasSLConvPf[nPho]/I");
  if (develop_) {
    tree_->Branch("MustacheEin", MustacheEin_, "MustacheEin[nPho]/F");
    tree_->Branch("MustacheEOut", MustacheEOut_, "MustacheEOut[nPho]/F");
    tree_->Branch("MustacheEtOut", MustacheEtOut_, "MustacheEtOut[nPho]/F");
    tree_->Branch("PFRecoMatch", PFRecoMatch_, "PFRecoMatch[nPho]/I");
    tree_->Branch("PFEleMatch", PFEleMatch_, "PFEleMatch[nPho]/I");
    tree_->Branch("PFEleVeto", PFEleVeto_, "PFEleVeto[nPho]/I");    
    tree_->Branch("PFLowestClustE",PFLowestClustE_, "PFLowestClustE[nPho]/F");
    tree_->Branch("PFClustdEta",PFClustdEta_, "PFClustdEta[nPho]/F");
    tree_->Branch("PFClustdPhi",PFClustdPhi_, "PFClustdPhi[nPho]/F");
    tree_->Branch("PFClustRMSPhi",PFClustRMSPhi_, "PFClustRMSPhi[nPho]/F");
    tree_->Branch("PFClustRMSPhiMust",PFClustRMSPhiMust_, "PFClustRMSPhiMust[nPho]/F");    
    tree_->Branch("PFPreShowerE1", PFPreShowerE1_, "PFPreShowerE1[nPho]/F");
    tree_->Branch("PFPreShowerE2", PFPreShowerE2_, "PFPreShowerE2[nPho]/F");
  }
  tree_->Branch("pho_pfconvVtxZ", pho_pfconvVtxZ_, "pho_pfconvVtxZ[nPho]/F");
  tree_->Branch("pho_pfconvVtxZErr", pho_pfconvVtxZErr_, "pho_pfconvVtxZErr[nPho]/F");
  tree_->Branch("pho_nSLConv",pho_nSLConv_, "pho_nSLConv[nPho]/I");
  tree_->Branch("pho_pfSLConvPos", pho_pfSLConvPos_, "pho_pfSLConvPos[nPho]]20][3]/F");
  tree_->Branch("pho_pfSLConvVtxZ", pho_pfSLConvVtxZ_, "pho_pfSLConvVtxZ[nPho][20]/F");
  //Conversion Branches
  tree_->Branch("phoIsConv",               phoIsConv_,               "phoIsConv[nPho]/I");
  tree_->Branch("phoNConv",                phoNConv_,                "phoNConv[nPho]/I"); 
  tree_->Branch("phoConvInvMass",          phoConvInvMass_,          "phoConvInvMass[nPho]/F");
  tree_->Branch("phoConvCotTheta",         phoConvCotTheta_,         "phoConvCotTheta[nPho]/F");
  tree_->Branch("phoConvEoverP",           phoConvEoverP_,           "phoConvEoverP[nPho]/F");
  tree_->Branch("phoConvZofPVfromTrks",    phoConvZofPVfromTrks_,    "phoConvZofPVfromTrks[nPho] /F");
  tree_->Branch("phoConvMinDist",          phoConvMinDist_,          "phoConvMinDist[nPho]/F");
  tree_->Branch("phoConvdPhiAtVtx",        phoConvdPhiAtVtx_,        "phoConvdPhiAtVtx[nPho]/F");
  tree_->Branch("phoConvdPhiAtCalo",       phoConvdPhiAtCalo_,       "phoConvdPhiAtCalo[nPho]/F");
  tree_->Branch("phoConvdEtaAtCalo",       phoConvdEtaAtCalo_,       "phoConvdEtaAtCalo[nPho]/F");
  tree_->Branch("phoConvTrkd0",            phoConvTrkd0_,            "phoConvTrkd0[nPho][2]/F");
  tree_->Branch("phoConvTrkPin",           phoConvTrkPin_,           "phoConvTrkPin[nPho][2]/F");
  tree_->Branch("phoConvTrkPout",          phoConvTrkPout_,          "phoConvTrkPout[nPho][2]/F");
  tree_->Branch("phoConvTrkdz",            phoConvTrkdz_,            "phoConvTrkdz[nPho][2]/F");
  tree_->Branch("phoConvTrkdzErr",         phoConvTrkdzErr_,         "phoConvTrkdzErr[nPho][2]/F");
  tree_->Branch("phoConvChi2",             phoConvChi2_,             "phoConvChi2[nPho]/F");
  tree_->Branch("phoConvChi2Prob",         phoConvChi2Prob_,         "phoConvChi2Prob[nPho]/F");
  tree_->Branch("phoConvNTrks",            phoConvNTrks_,            "phoConvNTrks[nPho]/I");
  tree_->Branch("phoConvCharge",           phoConvCharge_,           "phoConvCharge[nPho][2]/F");
  tree_->Branch("phoConvValidVtx",         phoConvValidVtx_,         "phoConvValidVtx[nPho]/F");
  tree_->Branch("phoConvLikeLihood",       phoConvLikeLihood_,       "phoConvLikeLihood[nPho]/F");
  tree_->Branch("phoConvP4",               phoConvP4_,               "phoConvP4[nPho][4]/F");
  tree_->Branch("phoConvVtx",              phoConvVtx_,              "phoConvVtx[nPho][3]/F");
  tree_->Branch("phoConvVtxErr",           phoConvVtxErr_,           "phoConvVtxErr[nPho][3]/F");
  tree_->Branch("phoConvPairMomentum",     phoConvPairMomentum_,     "phoConvPairMomentum[nPho][3]/F");
  tree_->Branch("phoConvRefittedMomentum", phoConvRefittedMomentum_, "phoConvRefittedMomentum[nPho][3]/F");
  tree_->Branch("SingleLegConv", SingleLegConv_,"SingleLegConv[nPho]/I");
  tree_->Branch("phoPFConvVtx", phoPFConvVtx_,"phoPFConvVtx[nPho][3]/F");
  tree_->Branch("phoPFConvMom", phoPFConvMom_,"phoPFConvMom[nPho][3]/F");
  tree_->Branch("phoESEffSigmaRR", phoESEffSigmaRR_, "phoESEffSigmaRR_[nPho][3]/F");
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
  tree_->Branch("muTrg", muTrg_, "muTrg[nMu][10]/I");
  tree_->Branch("muEta", muEta_, "muEta[nMu]/F");
  tree_->Branch("muPhi", muPhi_, "muPhi[nMu]/F");
  tree_->Branch("muCharge", muCharge_, "muCharge[nMu]/I");
  tree_->Branch("muPt", muPt_, "muPt[nMu]/F");
  tree_->Branch("muPz", muPz_, "muPz[nMu]/F");
  tree_->Branch("muVtx", muVtx_, "muVtx[nMu][3]/F");
  tree_->Branch("muVtxGlb", muVtxGlb_, "muVtxGlb[nMu][3]/F");
  if (doGenParticles_)
    tree_->Branch("muGenIndex", muGenIndex_, "muGenIndex[nMu]/I");
  tree_->Branch("mucktPt", mucktPt_, "mucktPt[nMu]/F");
  tree_->Branch("mucktPtErr", mucktPtErr_, "mucktPtErr[nMu]/F");
  tree_->Branch("mucktEta", mucktEta_, "mucktEta[nMu]/F");
  tree_->Branch("mucktPhi", mucktPhi_, "mucktPhi[nMu]/F");
  tree_->Branch("mucktdxy", mucktdxy_, "mucktdxy[nMu]/F");
  tree_->Branch("mucktdz", mucktdz_, "mucktdz[nMu]/F");
  tree_->Branch("muIsoTrk", muIsoTrk_, "muIsoTrk[nMu]/F");
  tree_->Branch("muIsoCalo", muIsoCalo_, "muIsoCalo[nMu]/F");
  tree_->Branch("muIsoEcal", muIsoEcal_, "muIsoEcal[nMu]/F");
  tree_->Branch("muIsoHcal", muIsoHcal_, "muIsoHcal[nMu]/F");
  tree_->Branch("muChi2NDF", muChi2NDF_, "muChi2NDF[nMu]/F");
  tree_->Branch("muInnerChi2NDF", muInnerChi2NDF_, "muInnerChi2NDF[nMu]/F");
  tree_->Branch("muPFIsoR04_CH",muPFIsoR04_CH_,"muPFIsoR04_CH[nMu]/F");
  tree_->Branch("muPFIsoR04_NH",muPFIsoR04_NH_,"muPFIsoR04_NH[nMu]/F");
  tree_->Branch("muPFIsoR04_Pho",muPFIsoR04_Pho_,"muPFIsoR04_Pho[nMu]/F");
  tree_->Branch("muPFIsoR04_PU",muPFIsoR04_PU_,"muPFIsoR04_PU[nMu]/F");
  tree_->Branch("muPFIsoR04_CPart",muPFIsoR04_CPart_,"muPFIsoR04_CPart[nMu]/F");
  tree_->Branch("muPFIsoR04_NHHT",muPFIsoR04_NHHT_,"muPFIsoR04_NHHT[nMu]/F");
  tree_->Branch("muPFIsoR04_PhoHT",muPFIsoR04_PhoHT_,"muPFIsoR04_PhoHT[nMu]/F");
  tree_->Branch("muPFIsoR03_CH",muPFIsoR03_CH_,"muPFIsoR03_CH[nMu]/F");
  tree_->Branch("muPFIsoR03_NH",muPFIsoR03_NH_,"muPFIsoR03_NH[nMu]/F");
  tree_->Branch("muPFIsoR03_Pho",muPFIsoR03_Pho_,"muPFIsoR03_Pho[nMu]/F");
  tree_->Branch("muPFIsoR03_PU",muPFIsoR03_PU_,"muPFIsoR03_PU[nMu]/F");
  tree_->Branch("muPFIsoR03_CPart",muPFIsoR03_CPart_,"muPFIsoR03_CPart[nMu]/F");
  tree_->Branch("muPFIsoR03_NHHT",muPFIsoR03_NHHT_,"muPFIsoR03_NHHT[nMu]/F");
  tree_->Branch("muPFIsoR03_PhoHT",muPFIsoR03_PhoHT_,"muPFIsoR03_PhoHT[nMu]/F");
  tree_->Branch("muType", muType_, "muType[nMu]/I");
  tree_->Branch("muD0", muD0_, "muD0[nMu]/F");
  tree_->Branch("muDz", muDz_, "muDz[nMu]/F");
  tree_->Branch("muD0GV", muD0GV_, "muD0GV[nMu]/F");
  tree_->Branch("muDzGV", muDzGV_, "muDzGV[nMu]/F");
  tree_->Branch("muD0Vtx", muD0Vtx_, "muD0Vtx[nMu][100]/F");
  tree_->Branch("muDzVtx", muDzVtx_, "muDzVtx[nMu][100]/F");   
  tree_->Branch("muInnerD0", muInnerD0_, "muInnerD0[nMu]/F");
  tree_->Branch("muInnerDz", muInnerDz_, "muInnerDz[nMu]/F");
  tree_->Branch("muInnerD0GV", muInnerD0GV_, "muInnerD0GV[nMu]/F");
  tree_->Branch("muInnerDzGV", muInnerDzGV_, "muInnerDzGV[nMu]/F");
  tree_->Branch("muInnerPt", muInnerPt_, "muInnerPt[nMu]/F");
  tree_->Branch("muInnerPtErr", muInnerPtErr_, "muInnerPtErr[nMu]/F");
  tree_->Branch("muNumberOfValidTrkLayers", muNumberOfValidTrkLayers_, "muNumberOfValidTrkLayers[nMu]/I");
  tree_->Branch("muNumberOfValidTrkHits", muNumberOfValidTrkHits_, "muNumberOfValidTrkHits[nMu]/I");
  tree_->Branch("muNumberOfValidPixelLayers", muNumberOfValidPixelLayers_, "muNumberOfValidPixelLayers[nMu]/I");
  tree_->Branch("muNumberOfValidPixelHits", muNumberOfValidPixelHits_, "muNumberOfValidPixelHits[nMu]/I");
  tree_->Branch("muNumberOfValidMuonHits", muNumberOfValidMuonHits_, "muNumberOfValidMuonHits[nMu]/I");
  tree_->Branch("muStations", muStations_, "muStations[nMu]/I");
  tree_->Branch("muChambers", muChambers_, "muChambers[nMu]/I");
  tree_->Branch("muIP3D", muIP3D_, "muIP3D[nMu]/F");
  tree_->Branch("muIP3DErr", muIP3DErr_, "muIP3DErr[nMu]/F");
  tree_->Branch("nPFPho", &nPFPho_, "nPFPho_/I");    
  tree_->Branch("PFPhoEt", PFPhoEt_, "PFPhoEt_[nPFPho_]/F");  
  tree_->Branch("PFPhoEta", PFPhoEta_, "PFPhoEta_[nPFPho_]/F");  
  tree_->Branch("PFPhoPhi", PFPhoPhi_, "PFPhoPhi_[nPFPho_]/F");    
  tree_->Branch("PFPhoType", PFPhoType_, "PFPhoType_[nPFPho_]/I");
  tree_->Branch("PFPhoIso", PFPhoIso_, "PFPhoIso_[nPFPho_]/F");
  if (develop_) {
    tree_->Branch("PFPhoE", PFPhoE_, "PFPhoE_[nPFPho_]/F");  
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
  // Jet
  if (dumpJets_) {
    tree_->Branch("nJet", &nJet_, "nJet/I");
    tree_->Branch("jetTrg", jetTrg_, "jetTrg[nJet][14]/I");
    tree_->Branch("jetEn", jetEn_, "jetEn[nJet]/F");
    tree_->Branch("jetPt", jetPt_, "jetPt[nJet]/F");
    tree_->Branch("jetEta", jetEta_, "jetEta[nJet]/F");
    tree_->Branch("jetPhi", jetPhi_, "jetPhi[nJet]/F");
    tree_->Branch("jetCharge", jetCharge_, "jetCharge[nJet]/F");
    tree_->Branch("jetEt", jetEt_, "jetEt[nJet]/F");
    tree_->Branch("jetRawPt", jetRawPt_, "jetRawPt[nJet]/F");
    tree_->Branch("jetRawEn", jetRawEn_, "jetRawEn[nJet]/F");
    tree_->Branch("jetArea", jetArea_, "jetArea[nJet]/F");
    tree_->Branch("jetCHF", jetCHF_, "jetCHF[nJet]/F");
    tree_->Branch("jetNHF", jetNHF_, "jetNHF[nJet]/F");
    tree_->Branch("jetCEF", jetCEF_, "jetCEF[nJet]/F");
    tree_->Branch("jetNEF", jetNEF_, "jetNEF[nJet]/F");
    tree_->Branch("jetNCH", jetNCH_, "jetNCH[nJet]/I");
    tree_->Branch("jetHFHAE", jetHFHAE_, "jetHFHAE[nJet]/F");
    tree_->Branch("jetHFEME", jetHFEME_, "jetHFEME[nJet]/F");
    tree_->Branch("jetNConstituents", jetNConstituents_, "jetNConstituents[nJet]/I");
    tree_->Branch("jetCombinedSecondaryVtxBJetTags", jetCombinedSecondaryVtxBJetTags_, "jetCombinedSecondaryVtxBJetTags[nJet]/F");  // rob !!!
    tree_->Branch("jetCombinedSecondaryVtxMVABJetTags", jetCombinedSecondaryVtxMVABJetTags_, "jetCombinedSecondaryVtxMVABJetTags[nJet]/F");  // rob !!!
    tree_->Branch("jetJetProbabilityBJetTags", jetJetProbabilityBJetTags_, "jetJetProbabilityBJetTags[nJet]/F");  // rob !!!
    tree_->Branch("jetJetBProbabilityBJetTags", jetJetBProbabilityBJetTags_, "jetJetBProbabilityBJetTags[nJet]/F");  // rob !!!
    tree_->Branch("jetTrackCountingHighPurBJetTags", jetTrackCountingHighPurBJetTags_, "jetTrackCountingHighPurBJetTags[nJet]/F");  // rob !!!
    tree_->Branch("jetBetaStar", jetBetaStar_, "jetBetaStar[nJet][100]/F");
    // CMG Jet Id Variables
    tree_->Branch("jetPFLooseId", jetPFLooseId_, "jetPFLooseId[nJet]/O");
    tree_->Branch("jetDRMean", jetDRMean_, "jetDRMean[nJet]/F");
    tree_->Branch("jetDR2Mean", jetDR2Mean_, "jetDR2Mean[nJet]/F");
    tree_->Branch("jetDZ", jetDZ_, "jetDZ[nJet]/F");
    tree_->Branch("jetFrac01", jetFrac01_, "jetFrac01[nJet]/F");
    tree_->Branch("jetFrac02", jetFrac02_, "jetFrac02[nJet]/F");
    tree_->Branch("jetFrac03", jetFrac03_, "jetFrac03[nJet]/F");
    tree_->Branch("jetFrac04", jetFrac04_, "jetFrac04[nJet]/F");
    tree_->Branch("jetFrac05", jetFrac05_, "jetFrac05[nJet]/F");
    tree_->Branch("jetFrac06", jetFrac06_, "jetFrac06[nJet]/F");
    tree_->Branch("jetFrac07", jetFrac07_, "jetFrac07[nJet]/F");
    tree_->Branch("jetBeta", jetBeta_, "jetBeta[nJet]/F");
    tree_->Branch("jetBetaStarCMG", jetBetaStarCMG_, "jetBetaStarCMG[nJet]/F");
    tree_->Branch("jetBetaStarClassic", jetBetaStarClassic_, "jetBetaStarClassic[nJet]/F");
    tree_->Branch("jetBetaExt", jetBetaExt_, "jetBetaExt[nJet][100]/F");
    tree_->Branch("jetBetaStarCMGExt", jetBetaStarCMGExt_, "jetBetaStarCMGExt[nJet][100]/F");
    tree_->Branch("jetBetaStarClassicExt", jetBetaStarClassicExt_, "jetBetaStarClassicExt[nJet][100]/F");
    tree_->Branch("jetNNeutrals", jetNNeutrals_, "jetNNeutrals[nJet]/F");
    tree_->Branch("jetNCharged", jetNCharged_, "jetNCharged[nJet]/F");
    tree_->Branch("jetMVAs", jetMVAs_, "jetMVAs[nJet][4]/F");             // 0: simple, 1: full, 2: cut based, 3: philv1
    tree_->Branch("jetWPLevels", jetWPLevels_, "jetWPLevels[nJet][4]/I"); // 0: simple, 1: full, 2: cut based, 3: philv1
    tree_->Branch("jetMVAsExt", jetMVAsExt_, "jetMVAsExt[nJet][4][100]/F");             // 0: simple, 1: full, 2: cut based, 3: philv1; vtxbs
    tree_->Branch("jetWPLevelsExt", jetWPLevelsExt_, "jetWPLevelsExt[nJet][4][100]/I"); // 0: simple, 1: full, 2: cut based, 3: philv1; vtxbs
    //b-jet regression variables
    tree_->Branch("jetMt", jetMt_, "jetMt[nJet]/F");
    tree_->Branch("jetJECUnc", jetJECUnc_, "jetJECUnc[nJet]/F");
    tree_->Branch("jetLeadTrackPt", jetLeadTrackPt_, "jetLeadTrackPt[nJet]/F");
    tree_->Branch("jetVtxPt", jetVtxPt_, "jetVtxPt[nJet]/F");
    tree_->Branch("jetVtxMass", jetVtxMass_, "jetVtxMass[nJet]/F");
    tree_->Branch("jetVtx3dL", jetVtx3dL_, "jetVtx3dL[nJet]/F");
    tree_->Branch("jetVtx3deL", jetVtx3deL_, "jetVtx3deL[nJet]/F");
    tree_->Branch("jetSoftLeptPt", jetSoftLeptPt_, "jetSoftLeptPt[nJet]/F");
    tree_->Branch("jetSoftLeptPtRel", jetSoftLeptPtRel_, "jetSoftLeptPtRel[nJet]/F");
    tree_->Branch("jetSoftLeptdR", jetSoftLeptdR_, "jetSoftLeptdR[nJet]/F");
    tree_->Branch("jetSoftLeptIdlooseMu", jetSoftLeptIdlooseMu_, "jetSoftLeptIdlooseMu[nJet]/F");
    tree_->Branch("jetSoftLeptIdEle95", jetSoftLeptIdEle95_, "jetSoftLeptIdEle95[nJet]/F");
    tree_->Branch("jetDPhiMETJet", jetDPhiMETJet_, "jetDPhiMETJet[nJet]/F");
    tree_->Branch("jetPuJetIdL", jetPuJetIdL_, "jetPuJetIdL[nJet]/F");
    tree_->Branch("jetPuJetIdM", jetPuJetIdM_, "jetPuJetIdM[nJet]/F");
    tree_->Branch("jetPuJetIdT", jetPuJetIdT_, "jetPuJetIdT[nJet]/F");
    if (doGenParticles_) {
      tree_->Branch("jetPartonID", jetPartonID_, "jetPartonID[nJet]/I");
      tree_->Branch("jetGenJetIndex", jetGenJetIndex_, "jetGenJetIndex[nJet]/I");
      tree_->Branch("jetGenJetEn", jetGenJetEn_, "jetGenJetEn[nJet]/F");
      tree_->Branch("jetGenJetPt", jetGenJetPt_, "jetGenJetPt[nJet]/F");
      tree_->Branch("jetGenJetEta", jetGenJetEta_, "jetGenJetEta[nJet]/F");
      tree_->Branch("jetGenJetPhi", jetGenJetPhi_, "jetGenJetPhi[nJet]/F");
      tree_->Branch("jetGenPartonID", jetGenPartonID_, "jetGenPartonID[nJet]/I");
      tree_->Branch("jetGenEn", jetGenEn_, "jetGenEn[nJet]/F");
      tree_->Branch("jetGenPt", jetGenPt_, "jetGenPt[nJet]/F");
      tree_->Branch("jetGenEta", jetGenEta_, "jetGenEta[nJet]/F");
      tree_->Branch("jetGenPhi", jetGenPhi_, "jetGenPhi[nJet]/F");
    }
    // Low Pt Jets
    if (dumpTrks_) {
      tree_->Branch("nLowPtJet", &nLowPtJet_, "nLowPtJet/I");
      tree_->Branch("jetLowPtEn", jetLowPtEn_, "jetLowPtEn[nLowPtJet]/F");
      tree_->Branch("jetLowPtPt", jetLowPtPt_, "jetLowPtPt[nLowPtJet]/F");
      tree_->Branch("jetLowPtEta", jetLowPtEta_, "jetLowPtEta[nLowPtJet]/F");
      tree_->Branch("jetLowPtPhi", jetLowPtPhi_, "jetLowPtPhi[nLowPtJet]/F");
      tree_->Branch("jetLowPtCharge", jetLowPtCharge_, "jetLowPtCharge[nLowPtJet]/F");
      tree_->Branch("jetLowPtEt", jetLowPtEt_, "jetLowPtEt[nLowPtJet]/F");
      tree_->Branch("jetLowPtRawPt", jetLowPtRawPt_, "jetLowPtRawPt[nLowPtJet]/F");
      tree_->Branch("jetLowPtRawEn", jetLowPtRawEn_, "jetLowPtRawEn[nLowPtJet]/F");
      tree_->Branch("jetLowPtArea", jetLowPtArea_, "jetLowPtArea[nLowPtJet]/F");
      if (doGenParticles_) {
	tree_->Branch("jetLowPtPartonID", jetLowPtPartonID_, "jetLowPtPartonID[nLowPtJet]/I");
	tree_->Branch("jetLowPtGenJetEn", jetLowPtGenJetEn_, "jetLowPtGenJetEn[nLowPtJet]/F");
	tree_->Branch("jetLowPtGenJetPt", jetLowPtGenJetPt_, "jetLowPtGenJetPt[nLowPtJet]/F");
	tree_->Branch("jetLowPtGenJetEta", jetLowPtGenJetEta_, "jetLowPtGenJetEta[nLowPtJet]/F");
	tree_->Branch("jetLowPtGenJetPhi", jetLowPtGenJetPhi_, "jetLowPtGenJetPhi[nLowPtJet]/F");
	tree_->Branch("jetLowPtGenPartonID", jetLowPtGenPartonID_, "jetLowPtGenPartonID[nLowPtJet]/I");
	tree_->Branch("jetLowPtGenEn", jetLowPtGenEn_, "jetLowPtGenEn[nLowPtJet]/F");
	tree_->Branch("jetLowPtGenPt", jetLowPtGenPt_, "jetLowPtGenPt[nLowPtJet]/F");
	tree_->Branch("jetLowPtGenEta", jetLowPtGenEta_, "jetLowPtGenEta[nLowPtJet]/F");
	tree_->Branch("jetLowPtGenPhi", jetLowPtGenPhi_, "jetLowPtGenPhi[nLowPtJet]/F");
      }
    }
  }
  // Converted Photon Collection
  tree_->Branch("nConv", &nConv_, "nConv/I");
  tree_->Branch("convP4", convP4_, "convP4[nConv][4]/F");
  tree_->Branch("convVtx", convVtx_, "convVtx[nConv][3]/F");
  tree_->Branch("convVtxErr", convVtxErr_, "convVtxErr[nConv][3]/F");
  tree_->Branch("convPairMomentum", convPairMomentum_, "convPairMomentum[nConv][3]/F");
  tree_->Branch("convRefittedMomentum", convRefittedMomentum_, "convRefittedMomentum[nConv][3]/F");
  tree_->Branch("convNTracks", convNTracks_, "convNTracks[nConv]/I");
  tree_->Branch("convPairInvMass", convPairInvMass_, "convPairInvMass[nConv]/F");
  tree_->Branch("convPairCotThetaSep", convPairCotThetaSep_,"convPairCotThetaSep[nConv]/F");
  tree_->Branch("convEoverP", convEoverP_, "convEoverP[nConv]/F");
  tree_->Branch("convDistOfMinApproach", convDistOfMinApproach_, "convDistOfMinApproach[nConv]/F");
  tree_->Branch("convDPhiTrksAtVtx", convDPhiTrksAtVtx_, "convDphiTrksAtVtx[nConv]/F");
  tree_->Branch("convDPhiTrksAtEcal", convDPhiTrksAtEcal_, "convDPhiTrksAtEcal[nConv]/F");
  tree_->Branch("convDEtaTrksAtEcal", convDEtaTrksAtEcal_, "convDEtaTrksAtEcal[nConv]/F");
  tree_->Branch("convDxy", convDxy_, "convDxy[nConv]/F"); 
  tree_->Branch("convDz", convDz_, "convDz[nConv]/F");    
  tree_->Branch("convLxy", convLxy_, "convLxy[nConv]/F"); 
  tree_->Branch("convLz", convLz_, "convLz[nConv]/F");    
  tree_->Branch("convZofPrimVtxFromTrks", convZofPrimVtxFromTrks_, "convZofPrimVtxFromTrks[nConv]/F");
  tree_->Branch("convNHitsBeforeVtx", convNHitsBeforeVtx_, "convNHitsBeforeVtx[nConv][2]/I");
  tree_->Branch("convNSharedHits", convNSharedHits_, "convNSharedHits[nConv]/I");
  tree_->Branch("convValidVtx", convValidVtx_, "convValidVtx[nConv]/I");
  tree_->Branch("convMVALikelihood", convMVALikelihood_, "convMVALikelihood[nConv]/F");
  tree_->Branch("convChi2", convChi2_, "convChi2[nConv]/F");
  tree_->Branch("convChi2Probability", convChi2Probability_, "convChi2Probability[nConv]/F");
  // per track quantities
  tree_->Branch("convTk1Dz",    convTk1Dz_,    "convTk1Dz[nConv]/F");
  tree_->Branch("convTk2Dz",    convTk2Dz_,    "convTk2Dz[nConv]/F");
  tree_->Branch("convTk1DzErr", convTk1DzErr_, "convTk1DzErr[nConv]/F");
  tree_->Branch("convTk2DzErr", convTk2DzErr_, "convTk2DzErr[nConv]/F");
  //tree_->Branch("convTk1Nh",    convTk1Nh_,    "convTk1Nh[nConv]/I");
  //tree_->Branch("convTk2Nh",    convTk2Nh_,    "convTk2Nh[nConv]/I");
  tree_->Branch("convCh1Ch2",   convCh1Ch2_,   "convCh1Ch2[nConv]/I");
  tree_->Branch("convTk1D0",    convTk1D0_,    "convTk1D0[nConv]/F");
  tree_->Branch("convTk1Pout",  convTk1Pout_,  "convTk1Pout[nConv]/F");
  tree_->Branch("convTk1Pin",   convTk1Pin_,   "convTk1Pin[nConv]/F");
  tree_->Branch("convTk2D0",    convTk2D0_,    "convTk2D0[nConv]/F");
  tree_->Branch("convTk2Pout",  convTk2Pout_,  "convTk2Pout[nConv]/F");
  tree_->Branch("convTk2Pin",   convTk2Pin_,   "convTk2Pin[nConv]/F");

  vtxbsTkIndex_  = new std::vector<std::vector<Int_t> >;   vtxbsTkIndex_->clear();
  vtxbsTkWeight_ = new std::vector<std::vector<Float_t> >; vtxbsTkWeight_->clear();

}

ggNtuplizer::~ggNtuplizer() {

  delete cicPhotonId_;
  delete trackMET_;
  delete vtxbsTkIndex_;
  delete vtxbsTkWeight_;

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
			     edm::Handle<reco::PFMETCollection> &       recoPfMETHandle,
			     edm::Handle<edm::View<pat::Electron> > &   electronHandle,
			     edm::Handle<edm::View<pat::Photon> > &     photonHandle,
			     edm::Handle<reco::PhotonCollection> &      recoPhotonHandle,
			     edm::Handle<edm::View<pat::Muon> > &       muonHandle,
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
  event.getByLabel(InputTag("offlinePrimaryVerticesWithBS",""), recVtxsBS);
  event.getByLabel(trgResults_            , trgResultsHandle_);
  event.getByLabel(trgEvent_              , triggerEvent);
  event.getByLabel(tracklabel_            , tracksHandle);
  event.getByLabel(gsfElectronlabel_      , gsfElectronHandle);
  event.getByLabel(METCollection_         , METHandle);
  event.getByLabel(pfMETlabel_            , pfMETHandle);
  event.getByLabel(recoPfMETlabel_        , recoPfMETHandle);
  event.getByLabel(electronCollection_    , electronHandle);
  event.getByLabel(photonCollection_      , photonHandle);
  event.getByLabel(recophotonCollection_  , recoPhotonHandle);
  event.getByLabel(muonCollection_        , muonHandle);
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
		   recoPfMETHandle_,
		   electronHandle_,
		   photonHandle_,
		   recoPhotonHandle_,
		   muonHandle_,
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
  SuperClusterFootprintRemoval scRemover(e, edm::ParameterSet(), es);

  // rho collection
  Handle<double> rhoHandle25;
  e.getByLabel(rhoCollection25_, rhoHandle25);
  rho25_ = *(rhoHandle25.product());

  Handle<double> rhoHandle25_neu;
  e.getByLabel(rhoCollection25_neu_, rhoHandle25_neu);
  rho25_neu_ = *(rhoHandle25_neu.product()) ;

  Handle<double> rhoHandle_elePFiso;
  e.getByLabel(edm::InputTag("kt6PFJetsCentralChargedPileUp","rho"), rhoHandle_elePFiso);
  rho25_elePFiso_ = *(rhoHandle_elePFiso.product());

  Handle<double> rhoHandle_muPFiso;
  e.getByLabel(rhoLepPFisoCollection_ ,rhoHandle_muPFiso);
  rho25_muPFiso_  = *(rhoHandle_muPFiso.product());

  edm::Handle<double> rho2011Handle;
  e.getByLabel(rho2011Label_, rho2011Handle);
  rho2011_ = *(rho2011Handle.product());

  Handle<double> rhoHandle_enReg;
  e.getByLabel(edm::InputTag("kt6PFJets","rho"), rhoHandle_enReg); 
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
	vtx_[nVtx_][0] = (*recVtxs_)[i].x();
	vtx_[nVtx_][1] = (*recVtxs_)[i].y();
	vtx_[nVtx_][2] = (*recVtxs_)[i].z();
	vtxNTrk_[nVtx_] = (*recVtxs_)[i].tracksSize();
	vtxNDF_[nVtx_] = (*recVtxs_)[i].ndof();
	vtxD0_[nVtx_] = (*recVtxs_)[i].position().rho();
	
	if (vtxNDF_[nVtx_] > 4 && fabs(vtx_[nVtx_][2]) <= 24 && fabs(vtxD0_[nVtx_]) <= 2) {
	  nGoodVtx_++;
	  if (nGoodVtx_ == 1) firstGoodVtx = i;
	}
	nVtx_++;
      }
    }
  }
  if (nGoodVtx_ > 0) IsVtxGood_ = firstGoodVtx;

  // Set PV and first good vertex
  math::XYZPoint pv(vtx_[0][0], vtx_[0][1], vtx_[0][2]);
  if (firstGoodVtx < 0) firstGoodVtx = 0; 
  math::XYZPoint gv(vtx_[firstGoodVtx][0], vtx_[firstGoodVtx][1], vtx_[firstGoodVtx][2]);

  // vertex with beamspot
  vtxbsTkIndex_->clear();
  vtxbsTkWeight_->clear();

  nVtxBS_ = 0;
  if (recVtxsBS_.isValid()) {
    for (size_t i=0; i<recVtxsBS_->size(); ++i) {
      
      reco::VertexRef vtxRef(recVtxsBS_, i);
      
      vtxbsPtMod_[nVtxBS_] = 0.;
      vtxbsSumPt2_[nVtxBS_] = 0.;
      
      vtxbs_[nVtxBS_][0] = (*recVtxsBS_)[i].x();
      vtxbs_[nVtxBS_][1] = (*recVtxsBS_)[i].y();
      vtxbs_[nVtxBS_][2] = (*recVtxsBS_)[i].z();
      
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
	    vtxbsSumPt2_[nVtxBS_] += tkPtVec.Mod2();

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

      vtxbsPtMod_[nVtxBS_] = vtxbsP_.XYvector().Mod();
      
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
	trkP_[nTrk_][0]    = aTrk->px(); 
	trkP_[nTrk_][1]    = aTrk->py(); 
	trkP_[nTrk_][2]    = aTrk->pz(); 
	trkVtx_[nTrk_][0]  = aTrk->vx();
	trkVtx_[nTrk_][1]  = aTrk->vy();
	trkVtx_[nTrk_][2]  = aTrk->vz();
	trkd0_[nTrk_]      = aTrk->d0();
	trkd0Err_[nTrk_]   = aTrk->d0Error();
	trkdz_[nTrk_]      = aTrk->dz();
	trkdzErr_[nTrk_]   = aTrk->dzError();
	trkPtErr_[nTrk_]   = aTrk->ptError();
	trkQuality_[nTrk_] = aTrk->qualityMask();
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

    for (int i=0; i<7; ++i) pdf_[i] = -99;
    pthat_ = -99;
    processID_ = -99;

    Handle<GenEventInfoProduct> pdfInfoHandle;
    if (e.getByLabel("generator", pdfInfoHandle)) {
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
    if (e.getByLabel("generator", genEventScale)) {
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
	nPU_[nPUInfo_]    = PVI->getPU_NumInteractions();
	puTrue_[nPUInfo_] = PVI->getTrueNumInteractions();
	puBX_[nPUInfo_]   = PVI->getBunchCrossing();
	nPUInfo_++;
      }
    }
  }

  // GenParticle
  if (!isData_ && genParticlesHandle_.isValid() ) {
    
    nMC_ = 0;
    int genIndex = 0;
    
    for (vector<GenParticle>::const_iterator ip = genParticlesHandle_->begin(); ip != genParticlesHandle_->end(); ++ip) {
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
	(ip->pdgId() == 23 || abs(ip->pdgId()) == 24 || ip->pdgId() == 25 || abs(ip->pdgId()) == 6 || abs(ip->pdgId()) == 5);
      
      if ( stableFinalStateParticle || heavyParticle || photonOrLepton ) {
	const Candidate *p = (const Candidate*)&(*ip);
	if (!runOnParticleGun_ && !p->mother()) continue;
	
	reco::GenParticleRef partRef = reco::GenParticleRef(genParticlesHandle_,ip-genParticlesHandle_->begin());
	genpartparentage::GenParticleParentage particleHistory(partRef);
	
	mcPID[nMC_]       = p->pdgId();
	mcVtx[nMC_][0]    = p->vx();
	mcVtx[nMC_][1]    = p->vy();
	mcVtx[nMC_][2]    = p->vz();
	mcPt[nMC_]        = p->pt();
	mcMass[nMC_]      = p->mass();
	mcEta[nMC_]       = p->eta();
	mcPhi[nMC_]       = p->phi();
	mcE[nMC_]         = p->energy();
	mcEt[nMC_]        = p->et();
	mcGMomPID[nMC_]   = -999;
	mcMomPID[nMC_]    = -999;
	mcMomPt[nMC_]     = -999;
	mcMomMass[nMC_]   = -999;
	mcMomEta[nMC_]    = -999;
	mcMomPhi[nMC_]    = -999;
	mcDecayType[nMC_] = -999;
	mcParentage[nMC_] = 
	  particleHistory.hasLeptonParent()*16   + 
	  particleHistory.hasBosonParent()*8     + 
	  particleHistory.hasNonPromptParent()*4 +
	  particleHistory.hasQCDParent()*2       +
	  particleHistory.hasExoticParent();
	mcStatus[nMC_]    = p->status();

	// if genParticle is W or Z, check its decay type
	if ( ip->pdgId() == 23 || abs(ip->pdgId()) == 24 ) {
	  for (size_t k=0; k < p->numberOfDaughters(); ++k) {
	    const Candidate *dp = p->daughter(k);
	    if (abs(dp->pdgId())<=6)
	      mcDecayType[nMC_] = 1;
	    else if (abs(dp->pdgId())==11 || abs(dp->pdgId())==12)
	      mcDecayType[nMC_] = 2;
	    else if (abs(dp->pdgId())==13 || abs(dp->pdgId())==14)
	      mcDecayType[nMC_] = 3;
	    else if (abs(dp->pdgId())==15 || abs(dp->pdgId())==16)
	      mcDecayType[nMC_] = 4;
	  }
	}

	if ( particleHistory.hasRealParent() ) {
	  reco::GenParticleRef momRef = particleHistory.parent();
	  if ( momRef.isNonnull() && momRef.isAvailable() ) {
	    mcMomPID[nMC_]  = momRef->pdgId();
	    mcMomPt[nMC_]   = momRef->pt();
	    mcMomMass[nMC_] = momRef->mass();
	    mcMomEta[nMC_]  = momRef->eta();
	    mcMomPhi[nMC_]  = momRef->phi();
	    
	    // get Granny
	    genpartparentage::GenParticleParentage motherParticle(momRef);
	    if ( motherParticle.hasRealParent() ) {
	      reco::GenParticleRef granny = motherParticle.parent();
	      mcGMomPID[nMC_] = granny->pdgId();
	    }
	  }
	}
	
	mcIndex[nMC_] = genIndex-1;
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

    //cout<<i<<" "<<HLT_[i]<<" "<<hltConfigProvider_.prescaleValue(e, es, hlNames[i])<<" "<<hlNames[i]<<endl;

    //const vector<string> filterLabels = hltConfigProvider_.moduleLabels(hlNames[i]);
    //for (vector<string>::const_iterator labelIter = filterLabels.begin(); labelIter != filterLabels.end(); ++labelIter) {
    //edm::InputTag testTag(*labelIter, "", "HLT");
    //cout<<(*labelIter)<<" "<<hltConfigProvider_.saveTags(*labelIter)<<endl;
    //}
    
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
    trkMETx_[v]   = trkMET[v*4];
    trkMETy_[v]   = trkMET[v*4+1];
    trkMETPhi_[v] = trkMET[v*4+2];
    trkMET_[v]    = trkMET[v*4+3];
  }

  for (int i=0; i<10; ++i) metFilters_[i] = -1;
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
    if (e.getByLabel(InputTag("HBHENoiseFilterResultProducer","HBHENoiseFilterResult"), HBHENoiseFilter)) metFilters_[1] = (int) (*HBHENoiseFilter);
    
    Handle<bool> HcalLaserFilter;
    if (e.getByLabel("hcalLaserEventFilter", HcalLaserFilter)) metFilters_[2] = (int) (*HcalLaserFilter);
    
    Handle<bool> EcalDeadCellFilter;
    if (e.getByLabel("EcalDeadCellTriggerPrimitiveFilter", EcalDeadCellFilter)) metFilters_[3] = (int) (*EcalDeadCellFilter);
    
    Handle<bool> TrackingFailureFilter;
    if (e.getByLabel("trackingFailureFilter", TrackingFailureFilter)) metFilters_[4] = (int) (*TrackingFailureFilter);
    
    Handle<bool> EEBadScFilter;
    if (e.getByLabel("eeBadScFilter", EEBadScFilter)) metFilters_[5] = (int) (*EEBadScFilter);
    
    Handle<bool> EcalLaserFilter;
    if (e.getByLabel("ecalLaserCorrFilter", EcalLaserFilter)) metFilters_[6] = (int) (*EcalLaserFilter);
    
    Handle<bool> Manystripclus53X;
    if (e.getByLabel("manystripclus53X", Manystripclus53X)) metFilters_[7] = (int) (*Manystripclus53X);
    
    Handle<bool> Toomanystripclus53X;
    if (e.getByLabel("toomanystripclus53X", Toomanystripclus53X)) metFilters_[8] = (int) (*Toomanystripclus53X);
    
    Handle<bool> LogErrorTooManyClusters;
    if (e.getByLabel("logErrorTooManyClusters", LogErrorTooManyClusters)) metFilters_[9] = (int) (*LogErrorTooManyClusters);
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
      eleTrg_[nEle_][0] = (eleTrigRef1.isAvailable()) ? 1 : -99;
      eleTrg_[nEle_][1] = (eleTrigRef2.isAvailable()) ? 1 : -99;
      eleTrg_[nEle_][2] = (eleTrigRef3.isAvailable()) ? 1 : -99;
      eleTrg_[nEle_][3] = (eleTrigRef4.isAvailable()) ? 1 : -99;
      eleTrg_[nEle_][4] = (eleTrigRef5.isAvailable()) ? 1 : -99;
      eleTrg_[nEle_][5] = (eleTrigRef6.isAvailable()) ? 1 : -99;
      eleTrg_[nEle_][6] = (eleTrigRef7.isAvailable()) ? 1 : -99;
      eleTrg_[nEle_][7] = (eleTrigRef8.isAvailable()) ? 1 : -99;
      eleTrg_[nEle_][8] = (eleTrigRef9.isAvailable()) ? 1 : -99;
      eleTrg_[nEle_][9] = (eleTrigRef10.isAvailable()) ? 1 : -99;
      eleTrg_[nEle_][10] = (eleTrigRef11.isAvailable()) ? 1 : -99;
      eleTrg_[nEle_][11] = (eleTrigRef12.isAvailable()) ? 1 : -99;
      eleTrg_[nEle_][12] = (eleTrigRef13.isAvailable()) ? 1 : -99;
      eleTrg_[nEle_][13] = (eleTrigRef14.isAvailable()) ? 1 : -99;
      eleTrg_[nEle_][14] = (eleTrigRef15.isAvailable()) ? 1 : -99;
      eleTrg_[nEle_][15] = (eleTrigRef16.isAvailable()) ? 1 : -99;

      edm::Ptr<reco::Candidate> recoEleRef = iEle->originalObjectRef();                                                                                  
      const reco::GsfElectron *recoElectron = dynamic_cast<const reco::GsfElectron *>(recoEleRef.get());        

      eleClass_[nEle_]            = iEle->classification();
      eleIsEcalDriven_[nEle_]     = iEle->ecalDrivenSeed();
      eleCharge_[nEle_]           = iEle->charge();
      eleChargeConsistent_[nEle_] = (Int_t)iEle->isGsfCtfScPixChargeConsistent();

      eleEn_[nEle_]           = recoElectron->energy(); // iEle->energy();
      eleEcalEn_[nEle_]       = iEle->ecalEnergy();
      elePt_[nEle_]           = recoElectron->pt(); // iEle->pt();
      eleEta_[nEle_]          = iEle->eta();
      elePhi_[nEle_]          = iEle->phi();
      eleHoverE_[nEle_]       = iEle->hadronicOverEm();
      eleHoverE12_[nEle_]     = iEle->hcalOverEcalBc();
      eleEoverP_[nEle_]       = iEle->eSuperClusterOverP();

      elePin_[nEle_]       = iEle->trackMomentumAtVtx().R();
      elePout_[nEle_]      = iEle->trackMomentumOut().R();
      eleTrkMomErr_[nEle_] = iEle->trackMomentumError();
      eleBrem_[nEle_]      = iEle->fbrem();
      
      eleD0_[nEle_]   = iEle->gsfTrack()->dxy(pv);
      eleDz_[nEle_]   = iEle->gsfTrack()->dz(pv);
      eleD0GV_[nEle_] = iEle->gsfTrack()->dxy(gv);
      eleDzGV_[nEle_] = iEle->gsfTrack()->dz(gv);

      eledEtaAtVtx_[nEle_] = iEle->deltaEtaSuperClusterTrackAtVtx();
      eledPhiAtVtx_[nEle_] = iEle->deltaPhiSuperClusterTrackAtVtx();
      
      eleMissHits_[nEle_]   = iEle->gsfTrack()->trackerExpectedHitsInner().numberOfHits();
      eleConvDist_[nEle_]   = iEle->convDist();
      eleConvDcot_[nEle_]   = iEle->convDcot();

      eleConvVtxFit_[nEle_] = (int) ConversionTools::hasMatchedConversion(*recoElectron, convH_, beamSpotHandle_->position(), true, 2.0, 1e-6, 0);      

      // IP3D
      eleIP3D_[nEle_] = -999.;
      eleIP3DErr_[nEle_] = -99.;
      if (iEle->gsfTrack().isNonnull()) {
	const double gsfsign = ((- iEle->gsfTrack()->dxy(pv)) >= 0) ? 1. : -1.;
	const reco::TransientTrack tt_ele = thebuilder.build(iEle->gsfTrack());
	const std::pair<bool,Measurement1D> ip3dpv_ele = IPTools::absoluteImpactParameter3D(tt_ele, *recVtxs_->begin());
	if (ip3dpv_ele.first) {
	  eleIP3D_[nEle_]    = gsfsign*ip3dpv_ele.second.value();
	  eleIP3DErr_[nEle_] = ip3dpv_ele.second.error();
	}
      }

      // 2012 MVA eID
      eleIDMVANonTrig_[nEle_] = iEle->electronID("mvaNonTrigV0");
      eleIDMVATrig_[nEle_]    = iEle->electronID("mvaTrigV0");

      // Access PF isolation
      elePFChIso03_[nEle_]  = (*(*electronIsoVals03)[0])[recoEleRef];
      elePFPhoIso03_[nEle_] = (*(*electronIsoVals03)[1])[recoEleRef];
      elePFNeuIso03_[nEle_] = (*(*electronIsoVals03)[2])[recoEleRef];

      elePFChIso04_[nEle_]  = (*(*electronIsoVals04)[0])[recoEleRef];
      elePFPhoIso04_[nEle_] = (*(*electronIsoVals04)[1])[recoEleRef];
      elePFNeuIso04_[nEle_] = (*(*electronIsoVals04)[2])[recoEleRef];

      // Access super cluster
      eleSCEta_[nEle_]      = iEle->superCluster()->eta();
      eleSCPhi_[nEle_]      = iEle->superCluster()->phi();
      eleSCRawEn_[nEle_]    = iEle->superCluster()->rawEnergy();
      eleSCEn_[nEle_]       = iEle->superCluster()->energy();
      eleESEn_[nEle_]       = iEle->superCluster()->preshowerEnergy();
      eleSCEtaWidth_[nEle_] = iEle->superCluster()->etaWidth();
      eleSCPhiWidth_[nEle_] = iEle->superCluster()->phiWidth();

      eleVtx_[nEle_][0] = iEle->trackPositionAtVtx().x();
      eleVtx_[nEle_][1] = iEle->trackPositionAtVtx().y();
      eleVtx_[nEle_][2] = iEle->trackPositionAtVtx().z();
      
      // Gen Particle
      eleGenIndex_[nEle_] = -1;
      int EleGenIndex = 0;
      if (!isData_) {
        if ((*iEle).genLepton() && genParticlesHandle_.isValid() ) {
	  
          for (vector<GenParticle>::const_iterator iGen = genParticlesHandle_->begin(); iGen != genParticlesHandle_->end(); ++iGen) {
	    
            if (iGen->p4() == (*iEle).genLepton()->p4() && iGen->pdgId() == (*iEle).genLepton()->pdgId() && iGen->status() == (*iEle).genLepton()->status()) {
	      
	      eleGenIndex_[nEle_] = EleGenIndex;
	      
	      const Candidate *elep = (const Candidate*)&(*iGen);
	      
	      for (size_t j=0; j<elep->numberOfMothers(); ++j) {
		
		elemom = elep->mother(j);
		eleGenMomPID_[nEle_] = elemom->pdgId();
		eleGenMomPt_[nEle_] = elemom->pt();
		if (elemom->mother()) eleGenGMomPID_[nEle_] = elemom->mother()->pdgId();
	      }
            }
            EleGenIndex++;
          }
        }
      }
      
      eleIsoTrkDR03_[nEle_]       = iEle->dr03TkSumPt();
      eleIsoEcalDR03_[nEle_]      = iEle->dr03EcalRecHitSumEt();
      eleIsoHcalDR03_[nEle_]      = iEle->dr03HcalTowerSumEt();
      eleIsoHcalDR0312_[nEle_]    = iEle->dr03HcalDepth1TowerSumEt();
      eleIsoTrkDR04_[nEle_]       = iEle->dr04TkSumPt();
      eleIsoEcalDR04_[nEle_]      = iEle->dr04EcalRecHitSumEt();
      eleIsoHcalDR04_[nEle_]      = iEle->dr04HcalTowerSumEt();
      eleIsoHcalDR0412_[nEle_]    = iEle->dr04HcalDepth1TowerSumEt();

      eleModIsoTrk_[nEle_]        = iEle->userIso(0);
      eleModIsoEcal_[nEle_]       = iEle->userIso(1);
      eleModIsoHcal_[nEle_]       = iEle->userIso(2);
      
      const reco::CaloClusterPtr eleSeed = (*iEle).superCluster()->seed();
      
      vector<float> eleCov;
      eleCov = lazyTool->localCovariances(*eleSeed);
      eleSigmaIEtaIEta_[nEle_] = iEle->sigmaIetaIeta();
      eleSigmaIEtaIPhi_[nEle_] = eleCov[1];
      eleSigmaIPhiIPhi_[nEle_] = eleCov[2];
      
      eleEmax_[nEle_]     = lazyTool->eMax(*eleSeed);
      eleE1x5_[nEle_]     = lazyTool->e1x5(*eleSeed);
      eleE3x3_[nEle_]     = lazyTool->e3x3(*eleSeed);
      eleE5x5_[nEle_]     = iEle->e5x5();
      eleE2x5Max_[nEle_]  = iEle->e2x5Max();

      // Energy Regression Correction
      if (!egCorrEle_.IsInitialized()) {   
	//std::string filenameEle = "http://homepages.spa.umn.edu/~rekovic/cms/regweights52xV3/gbrv3ele_52x.root";
	std::string filenameEle = "gbrv3ele_52x.root";
	egCorrEle_.Initialize(es, filenameEle);
      }

      eleRegrE_[nEle_]    = iEle->p4(reco::GsfElectron::P4_COMBINATION).t();
      eleRegrEerr_[nEle_] = iEle->p4Error(reco::GsfElectron::P4_COMBINATION);

      std::pair<double,double> regrCorEle = egCorrEle_.CorrectedEnergyWithErrorV3(*eleRef, *(recVtxsBS_.product()), *rhoHandle_enReg, *lazyTool, es);
      elePhoRegrE_[nEle_]    = regrCorEle.first;
      elePhoRegrEerr_[nEle_] = regrCorEle.second;
      
      eleSeedTime_[nEle_] = -999.;
      eleRecoFlag_[nEle_] = -999.;

      // where is electron ? (0: EB, 1: EE, 2: EBGap, 3: EEGap, 4: EBEEGap)
      elePos_[nEle_] = -1;
      if (iEle->isEB() == true) elePos_[nEle_] = 0;
      if (iEle->isEE() == true) elePos_[nEle_] = 1;
      if (iEle->isEBGap() == true) elePos_[nEle_] = 2;
      if (iEle->isEEGap() == true) elePos_[nEle_] = 3;
      if (iEle->isEBEEGap() == true) elePos_[nEle_] = 4;

      DetId eleSeedDetId = lazyTool->getMaximum(*eleSeed).first;
      
      if ( iEle->isEB() && EBReducedRecHits_.isValid() ) {
        EcalRecHitCollection::const_iterator eleebrhit = EBReducedRecHits_->find(eleSeedDetId);
        if ( eleebrhit != EBReducedRecHits_->end() ) { 
	  eleSeedTime_[nEle_] = eleebrhit->time(); 
	  eleRecoFlag_[nEle_] = eleebrhit->recoFlag();
	}
      } else if ( EEReducedRecHits_.isValid() ) {
        EcalRecHitCollection::const_iterator eleeerhit = EEReducedRecHits_->find(eleSeedDetId);
        if ( eleeerhit != EEReducedRecHits_->end() ) { 
	  eleSeedTime_[nEle_] = eleeerhit->time(); 
	  eleRecoFlag_[nEle_] = eleeerhit->recoFlag();
	}
      }
      
      for (int i=0; i<2; i++) eleESDetId_[nEle_][i] = 0;
      for (int i=0; i<3; ++i) 
	for (int j=0; j<62; ++j)
	  eleESHits_[nEle_][i][j] = 0;
      for (int i=0; i<3; ++i) eleESEffSigmaRR_[nEle_][i] = 0.;
      if (dumpESClusterInfo_) {
	for (int i=0; i<2; i++) {
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
        for (int i=0; i<62; ++i) {
          eleESHits_[nEle_][0][i] = eleESHits0[i];
          eleESHits_[nEle_][1][i] = eleESHits1[i];
          eleESHits_[nEle_][2][i] = eleESHits2[i];
        }

        vector<float> eleESShape = getESEffSigmaRR(eleESHits0);
        eleESEffSigmaRR_[nEle_][0] = eleESShape[0];
        eleESEffSigmaRR_[nEle_][1] = eleESShape[1];
        eleESEffSigmaRR_[nEle_][2] = eleESShape[2];
	
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
      
      //------------ Applying corrections ----------
      if (develop_) {
	eleBrLinear_[nEle_]   = eleSCPhiWidth_[nEle_]/eleSCEtaWidth_[nEle_];
	eleNBC_[nEle_]        = iEle->superCluster()->clustersSize();
	
	float absEta = fabs(eleSCEta_[nEle_]);
	float energyWithEtaCorrection;
	float pTWithEtaCorrection;
	
	// C(eta) for EB only
	if (iEle->isEB()){
	  energyWithEtaCorrection = eleSCRawEn_[nEle_]/fcorrs::f5x5((int)(absEta*(5/0.087)));
	  pTWithEtaCorrection     = energyWithEtaCorrection / cosh(eleSCEta_[nEle_]);
	  
	  eleCetaCorrE_[nEle_]  = energyWithEtaCorrection;
	  eleCetaCorrEt_[nEle_] = pTWithEtaCorrection;
	  
	  // f(brem)-corrected energies
	  eleBremCorrE_[nEle_]  = fcorrs::fBrem(eleBrLinear_[nEle_], eleCetaCorrE_[nEle_]);
	  eleBremCorrEt_[nEle_] = eleBremCorrE_[nEle_] / cosh(eleSCEta_[nEle_]);
	  
	  // fully-corrected energies
	  eleFullCorrEt_[nEle_] = fcorrs::fullCorr(eleBremCorrEt_[nEle_], absEta);
	  eleFullCorrE_[nEle_]  = eleFullCorrEt_[nEle_] * cosh(eleSCEta_[nEle_]); 
	}
	else{
	  energyWithEtaCorrection = eleSCRawEn_[nEle_] + (*iEle).superCluster()->preshowerEnergy();
	  pTWithEtaCorrection     = energyWithEtaCorrection / cosh(eleSCEta_[nEle_]);
	  
	  eleCetaCorrE_[nEle_]  = energyWithEtaCorrection;
	  eleCetaCorrEt_[nEle_] = pTWithEtaCorrection;
	  
	  // f(brem)-corrected energies
	  eleBremCorrE_[nEle_]  = fcorrs::fBrem_ee(eleBrLinear_[nEle_], eleCetaCorrE_[nEle_]);
	  eleBremCorrEt_[nEle_] = eleBremCorrE_[nEle_] / cosh(eleSCEta_[nEle_]);
	  
	  // fully-corrected energies
	  eleFullCorrEt_[nEle_] = fcorrs::fullCorr_ee(eleBremCorrEt_[nEle_], absEta);
	  eleFullCorrE_[nEle_]  = eleFullCorrEt_[nEle_] * cosh(eleSCEta_[nEle_]);
	}
      }

      math::XYZPoint vtxPoint[100];
      for (size_t iv=0; iv<recVtxsBS_->size(); ++iv) { 
	vtxPoint[iv] = math::XYZPoint(vtxbs_[iv][0], vtxbs_[iv][1], vtxbs_[iv][2]);
	TVector3 *v3 = new TVector3(iEle->caloPosition().x()-(*recVtxsBS_)[iv].x(), iEle->caloPosition().y()-(*recVtxsBS_)[iv].y(), iEle->caloPosition().z()-(*recVtxsBS_)[iv].z());
	eleEtaVtx_[nEle_][iv] = v3->Eta();
	elePhiVtx_[nEle_][iv] = v3->Phi();
	eleEtVtx_[nEle_][iv]  = iEle->energy() * TMath::Sin(2*TMath::ATan(TMath::Exp( - v3->Eta() )));
	eleD0Vtx_[nEle_][iv]  = iEle->gsfTrack()->dxy(vtxPoint[iv]);
	eleDzVtx_[nEle_][iv]  = iEle->gsfTrack()->dz(vtxPoint[iv]);
	delete v3;
      }
      
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

      phoTrg_[nPho_][0] = (phoTrigRef1.isAvailable()) ? 1 : -99;
      phoTrg_[nPho_][1] = (phoTrigRef2.isAvailable()) ? 1 : -99;
      phoTrg_[nPho_][2] = (phoTrigRef3.isAvailable()) ? 1 : -99;
      phoTrg_[nPho_][3] = (phoTrigRef4.isAvailable()) ? 1 : -99;
      phoTrg_[nPho_][4] = (phoTrigRef5.isAvailable()) ? 1 : -99;
      phoTrg_[nPho_][5] = (phoTrigRef6.isAvailable()) ? 1 : -99;
      phoTrg_[nPho_][6] = (phoTrigRef7.isAvailable()) ? 1 : -99;
      for (int i=7; i<8; ++i) phoTrg_[nPho_][i] = -99;

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

      phoTrgFilter_[nPho_][0]  = (eleTrigFRef1.isAvailable()) ? 1 : -99;
      phoTrgFilter_[nPho_][1]  = (eleTrigFRef2.isAvailable()) ? 1 : -99;
      phoTrgFilter_[nPho_][2]  = (eleTrigFRef3.isAvailable()) ? 1 : -99;
      phoTrgFilter_[nPho_][3]  = (eleTrigFRef4.isAvailable()) ? 1 : -99;
      phoTrgFilter_[nPho_][4]  = (eleTrigFRef5.isAvailable()) ? 1 : -99;
      phoTrgFilter_[nPho_][5]  = (eleTrigFRef6.isAvailable()) ? 1 : -99;
      phoTrgFilter_[nPho_][6]  = -99;
      phoTrgFilter_[nPho_][7]  = -99;
      phoTrgFilter_[nPho_][8]  = -99;
      phoTrgFilter_[nPho_][9]  = -99;
      phoTrgFilter_[nPho_][10] = (phoTrigFRef1.isAvailable()) ? 1 : -99;
      phoTrgFilter_[nPho_][11] = (phoTrigFRef2.isAvailable()) ? 1 : -99;
      phoTrgFilter_[nPho_][12] = (phoTrigFRef3.isAvailable()) ? 1 : -99;
      phoTrgFilter_[nPho_][13] = (phoTrigFRef4.isAvailable()) ? 1 : -99;
      phoTrgFilter_[nPho_][14] = (phoTrigFRef5.isAvailable()) ? 1 : -99;
      phoTrgFilter_[nPho_][15] = (phoTrigFRef6.isAvailable()) ? 1 : -99;
      phoTrgFilter_[nPho_][16] = (phoTrigFRef7.isAvailable()) ? 1 : -99;
      phoTrgFilter_[nPho_][17] = (phoTrigFRef8.isAvailable()) ? 1 : -99;
      phoTrgFilter_[nPho_][18] = (phoTrigFRef9.isAvailable()) ? 1 : -99;
      phoTrgFilter_[nPho_][19] = (phoTrigFRef10.isAvailable()) ? 1 : -99;
      phoTrgFilter_[nPho_][20] = (phoTrigFRef11.isAvailable()) ? 1 : -99;
      phoTrgFilter_[nPho_][21] = (phoTrigFRef12.isAvailable()) ? 1 : -99;
      phoTrgFilter_[nPho_][22] = (phoTrigFRef13.isAvailable()) ? 1 : -99;
      phoTrgFilter_[nPho_][23] = (phoTrigFRef14.isAvailable()) ? 1 : -99;
      phoTrgFilter_[nPho_][24] = (phoTrigFRef15.isAvailable()) ? 1 : -99;
      phoTrgFilter_[nPho_][25] = (phoTrigFRef16.isAvailable()) ? 1 : -99;
      phoTrgFilter_[nPho_][26] = (phoTrigFRef17.isAvailable()) ? 1 : -99;
      phoTrgFilter_[nPho_][27] = (phoTrigFRef18.isAvailable()) ? 1 : -99;
      phoTrgFilter_[nPho_][28] = (phoTrigFRef19.isAvailable()) ? 1 : -99;
      phoTrgFilter_[nPho_][28] = (phoTrigFRef20.isAvailable()) ? 1 : -99;
      phoTrgFilter_[nPho_][30] = (phoTrigFRef21.isAvailable()) ? 1 : -99;
      phoTrgFilter_[nPho_][31] = (phoTrigFRef22.isAvailable()) ? 1 : -99;
      phoTrgFilter_[nPho_][32] = (phoTrigFRef23.isAvailable()) ? 1 : -99;
      phoTrgFilter_[nPho_][33] = (phoTrigFRef24.isAvailable()) ? 1 : -99;
      phoTrgFilter_[nPho_][34] = (phoTrigFRef25.isAvailable()) ? 1 : -99;
      phoTrgFilter_[nPho_][35] = (phoTrigFRef26.isAvailable()) ? 1 : -99;
      phoTrgFilter_[nPho_][36] = (phoTrigFRef27.isAvailable()) ? 1 : -99;
      phoTrgFilter_[nPho_][37] = (phoTrigFRef28.isAvailable()) ? 1 : -99;
      phoTrgFilter_[nPho_][38] = (phoTrigFRef29.isAvailable()) ? 1 : -99;
      phoTrgFilter_[nPho_][39] = (phoTrigFRef30.isAvailable()) ? 1 : -99;
      phoTrgFilter_[nPho_][40] = (phoTrigFRef31.isAvailable()) ? 1 : -99;
      for (int i=41; i<50; ++i) phoTrgFilter_[nPho_][i] = -99;

      math::XYZPoint photonXYZ(iPho->caloPosition().x(),iPho->caloPosition().y(),iPho->caloPosition().z());
     
      edm::Ptr<reco::Candidate> recoPhoRef = iPho->originalObjectRef();
      const reco::Photon *recoPhoton = dynamic_cast<const reco::Photon *>(recoPhoRef.get());
      //isolation.mvaID(pfCandidates.product(),recoPhoton,recVtxs_);

      // Access PF isolation
      //phoPFChIso_[nPho_]  = (*(*photonIsoVals)[0])[recoPhoRef];
      //phoPFPhoIso_[nPho_] = (*(*photonIsoVals)[1])[recoPhoRef];
      //phoPFNeuIso_[nPho_] = (*(*photonIsoVals)[2])[recoPhoRef];

      // PF isolation from Alternate code
      isolator.fGetIsolation(recoPhoton, &thePfColl, myVtxRef, recVtxs_);
      phoPFChIso_[nPho_]  = isolator.getIsolationCharged();
      phoPFNeuIso_[nPho_] = isolator.getIsolationNeutral();
      phoPFPhoIso_[nPho_] = isolator.getIsolationPhoton();

      phoIsPhoton_[nPho_] = iPho->isPhoton();
      phoE_[nPho_]   = iPho->energy();
      phoEt_[nPho_]  = iPho->et();
      phoEta_[nPho_] = iPho->eta();
      phoPhi_[nPho_] = iPho->phi();

      phoVtx_[nPho_][0] = iPho->vx();
      phoVtx_[nPho_][1] = iPho->vy();
      phoVtx_[nPho_][2] = iPho->vz();

      phoSCPos_[nPho_][0] = (*iPho).superCluster()->x();
      phoSCPos_[nPho_][1] = (*iPho).superCluster()->y();
      phoSCPos_[nPho_][2] = (*iPho).superCluster()->z();
      phoCaloPos_[nPho_][0] = iPho->caloPosition().x();
      phoCaloPos_[nPho_][1] = iPho->caloPosition().y();
      phoCaloPos_[nPho_][2] = iPho->caloPosition().z();

      phoR9_[nPho_]               = iPho->r9();
      phoTrkIsoHollowDR03_[nPho_] = iPho->trkSumPtHollowConeDR03();
      phoEcalIsoDR03_[nPho_]      = iPho->ecalRecHitSumEtConeDR03();
      phoHcalIsoDR03_[nPho_]      = iPho->hcalTowerSumEtConeDR03();
      phoHcalIsoDR0312_[nPho_]    = iPho->hcalTowerSumEtConeDR03() + (iPho->hadronicOverEm() - iPho->hadTowOverEm())*(*iPho).superCluster()->energy()/cosh((*iPho).superCluster()->eta());

      phoTrkIsoHollowDR04_[nPho_] = iPho->trkSumPtHollowConeDR04();
      phoEcalIsoDR04_[nPho_]      = iPho->ecalRecHitSumEtConeDR04();
      phoHcalIsoDR04_[nPho_]      = iPho->hcalTowerSumEtConeDR04();
      phoHcalIsoDR0412_[nPho_]    = iPho->hcalTowerSumEtConeDR04() + (iPho->hadronicOverEm() - iPho->hadTowOverEm())*(*iPho).superCluster()->energy()/cosh((*iPho).superCluster()->eta());

      phoHoverE_[nPho_]           = iPho->hadronicOverEm();     
      phoHoverE12_[nPho_]         = iPho->hadTowOverEm();

      phoOverlap_[nPho_]          = (int) iPho->hasOverlaps("electrons");
      phohasPixelSeed_[nPho_]     = (int) iPho->hasPixelSeed();
      phoEleVeto_[nPho_]          = (int) ConversionTools::hasMatchedPromptElectron(recoPhoton->superCluster(), gsfElectronHandle_, convH_, beamSpotHandle_->position());

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
	phoCiCPF4phopfIso005_[nPho_] = val_pfiso_photon005;
	phoCiCPF4phopfIso01_[nPho_]  = val_pfiso_photon01;
	phoCiCPF4phopfIso02_[nPho_]  = val_pfiso_photon02;
      }
      float val_pfiso_photon03 = cicPhotonId_->pfEcalIso(recophoRef, 0.3, 0.0, 0.070, 0.015, 0.0, 0.0, 0.0, reco::PFCandidate::gamma);
      float val_pfiso_photon04 = cicPhotonId_->pfEcalIso(recophoRef, 0.4, 0.0, 0.070, 0.015, 0.0, 0.0, 0.0, reco::PFCandidate::gamma);
      phoCiCPF4phopfIso03_[nPho_] = val_pfiso_photon03;
      phoCiCPF4phopfIso04_[nPho_] = val_pfiso_photon04;
      if (dumpESClusterInfo_) {
	float val_pfiso_photon05 = cicPhotonId_->pfEcalIso(recophoRef, 0.5, 0.0, 0.070, 0.015, 0.0, 0.0, 0.0, reco::PFCandidate::gamma);
	float val_pfiso_photon06 = cicPhotonId_->pfEcalIso(recophoRef, 0.6, 0.0, 0.070, 0.015, 0.0, 0.0, 0.0, reco::PFCandidate::gamma);
	float val_pfiso_photon07 = cicPhotonId_->pfEcalIso(recophoRef, 0.7, 0.0, 0.070, 0.015, 0.0, 0.0, 0.0, reco::PFCandidate::gamma);
	float val_pfiso_photon08 = cicPhotonId_->pfEcalIso(recophoRef, 0.8, 0.0, 0.070, 0.015, 0.0, 0.0, 0.0, reco::PFCandidate::gamma);
	phoCiCPF4phopfIso05_[nPho_] = val_pfiso_photon05;
	phoCiCPF4phopfIso06_[nPho_] = val_pfiso_photon06;
	phoCiCPF4phopfIso07_[nPho_] = val_pfiso_photon07;
	phoCiCPF4phopfIso08_[nPho_] = val_pfiso_photon08;
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
	phoCiCPF4phopfIsoNoVETO005_[nPho_] = val_pfiso_photon005NoV;
	phoCiCPF4phopfIsoNoVETO01_[nPho_]  = val_pfiso_photon01NoV;
	phoCiCPF4phopfIsoNoVETO02_[nPho_]  = val_pfiso_photon02NoV;
	phoCiCPF4phopfIsoNoVETO03_[nPho_]  = val_pfiso_photon03NoV;
	phoCiCPF4phopfIsoNoVETO04_[nPho_]  = val_pfiso_photon04NoV;
	phoCiCPF4phopfIsoNoVETO05_[nPho_]  = val_pfiso_photon05NoV;
	phoCiCPF4phopfIsoNoVETO06_[nPho_]  = val_pfiso_photon06NoV;
	phoCiCPF4phopfIsoNoVETO07_[nPho_]  = val_pfiso_photon07NoV;
	phoCiCPF4phopfIsoNoVETO08_[nPho_]  = val_pfiso_photon08NoV;
      }
      temp.clear();
      temp.push_back(reco::PFCandidate::h); 
      if (dumpESClusterInfo_) {
	std::vector<float> vtxIsolations005   = cicPhotonId_->pfTkIsoWithVertex(recophoRef, 0.05, 0.02, 0.02, 0.0, 0.2, 0.1, reco::PFCandidate::h);
	std::vector<float> vtxIsolations01    = cicPhotonId_->pfTkIsoWithVertex(recophoRef, 0.1, 0.02, 0.02, 0.0, 0.2, 0.1, reco::PFCandidate::h);
	for (unsigned int v=0; v<vtxIsolations01.size(); ++v) {
	  phoCiCPF4chgpfIso005_[nPho_][v]   = vtxIsolations005[v];
	  phoCiCPF4chgpfIso01_[nPho_][v]    = vtxIsolations01[v];
	}
      }
      std::vector<float> vtxIsolations02    = cicPhotonId_->pfTkIsoWithVertex(recophoRef, 0.2, 0.02, 0.02, 0.0, 0.2, 0.1, reco::PFCandidate::h);
      std::vector<float> vtxIsolations03    = cicPhotonId_->pfTkIsoWithVertex(recophoRef, 0.3, 0.02, 0.02, 0.0, 0.2, 0.1, reco::PFCandidate::h);
      std::vector<float> vtxIsolations04    = cicPhotonId_->pfTkIsoWithVertex(recophoRef, 0.4, 0.02, 0.02, 0.0, 0.2, 0.1, reco::PFCandidate::h);
      for (unsigned int v=0; v<vtxIsolations03.size(); ++v) {
	phoCiCPF4chgpfIso02_[nPho_][v]    = vtxIsolations02[v];
	phoCiCPF4chgpfIso03_[nPho_][v]    = vtxIsolations03[v];
	phoCiCPF4chgpfIso04_[nPho_][v]    = vtxIsolations04[v];
	//cout<<"trk : "<<v<<" "<<phoCiCPF4chgpfIso02_[nPho_][v]<<" "<<phoCiCPF4chgpfIso03_[nPho_][v]<<" "<<phoCiCPF4chgpfIso04_[nPho_][v]<<endl;
      }
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

	for (unsigned int v=0; v<vtxIsolations05.size(); ++v) {
	  phoCiCPF4chgpfIso05_[nPho_][v]    = vtxIsolations05[v];
	  phoCiCPF4chgpfIso06_[nPho_][v]    = vtxIsolations06[v];
	  phoCiCPF4chgpfIso07_[nPho_][v]    = vtxIsolations07[v];
	  phoCiCPF4chgpfIso08_[nPho_][v]    = vtxIsolations08[v];
	  
	  phoCiCPF4chgpfIsoNoVETO005_[nPho_][v]   = vtxIsolations005NoV[v];
	  phoCiCPF4chgpfIsoNoVETO01_[nPho_][v]    = vtxIsolations01NoV[v];
	  phoCiCPF4chgpfIsoNoVETO02_[nPho_][v]    = vtxIsolations02NoV[v];
	  phoCiCPF4chgpfIsoNoVETO03_[nPho_][v]    = vtxIsolations03NoV[v];
	  phoCiCPF4chgpfIsoNoVETO04_[nPho_][v]    = vtxIsolations04NoV[v];
	  phoCiCPF4chgpfIsoNoVETO05_[nPho_][v]    = vtxIsolations05NoV[v];
	  phoCiCPF4chgpfIsoNoVETO06_[nPho_][v]    = vtxIsolations06NoV[v];
	  phoCiCPF4chgpfIsoNoVETO07_[nPho_][v]    = vtxIsolations07NoV[v];
	  phoCiCPF4chgpfIsoNoVETO08_[nPho_][v]    = vtxIsolations08NoV[v];
	} 
      }
      
      for (size_t iv=0; iv<recVtxsBS_->size(); ++iv) {
        TVector3 *v3 = new TVector3(iPho->caloPosition().x()-(*recVtxsBS_)[iv].x(), iPho->caloPosition().y()-(*recVtxsBS_)[iv].y(), iPho->caloPosition().z()-(*recVtxsBS_)[iv].z());
        phoEtaVtx_[nPho_][iv] = v3->Eta();
        phoPhiVtx_[nPho_][iv] = v3->Phi();
        phoEtVtx_[nPho_][iv]  = iPho->energy() * TMath::Sin(2*TMath::ATan(TMath::Exp( - v3->Eta() )));

        phoCiCTrkIsoDR03_[nPho_][iv] = getPhotonTrkIso((*recVtxsBS_)[iv].z(), (*recVtxsBS_)[iv].x(), (*recVtxsBS_)[iv].y(), (*recVtxsBS_)[iv].z(), v3->Eta(), v3->Phi(), tracksHandle_, 0.3, 0.02, 0.0, 1., 0.1, 1);
        phoCiCTrkIsoDR04_[nPho_][iv] = getPhotonTrkIso((*recVtxsBS_)[iv].z(), (*recVtxsBS_)[iv].x(), (*recVtxsBS_)[iv].y(), (*recVtxsBS_)[iv].z(), v3->Eta(), v3->Phi(), tracksHandle_, 0.4, 0.02, 0.0, 1., 0.1, 1);

	delete v3;
      }

      phoCiCdRtoTrk_[nPho_] = getPhotondRtoTrk(gsfElectronHandle_, *iPho, 2.5, 0);

      // where is photon ? (0: EB, 1: EE, 2: EBGap, 3: EEGap, 4: EBEEGap)
      phoPos_[nPho_] = -1;
      if (iPho->isEB() == true) phoPos_[nPho_] = 0;
      if (iPho->isEE() == true) phoPos_[nPho_] = 1;
      if (iPho->isEBGap() == true) phoPos_[nPho_] = 2;
      if (iPho->isEEGap() == true) phoPos_[nPho_] = 3;
      if (iPho->isEBEEGap() == true) phoPos_[nPho_] = 4;

      phoSeedTime_[nPho_] = -999.;
      phoRecoFlag_[nPho_] = -999.;
      phoSeedDetId1_[nPho_] = -999.;
      phoSeedDetId2_[nPho_] = -999.;
      const reco::CaloClusterPtr phoSeed = (*iPho).superCluster()->seed();
      DetId phoSeedDetId = lazyTool->getMaximum(*phoSeed).first;
   
      vector<float> phoCov;
      phoCov = lazyTool->localCovariances(*phoSeed);
      phoSigmaIEtaIEta_[nPho_] = iPho->sigmaIetaIeta();
      phoSigmaIEtaIPhi_[nPho_] = phoCov[1];
      phoSigmaIPhiIPhi_[nPho_] = phoCov[2];

      if ( iPho->isEB() && EBReducedRecHits_.isValid() ) {
        EcalRecHitCollection::const_iterator ebrhit = EBReducedRecHits_->find(phoSeedDetId);
        if ( ebrhit != EBReducedRecHits_->end() ) { 
	   phoSeedTime_[nPho_] = ebrhit->time();
           phoRecoFlag_[nPho_] = ebrhit->recoFlag();
	   EBDetId phoEBid = EBDetId(ebrhit->id());
	   phoSeedDetId1_[nPho_] = phoEBid.ieta();
	   phoSeedDetId2_[nPho_] = phoEBid.iphi();
        }
      } else if ( EEReducedRecHits_.isValid() ) {
        EcalRecHitCollection::const_iterator eerhit = EEReducedRecHits_->find(phoSeedDetId);
        if ( eerhit != EEReducedRecHits_->end() ) { 
	   phoSeedTime_[nPho_] = eerhit->time(); 
           phoRecoFlag_[nPho_] = eerhit->recoFlag();
	   EEDetId phoEEid = EEDetId(eerhit->id());
	   phoSeedDetId1_[nPho_] = phoEEid.ix();
	   phoSeedDetId2_[nPho_] = phoEEid.iy();
	}
      }

      // computing ECAL timing for spikes
      phoLICTD_[nPho_] = 99.;
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
	    if (fabs(thishit->time() - phoSeedTime_[nPho_]) > LICTD) 
	      LICTD = fabs(thishit->time() - phoSeedTime_[nPho_]);   
	  }
	} else if (((*detitr).first).det() == DetId::Ecal && ((*detitr).first).subdetId() == EcalEndcap) {
	  if (((*detitr).first).det() == DetId::Ecal && ((*detitr).first).subdetId() == EcalEndcap) {
	    EcalRecHitCollection::const_iterator eehit = EEReducedRecHits_->find(((*detitr).first));
	    EcalRecHitCollection::const_iterator thishit;
	    if (eehit != EEReducedRecHits_->end()) thishit = eehit;
	    if (eehit == EEReducedRecHits_->end()) continue;
	    if (thishit->energy() > 1.) {
	      if (fabs(thishit->time() - phoSeedTime_[nPho_]) > LICTD)
		LICTD = fabs(thishit->time() - phoSeedTime_[nPho_]);
	    }
	  }
	}
      }
      phoLICTD_[nPho_] = LICTD;
      
      phoEmax_[nPho_]       = iPho->maxEnergyXtal();
      phoE3x3_[nPho_]       = iPho->e3x3();
      phoE5x5_[nPho_]       = iPho->e5x5();
      phoE1x5_[nPho_]       = iPho->e1x5();
      phoE3x1_[nPho_]       = lazyTool->e3x1(*phoSeed);
      phoE1x3_[nPho_]       = lazyTool->e1x3(*phoSeed);
      phoE2x2_[nPho_]       = lazyTool->e2x2(*phoSeed);
      phoE2x5Max_[nPho_]    = iPho->e2x5();

      // Regression Correction
      if (!egCorrPho_.IsInitialized()) {
	//std::string filenamePho = "http://homepages.spa.umn.edu/~rekovic/cms/regweights52xV3/gbrv3ph_52x.root";
	std::string filenamePho = "gbrv3ph_52x.root";
	egCorrPho_.Initialize(es, filenamePho);
      }

      std::pair<double,double> regrCorPho = egCorrPho_.CorrectedEnergyWithErrorV3(*phoRef, *(recVtxsBS_.product()), *rhoHandle_enReg, *lazyTool, es);
      phoRegrE_[nPho_]    = regrCorPho.first;
      phoRegrEerr_[nPho_] = regrCorPho.second;

      // Gen Particle
      phoGenIndex_[nPho_]  = -999;
      phoGenMomPID[nPho_]  = -999;
      phoGenMomPt[nPho_]   = -999;
      phoGenGMomPID[nPho_] = -999;
      int phoGenIndex = 0;
      if ( !isData_ && genParticlesHandle_.isValid() ) {
        if ((*iPho).genPhoton()) {
          for (vector<GenParticle>::const_iterator iGen = genParticlesHandle_->begin(); iGen != genParticlesHandle_->end(); ++iGen) {

            if (iGen->p4() == (*iPho).genPhoton()->p4() && iGen->pdgId() == (*iPho).genPhoton()->pdgId() && iGen->status() == (*iPho).genPhoton()->status()) {

              phoGenIndex_[nPho_] = phoGenIndex;

              const Candidate *phop = (const Candidate*)&(*iGen);

              for (size_t j=0; j<phop->numberOfMothers(); ++j) {
                phomom = phop->mother(j);
                phoGenMomPID[nPho_] = phomom->pdgId();
                phoGenMomPt[nPho_] = phomom->pt();
                if (phomom->mother()) phoGenGMomPID[nPho_] = phomom->mother()->pdgId();
              }
            }

            phoGenIndex++;
          }
        }
      }

      // Super Cluster
      phoSCE_[nPho_]        = (*iPho).superCluster()->energy();
      phoSCRawE_[nPho_]     = (*iPho).superCluster()->rawEnergy();
      phoESEn_[nPho_]       = (*iPho).superCluster()->preshowerEnergy();
      phoSCEta_[nPho_]      = (*iPho).superCluster()->eta();
      phoSCPhi_[nPho_]      = (*iPho).superCluster()->phi();
      phoSCEt_[nPho_]       = (*iPho).superCluster()->energy()/cosh(phoSCEta_[nPho_]);
      phoSCEtaWidth_[nPho_] = (*iPho).superCluster()->etaWidth();
      phoSCPhiWidth_[nPho_] = (*iPho).superCluster()->phiWidth();
      phoSCBrem_[nPho_]     = phoSCPhiWidth_[nPho_]/phoSCEtaWidth_[nPho_]; 
      
      // supercluster removal PF isolations
      phoSCRChIso_[nPho_]   = scRemover.PFIsolation("charged", (*iPho).superCluster(), 0);
      phoSCRPhoIso_[nPho_]  = scRemover.PFIsolation("photon",  (*iPho).superCluster());
      phoSCRNeuIso_[nPho_]  = scRemover.PFIsolation("neutral", (*iPho).superCluster());
    
      //cout<<"==="<<endl;
      //cout<<phoPFChIso_[nPho_]<<" "<<phoPFPhoIso_[nPho_]<<" "<<phoPFNeuIso_[nPho_]<<endl;
      //cout<<phoSCRChIso_[nPho_]<<" "<<phoSCRPhoIso_[nPho_]<<" "<<phoSCRNeuIso_[nPho_]<<endl;
      
      // skim on good photons
      if (phoHoverE_[nPho_] < 0.15 && fabs(phoSCEta_[nPho_]) < 2.5) {
	if (fabs(phoSCEta_[nPho_]) < 1.4442) {
	  if (phoSigmaIEtaIEta_[nPho_] < 0.017) nGoodPho++;
	} else {
	  if (phoSigmaIEtaIEta_[nPho_] < 0.04)  nGoodPho++;
	}
      }

      // Conversion
      phoIsConv_[nPho_] = iPho->hasConversionTracks();
      phoNConv_[nPho_]  = iPho->conversions().size();
 
      phoConvInvMass_[nPho_]       = 0;
      phoConvCotTheta_[nPho_]      = 0;
      phoConvEoverP_[nPho_]        = 0;
      phoConvZofPVfromTrks_[nPho_] = 0;
      phoConvMinDist_[nPho_]       = 0;
      phoConvdPhiAtVtx_[nPho_]     = 0;
      phoConvdPhiAtCalo_[nPho_]    = 0;
      phoConvdEtaAtCalo_[nPho_]    = 0;
      phoConvTrkd0_[nPho_][0]      = 0;
      phoConvTrkd0_[nPho_][1]      = 0;
      phoConvTrkPin_[nPho_][0]     = 0;
      phoConvTrkPin_[nPho_][1]     = 0;
      phoConvTrkPout_[nPho_][0]    = 0;
      phoConvTrkPout_[nPho_][1]    = 0;
      phoConvTrkdz_[nPho_][0]      = 0;
      phoConvTrkdz_[nPho_][1]      = 0;
      phoConvTrkdzErr_[nPho_][0]   = 0;
      phoConvTrkdzErr_[nPho_][1]   = 0;
      phoConvChi2_[nPho_]          = 0;
      phoConvChi2Prob_[nPho_]      = 0;
      phoConvNTrks_[nPho_]         = 0;
      phoConvCharge_[nPho_][0]     = 0;
      phoConvCharge_[nPho_][1]     = 0;
      phoConvValidVtx_[nPho_]      = 0;
      phoConvLikeLihood_[nPho_]    = 0;
      for (int ind=0; ind<3; ind++) {
        phoConvP4_[nPho_][ind] = 0.;
        phoConvVtx_[nPho_][ind] = 0.;
        phoConvVtxErr_[nPho_][ind] = 0.;
        phoConvPairMomentum_[nPho_][ind] = 0.;
        phoConvRefittedMomentum_[nPho_][ind] = 0.;
      }
      phoConvP4_[nPho_][3] = 0.;

      if (iPho->hasConversionTracks()) {
 
        reco::ConversionRefVector conversions = iPho->conversions();
        reco::ConversionRef conv = conversions[0];
 
        phoConvValidVtx_[nPho_] = conv->conversionVertex().isValid();
 
        for (unsigned int j=0; j<conversions.size(); ++j) {
          conv = conversions[j];
          reco::Vertex vtx=conv->conversionVertex();

          phoConvP4_[nPho_][0] = conv->refittedPair4Momentum().px();
          phoConvP4_[nPho_][1] = conv->refittedPair4Momentum().py();
          phoConvP4_[nPho_][2] = conv->refittedPair4Momentum().pz();
          phoConvP4_[nPho_][3] = conv->refittedPair4Momentum().energy();

          phoConvVtx_[nPho_][0] = vtx.x();
          phoConvVtx_[nPho_][1] = vtx.y();
          phoConvVtx_[nPho_][2] = vtx.z();
          phoConvVtxErr_[nPho_][0]= vtx.xError();
          phoConvVtxErr_[nPho_][1]= vtx.yError();
          phoConvVtxErr_[nPho_][2]= vtx.zError();
          phoConvPairMomentum_[nPho_][0] = conv->pairMomentum().x();
          phoConvPairMomentum_[nPho_][1] = conv->pairMomentum().y();
          phoConvPairMomentum_[nPho_][2] = conv->pairMomentum().z();
          phoConvRefittedMomentum_[nPho_][0] = conv->refittedPairMomentum().x();
          phoConvRefittedMomentum_[nPho_][1] = conv->refittedPairMomentum().y();
          phoConvRefittedMomentum_[nPho_][2] = conv->refittedPairMomentum().z();
 
          phoConvChi2_[nPho_]       = vtx.chi2();
          phoConvChi2Prob_[nPho_]   = ChiSquaredProbability(vtx.chi2(), vtx.ndof());
          phoConvNTrks_[nPho_]      = conv->nTracks();
          phoConvLikeLihood_[nPho_] = conv->MVAout();

          const std::vector<edm::RefToBase<reco::Track> > tracks = conv->tracks();
          for (unsigned int k=0; k<tracks.size(); ++k) {
            if (k<=1) {
              phoConvTrkdz_[nPho_][k] = tracks[k]->dz();
              phoConvTrkdzErr_[nPho_][k] = tracks[k]->dzError();
              phoConvCharge_[nPho_][k] = tracks[k]->charge();
            }
          }
 
        }
 
        phoConvInvMass_[nPho_]       = conv->pairInvariantMass();
        phoConvCotTheta_[nPho_]      = conv->pairCotThetaSeparation();
        phoConvEoverP_[nPho_]        = conv->EoverPrefittedTracks();
        phoConvZofPVfromTrks_[nPho_] = conv->zOfPrimaryVertexFromTracks();
        phoConvMinDist_[nPho_]       = conv->distOfMinimumApproach();
        phoConvdPhiAtVtx_[nPho_]     = conv->dPhiTracksAtVtx();
	phoConvdPhiAtCalo_[nPho_]    = conv->dPhiTracksAtEcal();
	phoConvdEtaAtCalo_[nPho_]    = conv->dEtaTracksAtEcal();

        if (conv->tracks().size() > 0) {
          phoConvTrkd0_[nPho_][0]   = conv->tracksSigned_d0()[0];
          phoConvTrkPin_[nPho_][0]  = sqrt(conv->tracksPin()[0].Mag2());
          phoConvTrkPout_[nPho_][0] = sqrt(conv->tracksPout()[0].Mag2());
        }
 
        if (conv->tracks().size() > 1) {
          phoConvTrkd0_[nPho_][1]   = conv->tracksSigned_d0()[1];
          phoConvTrkPin_[nPho_][1]  = sqrt(conv->tracksPin()[1].Mag2());
          phoConvTrkPout_[nPho_][1] = sqrt(conv->tracksPout()[1].Mag2());
        }
 
      }
      //do the same for Single Legs:
      SingleLegConv_[nPho_]=0;
      reco::PhotonCollection::const_iterator pfphot=pfPhoTranslator_->begin();
      for(;pfphot!=pfPhoTranslator_->end();++pfphot){
	
	if(recophoRef->superCluster()!=pfphot->superCluster())continue;
	//	cout<<"MATCHED PHOTON "<<endl;
	reco::ConversionRefVector SLconversions=pfphot->conversionsOneLeg();
	SingleLegConv_[nPho_]=SLconversions.size();
	for(unsigned int SL=0; SL<SLconversions.size(); ++SL){
	  const std::vector<edm::RefToBase<reco::Track> > tracks = SLconversions[SL]->tracks();
	  phoPFConvVtx_[SL][0]=SLconversions[SL]->conversionVertex().x();
	  phoPFConvVtx_[SL][1]=SLconversions[SL]->conversionVertex().y();
	  phoPFConvVtx_[SL][2]=SLconversions[SL]->conversionVertex().z();
	  /*
	  phoPFConvRefitMom_[SL][0]=SLconversions[SL]->refittedPair4Momentum.px();
	  phoPFConvRefitMom_[SL][1]=SLconversions[SL]->refittedPair4Momentum.py();
	  phoPFConvRefitMom_[SL][2]=SLconversions[SL]->refittedPair4Momentum.pz();
	  */
	  phoPFConvMom_[SL][0]=tracks[0]->px();
	  phoPFConvMom_[SL][0]=tracks[0]->py();
	  phoPFConvMom_[SL][0]=tracks[0]->pz();
	}
      }
      
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
      
      MustacheEin_[nPho_]=-9999;
      MustacheEOut_[nPho_]=-9999;
      MustacheEtOut_[nPho_]=-9999;
      PFLowestClustE_[nPho_]=-9999;
      PFClustdEta_[nPho_]=-9999;
      PFClustdPhi_[nPho_]=-9999;
      PFClustRMSPhi_[nPho_]=-9999;
      PFClustRMSPhiMust_[nPho_]=-9999;
      pho_pfconvVtxZ_[nPho_]=-9999;
      pho_pfconvVtxZErr_[nPho_]=-9999;
      PFClustEneCorr_[nPho_]=-9999;
      pho_hasConvPf_[nPho_] = 0;
      pho_hasSLConvPf_[nPho_] = 0;
      PFEleMatch_[nPho_]=0;
      PFEleVeto_[nPho_]=0;
      pho_nSLConv_[nPho_]=0;

      if(ggPFPhoton.MatchPFReco()){
	PFRecoMatch_[nPho_]=1;

	//ggPFPhoton.fillPFClusters();

	if(ggPFPhoton.isConv()){
	  pho_hasConvPf_[nPho_] = 1;
	  
	  std::pair<float, float>VertexZ=ggPFPhoton.SLPoint();
	  pho_pfconvVtxZ_[nPho_] = VertexZ.first;
	  pho_pfconvVtxZErr_[nPho_] = VertexZ.second;
	}
	
	if(ggPFPhoton.hasSLConv()){
	  pho_hasSLConvPf_[nPho_] = 1;
	  reco::ConversionRefVector SLconversions;
	  std::vector<float> Zint;

	  if(SLPoint(*recophoRef, SLconversions, Zint)){
	    for(unsigned int z=0; z<SLconversions.size(); ++z){
	      pho_nSLConv_[nPho_]=z+1;
	      pho_pfSLConvVtxZ_[nPho_][z]=Zint[z];
	      pho_pfSLConvPos_[nPho_][z][0]=SLconversions[z]->conversionVertex().x();
	      pho_pfSLConvPos_[nPho_][z][1]=SLconversions[z]->conversionVertex().y();
	      pho_pfSLConvPos_[nPho_][z][2]=SLconversions[z]->conversionVertex().z();
	      
	    }
	  }
	}
	
	if(ggPFPhoton.isPFEle())PFEleMatch_[nPho_]=1;
	
	if(ggPFPhoton.PFElectronVeto(convH_, gsfElectronHandle_))PFEleVeto_[nPho_]=1;
	MustacheEin_[nPho_]=ggPFPhoton.MustE();
	MustacheEOut_[nPho_]=ggPFPhoton.MustEOut();
	MustacheEtOut_[nPho_]=ggPFPhoton.MustEtOut();
	PFLowestClustE_[nPho_]=ggPFPhoton.PFLowE();
	PFClustdEta_[nPho_]=ggPFPhoton.PFdEta();
	PFClustdPhi_[nPho_]=ggPFPhoton.PFdPhi();
	PFClustRMSPhi_[nPho_]=ggPFPhoton.PFClusRMSTot();
	PFClustRMSPhiMust_[nPho_]=ggPFPhoton.PFClusRMSMust();	
	std::vector<reco::CaloCluster>PFC=ggPFPhoton.PFClusters();
	//PFClustEneCorr_[nPho_]=ggPFPhoton.getPFPhoECorr(PFC, ReaderLCEB_, ReaderLCEE_);
      }
      else{	
	PFRecoMatch_[nPho_]=0;	
	ggPFPhoton.recoPhotonClusterLink(*iPho, pfCandidates);  
	MustacheEin_[nPho_]=ggPFPhoton.MustE();
	MustacheEOut_[nPho_]=ggPFPhoton.MustEOut();
	MustacheEtOut_[nPho_]=ggPFPhoton.MustEtOut();
	PFLowestClustE_[nPho_]=ggPFPhoton.PFLowE();
	PFClustdEta_[nPho_]=ggPFPhoton.PFdEta();
	PFClustdPhi_[nPho_]=ggPFPhoton.PFdPhi();
	PFClustRMSPhi_[nPho_]=ggPFPhoton.PFClusRMSTot();
	PFClustRMSPhiMust_[nPho_]=ggPFPhoton.PFClusRMSMust();	
      }

      //------------ Applying corrections ----------      
      if (develop_) {
	float absEta = fabs(phoSCEta_[nPho_]);
	float energyWithEtaCorrection;
	float pTWithEtaCorrection;
	
	// C(eta) for EB only
	if (iPho->isEB()){
	  energyWithEtaCorrection = phoSCRawE_[nPho_]/fcorrs::f5x5((int)(absEta*(5/0.087)));
	  pTWithEtaCorrection     = energyWithEtaCorrection / cosh(phoSCEta_[nPho_]);
	  
	  phoCetaCorrE_[nPho_]  = energyWithEtaCorrection;
	  phoCetaCorrEt_[nPho_] = pTWithEtaCorrection;
	  
	  // f(brem)-corrected energies
	  phoBremCorrE_[nPho_]  = fcorrs::fBrem(phoSCBrem_[nPho_], phoCetaCorrE_[nPho_]);
	  phoBremCorrEt_[nPho_] = phoBremCorrE_[nPho_] / cosh(phoSCEta_[nPho_]);
	  
	  // fully-corrected energies
	  phoFullCorrEt_[nPho_] = fcorrs::fullCorr(phoBremCorrEt_[nPho_], absEta);
	  phoFullCorrE_[nPho_]  = phoFullCorrEt_[nPho_] * cosh(phoSCEta_[nPho_]); 
	}
	else{
	  energyWithEtaCorrection = phoSCRawE_[nPho_] + (*iPho).superCluster()->preshowerEnergy();
	  pTWithEtaCorrection     = energyWithEtaCorrection / cosh(phoSCEta_[nPho_]);
	  
	  phoCetaCorrE_[nPho_]  = energyWithEtaCorrection;
	  phoCetaCorrEt_[nPho_] = pTWithEtaCorrection;
	  
	  // f(brem)-corrected energies
	  phoBremCorrE_[nPho_]  = fcorrs::fBrem_ee(phoSCBrem_[nPho_], phoCetaCorrE_[nPho_]);
	  phoBremCorrEt_[nPho_] = phoBremCorrE_[nPho_] / cosh(phoSCEta_[nPho_]);
	  
	  // fully-corrected energies
	  phoFullCorrEt_[nPho_] = fcorrs::fullCorr_ee(phoBremCorrEt_[nPho_], absEta);
	  phoFullCorrE_[nPho_]  = phoFullCorrEt_[nPho_] * cosh(phoSCEta_[nPho_]);
	}
      }
      //--------------------------------------------
      
      for (int i=0; i<2; i++) phoESDetId_[nPho_][i] = 0;
      for (int i=0; i<3; ++i) 
	for (int j=0; j<62; ++j) 
	  phoESHits_[nPho_][i][j] = 0;
      for (int i=0; i<3; ++i) phoESEffSigmaRR_[nPho_][i] = 0.;
      if (dumpESClusterInfo_) {
	for (int i=0; i<2; i++) {
	  phoESE1_[nPho_][i] = 0.;
	  phoESE3_[nPho_][i] = 0.;
	  phoESE5_[nPho_][i] = 0.;
	  phoESE7_[nPho_][i] = 0.;
	  phoESE11_[nPho_][i] = 0.;
	  phoESE21_[nPho_][i] = 0.;
	}
      }

      if (ESRecHits_.isValid() && (fabs(phoSCEta_[nPho_]) > 1.6 && fabs(phoSCEta_[nPho_]) < 3)) {
        const GlobalPoint scPoint((*iPho).superCluster()->x(), (*iPho).superCluster()->y(), (*iPho).superCluster()->z());
        DetId extrapolatedESId1 = (dynamic_cast<const EcalPreshowerGeometry*>(geometry_p))->getClosestCellInPlane(scPoint, 1);
        DetId extrapolatedESId2 = (dynamic_cast<const EcalPreshowerGeometry*>(geometry_p))->getClosestCellInPlane(scPoint, 2);
        phoESDetId_[nPho_][0] = extrapolatedESId1.rawId();
        phoESDetId_[nPho_][1] = extrapolatedESId2.rawId();

        vector<float> phoESHits0 = getESHits((*iPho).superCluster()->x(), (*iPho).superCluster()->y(), (*iPho).superCluster()->z(), rechits_map_, geometry_p, topology_p, 0);
        vector<float> phoESHits1 = getESHits((*iPho).superCluster()->x(), (*iPho).superCluster()->y(), (*iPho).superCluster()->z(), rechits_map_, geometry_p, topology_p, 1);
        vector<float> phoESHits2 = getESHits((*iPho).superCluster()->x(), (*iPho).superCluster()->y(), (*iPho).superCluster()->z(), rechits_map_, geometry_p, topology_p, -1);
        for (int i=0; i<62; ++i) {
          phoESHits_[nPho_][0][i] = phoESHits0[i];
          phoESHits_[nPho_][1][i] = phoESHits1[i];
          phoESHits_[nPho_][2][i] = phoESHits2[i];
        }

        vector<float> phoESShape = getESEffSigmaRR(phoESHits0);
        phoESEffSigmaRR_[nPho_][0] = phoESShape[0];
        phoESEffSigmaRR_[nPho_][1] = phoESShape[1];
        phoESEffSigmaRR_[nPho_][2] = phoESShape[2];

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
      muTrg_[nMu_][0] = (muTrigRef1.isAvailable()) ? 1 : -99;
      muTrg_[nMu_][1] = (muTrigRef2.isAvailable()) ? 1 : -99;
      muTrg_[nMu_][2] = (muTrigRef3.isAvailable()) ? 1 : -99;
      muTrg_[nMu_][3] = (muTrigRef4.isAvailable()) ? 1 : -99;
      muTrg_[nMu_][4] = (muTrigRef5.isAvailable()) ? 1 : -99;
      muTrg_[nMu_][5] = (muTrigRef6.isAvailable()) ? 1 : -99;
      muTrg_[nMu_][6] = (muTrigRef7.isAvailable()) ? 1 : -99;
      muTrg_[nMu_][7] = (muTrigRef8.isAvailable()) ? 1 : -99;
      muTrg_[nMu_][8] = (muTrigRef9.isAvailable()) ? 1 : -99;
      muTrg_[nMu_][9] = (muTrigRef10.isAvailable()) ? 1 : -99;

      Int_t goodMuonTrack = 1;
      const reco::TrackRef trkr   = iMu->globalTrack();
      const reco::TrackRef inntrk = iMu->innerTrack();
      if (trkr.isNull()) {
	goodMuonTrack = 0;
        muD0_[nMu_]   = -99.;
	muDz_[nMu_]   = -99.;
        muD0GV_[nMu_] = -99.;
	muDzGV_[nMu_] = -99.;
	muChi2NDF_[nMu_] = -99.;
	muNumberOfValidTrkHits_[nMu_]     = -99;
	muNumberOfValidMuonHits_[nMu_]    = -99;
	muVtxGlb_[nMu_][0] = -99;
	muVtxGlb_[nMu_][1] = -99;
	muVtxGlb_[nMu_][2] = -99;
      } else {
        muD0_[nMu_]   = trkr->dxy(pv);
        muDz_[nMu_]   = trkr->dz(pv);
        muD0GV_[nMu_] = trkr->dxy(gv);
        muDzGV_[nMu_] = trkr->dz(gv);
  	muChi2NDF_[nMu_] = trkr->normalizedChi2();
        muNumberOfValidTrkHits_[nMu_]     = trkr->hitPattern().numberOfValidTrackerHits();
        muNumberOfValidMuonHits_[nMu_]    = trkr->hitPattern().numberOfValidMuonHits();
	muVtxGlb_[nMu_][0] = trkr->vx();
	muVtxGlb_[nMu_][1] = trkr->vy();
	muVtxGlb_[nMu_][2] = trkr->vz();
      }

      if (inntrk.isNull()) {
	goodMuonTrack      = 0;
        muInnerD0_[nMu_]   = -99.;
	muInnerDz_[nMu_]   = -99.;
        muInnerD0GV_[nMu_] = -99.;
	muInnerDzGV_[nMu_] = -99.;
	muInnerChi2NDF_[nMu_] = -99.;
	muInnerPt_[nMu_]      = -99.;
	muInnerPtErr_[nMu_]   = -99.;
	muNumberOfValidTrkLayers_[nMu_]   = -99;
	muNumberOfValidPixelLayers_[nMu_] = -99;
	muNumberOfValidPixelHits_[nMu_]   = -99;
      } else {
        muInnerD0_[nMu_]   = inntrk->dxy(pv);
        muInnerDz_[nMu_]   = inntrk->dz(pv);
        muInnerD0GV_[nMu_] = inntrk->dxy(gv);
        muInnerDzGV_[nMu_] = inntrk->dz(gv);
  	muInnerChi2NDF_[nMu_] = inntrk->normalizedChi2();
	muInnerPt_[nMu_]      = inntrk->pt();
	muInnerPtErr_[nMu_]   = inntrk->ptError();
	muNumberOfValidTrkLayers_[nMu_]   = inntrk->hitPattern().trackerLayersWithMeasurement();
	muNumberOfValidPixelLayers_[nMu_] = inntrk->hitPattern().pixelLayersWithMeasurement();
        muNumberOfValidPixelHits_[nMu_]   = inntrk->hitPattern().numberOfValidPixelHits();
	math::XYZPoint vtxPoint[100];
	for (size_t iv=0; iv<recVtxsBS_->size(); ++iv) {   
	  vtxPoint[iv] = math::XYZPoint(vtxbs_[iv][0], vtxbs_[iv][1], vtxbs_[iv][2]);	 
	  muD0Vtx_[nMu_][iv] = inntrk->dxy(vtxPoint[iv]);
	  muDzVtx_[nMu_][iv] = inntrk->dz(vtxPoint[iv]);	  
	}
      }

      muStations_[nMu_] = iMu->numberOfMatchedStations();;
      muChambers_[nMu_] = iMu->numberOfMatches();

      // IP3D
      muIP3D_[nMu_]    = -999.;
      muIP3DErr_[nMu_] = -99.;
      muVtx_[nMu_][0]  = -99; 
      muVtx_[nMu_][1]  = -99; 
      muVtx_[nMu_][2]  = -99; 
      if (! iMu->track().isNull()) {
	const double musign = ((- iMu->track()->dxy(pv)) >= 0) ? 1. : -1.;
	const reco::TransientTrack tt_mu = thebuilder.build(iMu->track());
	const std::pair<bool,Measurement1D> ip3dpv_mu = IPTools::absoluteImpactParameter3D(tt_mu, *recVtxs_->begin());
	if (ip3dpv_mu.first) {
	  muIP3D_[nMu_]    = musign*ip3dpv_mu.second.value();
	  muIP3DErr_[nMu_] = ip3dpv_mu.second.error();
	}
	muVtx_[nMu_][0] = iMu->track()->vx();
	muVtx_[nMu_][1] = iMu->track()->vy();
	muVtx_[nMu_][2] = iMu->track()->vz();
      }

      muType_[nMu_]   = iMu->type();
      muEta_[nMu_]    = iMu->eta();
      muPhi_[nMu_]    = iMu->phi();
      muCharge_[nMu_] = iMu->charge();
      muPt_[nMu_]     = iMu->pt();
      muPz_[nMu_]     = iMu->pz();

      // muon cocktail
      edm::Ptr<reco::Candidate> recoMuRef = iMu->originalObjectRef();
      const reco::Muon *recoMu = dynamic_cast<const reco::Muon *>(recoMuRef.get());

      if (iMu->tpfmsTrack().isNull()) goodMuonTrack = 0;
      if (iMu->pickyTrack().isNull()) goodMuonTrack = 0;
      reco::TrackRef cktTrack;
      if (goodMuonTrack == 1) cktTrack = (muon::tevOptimized(*recoMu, 200, 17., 40., 0.25)).first;
      mucktPt_[nMu_]    = (! cktTrack.isNull()) ? cktTrack->pt() : -99;
      mucktPtErr_[nMu_] = (! cktTrack.isNull()) ? cktTrack->ptError() : -99;
      mucktEta_[nMu_]   = (! cktTrack.isNull()) ? cktTrack->eta() : -99;
      mucktPhi_[nMu_]   = (! cktTrack.isNull()) ? cktTrack->phi() : -99;
      mucktdxy_[nMu_]   = (! cktTrack.isNull()) ? cktTrack->dxy(pv) : -99;
      mucktdz_[nMu_]    = (! cktTrack.isNull()) ? cktTrack->dz(pv) : -99;

      muIsoTrk_[nMu_]  = iMu->trackIso();
      muIsoCalo_[nMu_] = iMu->caloIso();
      muIsoEcal_[nMu_] = iMu->ecalIso();
      muIsoHcal_[nMu_] = iMu->hcalIso();

      muPFIsoR03_CH_[nMu_]    = iMu->pfIsolationR03().sumChargedHadronPt;
      muPFIsoR03_NH_[nMu_]    = iMu->pfIsolationR03().sumNeutralHadronEt;
      muPFIsoR03_Pho_[nMu_]   = iMu->pfIsolationR03().sumPhotonEt;
      
      muPFIsoR03_PU_[nMu_]    = iMu->pfIsolationR03().sumPUPt;
      muPFIsoR03_CPart_[nMu_] = iMu->pfIsolationR03().sumChargedParticlePt;
      muPFIsoR03_NHHT_[nMu_]  = iMu->pfIsolationR03().sumNeutralHadronEtHighThreshold;
      muPFIsoR03_PhoHT_[nMu_] = iMu->pfIsolationR03().sumPhotonEtHighThreshold; 

      muPFIsoR04_CH_[nMu_]    = iMu->pfIsolationR04().sumChargedHadronPt;
      muPFIsoR04_NH_[nMu_]    = iMu->pfIsolationR04().sumNeutralHadronEt;
      muPFIsoR04_Pho_[nMu_]   = iMu->pfIsolationR04().sumPhotonEt;
      
      muPFIsoR04_PU_[nMu_]    = iMu->pfIsolationR04().sumPUPt;
      muPFIsoR04_CPart_[nMu_] = iMu->pfIsolationR04().sumChargedParticlePt;
      muPFIsoR04_NHHT_[nMu_]  = iMu->pfIsolationR04().sumNeutralHadronEtHighThreshold;
      muPFIsoR04_PhoHT_[nMu_] = iMu->pfIsolationR04().sumPhotonEtHighThreshold; 

      muGenIndex_[nMu_] = -1;
      int MuGenIndex = 0;
      if (!isData_) {
        if ( (*iMu).genLepton() && genParticlesHandle_.isValid() ) {
          if (fabs((*iMu).genLepton()->pdgId())==13) {
            for (vector<GenParticle>::const_iterator iGen = genParticlesHandle_->begin(); iGen != genParticlesHandle_->end(); ++iGen) {

              if (iGen->p4() == (*iMu).genLepton()->p4() && iGen->pdgId() == (*iMu).genLepton()->pdgId() && iGen->status() == (*iMu).genLepton()->status()) 
	        muGenIndex_[nMu_] = MuGenIndex;

              MuGenIndex++;
            }
          }
        }
      }
      nMu_++;
    }  // iMu loop
  } // muonHandle is valid
 
  nPFPho_  = 0;
  nPFEle_  = 0;
  nPFchad_ = 0;
  nPFnhad_ = 0;    

  if (!develop_) { // loop all PF candidates

    for (unsigned iCand=0; iCand< pfAllCandidates->size(); ++iCand) {
      
      const reco::PFCandidate  & pfParticleT((*pfAllCandidates)[iCand]);
      const reco::PFCandidate  *pfParticle = (&pfParticleT);
      
      if (pfParticle->pdgId() == 22 && pfParticle->pt() > 2 && fabs(pfParticle->eta()) < 2.5) {
	
	PFPhoE_[nPFPho_]    = pfParticle->energy();
	PFPhoEt_[nPFPho_]   = pfParticle->pt();
	PFPhoEta_[nPFPho_]  = pfParticle->eta();
	PFPhoPhi_[nPFPho_]  = pfParticle->phi();
	PFPhoIso_[nPFPho_]  = fsrPhotonIso03(pfAllCandidates.product(), pfCandidates.product(), *pfParticle);
	PFPhoType_[nPFPho_] = 1;

	nPFPho_++;
      }
      
      if (abs(pfParticle->pdgId()) == 13 && pfParticle->muonRef().isNonnull() && pfParticle->ecalEnergy() > 2) {
	if ((pfParticle->ecalEnergy() * pfParticle->pt() / pfParticle-> energy()) > 2) {
	  
	  PFPhoE_[nPFPho_]    = pfParticle->ecalEnergy();
	  PFPhoEt_[nPFPho_]   = pfParticle->ecalEnergy() * pfParticle->pt() / pfParticle->energy();
	  PFPhoEta_[nPFPho_]  = pfParticle->eta();
	  PFPhoPhi_[nPFPho_]  = pfParticle->phi();
	  PFPhoIso_[nPFPho_]  = fsrPhotonIso03(pfAllCandidates.product(), pfCandidates.product(), *pfParticle);
	  PFPhoType_[nPFPho_] = 2;
	  
	  nPFPho_++;
	}
      }
    }

    /*
    for (int ii=0; ii<nMu_; ++ii) {
      cout<<"iMu : "<<ii<<" ==="<<endl;
      double myPhoIso04 = 0;
      double myPhoIso03 = 0;
      for (int jj=0; jj<nPFPho_; ++jj) {
	if (PFPhoType_[jj] != 1) continue;
	if (deltaR(PFPhoEta_[jj], PFPhoPhi_[jj], muEta_[ii], muPhi_[ii]) < 0.4) {
	  if (deltaR(PFPhoEta_[jj], PFPhoPhi_[jj], muEta_[ii], muPhi_[ii]) > 0.01) {
	    if (PFPhoEt_[jj] > 0.5) {
	      myPhoIso04 += PFPhoEt_[jj];
	      if (deltaR(PFPhoEta_[jj], PFPhoPhi_[jj], muEta_[ii], muPhi_[ii]) < 0.3) myPhoIso03 += PFPhoEt_[jj];
	    }
	  }
	}
      }
      cout<<myPhoIso04<<" "<<muPFIsoR04_Pho_[ii]<<" "<<myPhoIso03<<" "<<muPFIsoR03_Pho_[ii]<<endl;
    }
    */

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
	PFPhoE_[nPFPho_]    = pfParticle->energy();
	PFPhoEt_[nPFPho_]   = pfParticle->pt();
	PFPhoEta_[nPFPho_]  = pfParticle->eta();
	PFPhoPhi_[nPFPho_]  = pfParticle->phi();
	PFPhoType_[nPFPho_] = 1;
      
	PFCandidate::ElementsInBlocks eleInBlocks = pfParticle->elementsInBlocks();
	nPFPho_tks=0;
	nPFPhoTks_[nPFPho_]=nPFPho_tks;
	nPFPho_clust=0;
	nPFPho_ES1clust=0;
	nPFPho_ES2clust=0;
	fillBlockInfo(*pfParticle, true, false,geomBar, geomEnd, geometry_p); 
	
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
      
      PFPhoE_[nPFPho_]=matchPFEle[pfe].energy();
      PFPhoEt_[nPFPho_]=matchPFEle[pfe].pt();
      PFPhoEta_[nPFPho_]=matchPFEle[pfe].eta();
      PFPhoPhi_[nPFPho_]=matchPFEle[pfe].phi();
      
      PFCandidate::ElementsInBlocks eleInBlocks = matchPFEle[pfe].elementsInBlocks();
      nPFPho_tks=0;
      nPFPhoTks_[nPFPho_]=nPFPho_tks;
      nPFPho_clust=0;
      nPFPho_ES1clust=0;
      nPFPho_ES2clust=0;
      //cout << endl << "Called by photons (fill)";
      if (develop_) fillBlockInfo(matchPFEle[pfe], true,true, geomBar, geomEnd, geometry_p);
      
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
    
    for (PFCandidateCollection::const_iterator pfParticle =pfCandidates->begin(); pfParticle!=pfCandidates->end(); pfParticle++) {
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
        jetLowPtEn_[nLowPtJet_]       = iJet->energy();
        jetLowPtPt_[nLowPtJet_]       = iJet->pt();
        jetLowPtEta_[nLowPtJet_]      = iJet->eta();
        jetLowPtPhi_[nLowPtJet_]      = iJet->phi();
        jetLowPtCharge_[nLowPtJet_]   = iJet->jetCharge();
        jetLowPtEt_[nLowPtJet_]       = iJet->et();
        jetLowPtRawPt_[nLowPtJet_]    = (*iJet).correctedJet("Uncorrected").pt();
        jetLowPtRawEn_[nLowPtJet_]    = (*iJet).correctedJet("Uncorrected").energy();
        jetLowPtArea_[nLowPtJet_]     = iJet->jetArea();
        jetLowPtPartonID_[nLowPtJet_] = iJet->partonFlavour();
  
        jetLowPtGenPartonID_[nLowPtJet_]    = -99;
        if (!isData_ && genParticlesHandle_.isValid() ) {
          if ((*iJet).genParton()) {
            jetLowPtGenPartonID_[nLowPtJet_] = (*iJet).genParton()->pdgId();
            jetLowPtGenEn_[nLowPtJet_]       = (*iJet).genParton()->energy();
            jetLowPtGenPt_[nLowPtJet_]       = (*iJet).genParton()->pt();
            jetLowPtGenEta_[nLowPtJet_]      = (*iJet).genParton()->eta();
            jetLowPtGenPhi_[nLowPtJet_]      = (*iJet).genParton()->phi();
          }
        }

        jetLowPtGenJetEn_[nJet_]    = -1; 
        jetLowPtGenJetPt_[nJet_]    = -999;
        jetLowPtGenJetEta_[nJet_]   = -999;
        jetLowPtGenJetPhi_[nJet_]   = -999;
        if (!isData_ && genParticlesHandle_.isValid() ) {
          if ((*iJet).genJet()) {
            jetLowPtGenJetEn_[nJet_]  = (*iJet).genJet()->energy();
            jetLowPtGenJetPt_[nJet_]  = (*iJet).genJet()->pt();
            jetLowPtGenJetEta_[nJet_] = (*iJet).genJet()->eta();
            jetLowPtGenJetPhi_[nJet_] = (*iJet).genJet()->phi();
          }
        }
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

      jetTrg_[nJet_][0] = (jetTrigRef1.isAvailable()) ? 1 : -99;
      jetTrg_[nJet_][1] = (jetTrigRef2.isAvailable()) ? 1 : -99;
      jetTrg_[nJet_][2] = (jetTrigRef3.isAvailable()) ? 1 : -99;
      jetTrg_[nJet_][3] = (jetTrigRef4.isAvailable()) ? 1 : -99;
      jetTrg_[nJet_][4] = (jetTrigRef5.isAvailable()) ? 1 : -99;
      jetTrg_[nJet_][5] = (jetTrigRef6.isAvailable()) ? 1 : -99;
      jetTrg_[nJet_][6] = (jetTrigRef7.isAvailable()) ? 1 : -99;
      jetTrg_[nJet_][7] = (jetTrigRef8.isAvailable()) ? 1 : -99;
      jetTrg_[nJet_][8] = (jetTrigRef9.isAvailable()) ? 1 : -99;
      for (int i=9; i<14; ++i) jetTrg_[nJet_][i] = -99;

      jetEn_[nJet_]     = iJet->energy();
      jetPt_[nJet_]     = iJet->pt();
      jetEta_[nJet_]    = iJet->eta();
      jetPhi_[nJet_]    = iJet->phi();
      jetCharge_[nJet_] = iJet->jetCharge();
      jetEt_[nJet_]     = iJet->et();
      jetMt_[nJet_]     = iJet->mt();
      jetRawPt_[nJet_]  = (*iJet).correctedJet("Uncorrected").pt();
      jetRawEn_[nJet_]  = (*iJet).correctedJet("Uncorrected").energy();

      jetArea_[nJet_]   = iJet->jetArea();
      jetCEF_[nJet_]    = iJet->chargedEmEnergyFraction();
      jetNEF_[nJet_]    = iJet->neutralEmEnergyFraction();
      jetCHF_[nJet_]    = iJet->chargedHadronEnergyFraction();
      jetNHF_[nJet_]    = iJet->neutralHadronEnergyFraction();
      jetHFHAE_[nJet_]  = iJet->HFHadronEnergy();
      jetHFEME_[nJet_]  = iJet->HFEMEnergy();
      jetNCH_[nJet_]    = iJet->chargedMultiplicity();
      jetNConstituents_[nJet_] = iJet->getPFConstituents().size();

      if (fabs(iJet->eta()) < 5.2) {
	jecUnc->setJetEta(iJet->eta());
	jecUnc->setJetPt(iJet->pt()); // here you must use the CORRECTED jet pt
	jetJECUnc_[nJet_] = jecUnc->getUncertainty(true);
      } else {
	jetJECUnc_[nJet_] = -1.; 
      }

      float PI = 3.1415927;
      float dphiTemp=fabs(pfMETPhi_-iJet->phi());
      if(dphiTemp>PI) dphiTemp = 2*PI-dphiTemp;
      jetDPhiMETJet_[nJet_] = dphiTemp;

      jetLeadTrackPt_[nJet_]=-9999.;
      if (iJet->isPFJet() == true) {
	std::vector <reco::PFCandidatePtr> constituents = iJet->getPFConstituents ();
	for (unsigned ic = 0; ic < constituents.size (); ++ic) {
	  if ( constituents[ic]->particleId() > 3 ) continue;
	  reco::TrackRef trackRef = constituents[ic]->trackRef();
	  if ( trackRef.isNonnull() ) { if(trackRef->pt() > jetLeadTrackPt_[nJet_]) jetLeadTrackPt_[nJet_]=trackRef->pt(); }
	}
      }

      jetVtxPt_[nJet_]=jetVtxMass_[nJet_]=jetVtx3dL_[nJet_]=jetVtx3deL_[nJet_]=-99;
      const reco::SecondaryVertexTagInfo * tf = iJet->tagInfoSecondaryVertex("secondaryVertex");
      if (tf){
      	math::XYZTLorentzVectorD vertexSum;
      	for(size_t vi=0;vi< tf->nVertices();vi++)
      	  {
      	    vertexSum+=tf->secondaryVertex(vi).p4();
      	  }
      	jetVtxPt_[nJet_] = vertexSum.Pt();
      	if (tf->nVertices() >0){
      	  jetVtxMass_[nJet_] =  tf->secondaryVertex(0).p4().mass();
      	  Measurement1D m = tf->flightDistance(0);
      	  jetVtx3dL_[nJet_] = m.value();
      	  jetVtx3deL_[nJet_] = m.error();
      	}
      }

      int isSemiLept=0;
      jetSoftLeptdR_[nJet_]=jetSoftLeptPt_[nJet_]=jetSoftLeptPtRel_[nJet_]=jetSoftLeptIdlooseMu_[nJet_]=jetSoftLeptIdEle95_[nJet_]=-99;
      edm::Handle<edm::View<reco::Candidate> > muonNoCutsHandle;
      e.getByLabel("muons",muonNoCutsHandle);
      edm::View<reco::Candidate> muonsNoCuts = *muonNoCutsHandle;
      for(edm::View<reco::Candidate>::const_iterator mu = muonsNoCuts.begin(); mu!=muonsNoCuts.end() && isSemiLept!=1; ++mu){
	const pat::Muon& m = static_cast <const pat::Muon&> (*mu);
	float Smpt = m.pt();
	float Smeta = m.eta();
	float Smphi = m.phi();
	float SmJdR = deltaR(Smeta, Smphi, iJet->eta(), iJet->phi());
	if   ( Smpt > 5 && SmJdR < 0.5) { //lep_ptCutForBjets_ = 5 GeV
	  isSemiLept=1;
	  jetSoftLeptdR_[nJet_]= SmJdR;
	  jetSoftLeptPt_[nJet_]=Smpt;
	  TLorentzVector jvec (iJet->p4().X(), iJet->p4().Y(), iJet->p4().Z(), iJet->p4().T());
	  TVector3 mvec ( m.p4().Vect().X(), m.p4().Vect().Y(), m.p4().Vect().Z()  );
	  jetSoftLeptPtRel_[nJet_]= jvec.Perp(  mvec );
	  jetSoftLeptIdlooseMu_[nJet_]=m.muonID("TMLastStationLoose");
	}
      }

      edm::Handle<edm::View<reco::Candidate> > eleNoCutsHandle;
      e.getByLabel("gsfElectrons",eleNoCutsHandle);
      edm::View<reco::Candidate> elesNoCuts = *eleNoCutsHandle;
      for(edm::View<reco::Candidate>::const_iterator ele = elesNoCuts.begin(); ele!=elesNoCuts.end() && isSemiLept!=1; ++ele){
        const pat::Electron& e = static_cast <const pat::Electron&> (*ele);
        float Smpt = e.pt();
        float Smeta = e.eta();
        float Smphi = e.phi();
        float SmJdR = deltaR(Smeta, Smphi, iJet->eta(), iJet->phi());
        if   ( Smpt> 5 && SmJdR <0.5) { //lep_ptCutForBjets_ = 5 GeV
	  isSemiLept=1;
	  jetSoftLeptdR_[nJet_]= SmJdR;
	  jetSoftLeptPt_[nJet_]=Smpt;
	  TLorentzVector jvec (iJet->p4().X(), iJet->p4().Y(), iJet->p4().Z(), iJet->p4().T());
	  TVector3 mvec ( e.p4().Vect().X(), e.p4().Vect().Y(), e.p4().Vect().Z()  );
	  jetSoftLeptPtRel_[nJet_]= jvec.Perp(  mvec );
	  //jetSoftLeptIdEle95_[nJet_]=e.electronID("eidVBTFCom95");
	  if (
	      ( fabs(Smeta)<2.5 && !( fabs(Smeta)>1.4442 && fabs(Smeta)<1.566))  &&
	      (( fabs(Smeta)>1.566  && e.sigmaIetaIeta()<0.01 && fabs(e.deltaPhiSuperClusterTrackAtVtx())<0.8 && fabs(e.deltaEtaSuperClusterTrackAtVtx())<0.007 ) ||
	       ( fabs(Smeta)<1.4442 && e.sigmaIetaIeta()<0.03 && fabs(e.deltaPhiSuperClusterTrackAtVtx())<0.7 && fabs(e.deltaEtaSuperClusterTrackAtVtx())<0.01 ))     )
	    jetSoftLeptIdEle95_[nJet_]=1;
	}
      }

      // b-tagging
      jetCombinedSecondaryVtxBJetTags_[nJet_]    = iJet->bDiscriminator("combinedSecondaryVertexBJetTags");   //   rob !!!!
      jetCombinedSecondaryVtxMVABJetTags_[nJet_] = iJet->bDiscriminator("combinedSecondaryVertexMVABJetTags");   //   rob !!!!
      jetJetProbabilityBJetTags_[nJet_]          = iJet->bDiscriminator("jetProbabilityBJetTags"); 
      jetJetBProbabilityBJetTags_[nJet_]         = iJet->bDiscriminator("jetBProbabilityBJetTags"); 
      jetTrackCountingHighPurBJetTags_[nJet_]    = iJet->bDiscriminator("trackCountingHighPurBJetTags"); 

      // betastar
      for (int iVTX = 0; iVTX < 100; iVTX++) jetBetaStar_[nJet_][iVTX] = 0.;
      for (int iVTX = 0; iVTX < nVtx_; iVTX++) {
        double tracks_x = 0.;
        double tracks_y = 0.;
        double tracks_x_tot = 0.;
        double tracks_y_tot = 0.;
        for (unsigned i = 0;  i <  iJet->numberOfDaughters (); i++) {	     
          const reco::PFCandidatePtr pfcand = iJet->getPFConstituent(i);
          reco::TrackRef trackref = pfcand->trackRef();
          if( trackref.isNonnull()) {
            // track_all
            tracks_x_tot += (trackref)->px();
            tracks_y_tot += (trackref)->py();
            for (int jVTX = 0; jVTX < nVtx_; jVTX++) {
              if (jVTX == iVTX) continue;
 
              if (fabs((trackref)->vz()-vtx_[jVTX][2]) < 0.1) {        	  
                tracks_x += (trackref)->px();
                tracks_y += (trackref)->py();	
                break;	
              } // track_PU
            } // non-PV loop (assume that iVTX is PV)
          }
        } // jet associated tracks

        if (tracks_x_tot!=0. || tracks_y_tot!=0.) jetBetaStar_[nJet_][iVTX] = sqrt(tracks_x*tracks_x+tracks_y*tracks_y)/sqrt(tracks_x_tot*tracks_x_tot+tracks_y_tot*tracks_y_tot);
      } // vtx loop

      //jet PF Loose ID
      pat::strbitset retjet = pfLooseId_.getBitTemplate();
      jetPFLooseId_[nJet_]  = pfLooseId_(*iJet, retjet);

      //PileupJet ID variables
      if (pujetIDalgos_.size()>0) {
        float jecSetup = 0.0;
        const reco::VertexCollection vertexCollection = *(recVtxsBS_.product());
        const reco::Vertex* selectedVtx  = &(*vertexCollection.begin());;
        const pat::Jet* thisjet = &(*iJet);
  
        PileupJetIdentifier jetMVAinputs = jetMVACalculator->computeIdVariables( thisjet, jecSetup, selectedVtx, vertexCollection);
  
        jetDRMean_[nJet_]  = jetMVAinputs.dRMean();
        jetDR2Mean_[nJet_] = jetMVAinputs.dR2Mean();
        jetDZ_[nJet_]      = jetMVAinputs.dZ();
        jetFrac01_[nJet_]  = jetMVAinputs.frac01();
        jetFrac02_[nJet_]  = jetMVAinputs.frac02();
        jetFrac03_[nJet_]  = jetMVAinputs.frac03();
        jetFrac04_[nJet_]  = jetMVAinputs.frac04();
        jetFrac05_[nJet_]  = jetMVAinputs.frac05();
        jetFrac06_[nJet_]  = jetMVAinputs.frac06();
        jetFrac07_[nJet_]  = jetMVAinputs.frac07();
        jetBeta_[nJet_]    = jetMVAinputs.beta();
        jetBetaStarCMG_[nJet_] = jetMVAinputs.betaStar();
        jetBetaStarClassic_[nJet_] = jetMVAinputs.betaStarClassic();
        jetNNeutrals_[nJet_] = jetMVAinputs.nNeutrals();
        jetNCharged_[nJet_]  = jetMVAinputs.nCharged();
  
	jetPuJetIdL_[nJet_]=jetPuJetIdM_[nJet_]=jetPuJetIdT_[nJet_]=0;
	for(unsigned int imva=0; imva<jetMVAAlgos_.size(); imva++){
          PileupJetIdAlgo* ialgo = (pujetIDalgos_[imva]);
          ialgo->set(jetMVAinputs);
          PileupJetIdentifier id = ialgo->computeMva();
          jetMVAs_[nJet_][imva] = id.mva();
          jetWPLevels_[nJet_][imva] = id.idFlag();

	  if(imva==1){
	    int    idflag = id.idFlag();
	    if( PileupJetIdentifier::passJetId( idflag, PileupJetIdentifier::kLoose )) {
	      jetPuJetIdL_[nJet_] =1;
	    }
	    if( PileupJetIdentifier::passJetId( idflag, PileupJetIdentifier::kMedium )) {
	      jetPuJetIdM_[nJet_] =1;
	    }
	    if( PileupJetIdentifier::passJetId( idflag, PileupJetIdentifier::kTight )) {
	      jetPuJetIdT_[nJet_] =1;
	    }
	  }
        }

        for (int iVtx = 0; iVtx < 100; iVtx++) {
          jetBetaExt_[nJet_][iVtx]            = -1.;
          jetBetaStarCMGExt_[nJet_][iVtx]     = -1.;
          jetBetaStarClassicExt_[nJet_][iVtx] = -1.;

          for(unsigned int imva=0; imva<jetMVAAlgos_.size(); imva++){
            jetMVAsExt_[nJet_][imva][iVtx]  = -3.;
            jetWPLevelsExt_[nJet_][imva][iVtx] = -1;
          }
        }
        for(size_t iVtx=0; iVtx<vertexCollection.size(); ++iVtx) {
          PileupJetIdentifier jetMVAinputsExt = jetMVACalculator->computeIdVariables( thisjet, jecSetup, &vertexCollection[iVtx], vertexCollection);
          jetBetaExt_[nJet_][iVtx]            = jetMVAinputsExt.beta();
          jetBetaStarCMGExt_[nJet_][iVtx]     = jetMVAinputsExt.betaStar();
          jetBetaStarClassicExt_[nJet_][iVtx] = jetMVAinputsExt.betaStarClassic();
          
          for(unsigned int imva=0; imva<jetMVAAlgos_.size(); imva++){
            PileupJetIdAlgo* ialgo = (pujetIDalgos_[imva]);
            ialgo->set(jetMVAinputsExt);
            PileupJetIdentifier id = ialgo->computeMva();
            jetMVAsExt_[nJet_][imva][iVtx]  = id.mva() ;
            jetWPLevelsExt_[nJet_][imva][iVtx] = id.idFlag();
          }
        }
      }

      // gen jet and parton
      jetPartonID_[nJet_] = iJet->partonFlavour();

      jetGenPartonID_[nJet_]    = -99;
      if (!isData_ && genParticlesHandle_.isValid() ) {
        if ((*iJet).genParton()) {
          jetGenPartonID_[nJet_] = (*iJet).genParton()->pdgId();
          jetGenEn_[nJet_]       = (*iJet).genParton()->energy();
          jetGenPt_[nJet_]       = (*iJet).genParton()->pt();
          jetGenEta_[nJet_]      = (*iJet).genParton()->eta();
          jetGenPhi_[nJet_]      = (*iJet).genParton()->phi();
	}
      }

      jetGenJetIndex_[nJet_] = -1;
      jetGenJetEn_[nJet_]    = -1;
      jetGenJetPt_[nJet_]    = -999;
      jetGenJetEta_[nJet_]   = -999;
      jetGenJetPhi_[nJet_]   = -999;
      
      if (!isData_ && genParticlesHandle_.isValid() ) {
	if ((*iJet).genJet()) {
	  jetGenJetIndex_[nJet_] = 1;
	  jetGenJetEn_[nJet_] = (*iJet).genJet()->energy();
	  jetGenJetPt_[nJet_] = (*iJet).genJet()->pt();
	  jetGenJetEta_[nJet_] = (*iJet).genJet()->eta();
	  jetGenJetPhi_[nJet_] = (*iJet).genJet()->phi();
	}
      }
      nJet_++;
    }

  nConv_=0;

  if (convH_.isValid()) {
    for( reco::ConversionCollection::const_iterator  iConv = convH_->begin(); iConv != convH_->end(); iConv++) {
      if (nConv_ >= 500) {
        cout << "WARNING :TOO MANY CONVERTED CANDIDATES: " << convH_->size() << endl;
        break;
      }
      reco::Conversion localConv = reco::Conversion(*iConv);

      convNTracks_[nConv_]=0;
      convPairInvMass_[nConv_]=-999.;
      convPairCotThetaSep_[nConv_]=-999.;
      convEoverP_[nConv_]=-999.;
      convDistOfMinApproach_[nConv_]=-999.;
      convDPhiTrksAtVtx_[nConv_]=-999.;
      convDPhiTrksAtEcal_[nConv_]=-999.;
      convDEtaTrksAtEcal_[nConv_]=-999.;
      convDxy_[nConv_]=-999.;
      convDz_[nConv_]=-999.;
      convLxy_[nConv_]=-999.;
      convLz_[nConv_]=-999.;
      convNSharedHits_[nConv_]=0;
      convZofPrimVtxFromTrks_[nConv_]=-999.;
      convTk1D0_[nConv_]=-999.;
      convTk1Pout_[nConv_]=-999.;
      convTk1Pin_[nConv_]=-999.;
      convTk2D0_[nConv_]=-999.;
      convTk2Pout_[nConv_]=-999.;
      convTk2Pin_[nConv_]=-999.;
      convValidVtx_[nConv_]=0;
      convMVALikelihood_[nConv_]=0.;
      convCh1Ch2_[nConv_] = -999.;
      convTk1Dz_[nConv_]=-999.;
      convTk1DzErr_[nConv_]=-999.;
      convTk2Dz_[nConv_]=-999.;
      convTk2DzErr_[nConv_]=-999.;
      for (int ind=0; ind<3; ind++) {
        convP4_[nConv_][ind] = -999.;
        convVtx_[nConv_][ind] = -999.;
        convVtxErr_[nConv_][ind]=-999.;
        convPairMomentum_[nConv_][ind] = -999.;
        convRefittedMomentum_[nConv_][ind] = -999.;
      }
      convP4_[nConv_][3]  = -999.;
      convNHitsBeforeVtx_[nConv_][0] = -999;
      convNHitsBeforeVtx_[nConv_][1] = -999;

      convP4_[nConv_][0] = localConv.refittedPair4Momentum().px();
      convP4_[nConv_][1] = localConv.refittedPair4Momentum().py();
      convP4_[nConv_][2] = localConv.refittedPair4Momentum().pz();
      convP4_[nConv_][3] = localConv.refittedPair4Momentum().energy();

      convValidVtx_[nConv_]=localConv.conversionVertex().isValid();
      if ( !localConv.conversionVertex().isValid() ) continue;
      reco::Vertex vtx=localConv.conversionVertex();

      convVtx_[nConv_][0] = vtx.x();
      convVtx_[nConv_][1] = vtx.y();
      convVtx_[nConv_][2] = vtx.z();
      convVtxErr_[nConv_][0]= vtx.xError();
      convVtxErr_[nConv_][1]= vtx.yError();
      convVtxErr_[nConv_][2]= vtx.zError();
      convPairMomentum_[nConv_][0] = localConv.pairMomentum().x();
      convPairMomentum_[nConv_][1] = localConv.pairMomentum().y();
      convPairMomentum_[nConv_][2] = localConv.pairMomentum().z();
      convRefittedMomentum_[nConv_][0] = localConv.refittedPairMomentum().x();
      convRefittedMomentum_[nConv_][1] = localConv.refittedPairMomentum().y();
      convRefittedMomentum_[nConv_][2] = localConv.refittedPairMomentum().z();

      convChi2_[nConv_]=vtx.chi2();
      convChi2Probability_[nConv_]=ChiSquaredProbability(vtx.chi2(), vtx.ndof());
      convNTracks_[nConv_]=localConv.nTracks();
      convMVALikelihood_[nConv_]=localConv.MVAout();

      if( localConv.nTracks()) {
        const std::vector<edm::RefToBase<reco::Track> > tracks = localConv.tracks();
        for (unsigned int i=0; i<tracks.size(); i++) {
          if(i==0) {
            convTk1Dz_[nConv_]=tracks[i]->dz();
            convTk1DzErr_[nConv_]=tracks[i]->dzError();
            convCh1Ch2_[nConv_]=tracks[i]->charge();
          }
          else if(i==1) {
            convTk2Dz_[nConv_]=tracks[i]->dz();
            convTk2DzErr_[nConv_]=tracks[i]->dzError();
            convCh1Ch2_[nConv_]*=tracks[i]->charge();
          }
        }
      }
  
      convPairInvMass_[nConv_]=localConv.pairInvariantMass();
      convPairCotThetaSep_[nConv_]=localConv.pairCotThetaSeparation();
      convEoverP_[nConv_]=localConv.EoverPrefittedTracks();
      //if (!runOn44X_){
      //convDPhiTrksAtEcal_[nConv_]=localConv.dPhiTracksAtEcal();
      //convDEtaTrksAtEcal_[nConv_]=localConv.dEtaTracksAtEcal();
      //}
      convZofPrimVtxFromTrks_[nConv_]=localConv.zOfPrimaryVertexFromTracks();
      convDistOfMinApproach_[nConv_]=localConv.distOfMinimumApproach();
      convDPhiTrksAtVtx_[nConv_]=localConv.dPhiTracksAtVtx();
      convDxy_[nConv_]=localConv.dxy();
      convDz_[nConv_]=localConv.dz();
      convLxy_[nConv_]=localConv.lxy();
      convLz_[nConv_]=localConv.lz();
  
      for (unsigned int i=0; i<localConv.nHitsBeforeVtx().size(); ++i) {
        if (i<=1) convNHitsBeforeVtx_[nConv_][i] = localConv.nHitsBeforeVtx()[i];
      }
      convNSharedHits_[nConv_] = localConv.nSharedHits();
  
      if(localConv.tracks().size() > 0) {
        convTk1D0_[nConv_]=localConv.tracksSigned_d0()[0];
        convTk1Pout_[nConv_]=sqrt(localConv.tracksPout()[0].Mag2());
        convTk1Pin_[nConv_]=sqrt(localConv.tracksPin()[0].Mag2());
      }
      
      if(localConv.tracks().size() > 1) {
          convTk2D0_[nConv_]=localConv.tracksSigned_d0()[1];
          convTk2Pout_[nConv_]=sqrt(localConv.tracksPout()[1].Mag2());
          convTk2Pin_[nConv_]=sqrt(localConv.tracksPin()[1].Mag2());
      }
      
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
    for (int i=0; i<31; ++i) esHits.push_back(0);
  } else {
    
    it = rechits_map.find(strip1);
    if (it->second.energy() > 1.0e-10 && it != rechits_map.end()) esHits.push_back(it->second.energy());  
    else esHits.push_back(0);
    //cout<<"center : "<<strip1<<" "<<it->second.energy()<<endl;      

    // east road 
    for (int i=0; i<15; ++i) {
      next = theESNav1.east();
      if (next != ESDetId(0)) {
	it = rechits_map.find(next);
	if (it->second.energy() > 1.0e-10 && it != rechits_map.end()) esHits.push_back(it->second.energy());  
	else esHits.push_back(0);
	//cout<<"east "<<i<<" : "<<next<<" "<<it->second.energy()<<endl;
      } else {
        for (int j=i; j<15; j++) esHits.push_back(0);
        break;
	//cout<<"east "<<i<<" : "<<next<<" "<<0<<endl;
      }
    }

    // west road 
    theESNav1.setHome(strip1);
    theESNav1.home();
    for (int i=0; i<15; ++i) {
      next = theESNav1.west();
      if (next != ESDetId(0)) {
	it = rechits_map.find(next);
	if (it->second.energy() > 1.0e-10 && it != rechits_map.end()) esHits.push_back(it->second.energy());  
	else esHits.push_back(0);
	//cout<<"west "<<i<<" : "<<next<<" "<<it->second.energy()<<endl;
      } else {
        for (int j=i; j<15; j++) esHits.push_back(0);
        break;
	//cout<<"west "<<i<<" : "<<next<<" "<<0<<endl;
      }
    }
  }

  if (strip2 == ESDetId(0)) {
    for (int i=0; i<31; ++i) esHits.push_back(0);
  } else {

    it = rechits_map.find(strip2);
    if (it->second.energy() > 1.0e-10 && it != rechits_map.end()) esHits.push_back(it->second.energy());  
    else esHits.push_back(0);
    //cout<<"center : "<<strip2<<" "<<it->second.energy()<<endl;      

    // north road 
    for (int i=0; i<15; ++i) {
      next = theESNav2.north();
      if (next != ESDetId(0)) {
	it = rechits_map.find(next);
	if (it->second.energy() > 1.0e-10 && it != rechits_map.end()) esHits.push_back(it->second.energy());
	else esHits.push_back(0);
	//cout<<"north "<<i<<" : "<<next<<" "<<it->second.energy()<<endl;  
      } else {
        for (int j=i; j<15; j++) esHits.push_back(0);
        break;
	//cout<<"north "<<i<<" : "<<next<<" "<<0<<endl;
      }
    }

    // south road 
    theESNav2.setHome(strip2);
    theESNav2.home();
    for (int i=0; i<15; ++i) {
      next = theESNav2.south();
      if (next != ESDetId(0)) {
	it = rechits_map.find(next);
	if (it->second.energy() > 1.0e-10 && it != rechits_map.end()) esHits.push_back(it->second.energy());  
	else esHits.push_back(0);
	//cout<<"south "<<i<<" : "<<next<<" "<<it->second.energy()<<endl;
      } else {
        for (int j=i; j<15; j++) esHits.push_back(0);
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

  for(int ibin=0; ibin<((nBIN+1)/2); ibin++) {
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
  for (int id_X=1; id_X<=21; id_X++) {
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
  for(int ibin=1; ibin<((nBIN+1)/2); ibin++) {
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
	  for (unsigned int t=0; t<tracks.size(); t++){
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
	  for (unsigned int t=0; t<tracks.size(); t++) {
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
