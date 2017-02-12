#include <TString.h>
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "ggAnalysis/ggNtuplizer/interface/GEDPhoIDTools.h"
#include "ggAnalysis/ggNtuplizer/interface/ggNtuplizer.h"

using namespace std;

Int_t          nPho_;
vector<float>  phoE_;
vector<float>  phoEt_;
vector<float>  phoEta_;
vector<float>  phoPhi_;
vector<float>  phoCalibE_;
vector<float>  phoCalibEt_;
vector<float>  phoSCE_;
vector<float>  phoSCRawE_;
vector<float>  phoESEn_;
vector<float>  phoESEnP1_;
vector<float>  phoESEnP2_;
vector<float>  phoSCEta_;
vector<float>  phoSCPhi_;
vector<float>  phoSCEtaWidth_;
vector<float>  phoSCPhiWidth_;
vector<float>  phoSCBrem_;
vector<int>    phohasPixelSeed_;
vector<int>    phoEleVeto_;
vector<float>  phoR9_;
vector<float>  phoHoverE_;
//vector<float>  phoSigmaIEtaIEta_;
//vector<float>  phoSigmaIEtaIPhi_;
//vector<float>  phoSigmaIPhiIPhi_;
vector<float>  phoE1x3_;
vector<float>  phoE1x5_;
vector<float>  phoE2x2_;
vector<float>  phoE2x5Max_;
vector<float>  phoE5x5_;
vector<float>  phoESEffSigmaRR_;
vector<float>  phoSigmaIEtaIEtaFull5x5_;
vector<float>  phoSigmaIEtaIPhiFull5x5_;
vector<float>  phoSigmaIPhiIPhiFull5x5_;
vector<float>  phoE1x3Full5x5_;
vector<float>  phoE1x5Full5x5_;
vector<float>  phoE2x2Full5x5_;
vector<float>  phoE2x5MaxFull5x5_;
vector<float>  phoE5x5Full5x5_;
vector<float>  phoR9Full5x5_;
vector<float>  phoPFChIso_;
vector<float>  phoPFPhoIso_;
vector<float>  phoPFNeuIso_;
vector<float>  phoPFChWorstIso_;
vector<float>  phoPFChIsoFrix1_;
vector<float>  phoPFChIsoFrix2_;
vector<float>  phoPFChIsoFrix3_;
vector<float>  phoPFChIsoFrix4_;
vector<float>  phoPFChIsoFrix5_;
vector<float>  phoPFChIsoFrix6_;
vector<float>  phoPFChIsoFrix7_;
vector<float>  phoPFChIsoFrix8_;
vector<float>  phoPFPhoIsoFrix1_;
vector<float>  phoPFPhoIsoFrix2_;
vector<float>  phoPFPhoIsoFrix3_;
vector<float>  phoPFPhoIsoFrix4_;
vector<float>  phoPFPhoIsoFrix5_;
vector<float>  phoPFPhoIsoFrix6_;
vector<float>  phoPFPhoIsoFrix7_;
vector<float>  phoPFPhoIsoFrix8_;
vector<float>  phoPFNeuIsoFrix1_;
vector<float>  phoPFNeuIsoFrix2_;
vector<float>  phoPFNeuIsoFrix3_;
vector<float>  phoPFNeuIsoFrix4_;
vector<float>  phoPFNeuIsoFrix5_;
vector<float>  phoPFNeuIsoFrix6_;
vector<float>  phoPFNeuIsoFrix7_;
vector<float>  phoPFNeuIsoFrix8_;
vector<float>  phoCITKChIso_;
vector<float>  phoCITKPhoIso_;
vector<float>  phoCITKNeuIso_;
//vector<float>  phoPUPPIChIso_;
//vector<float>  phoPUPPIPhoIso_;
//vector<float>  phoPUPPINeuIso_;
//vector<float>  phoSeedBCE_;
//vector<float>  phoSeedBCEta_;
vector<float>  phoIDMVA_;
vector<UInt_t> phoFiredSingleTrgs_;
vector<UInt_t> phoFiredDoubleTrgs_;
vector<UInt_t> phoFiredL1Trgs_;
//vector<float>  phoEcalRecHitSumEtConeDR03_;
//vector<float>  phohcalDepth1TowerSumEtConeDR03_;
//vector<float>  phohcalDepth2TowerSumEtConeDR03_;
//vector<float>  phohcalTowerSumEtConeDR03_;
//vector<float>  photrkSumPtHollowConeDR03_;
//vector<float>  photrkSumPtSolidConeDR03_;
vector<float>  phoSeedTime_;
vector<float>  phoSeedEnergy_;
//vector<float>  phoSeedTimeFull5x5_;
//vector<float>  phoMIPChi2_;
//vector<float>  phoMIPTotEnergy_;
//vector<float>  phoMIPSlope_;
//vector<float>  phoMIPIntercept_;
//vector<float>  phoMIPNhitCone_;
//vector<float>  phoMIPIsHalo_;

vector<UShort_t> phoxtalBits_;

vector<UShort_t> phoIDbit_;

//Necessary for the Photon Footprint removal
template <class T, class U>
bool isInFootprint(const T& thefootprint, const U& theCandidate) {
  for ( auto itr = thefootprint.begin(); itr != thefootprint.end(); ++itr ) {

    if( itr.key() == theCandidate.key() ) return true;
    
  }
  return false;
}

void ggNtuplizer::branchesPhotons(TTree* tree) {
  
  tree->Branch("nPho",                    &nPho_);
  tree->Branch("phoE",                    &phoE_);
  tree->Branch("phoEt",                   &phoEt_);
  tree->Branch("phoEta",                  &phoEta_);
  tree->Branch("phoPhi",                  &phoPhi_);
  tree->Branch("phoCalibE",               &phoCalibE_);
  tree->Branch("phoCalibEt",              &phoCalibEt_);
  tree->Branch("phoSCE",                  &phoSCE_);
  tree->Branch("phoSCRawE",               &phoSCRawE_);
  tree->Branch("phoESEn",                 &phoESEn_);
  tree->Branch("phoESEnP1",               &phoESEnP1_);
  tree->Branch("phoESEnP2",               &phoESEnP2_);
  tree->Branch("phoSCEta",                &phoSCEta_);
  tree->Branch("phoSCPhi",                &phoSCPhi_);
  tree->Branch("phoSCEtaWidth",           &phoSCEtaWidth_);
  tree->Branch("phoSCPhiWidth",           &phoSCPhiWidth_);
  tree->Branch("phoSCBrem",               &phoSCBrem_);
  tree->Branch("phohasPixelSeed",         &phohasPixelSeed_);
  tree->Branch("phoEleVeto",              &phoEleVeto_);
  tree->Branch("phoR9",                   &phoR9_);
  tree->Branch("phoHoverE",               &phoHoverE_);
  //tree->Branch("phoSigmaIEtaIEta",        &phoSigmaIEtaIEta_);
  //tree->Branch("phoSigmaIEtaIPhi",        &phoSigmaIEtaIPhi_);
  //tree->Branch("phoSigmaIPhiIPhi",        &phoSigmaIPhiIPhi_);
  tree->Branch("phoE1x3",                 &phoE1x3_);
  tree->Branch("phoE1x5",                 &phoE1x5_);
  tree->Branch("phoE2x2",                 &phoE2x2_);
  tree->Branch("phoE2x5Max",              &phoE2x5Max_);
  tree->Branch("phoE5x5",                 &phoE5x5_);
  tree->Branch("phoESEffSigmaRR",         &phoESEffSigmaRR_);
  tree->Branch("phoSigmaIEtaIEtaFull5x5", &phoSigmaIEtaIEtaFull5x5_);
  tree->Branch("phoSigmaIEtaIPhiFull5x5", &phoSigmaIEtaIPhiFull5x5_);
  tree->Branch("phoSigmaIPhiIPhiFull5x5", &phoSigmaIPhiIPhiFull5x5_);
  tree->Branch("phoE1x3Full5x5",          &phoE1x3Full5x5_);
  tree->Branch("phoE1x5Full5x5",          &phoE1x5Full5x5_);
  tree->Branch("phoE2x2Full5x5",          &phoE2x2Full5x5_);
  tree->Branch("phoE2x5MaxFull5x5",       &phoE2x5MaxFull5x5_);
  tree->Branch("phoE5x5Full5x5",          &phoE5x5Full5x5_);
  tree->Branch("phoR9Full5x5",            &phoR9Full5x5_);
  //tree->Branch("phoSeedBCE",              &phoSeedBCE_);
  //tree->Branch("phoSeedBCEta",            &phoSeedBCEta_);
  tree->Branch("phoPFChIso",              &phoPFChIso_);
  tree->Branch("phoPFPhoIso",             &phoPFPhoIso_);
  tree->Branch("phoPFNeuIso",             &phoPFNeuIso_);
  tree->Branch("phoPFChWorstIso",         &phoPFChWorstIso_);
  /*
  tree->Branch("phoPFChIsoFrix1",         &phoPFChIsoFrix1_);
  tree->Branch("phoPFChIsoFrix2",         &phoPFChIsoFrix2_);
  tree->Branch("phoPFChIsoFrix3",         &phoPFChIsoFrix3_);
  tree->Branch("phoPFChIsoFrix4",         &phoPFChIsoFrix4_);
  tree->Branch("phoPFChIsoFrix5",         &phoPFChIsoFrix5_);
  tree->Branch("phoPFChIsoFrix6",         &phoPFChIsoFrix6_);
  tree->Branch("phoPFChIsoFrix7",         &phoPFChIsoFrix7_);
  tree->Branch("phoPFChIsoFrix8",         &phoPFChIsoFrix8_);
  tree->Branch("phoPFPhoIsoFrix1",        &phoPFPhoIsoFrix1_);
  tree->Branch("phoPFPhoIsoFrix2",        &phoPFPhoIsoFrix2_);
  tree->Branch("phoPFPhoIsoFrix3",        &phoPFPhoIsoFrix3_);
  tree->Branch("phoPFPhoIsoFrix4",        &phoPFPhoIsoFrix4_);
  tree->Branch("phoPFPhoIsoFrix5",        &phoPFPhoIsoFrix5_);
  tree->Branch("phoPFPhoIsoFrix6",        &phoPFPhoIsoFrix6_);
  tree->Branch("phoPFPhoIsoFrix7",        &phoPFPhoIsoFrix7_);
  tree->Branch("phoPFPhoIsoFrix8",        &phoPFPhoIsoFrix8_);
  tree->Branch("phoPFNeuIsoFrix1",        &phoPFNeuIsoFrix1_);
  tree->Branch("phoPFNeuIsoFrix2",        &phoPFNeuIsoFrix2_);
  tree->Branch("phoPFNeuIsoFrix3",        &phoPFNeuIsoFrix3_);
  tree->Branch("phoPFNeuIsoFrix4",        &phoPFNeuIsoFrix4_);
  tree->Branch("phoPFNeuIsoFrix5",        &phoPFNeuIsoFrix5_);
  tree->Branch("phoPFNeuIsoFrix6",        &phoPFNeuIsoFrix6_);
  tree->Branch("phoPFNeuIsoFrix7",        &phoPFNeuIsoFrix7_);
  tree->Branch("phoPFNeuIsoFrix8",        &phoPFNeuIsoFrix8_);
  */
  tree->Branch("phoCITKChIso",            &phoCITKChIso_);
  tree->Branch("phoCITKPhoIso",           &phoCITKPhoIso_);
  tree->Branch("phoCITKNeuIso",           &phoCITKNeuIso_);
  //tree->Branch("phoPUPPIChIso",           &phoPUPPIChIso_);
  //tree->Branch("phoPUPPIPhoIso",          &phoPUPPIPhoIso_);
  //tree->Branch("phoPUPPINeuIso",          &phoPUPPINeuIso_);
  //tree->Branch("phoEcalRecHitSumEtConeDR03",      &phoEcalRecHitSumEtConeDR03_);
  //tree->Branch("phohcalDepth1TowerSumEtConeDR03", &phohcalDepth1TowerSumEtConeDR03_);
  //tree->Branch("phohcalDepth2TowerSumEtConeDR03", &phohcalDepth2TowerSumEtConeDR03_);
  //tree->Branch("phohcalTowerSumEtConeDR03",       &phohcalTowerSumEtConeDR03_);
  //tree->Branch("photrkSumPtHollowConeDR03",       &photrkSumPtHollowConeDR03_);
  //tree->Branch("photrkSumPtSolidConeDR03",        &photrkSumPtSolidConeDR03_);
  tree->Branch("phoIDMVA",                        &phoIDMVA_);
  tree->Branch("phoFiredSingleTrgs",              &phoFiredSingleTrgs_);
  tree->Branch("phoFiredDoubleTrgs",              &phoFiredDoubleTrgs_);
  tree->Branch("phoFiredL1Trgs",                  &phoFiredL1Trgs_);
  tree->Branch("phoSeedTime",                     &phoSeedTime_);
  tree->Branch("phoSeedEnergy",                   &phoSeedEnergy_);
  //tree->Branch("phoSeedTimeFull5x5",              &phoSeedTimeFull5x5_);
  //tree->Branch("phoMIPChi2",                      &phoMIPChi2_);
  //tree->Branch("phoMIPTotEnergy",                 &phoMIPTotEnergy_);
  //tree->Branch("phoMIPSlope",                     &phoMIPSlope_);
  //tree->Branch("phoMIPIntercept",                 &phoMIPIntercept_);
  //tree->Branch("phoMIPNhitCone",                  &phoMIPNhitCone_);
  //tree->Branch("phoMIPIsHalo",                    &phoMIPIsHalo_);

  tree->Branch("phoxtalBits", &phoxtalBits_);
  tree->Branch("phoIDbit",    &phoIDbit_);

}

void ggNtuplizer::fillPhotons(const edm::Event& e, const edm::EventSetup& es) {
  
  // cleanup from previous execution
  phoE_                 .clear();
  phoEt_                .clear();
  phoEta_               .clear();
  phoPhi_               .clear();
  phoCalibE_            .clear();
  phoCalibEt_           .clear();
  phoSCE_               .clear();
  phoSCRawE_            .clear();
  phoESEn_              .clear();
  phoESEnP1_            .clear();
  phoESEnP2_            .clear();
  phoSCEta_             .clear();
  phoSCPhi_             .clear();
  phoSCEtaWidth_        .clear();
  phoSCPhiWidth_        .clear();
  phoSCBrem_            .clear();
  phohasPixelSeed_      .clear();
  phoEleVeto_           .clear();
  phoR9_                .clear();
  phoHoverE_            .clear();
  //phoSigmaIEtaIEta_     .clear();
  //phoSigmaIEtaIPhi_     .clear();
  //phoSigmaIPhiIPhi_     .clear();
  phoE1x3_              .clear();
  phoE1x5_              .clear();
  phoE2x2_              .clear();
  phoE2x5Max_           .clear();
  phoE5x5_              .clear();
  phoESEffSigmaRR_      .clear();
  phoSigmaIEtaIEtaFull5x5_.clear();
  phoSigmaIEtaIPhiFull5x5_.clear();
  phoSigmaIPhiIPhiFull5x5_.clear();
  phoE1x3Full5x5_       .clear();
  phoE1x5Full5x5_       .clear();
  phoE2x2Full5x5_       .clear();
  phoE2x5MaxFull5x5_    .clear();
  phoE5x5Full5x5_       .clear();
  phoR9Full5x5_         .clear();
  phoPFChIso_           .clear();
  phoPFPhoIso_          .clear();
  phoPFNeuIso_          .clear();
  phoPFChWorstIso_      .clear();
  phoPFChIsoFrix1_      .clear();
  phoPFChIsoFrix2_      .clear();
  phoPFChIsoFrix3_      .clear();
  phoPFChIsoFrix4_      .clear();
  phoPFChIsoFrix5_      .clear();
  phoPFChIsoFrix6_      .clear();
  phoPFChIsoFrix7_      .clear();
  phoPFChIsoFrix8_      .clear();
  phoPFPhoIsoFrix1_     .clear();
  phoPFPhoIsoFrix2_     .clear();
  phoPFPhoIsoFrix3_     .clear();
  phoPFPhoIsoFrix4_     .clear();
  phoPFPhoIsoFrix5_     .clear();
  phoPFPhoIsoFrix6_     .clear();
  phoPFPhoIsoFrix7_     .clear();
  phoPFPhoIsoFrix8_     .clear();
  phoPFNeuIsoFrix1_     .clear();
  phoPFNeuIsoFrix2_     .clear();
  phoPFNeuIsoFrix3_     .clear();
  phoPFNeuIsoFrix4_     .clear();
  phoPFNeuIsoFrix5_     .clear();
  phoPFNeuIsoFrix6_     .clear();
  phoPFNeuIsoFrix7_     .clear();
  phoPFNeuIsoFrix8_     .clear();
  phoCITKChIso_         .clear();
  phoCITKPhoIso_        .clear();
  phoCITKNeuIso_        .clear();
  //phoPUPPIChIso_        .clear();
  //phoPUPPIPhoIso_       .clear();
  //phoPUPPINeuIso_       .clear();
  //phoSeedBCE_           .clear();
  //phoSeedBCEta_         .clear();
  phoIDMVA_             .clear();
  phoFiredSingleTrgs_   .clear();
  phoFiredDoubleTrgs_   .clear();
  phoFiredL1Trgs_       .clear();
  phoxtalBits_          .clear();
  phoSeedTime_          .clear();
  phoSeedEnergy_        .clear();
  /*
  phoSeedTimeFull5x5_   .clear();
  phoMIPChi2_           .clear();
  phoMIPTotEnergy_      .clear();
  phoMIPSlope_          .clear();
  phoMIPIntercept_      .clear();
  phoMIPNhitCone_       .clear();
  phoMIPIsHalo_         .clear();
  phoEcalRecHitSumEtConeDR03_     .clear();
  phohcalDepth1TowerSumEtConeDR03_.clear();
  phohcalDepth2TowerSumEtConeDR03_.clear();
  phohcalTowerSumEtConeDR03_      .clear();
  photrkSumPtHollowConeDR03_      .clear();
  photrkSumPtSolidConeDR03_       .clear();
  */

  phoIDbit_                       .clear();
  
  nPho_ = 0;

  edm::Handle<edm::View<pat::Photon> > photonHandle;
  e.getByToken(photonCollection_, photonHandle);

  edm::Handle<edm::View<pat::Photon> > calibphotonHandle;
  e.getByToken(calibphotonCollection_, calibphotonHandle);

  if (!photonHandle.isValid()) {
    edm::LogWarning("ggNtuplizer") << "no pat::Photons in event";
    return;
  }

  if (!calibphotonHandle.isValid()) {
    edm::LogWarning("ggNtuplizer") << "no calibrated pat::Photons in event";
    return;
  }

  edm::Handle<reco::PhotonCollection> recoPhotonHandle;
  e.getByToken(recophotonCollection_, recoPhotonHandle);

  ///Photon ID in VID framwork - 11th may, 2015
  // Get the photon ID data from the event stream.
  // Note: this implies that the VID ID modules have been run upstream.
  // If you need more info, check with the EGM group.
  edm::Handle<edm::ValueMap<bool> >  loose_id_decisions;
  edm::Handle<edm::ValueMap<bool> >  medium_id_decisions;
  edm::Handle<edm::ValueMap<bool> >  tight_id_decisions;
  edm::Handle<edm::ValueMap<float> > mvaValues;
  edm::Handle<edm::ValueMap<float> > phoChargedIsolationMap;
  edm::Handle<edm::ValueMap<float> > phoNeutralHadronIsolationMap;
  edm::Handle<edm::ValueMap<float> > phoPhotonIsolationMap;
  edm::Handle<edm::ValueMap<float> > phoWorstChargedIsolationMap;
  edm::Handle<std::vector<pat::PackedCandidate>> pfCndHandle; 
  edm::Handle<edm::View<reco::Candidate> > CndHandle; 
  
  e.getByToken(phoLooseIdMapToken_ ,  loose_id_decisions);
  e.getByToken(phoMediumIdMapToken_,  medium_id_decisions);
  e.getByToken(phoTightIdMapToken_ ,  tight_id_decisions);
  e.getByToken(phoMVAValuesMapToken_, mvaValues);
  
  e.getByToken(phoChargedIsolationToken_,       phoChargedIsolationMap);
  e.getByToken(phoNeutralHadronIsolationToken_, phoNeutralHadronIsolationMap);
  e.getByToken(phoPhotonIsolationToken_,        phoPhotonIsolationMap);
  e.getByToken(phoWorstChargedIsolationToken_,  phoWorstChargedIsolationMap);

  // Get the isolation maps for CITK
  edm::Handle<edm::ValueMap<float> > phoChargedIsolationMap_CITK;
  e.getByToken(phoChargedIsolationToken_CITK_, phoChargedIsolationMap_CITK);
  edm::Handle<edm::ValueMap<float> > phoPhotonIsolationMap_CITK;
  e.getByToken(phoPhotonIsolationToken_CITK_, phoPhotonIsolationMap_CITK);
  edm::Handle<edm::ValueMap<float> > phoNeutralHadronIsolationMap_CITK;
  e.getByToken(phoNeutralHadronIsolationToken_CITK_, phoNeutralHadronIsolationMap_CITK);

  // Get the isolation maps for PUPPI
  //edm::Handle<edm::ValueMap<float> > phoChargedIsolationMap_PUPPI;
  //e.getByToken(phoChargedIsolationToken_PUPPI_, phoChargedIsolationMap_PUPPI);
  //edm::Handle<edm::ValueMap<float> > phoPhotonIsolationMap_PUPPI;
  //e.getByToken(phoPhotonIsolationToken_PUPPI_, phoPhotonIsolationMap_PUPPI);
  //edm::Handle<edm::ValueMap<float> > phoNeutralHadronIsolationMap_PUPPI;
  //e.getByToken(phoNeutralHadronIsolationToken_PUPPI_, phoNeutralHadronIsolationMap_PUPPI);

  EcalClusterLazyTools       lazyTool    (e, es, ebReducedRecHitCollection_, eeReducedRecHitCollection_, esReducedRecHitCollection_);
  noZS::EcalClusterLazyTools lazyToolnoZS(e, es, ebReducedRecHitCollection_, eeReducedRecHitCollection_, esReducedRecHitCollection_);

  GEDPhoIDTools* GEDIdTool = NULL;
  if (isAOD_)
    GEDIdTool = new GEDPhoIDTools(e);

  edm::Handle<reco::VertexCollection> recVtxs;
  e.getByToken(vtxLabel_, recVtxs);

  edm::Handle<reco::VertexCollection> recVtxsBS;
  e.getByToken(vtxBSLabel_, recVtxsBS);

  edm::Handle<reco::TrackCollection> tracksHandle;
  e.getByToken(tracklabel_, tracksHandle);

  edm::Handle<reco::GsfElectronCollection> gsfElectronHandle;
  e.getByToken(gsfElectronlabel_, gsfElectronHandle);

  edm::Handle<double> rhoHandle;
  e.getByToken(rhoLabel_, rhoHandle);
  double rho    = *(rhoHandle.product());

  if (isAOD_) {
    edm::Handle<reco::PFCandidateCollection> pfAllCandidates;
    e.getByToken(pfAllParticles_, pfAllCandidates);

    //cicPhotonId_->configure(recVtxsBS, tracksHandle, gsfElectronHandle, pfAllCandidates, rho);
    cicPhotonId_->configure(recVtxs, tracksHandle, gsfElectronHandle, pfAllCandidates, rho);
  }

  for (edm::View<pat::Photon>::const_iterator iPho = photonHandle->begin(); iPho != photonHandle->end(); ++iPho) {

    Float_t corrPt = -1;
    Float_t corrEn = -1;
    for (edm::View<pat::Photon>::const_iterator iCPho = calibphotonHandle->begin(); iCPho != calibphotonHandle->end(); ++iCPho) {
      if (fabs(iPho->eta() - iCPho->eta()) < 0.001 && fabs(iPho->phi() - iCPho->phi()) < 0.001) {
	corrPt = iCPho->pt();
	corrEn = iCPho->energy();
      }
    }
    phoCalibEt_       .push_back(corrPt);
    phoCalibE_        .push_back(corrEn);

    phoE_             .push_back(iPho->energy());
    phoEt_            .push_back(iPho->et());
    phoEta_           .push_back(iPho->eta());
    phoPhi_           .push_back(iPho->phi());
    phoSCE_           .push_back((*iPho).superCluster()->energy());
    phoSCRawE_        .push_back((*iPho).superCluster()->rawEnergy());
    phoESEn_          .push_back((*iPho).superCluster()->preshowerEnergy());
    phoESEnP1_        .push_back((*iPho).superCluster()->preshowerEnergyPlane1());
    phoESEnP2_        .push_back((*iPho).superCluster()->preshowerEnergyPlane2());
    phoSCEta_         .push_back((*iPho).superCluster()->eta());
    phoSCPhi_         .push_back((*iPho).superCluster()->phi());
    phoSCEtaWidth_    .push_back((*iPho).superCluster()->etaWidth());
    phoSCPhiWidth_    .push_back((*iPho).superCluster()->phiWidth());
    phoSCBrem_        .push_back((*iPho).superCluster()->phiWidth()/(*iPho).superCluster()->etaWidth());
    phohasPixelSeed_  .push_back((Int_t)iPho->hasPixelSeed());
    phoEleVeto_       .push_back((Int_t)iPho->passElectronVeto());
    phoR9_            .push_back(iPho->r9());
    phoHoverE_        .push_back(iPho->hadTowOverEm());
    //phoSigmaIEtaIEta_ .push_back(iPho->see());
    //phoSigmaIEtaIPhi_ .push_back(iPho->sep());
    //phoSigmaIPhiIPhi_ .push_back(iPho->spp());
    phoE1x3_          .push_back(lazyTool.e1x3(*((*iPho).superCluster()->seed())));
    phoE1x5_          .push_back(iPho->e1x5());
    phoE2x2_          .push_back(lazyTool.e2x2(*((*iPho).superCluster()->seed())));
    phoE2x5Max_       .push_back(iPho->e2x5());
    phoE5x5_          .push_back(iPho->e5x5());
    phoESEffSigmaRR_  .push_back(lazyTool.eseffsirir(*((*iPho).superCluster())));
    //phoPFChIso_       .push_back(iPho->chargedHadronIso());
    //phoPFPhoIso_      .push_back(iPho->photonIso());
    //phoPFNeuIso_      .push_back(iPho->neutralHadronIso());


    ///////////////////////////////SATURATED/UNSATURATED ///from ggFlash////
    DetId seed = (iPho->superCluster()->seed()->hitsAndFractions())[0].first;
    bool isBarrel = seed.subdetId() == EcalBarrel;
    const EcalRecHitCollection * rechits = (isBarrel?lazyTool.getEcalEBRecHitCollection():lazyTool.getEcalEERecHitCollection());
            
    EcalRecHitCollection::const_iterator theSeedHit = rechits->find(seed);
    if (theSeedHit != rechits->end()) {
      //std::cout<<"(*theSeedHit).time()"<<(*theSeedHit).time()<<"seed energy: "<<(*theSeedHit).energy()<<std::endl;  

      phoSeedTime_  .push_back((*theSeedHit).time());
      phoSeedEnergy_.push_back((*theSeedHit).energy());
    } else{
      phoSeedTime_  .push_back(-99.);
      phoSeedEnergy_.push_back(-99.);
    }

    /// if( isBarrel ) {
    ///     EBDetId ebId(seed);
    ///     cout << "seed barrel " << ebId.ieta() << " " << ebId.iphi() << endl;
    /// } else {
    ///     EEDetId eeId(seed);
    ///     cout << "seed endpcas " << eeId.ix() << " " << eeId.iy() << endl;
    /// 
    /// }
    unsigned short nSaturated = 0, nLeRecovered = 0, nNeighRecovered = 0, nGain1 = 0, nGain6 = 0, nWeired = 0;

    int isSaturated = 0;
    int isSaturated_gain6 = 0;
    
    UShort_t tmpxtalbit = 0;

    auto matrix5x5 = lazyTool.matrixDetId(seed,-2,+2,-2,+2);
    for (auto & deId : matrix5x5 ) {
      /// cout << "matrix " << deId.rawId() << endl;
      auto rh = rechits->find(deId);
      if( rh != rechits->end() ) {
	nSaturated += rh->checkFlag( EcalRecHit::kSaturated );
	nLeRecovered += rh->checkFlag( EcalRecHit::kLeadingEdgeRecovered );
	nNeighRecovered += rh->checkFlag( EcalRecHit::kNeighboursRecovered );
	nGain1 += rh->checkFlag( EcalRecHit::kHasSwitchToGain1 );
	nGain6 += rh->checkFlag( EcalRecHit::kHasSwitchToGain6 );
	nWeired += rh->checkFlag( EcalRecHit::kWeird ) || rh->checkFlag( EcalRecHit::kDiWeird );
	
	if( rh->checkFlag( EcalRecHit::kHasSwitchToGain1 ) && rh->checkFlag( EcalRecHit::kSaturated ) && !isSaturated){ //this is to fill only once, i.e. only if xtal has this, no need to check for other xtals

	  setbit(tmpxtalbit, 0);
	  isSaturated = 1;
	  //break;
	}
	
	if( rh->checkFlag( EcalRecHit::kHasSwitchToGain6 ) && rh->checkFlag( EcalRecHit::kSaturated ) && !isSaturated_gain6){ //this is to fill only once, i.e. only if xtal has this, no need to check for other xtals

	  setbit(tmpxtalbit, 1);
	  isSaturated_gain6 = 1;
	  //break;
	}
	
      }//if( rh != rechits->end() ) 
       
      if (nWeired>0) setbit(tmpxtalbit,2);      
      if (nGain6>0) setbit(tmpxtalbit,3); 

    }//for(auto & deId : matrix5x5 )
  
    phoxtalBits_.push_back(tmpxtalbit);

    phoFiredSingleTrgs_     .push_back(matchSinglePhotonTriggerFilters(iPho->et(), iPho->eta(), iPho->phi()));
    phoFiredDoubleTrgs_     .push_back(matchDoublePhotonTriggerFilters(iPho->et(), iPho->eta(), iPho->phi()));
    phoFiredL1Trgs_         .push_back(matchL1TriggerFilters(iPho->et(), iPho->eta(), iPho->phi()));

    std::vector<float> vCov = lazyToolnoZS.localCovariances( *((*iPho).superCluster()->seed()) );
    //const float see = (isnan(vCov[0]) ? 0. : sqrt(vCov[0]));
    const float spp = (isnan(vCov[2]) ? 0. : sqrt(vCov[2]));
    const float sep = vCov[1];

    phoSigmaIEtaIEtaFull5x5_ .push_back(iPho->full5x5_sigmaIetaIeta());
    phoSigmaIEtaIPhiFull5x5_ .push_back(sep);
    phoSigmaIPhiIPhiFull5x5_ .push_back(spp);
    phoE1x3Full5x5_          .push_back(lazyToolnoZS.e1x3(*((*iPho).superCluster()->seed())));
    phoE1x5Full5x5_          .push_back(iPho->full5x5_e1x5());
    phoE2x2Full5x5_          .push_back(lazyToolnoZS.e2x2(*((*iPho).superCluster()->seed())));
    phoE2x5MaxFull5x5_       .push_back(iPho->full5x5_e2x5());
    phoE5x5Full5x5_          .push_back(iPho->full5x5_e5x5());
    phoR9Full5x5_            .push_back(iPho->full5x5_r9());

    if (isAOD_) {
      size_t rightRecoPho = -1;
      for (size_t iv = 0; iv < recoPhotonHandle->size(); ++iv) {
        reco::PhotonRef recophoRef2(recoPhotonHandle, iv);
        if (deltaR(iPho->eta(), iPho->phi(), recophoRef2->eta(), recophoRef2->phi()) < 0.01) rightRecoPho = iv;
      }

      reco::PhotonRef recophoRef(recoPhotonHandle, rightRecoPho);
      reco::Vertex pv = recVtxs->at(0);
      GEDIdTool->setPhotonP4(recophoRef, pv);

      //phoPFChIso_       .push_back(GEDIdTool->SolidConeIso(0.3, reco::PFCandidate::h));
      //phoPFPhoIso_      .push_back(GEDIdTool->SolidConeIso(0.3, reco::PFCandidate::gamma));
      //phoPFNeuIso_      .push_back(GEDIdTool->SolidConeIso(0.3, reco::PFCandidate::h0));

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
      //phoPFChWorstIso_  .push_back(*max_element(vtxIsolations03.begin(), vtxIsolations03.end()));
    }
      
    //
    //Filling Frix in case this is not AOD
    /*
    if (!isAOD_) {

      e.getByToken(pckPFCdsLabel_, pfCndHandle); 
      e.getByToken(recoCdsLabel_,  CndHandle); 
      auto_ptr<vector<pat::PackedCandidate> > CandColl( new vector<pat::PackedCandidate> (*pfCndHandle) );

      //register 8 variables for each type of Iso corresponding to each frix. ring 
      double IsoRNeu[8] = {0};
      double IsoRPho[8] = {0};
      double IsoRChg[8] = {0};
      
      //Quantities for charged hadron subtraction
      const float dxyMax     = 0.1; 
      const float dzMax      = 0.2; 
      reco::Vertex pv = recVtxs->at(0);   

      //Loop over PFCandidates in the event
      for(uint ika = 0; ika < CandColl->size() ;ika++){ //PFCand Loop
	pat::PackedCandidate & pfCand = (*CandColl)[ika];
	const auto& iCand = CndHandle->ptrAt(ika);
	
	double DRgamma_cand = deltaR(iPho->eta(),iPho->phi(),iCand->eta(),iCand->phi()); 
	
	if(DRgamma_cand <  0.8){
	  bool isF = isInFootprint(iPho->associatedPackedPFCandidates(),iCand); 
	  if(isF == 0){
	    int ir = DRgamma_cand*10;	    
	    if( pfCand.pdgId() == 22  ) IsoRPho[ir] += pfCand.pt();
	    if( pfCand.pdgId() == 130 ) IsoRNeu[ir] += pfCand.pt();
	    if( pfCand.pdgId() == 211 ){
	      float dz  = fabs(pfCand.pseudoTrack().dz(pv.position()));
	      float dxy = fabs(pfCand.pseudoTrack().dxy(pv.position()));
	      if( dxyMax > dxy ){
		if( dzMax > dz ){
		  IsoRChg[ir] += pfCand.pt();
		}
	      } 
	    }// Charge Hadron Iso
	  } // Is not in Photon Footprint
	} // DR photon cand check
      }//eof loop on pfCands 
      phoPFChIsoFrix1_.push_back(IsoRChg[0]);
      phoPFChIsoFrix2_.push_back(IsoRChg[1]);
      phoPFChIsoFrix3_.push_back(IsoRChg[2]);
      phoPFChIsoFrix4_.push_back(IsoRChg[3]);
      phoPFChIsoFrix5_.push_back(IsoRChg[4]);
      phoPFChIsoFrix6_.push_back(IsoRChg[5]);
      phoPFChIsoFrix7_.push_back(IsoRChg[6]);
      phoPFChIsoFrix8_.push_back(IsoRChg[7]);
   
      phoPFPhoIsoFrix1_.push_back(IsoRPho[0]);
      phoPFPhoIsoFrix2_.push_back(IsoRPho[1]);
      phoPFPhoIsoFrix3_.push_back(IsoRPho[2]);
      phoPFPhoIsoFrix4_.push_back(IsoRPho[3]);
      phoPFPhoIsoFrix5_.push_back(IsoRPho[4]);
      phoPFPhoIsoFrix6_.push_back(IsoRPho[5]);
      phoPFPhoIsoFrix7_.push_back(IsoRPho[6]);
      phoPFPhoIsoFrix8_.push_back(IsoRPho[7]);
   
      phoPFNeuIsoFrix1_.push_back(IsoRNeu[0]);
      phoPFNeuIsoFrix2_.push_back(IsoRNeu[1]);
      phoPFNeuIsoFrix3_.push_back(IsoRNeu[2]);
      phoPFNeuIsoFrix4_.push_back(IsoRNeu[3]);
      phoPFNeuIsoFrix5_.push_back(IsoRNeu[4]);
      phoPFNeuIsoFrix6_.push_back(IsoRNeu[5]);
      phoPFNeuIsoFrix7_.push_back(IsoRNeu[6]);
      phoPFNeuIsoFrix8_.push_back(IsoRNeu[7]);

    } // is not AOD for frix calculations 
    */

    //phoSeedBCE_        .push_back((*iPho).superCluster()->seed()->energy());
    //phoSeedBCEta_      .push_back((*iPho).superCluster()->seed()->eta());
    /*
    phoSeedTimeFull5x5_.push_back(lazyToolnoZS.SuperClusterSeedTime(*((*iPho).superCluster())));
    phoMIPChi2_        .push_back(iPho->mipChi2());
    phoMIPTotEnergy_   .push_back(iPho->mipTotEnergy());
    phoMIPSlope_       .push_back(iPho->mipSlope());
    phoMIPIntercept_   .push_back(iPho->mipIntercept());
    phoMIPNhitCone_    .push_back(iPho->mipNhitCone());
    phoMIPIsHalo_      .push_back(iPho->mipIsHalo());
    */
    ///SJ - isolation variables
    //phoEcalRecHitSumEtConeDR03_                   .push_back(iPho->ecalRecHitSumEtConeDR03());
    //phohcalDepth1TowerSumEtConeDR03_              .push_back(iPho->hcalDepth1TowerSumEtConeDR03());
    //phohcalDepth2TowerSumEtConeDR03_              .push_back(iPho->hcalDepth2TowerSumEtConeDR03());
    //phohcalTowerSumEtConeDR03_                    .push_back(iPho->hcalTowerSumEtConeDR03());
    //photrkSumPtHollowConeDR03_                    .push_back(iPho->trkSumPtHollowConeDR03());
    //photrkSumPtSolidConeDR03_                     .push_back(iPho->trkSumPtSolidConeDR03());

    const auto pho = photonHandle->ptrAt(nPho_);
    
    UShort_t tmpphoIDbit = 0;
    
    phoPFChIso_              .push_back((*phoChargedIsolationMap)[pho]);
    phoPFPhoIso_             .push_back((*phoPhotonIsolationMap)[pho]);
    phoPFNeuIso_             .push_back((*phoNeutralHadronIsolationMap)[pho]);
    phoPFChWorstIso_         .push_back((*phoWorstChargedIsolationMap)[pho]);
    
    phoCITKChIso_            .push_back((*phoChargedIsolationMap_CITK)[pho]);
    phoCITKPhoIso_           .push_back((*phoPhotonIsolationMap_CITK)[pho]);
    phoCITKNeuIso_           .push_back((*phoNeutralHadronIsolationMap_CITK)[pho]);
    //phoPUPPIChIso_           .push_back((*phoChargedIsolationMap_PUPPI)[pho]);
    //phoPUPPIPhoIso_          .push_back((*phoPhotonIsolationMap_PUPPI)[pho]);
    //phoPUPPINeuIso_          .push_back((*phoNeutralHadronIsolationMap_PUPPI)[pho]);
    
    bool isPassLoose  = (*loose_id_decisions)[pho];
    if(isPassLoose) setbit(tmpphoIDbit, 0);
    
    bool isPassMedium = (*medium_id_decisions)[pho];
    if(isPassMedium) setbit(tmpphoIDbit, 1);
    
    bool isPassTight  = (*tight_id_decisions)[pho];
    if(isPassTight) setbit(tmpphoIDbit, 2);
    
    phoIDMVA_.push_back((*mvaValues)[pho]);  
    phoIDbit_.push_back(tmpphoIDbit);      

    nPho_++;
  }

  if (GEDIdTool) delete GEDIdTool;
}

void ggNtuplizer::cleanupPhotons() {

}
