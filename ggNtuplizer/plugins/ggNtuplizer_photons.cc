#include <TString.h>
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "RecoCaloTools/Navigation/interface/CaloRectangle.h"
#include "ggAnalysis/ggNtuplizer/interface/ggNtuplizer.h"

using namespace std;

Int_t          nPho_;
vector<float>  phoE_;
vector<float>  phoSigmaE_;
vector<float>  phoEt_;
vector<float>  phoEta_;
vector<float>  phoPhi_;
vector<float>  phoCalibE_;
vector<float>  phoCalibEt_;
vector<float>  phoSCE_;
vector<float>  phoSCRawE_;
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
vector<float>  phoConeHoverE_;
vector<float>  phoESEffSigmaRR_;
vector<float>  phoSigmaIEtaIEtaFull5x5_;
vector<float>  phoSigmaIEtaIPhiFull5x5_;
vector<float>  phoSigmaIPhiIPhiFull5x5_;
vector<float>  phoE1x3Full5x5_;
vector<float>  phoE2x2Full5x5_;
vector<float>  phoE2x5Full5x5_;
vector<float>  phoE3x3Full5x5_;
vector<float>  phoE5x5Full5x5_;
vector<float>  phoEmax_;
vector<float>  phoE2nd_;
vector<float>  phoEtop_;
vector<float>  phoEleft_;
vector<float>  phoEright_;
vector<float>  phoEbottom_;
vector<float>  phoR9Full5x5_;
vector<float>  phoPFChIso_;
vector<float>  phoPFChPVIso_;
vector<float>  phoPFPhoIso_;
vector<float>  phoPFNeuIso_;
vector<float>  phoPFChWorstIso_;
vector<float>  phoPFChWorstVetoIso_;
vector<float>  phoTrkIsoHollowConeDR03_;
vector<float>  phoEcalPFClusterIso_;
vector<float>  phoHcalPFClusterIso_;
//vector<float>  phoSeedBCE_;
//vector<float>  phoSeedBCEta_;
vector<float>  phoIDMVA_;
vector<ULong64_t> phoFiredSingleTrgs_;
vector<ULong64_t> phoFiredDoubleTrgs_;
vector<ULong64_t> phoFiredTripleTrgs_;
vector<ULong64_t> phoFiredL1Trgs_;
vector<float>  phoSeedTime_;
vector<float>  phoSeedEnergy_;
vector<float>  phoMIPTotEnergy_;
//vector<float>  phoSeedTimeFull5x5_;
//vector<float>  phoMIPChi2_;
//vector<float>  phoMIPSlope_;
//vector<float>  phoMIPIntercept_;
//vector<float>  phoMIPNhitCone_;
//vector<float>  phoMIPIsHalo_;
vector<UShort_t> phoxtalBits_;
vector<UShort_t> phoIDbit_;
vector<float>    phoScale_stat_up_;
vector<float>    phoScale_stat_dn_;
vector<float>    phoScale_syst_up_;
vector<float>    phoScale_syst_dn_;
vector<float>    phoScale_gain_up_;
vector<float>    phoScale_gain_dn_;
vector<float>    phoResol_rho_up_;
vector<float>    phoResol_rho_dn_;
vector<float>    phoResol_phi_up_;
vector<float>    phoResol_phi_dn_;

//Necessary for the Photon Footprint removal
template <class T, class U>
bool isInFootprint(const T& thefootprint, const U& theCandidate) {
  for ( auto itr = thefootprint.begin(); itr != thefootprint.end(); ++itr ) {

    if( itr.key() == theCandidate.key() ) return true;
    
  }
  return false;
}

void ggNtuplizer::branchesPhotons(TTree* tree) {
  
  tree->Branch("nPho",                      &nPho_);
  tree->Branch("phoE",                      &phoE_);
  tree->Branch("phoSigmaE",                 &phoSigmaE_);
  tree->Branch("phoEt",                     &phoEt_);
  tree->Branch("phoEta",                    &phoEta_);
  tree->Branch("phoPhi",                    &phoPhi_);
  tree->Branch("phoCalibE",                 &phoCalibE_);
  tree->Branch("phoCalibEt",                &phoCalibEt_);
  tree->Branch("phoSCE",                    &phoSCE_);
  tree->Branch("phoSCRawE",                 &phoSCRawE_);
  tree->Branch("phoESEnP1",                 &phoESEnP1_);
  tree->Branch("phoESEnP2",                 &phoESEnP2_);
  tree->Branch("phoSCEta",                  &phoSCEta_);
  tree->Branch("phoSCPhi",                  &phoSCPhi_);
  tree->Branch("phoSCEtaWidth",             &phoSCEtaWidth_);
  tree->Branch("phoSCPhiWidth",             &phoSCPhiWidth_);
  tree->Branch("phoSCBrem",                 &phoSCBrem_);
  tree->Branch("phohasPixelSeed",           &phohasPixelSeed_);
  tree->Branch("phoEleVeto",                &phoEleVeto_);
  tree->Branch("phoR9",                     &phoR9_);
  tree->Branch("phoHoverE",                 &phoHoverE_);
  tree->Branch("phoConeHoverE",             &phoConeHoverE_);
  tree->Branch("phoESEffSigmaRR",           &phoESEffSigmaRR_);
  tree->Branch("phoSigmaIEtaIEtaFull5x5",   &phoSigmaIEtaIEtaFull5x5_);
  tree->Branch("phoSigmaIEtaIPhiFull5x5",   &phoSigmaIEtaIPhiFull5x5_);
  tree->Branch("phoSigmaIPhiIPhiFull5x5",   &phoSigmaIPhiIPhiFull5x5_);
  tree->Branch("phoE1x3Full5x5",            &phoE1x3Full5x5_);
  tree->Branch("phoE2x2Full5x5",            &phoE2x2Full5x5_);
  tree->Branch("phoE2x5Full5x5",            &phoE2x5Full5x5_);
  tree->Branch("phoE3x3Full5x5",            &phoE3x3Full5x5_);
  tree->Branch("phoE5x5Full5x5",            &phoE5x5Full5x5_);
  tree->Branch("phoEmax",                   &phoEmax_);
  tree->Branch("phoE2nd",                   &phoE2nd_);
  tree->Branch("phoEtop",                   &phoEtop_);
  tree->Branch("phoEleft",                  &phoEleft_);
  tree->Branch("phoEright",                 &phoEright_);
  tree->Branch("phoEbottom",                &phoEbottom_);
  tree->Branch("phoR9Full5x5",              &phoR9Full5x5_);
  //tree->Branch("phoSeedBCE",              &phoSeedBCE_);
  //tree->Branch("phoSeedBCEta",            &phoSeedBCEta_);
  tree->Branch("phoPFChIso",                &phoPFChIso_);
  tree->Branch("phoPFChPVIso",              &phoPFChPVIso_);
  tree->Branch("phoPFPhoIso",               &phoPFPhoIso_);
  tree->Branch("phoPFNeuIso",               &phoPFNeuIso_);
  tree->Branch("phoPFChWorstIso",           &phoPFChWorstIso_);
  tree->Branch("phoPFChWorstVetoIso",       &phoPFChWorstVetoIso_);
  tree->Branch("phoTrkIsoHollowConeDR03",   &phoTrkIsoHollowConeDR03_);
  tree->Branch("phoEcalPFClusterIso",       &phoEcalPFClusterIso_);
  tree->Branch("phoHcalPFClusterIso",       &phoHcalPFClusterIso_);
  tree->Branch("phoIDMVA",                  &phoIDMVA_);
  tree->Branch("phoFiredSingleTrgs",        &phoFiredSingleTrgs_);
  tree->Branch("phoFiredDoubleTrgs",        &phoFiredDoubleTrgs_);
  tree->Branch("phoFiredTripleTrgs",        &phoFiredTripleTrgs_);
  tree->Branch("phoFiredL1Trgs",            &phoFiredL1Trgs_);
  tree->Branch("phoSeedTime",               &phoSeedTime_);
  tree->Branch("phoSeedEnergy",             &phoSeedEnergy_);
  tree->Branch("phoMIPTotEnergy",           &phoMIPTotEnergy_);
  //tree->Branch("phoSeedTimeFull5x5",              &phoSeedTimeFull5x5_);
  //tree->Branch("phoMIPChi2",                      &phoMIPChi2_);
  //tree->Branch("phoMIPSlope",                     &phoMIPSlope_);
  //tree->Branch("phoMIPIntercept",                 &phoMIPIntercept_);
  //tree->Branch("phoMIPNhitCone",                  &phoMIPNhitCone_);
  //tree->Branch("phoMIPIsHalo",                    &phoMIPIsHalo_);

  tree->Branch("phoxtalBits",      &phoxtalBits_);
  tree->Branch("phoIDbit",         &phoIDbit_);
  tree->Branch("phoScale_stat_up", &phoScale_stat_up_);
  tree->Branch("phoScale_stat_dn", &phoScale_stat_dn_);
  tree->Branch("phoScale_syst_up", &phoScale_syst_up_);
  tree->Branch("phoScale_syst_dn", &phoScale_syst_dn_);
  tree->Branch("phoScale_gain_up", &phoScale_gain_up_);
  tree->Branch("phoScale_gain_dn", &phoScale_gain_dn_);
  tree->Branch("phoResol_rho_up",  &phoResol_rho_up_);
  tree->Branch("phoResol_rho_dn",  &phoResol_rho_dn_);
  tree->Branch("phoResol_phi_up",  &phoResol_phi_up_);
  tree->Branch("phoResol_phi_dn",  &phoResol_phi_dn_);

}

void ggNtuplizer::fillPhotons(const edm::Event& e, const edm::EventSetup& es) {
  
  // cleanup from previous execution
  phoE_                   .clear();
  phoSigmaE_              .clear();
  phoEt_                  .clear();
  phoEta_                 .clear();
  phoPhi_                 .clear();
  phoCalibE_              .clear();
  phoCalibEt_             .clear();
  phoSCE_                 .clear();
  phoSCRawE_              .clear();
  phoESEnP1_              .clear();
  phoESEnP2_              .clear();
  phoSCEta_               .clear();
  phoSCPhi_               .clear();
  phoSCEtaWidth_          .clear();
  phoSCPhiWidth_          .clear();
  phoSCBrem_              .clear();
  phohasPixelSeed_        .clear();
  phoEleVeto_             .clear();
  phoR9_                  .clear();
  phoHoverE_              .clear();
  phoConeHoverE_          .clear();
  phoESEffSigmaRR_        .clear();
  phoSigmaIEtaIEtaFull5x5_.clear();
  phoSigmaIEtaIPhiFull5x5_.clear();
  phoSigmaIPhiIPhiFull5x5_.clear();
  phoE1x3Full5x5_         .clear();
  phoE2x2Full5x5_         .clear();
  phoE2x5Full5x5_         .clear();
  phoE3x3Full5x5_         .clear();
  phoE5x5Full5x5_         .clear();
  phoEmax_                .clear();
  phoE2nd_                .clear();
  phoEtop_                .clear();
  phoEleft_               .clear();
  phoEright_              .clear();
  phoEbottom_             .clear();
  phoR9Full5x5_           .clear();
  phoPFChIso_             .clear();
  phoPFChPVIso_           .clear();
  phoPFPhoIso_            .clear();
  phoPFNeuIso_            .clear();
  phoPFChWorstIso_        .clear();
  phoPFChWorstVetoIso_    .clear();
  phoTrkIsoHollowConeDR03_.clear();
  phoEcalPFClusterIso_    .clear();
  phoHcalPFClusterIso_    .clear();
  //phoSeedBCE_           .clear();
  //phoSeedBCEta_         .clear();
  phoIDMVA_               .clear();
  phoFiredSingleTrgs_     .clear();
  phoFiredDoubleTrgs_     .clear();
  phoFiredTripleTrgs_     .clear();
  phoFiredL1Trgs_         .clear();
  phoxtalBits_            .clear();
  phoSeedTime_            .clear();
  phoSeedEnergy_          .clear();
  phoMIPTotEnergy_        .clear();
  /*
  phoSeedTimeFull5x5_   .clear();
  phoMIPChi2_           .clear();
  phoMIPSlope_          .clear();
  phoMIPIntercept_      .clear();
  phoMIPNhitCone_       .clear();
  phoMIPIsHalo_         .clear();
  */

  phoIDbit_        .clear();
  phoScale_stat_up_.clear();
  phoScale_stat_dn_.clear();
  phoScale_syst_up_.clear();
  phoScale_syst_dn_.clear();
  phoScale_gain_up_.clear();
  phoScale_gain_dn_.clear();
  phoResol_rho_up_ .clear();
  phoResol_rho_dn_ .clear();
  phoResol_phi_up_ .clear();
  phoResol_phi_dn_ .clear();  

  nPho_ = 0;

  edm::Handle<edm::View<pat::Photon> > photonHandle;
  e.getByToken(photonCollection_, photonHandle);

  if (!photonHandle.isValid()) {
    edm::LogWarning("ggNtuplizer") << "no pat::Photons in event";
    return;
  }

  edm::Handle<reco::PhotonCollection> recoPhotonHandle;
  e.getByToken(recophotonCollection_, recoPhotonHandle);

  EcalClusterLazyTools       lazyTool    (e, es, ebReducedRecHitCollection_, eeReducedRecHitCollection_, esReducedRecHitCollection_);
  noZS::EcalClusterLazyTools lazyToolnoZS(e, es, ebReducedRecHitCollection_, eeReducedRecHitCollection_, esReducedRecHitCollection_);

  for (edm::View<pat::Photon>::const_iterator iPho = photonHandle->begin(); iPho != photonHandle->end(); ++iPho) {

    phoE_                     .push_back(iPho->energy());
    phoCalibE_                .push_back(iPho->userFloat("ecalEnergyPostCorr"));
    phoEt_                    .push_back(iPho->et());
    phoCalibEt_               .push_back(iPho->et()*iPho->userFloat("ecalEnergyPostCorr")/iPho->energy());
    phoSigmaE_                .push_back(iPho->userFloat("ecalEnergyErrPostCorr"));
    phoEta_                   .push_back(iPho->eta());
    phoPhi_                   .push_back(iPho->phi());
    phoSCE_                   .push_back((*iPho).superCluster()->energy());
    phoSCRawE_                .push_back((*iPho).superCluster()->rawEnergy());
    phoESEnP1_                .push_back((*iPho).superCluster()->preshowerEnergyPlane1());
    phoESEnP2_                .push_back((*iPho).superCluster()->preshowerEnergyPlane2());
    phoSCEta_                 .push_back((*iPho).superCluster()->eta());
    phoSCPhi_                 .push_back((*iPho).superCluster()->phi());
    phoSCEtaWidth_            .push_back((*iPho).superCluster()->etaWidth());
    phoSCPhiWidth_            .push_back((*iPho).superCluster()->phiWidth());
    phoSCBrem_                .push_back((*iPho).superCluster()->phiWidth()/(*iPho).superCluster()->etaWidth());
    phohasPixelSeed_          .push_back((Int_t)iPho->hasPixelSeed());
    phoEleVeto_               .push_back((Int_t)iPho->passElectronVeto());
    phoR9_                    .push_back(iPho->r9());
    phoHoverE_                .push_back(iPho->hadTowOverEm());
    phoConeHoverE_            .push_back(iPho->hadronicOverEm());
    phoESEffSigmaRR_          .push_back(lazyTool.eseffsirir(*((*iPho).superCluster())));
    phoPFChIso_               .push_back(iPho->chargedHadronIso()); //charged hadron isolation with dxy,dz match to pv
    phoPFChPVIso_             .push_back(iPho->chargedHadronPFPVIso()); //only considers particles assigned to the primary vertex (PV) by particle flow, corresponds to <10_6 chargedHadronIso
    phoPFPhoIso_              .push_back(iPho->photonIso());
    phoPFNeuIso_              .push_back(iPho->neutralHadronIso());
    phoPFChWorstIso_          .push_back(iPho->chargedHadronWorstVtxIso()); //max charged hadron isolation when dxy/dz matching to given vtx
    phoPFChWorstVetoIso_      .push_back(iPho->chargedHadronWorstVtxGeomVetoIso()); //as chargedHadronWorstVtxIso but an additional geometry based veto cone
    phoTrkIsoHollowConeDR03_  .push_back(iPho->trkSumPtHollowConeDR03());
    phoEcalPFClusterIso_      .push_back(iPho->ecalPFClusterIso());
    phoHcalPFClusterIso_      .push_back(iPho->hcalPFClusterIso());
    phoIDMVA_                 .push_back(iPho->userFloat("PhotonMVAEstimatorRunIIFall17v2Values"));  

    // VID decisions     
    UShort_t tmpphoIDbit = 0;        
    bool isPassLoose  = iPho->photonID("cutBasedPhotonID-Fall17-94X-V2-loose");
    if (isPassLoose)  setbit(tmpphoIDbit, 0);   
    bool isPassMedium = iPho->photonID("cutBasedPhotonID-Fall17-94X-V2-medium");
    if (isPassMedium) setbit(tmpphoIDbit, 1);    
    bool isPassTight  = iPho->photonID("cutBasedPhotonID-Fall17-94X-V2-tight");
    if (isPassTight)  setbit(tmpphoIDbit, 2);
    
    phoIDbit_.push_back(tmpphoIDbit);      

    // systematics for energy scale and resolution
    phoScale_stat_up_.push_back(iPho->userFloat("energyScaleStatUp"));
    phoScale_stat_dn_.push_back(iPho->userFloat("energyScaleStatDown"));
    phoScale_syst_up_.push_back(iPho->userFloat("energyScaleSystUp"));
    phoScale_syst_dn_.push_back(iPho->userFloat("energyScaleSystDown"));
    phoScale_gain_up_.push_back(iPho->userFloat("energyScaleGainUp"));
    phoScale_gain_dn_.push_back(iPho->userFloat("energyScaleGainDown"));
    phoResol_rho_up_ .push_back(iPho->userFloat("energySigmaRhoUp"));
    phoResol_rho_dn_ .push_back(iPho->userFloat("energySigmaRhoDown"));
    phoResol_phi_up_ .push_back(iPho->userFloat("energySigmaPhiUp"));
    phoResol_phi_dn_ .push_back(iPho->userFloat("energySigmaPhiDown"));

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
    
    unsigned short nSaturated = 0, nLeRecovered = 0, nNeighRecovered = 0, nGain1 = 0, nGain6 = 0, nWeired = 0;
    int isSaturated       = 0;
    int isSaturated_gain6 = 0;
    
    UShort_t tmpxtalbit = 0;

    auto matrix5x5 = CaloRectangleRange(2, seed, *topology_);
    for (auto const& deId : matrix5x5 ) {

      auto rh = rechits->find(deId);
      if (rh != rechits->end()) {
	nSaturated      += rh->checkFlag( EcalRecHit::kSaturated );
	nLeRecovered    += rh->checkFlag( EcalRecHit::kLeadingEdgeRecovered );
	nNeighRecovered += rh->checkFlag( EcalRecHit::kNeighboursRecovered );
	nGain1          += rh->checkFlag( EcalRecHit::kHasSwitchToGain1 );
	nGain6          += rh->checkFlag( EcalRecHit::kHasSwitchToGain6 );
	nWeired         += rh->checkFlag( EcalRecHit::kWeird ) || rh->checkFlag( EcalRecHit::kDiWeird );
	
	if (rh->checkFlag( EcalRecHit::kHasSwitchToGain1 ) && rh->checkFlag( EcalRecHit::kSaturated ) && !isSaturated) { 
	  setbit(tmpxtalbit, 0);
	  isSaturated = 1;
	  //break;
	}
	
	if (rh->checkFlag( EcalRecHit::kHasSwitchToGain6 ) && rh->checkFlag( EcalRecHit::kSaturated ) && !isSaturated_gain6) {
	  setbit(tmpxtalbit, 1);
	  isSaturated_gain6 = 1;
	  //break;
	}
	
      }//if( rh != rechits->end() ) 
       
      if (nWeired>0) setbit(tmpxtalbit,2);      
      if (nGain6>0)  setbit(tmpxtalbit,3); 
    }
  
    phoxtalBits_.push_back(tmpxtalbit);

    phoFiredSingleTrgs_     .push_back(matchSinglePhotonTriggerFilters(iPho->et(), iPho->eta(), iPho->phi()));
    phoFiredDoubleTrgs_     .push_back(matchDoublePhotonTriggerFilters(iPho->et(), iPho->eta(), iPho->phi()));
    phoFiredTripleTrgs_     .push_back(matchTriplePhotonTriggerFilters(iPho->et(), iPho->eta(), iPho->phi()));
    phoFiredL1Trgs_         .push_back(matchL1TriggerFilters(iPho->et(), iPho->eta(), iPho->phi()));

    phoSigmaIEtaIEtaFull5x5_ .push_back(iPho->full5x5_sigmaIetaIeta());
    phoSigmaIEtaIPhiFull5x5_ .push_back(iPho->full5x5_showerShapeVariables().sigmaIetaIphi);
    phoSigmaIPhiIPhiFull5x5_ .push_back(iPho->full5x5_showerShapeVariables().sigmaIphiIphi);
    phoE1x3Full5x5_          .push_back(iPho->full5x5_showerShapeVariables().e1x3);
    phoE2x2Full5x5_          .push_back(iPho->full5x5_showerShapeVariables().e2x2);
    phoE2x5Full5x5_          .push_back(iPho->full5x5_e2x5());
    phoE3x3Full5x5_          .push_back(iPho->full5x5_e3x3());
    phoE5x5Full5x5_          .push_back(iPho->full5x5_e5x5());
    phoEmax_                 .push_back(iPho->full5x5_maxEnergyXtal());
    phoE2nd_                 .push_back(iPho->full5x5_showerShapeVariables().e2nd);
    phoEtop_                 .push_back(iPho->full5x5_showerShapeVariables().eTop); 
    phoEleft_                .push_back(iPho->full5x5_showerShapeVariables().eLeft); 
    phoEright_               .push_back(iPho->full5x5_showerShapeVariables().eRight); 
    phoEbottom_              .push_back(iPho->full5x5_showerShapeVariables().eBottom); 
    phoR9Full5x5_            .push_back(iPho->full5x5_r9());
    phoMIPTotEnergy_         .push_back(iPho->mipTotEnergy());

    //phoSeedBCE_        .push_back((*iPho).superCluster()->seed()->energy());
    //phoSeedBCEta_      .push_back((*iPho).superCluster()->seed()->eta());
    /*
    phoSeedTimeFull5x5_.push_back(lazyToolnoZS.SuperClusterSeedTime(*((*iPho).superCluster())));
    phoMIPChi2_        .push_back(iPho->mipChi2());
    phoMIPSlope_       .push_back(iPho->mipSlope());
    phoMIPIntercept_   .push_back(iPho->mipIntercept());
    phoMIPNhitCone_    .push_back(iPho->mipNhitCone());
    phoMIPIsHalo_      .push_back(iPho->mipIsHalo());
    */
    
    nPho_++;

    tester|=matchSinglePhotonTriggerFilters(iPho->et(), iPho->eta(), iPho->phi());
  }

}

void ggNtuplizer::cleanupPhotons() {

}
