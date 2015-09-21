#include <TString.h>
#include <TMVA/Reader.h>
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
vector<float>  phoSigmaIEtaIEta_;
vector<float>  phoSigmaIEtaIPhi_;
vector<float>  phoSigmaIPhiIPhi_;
vector<float>  phoE1x3_;
vector<float>  phoE2x2_;
vector<float>  phoE2x5Max_;
vector<float>  phoE5x5_;
vector<float>  phoESEffSigmaRR_;
vector<float>  phoSigmaIEtaIEtaFull5x5_;
vector<float>  phoSigmaIEtaIPhiFull5x5_;
vector<float>  phoSigmaIPhiIPhiFull5x5_;
vector<float>  phoE1x3Full5x5_;
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
vector<float>  phoSeedBCE_;
vector<float>  phoSeedBCEta_;
vector<float>  phoIDMVA_;
vector<Int_t>  phoFiredSingleTrgs_;
vector<Int_t>  phoFiredDoubleTrgs_;
vector<float>  phoEcalRecHitSumEtConeDR03_;
vector<float>  phohcalDepth1TowerSumEtConeDR03_;
vector<float>  phohcalDepth2TowerSumEtConeDR03_;
vector<float>  phohcalTowerSumEtConeDR03_;
vector<float>  photrkSumPtHollowConeDR03_;

vector<UShort_t> phoIDbit_;

// variables that will be containers on which TMVA Reader works
float varPhi_;
float varR9_;
float varSieie_;
float varSieip_;
float varE1x3overE5x5_;
float varE2x2overE5x5_;
float varE2x5overE5x5_;
float varSCEta_;
float varRawE_;
float varSCEtaWidth_;
float varSCPhiWidth_;
float varRho_;
float varPhoIsoRaw_;
float varChIsoRaw_;
float varWorstChRaw_;
float varESEnOverRawE_; // for endcap MVA only
float varESEffSigmaRR_; // for endcap MVA only
// The spectators
float varPt_;
float varEta_;

//Necessary for the Photon Footprint removal
template <class T, class U>
bool isInFootprint(const T& thefootprint, const U& theCandidate) {
  for ( auto itr = thefootprint.begin(); itr != thefootprint.end(); ++itr ) {

    if( itr.key() == theCandidate.key() ) return true;
    
  }
  return false;
}



// TMVA Reader for applying MVA
TMVA::Reader *tmvaReader_[2];
TString methodName_[2];

void ggNtuplizer::branchesPhotons(TTree* tree) {
  
  tree->Branch("nPho",                    &nPho_);
  tree->Branch("phoE",                    &phoE_);
  tree->Branch("phoEt",                   &phoEt_);
  tree->Branch("phoEta",                  &phoEta_);
  tree->Branch("phoPhi",                  &phoPhi_);
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
  tree->Branch("phoSigmaIEtaIEta",        &phoSigmaIEtaIEta_);
  tree->Branch("phoSigmaIEtaIPhi",        &phoSigmaIEtaIPhi_);
  tree->Branch("phoSigmaIPhiIPhi",        &phoSigmaIPhiIPhi_);
  tree->Branch("phoE1x3",                 &phoE1x3_);
  tree->Branch("phoE2x2",                 &phoE2x2_);
  tree->Branch("phoE2x5Max",              &phoE2x5Max_);
  tree->Branch("phoE5x5",                 &phoE5x5_);
  tree->Branch("phoESEffSigmaRR",         &phoESEffSigmaRR_);
  tree->Branch("phoSigmaIEtaIEtaFull5x5", &phoSigmaIEtaIEtaFull5x5_);
  tree->Branch("phoSigmaIEtaIPhiFull5x5", &phoSigmaIEtaIPhiFull5x5_);
  tree->Branch("phoSigmaIPhiIPhiFull5x5", &phoSigmaIPhiIPhiFull5x5_);
  tree->Branch("phoE1x3Full5x5",          &phoE1x3Full5x5_);
  tree->Branch("phoE2x2Full5x5",          &phoE2x2Full5x5_);
  tree->Branch("phoE2x5MaxFull5x5",       &phoE2x5MaxFull5x5_);
  tree->Branch("phoE5x5Full5x5",          &phoE5x5Full5x5_);
  tree->Branch("phoR9Full5x5",            &phoR9Full5x5_);
  tree->Branch("phoSeedBCE",              &phoSeedBCE_);
  tree->Branch("phoSeedBCEta",            &phoSeedBCEta_);
  tree->Branch("phoPFChIso",              &phoPFChIso_);
  tree->Branch("phoPFPhoIso",             &phoPFPhoIso_);
  tree->Branch("phoPFNeuIso",             &phoPFNeuIso_);
  tree->Branch("phoPFChWorstIso",         &phoPFChWorstIso_);
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
  tree->Branch("phoEcalRecHitSumEtConeDR03",      &phoEcalRecHitSumEtConeDR03_);
  tree->Branch("phohcalDepth1TowerSumEtConeDR03", &phohcalDepth1TowerSumEtConeDR03_);
  tree->Branch("phohcalDepth2TowerSumEtConeDR03", &phohcalDepth2TowerSumEtConeDR03_);
  tree->Branch("phohcalTowerSumEtConeDR03",       &phohcalTowerSumEtConeDR03_);
  tree->Branch("photrkSumPtHollowConeDR03",       &photrkSumPtHollowConeDR03_);
  tree->Branch("phoIDMVA",                        &phoIDMVA_);
  tree->Branch("phoFiredSingleTrgs",              &phoFiredSingleTrgs_);
  tree->Branch("phoFiredDoubleTrgs",              &phoFiredDoubleTrgs_);
  if (runphoIDVID_) tree->Branch("phoIDbit",      &phoIDbit_);

  if (isAOD_ && runphoMVAID_) {
    ////////////Prepare for photon ID MVA for Run II///////////
    //
    // Create and configure barrel MVA
    //
    tmvaReader_[0] = new TMVA::Reader( "!Color:!Silent:Error" );
    tmvaReader_[0]->SetVerbose(kFALSE);
    // Add all the vars, we take the string with variable name from the weights file (the Expression field)
    tmvaReader_[0]->AddVariable("recoPhi"   , &varPhi_);
    tmvaReader_[0]->AddVariable("r9"        , &varR9_);
    tmvaReader_[0]->AddVariable("sieieFull5x5", &varSieie_);
    tmvaReader_[0]->AddVariable("sieipFull5x5", &varSieip_);
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
    
  }//  if (isAOD_ && runphoMVAID_) 

}

void ggNtuplizer::fillPhotons(const edm::Event& e, const edm::EventSetup& es) {
  
  // cleanup from previous execution
  phoE_                 .clear();
  phoEt_                .clear();
  phoEta_               .clear();
  phoPhi_               .clear();
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
  phoSigmaIEtaIEta_     .clear();
  phoSigmaIEtaIPhi_     .clear();
  phoSigmaIPhiIPhi_     .clear();
  phoE1x3_              .clear();
  phoE2x2_              .clear();
  phoE2x5Max_           .clear();
  phoE5x5_              .clear();
  phoESEffSigmaRR_      .clear();
  phoSigmaIEtaIEtaFull5x5_.clear();
  phoSigmaIEtaIPhiFull5x5_.clear();
  phoSigmaIPhiIPhiFull5x5_.clear();
  phoE1x3Full5x5_       .clear();
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
  phoSeedBCE_           .clear();
  phoSeedBCEta_         .clear();
  phoIDMVA_             .clear();
  phoFiredSingleTrgs_   .clear();
  phoFiredDoubleTrgs_   .clear();

  phoEcalRecHitSumEtConeDR03_     .clear();
  phohcalDepth1TowerSumEtConeDR03_.clear();
  phohcalDepth2TowerSumEtConeDR03_.clear();
  phohcalTowerSumEtConeDR03_      .clear();
  photrkSumPtHollowConeDR03_      .clear();

  phoIDbit_                       .clear();
  
  nPho_ = 0;

  edm::Handle<edm::View<pat::Photon> > photonHandle;
  e.getByToken(photonCollection_, photonHandle);

  if (!photonHandle.isValid()) {
    edm::LogWarning("ggNtuplizer") << "no pat::Photons in event";
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
  
  if (runphoIDVID_) {
    e.getByToken(phoLooseIdMapToken_ ,  loose_id_decisions);
    e.getByToken(phoMediumIdMapToken_,  medium_id_decisions);
    e.getByToken(phoTightIdMapToken_ ,  tight_id_decisions);
    e.getByToken(phoMVAValuesMapToken_, mvaValues);

    e.getByToken(phoChargedIsolationToken_,       phoChargedIsolationMap);
    e.getByToken(phoNeutralHadronIsolationToken_, phoNeutralHadronIsolationMap);
    e.getByToken(phoPhotonIsolationToken_,        phoPhotonIsolationMap);
    e.getByToken(phoWorstChargedIsolationToken_,  phoWorstChargedIsolationMap);
  }

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
    phoSigmaIEtaIEta_ .push_back(iPho->see());
    phoSigmaIEtaIPhi_ .push_back(iPho->sep());
    phoSigmaIPhiIPhi_ .push_back(iPho->spp());
    phoE1x3_          .push_back(lazyTool.e1x3(*((*iPho).superCluster()->seed())));
    phoE2x2_          .push_back(lazyTool.e2x2(*((*iPho).superCluster()->seed())));
    phoE2x5Max_       .push_back(lazyTool.e2x5Max(*((*iPho).superCluster()->seed())));
    phoE5x5_          .push_back(lazyTool.e5x5(*((*iPho).superCluster()->seed())));
    phoESEffSigmaRR_  .push_back(lazyTool.eseffsirir(*((*iPho).superCluster())));
    //phoPFChIso_       .push_back(iPho->chargedHadronIso());
    //phoPFPhoIso_      .push_back(iPho->photonIso());
    //phoPFNeuIso_      .push_back(iPho->neutralHadronIso());

    phoFiredSingleTrgs_     .push_back(matchSinglePhotonTriggerFilters(iPho->et(), iPho->eta(), iPho->phi()));
    phoFiredDoubleTrgs_     .push_back(matchDoublePhotonTriggerFilters(iPho->et(), iPho->eta(), iPho->phi()));

    std::vector<float> vCov = lazyToolnoZS.localCovariances( *((*iPho).superCluster()->seed()) );
    //const float see = (isnan(vCov[0]) ? 0. : sqrt(vCov[0]));
    const float spp = (isnan(vCov[2]) ? 0. : sqrt(vCov[2]));
    const float sep = vCov[1];

    phoSigmaIEtaIEtaFull5x5_ .push_back(iPho->full5x5_sigmaIetaIeta());
    phoSigmaIEtaIPhiFull5x5_ .push_back(sep);
    phoSigmaIPhiIPhiFull5x5_ .push_back(spp);
    phoE1x3Full5x5_          .push_back(lazyToolnoZS.e1x3(*((*iPho).superCluster()->seed())));
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
    if (!isAOD_) {

      e.getByLabel(pckPFCdsLabel_,pfCndHandle); 
      e.getByLabel(pckPFCdsLabel_,CndHandle); 
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

    }// is not AOD for frix calculations 
     //

    phoSeedBCE_          .push_back((*iPho).superCluster()->seed()->energy());
    phoSeedBCEta_        .push_back((*iPho).superCluster()->seed()->eta());

    ///SJ - isolation variables
    phoEcalRecHitSumEtConeDR03_                   .push_back(iPho->ecalRecHitSumEtConeDR03());
    phohcalDepth1TowerSumEtConeDR03_              .push_back(iPho->hcalDepth1TowerSumEtConeDR03());
    phohcalDepth2TowerSumEtConeDR03_              .push_back(iPho->hcalDepth2TowerSumEtConeDR03());
    phohcalTowerSumEtConeDR03_                    .push_back(iPho->hcalTowerSumEtConeDR03());
    photrkSumPtHollowConeDR03_                    .push_back(iPho->trkSumPtHollowConeDR03());

    ////////9th April, 2015 - SJ
    ///MVA for photons in PHYS14
    // set MVA variables
  
    if (isAOD_ && runphoMVAID_ ) {
      varPhi_          = phoPhi_[nPho_];
      varR9_           = phoR9_[nPho_];
      varSieie_        = phoSigmaIEtaIEtaFull5x5_[nPho_];
      varSieip_        = phoSigmaIEtaIPhiFull5x5_[nPho_];
      varE1x3overE5x5_ = phoE1x3Full5x5_[nPho_]/phoE5x5Full5x5_[nPho_];
      varE2x2overE5x5_ = phoE2x2Full5x5_[nPho_]/phoE5x5Full5x5_[nPho_];
      varSCEta_        = phoSCEta_[nPho_];
      varRawE_         = phoSCRawE_[nPho_];
      varSCEtaWidth_   = phoSCEtaWidth_[nPho_];
      varSCPhiWidth_   = phoSCPhiWidth_[nPho_];
      varESEnOverRawE_ = phoESEn_[nPho_]/phoSCRawE_[nPho_];
      varESEffSigmaRR_ = phoESEffSigmaRR_[nPho_];
      varRho_          = rho;
      varPhoIsoRaw_    = phoPFPhoIso_[nPho_];
      varChIsoRaw_     = phoPFChIso_[nPho_];
      varWorstChRaw_   = phoPFNeuIso_[nPho_];

      // spectator vars
      varPt_           = phoEt_[nPho_];
      varEta_          = phoEta_[nPho_];
      
      // 0=ECAL barrel or 1=ECAL endcaps
      //int iBE = (fabs(phoSCEta_[nPho_]) < 1.479) ? 0 : 1;

      //phoIDMVA_.push_back(tmvaReader_[iBE]->EvaluateMVA("BDT"));
    }//if (isAOD_ && runphoMVAID_ )

    if (runphoIDVID_) {
      //Photon ID in VID framwork - 11th may, 2015
      // Look up and save the ID decisions
      // 
      const auto pho = photonHandle->ptrAt(nPho_);
      
      UShort_t tmpphoIDbit = 0;

      if (isAOD_) {
        phoPFChIso_              .push_back((*phoChargedIsolationMap)[pho->originalObjectRef()]);
        phoPFPhoIso_             .push_back((*phoPhotonIsolationMap)[pho->originalObjectRef()]);
        phoPFNeuIso_             .push_back((*phoNeutralHadronIsolationMap)[pho->originalObjectRef()]);
	phoPFChWorstIso_         .push_back((*phoWorstChargedIsolationMap)[pho->originalObjectRef()]);

        //cout<<"Photons "<<endl;
        bool isPassLoose  = (*loose_id_decisions)[pho->originalObjectRef()];
        if(isPassLoose) setbit(tmpphoIDbit, 0);
        //cout<<"isPassLoose "<<isPassLoose<<endl;

        bool isPassMedium = (*medium_id_decisions)[pho->originalObjectRef()];
        if(isPassMedium) setbit(tmpphoIDbit, 1);
        //cout<<"isPassMedium "<<isPassMedium<<endl;

        bool isPassTight  = (*tight_id_decisions)[pho->originalObjectRef()];
        if(isPassTight) setbit(tmpphoIDbit, 2);
        //cout<<"isPassTight "<<isPassTight<<endl;

        phoIDMVA_.push_back((*mvaValues)[pho->originalObjectRef()]);
      }

      if (!isAOD_) {
        phoPFChIso_              .push_back((*phoChargedIsolationMap)[pho]);
        phoPFPhoIso_             .push_back((*phoPhotonIsolationMap)[pho]);
        phoPFNeuIso_             .push_back((*phoNeutralHadronIsolationMap)[pho]);
	phoPFChWorstIso_         .push_back((*phoWorstChargedIsolationMap)[pho]);

        //cout<<"Photons "<<endl;
        bool isPassLoose  = (*loose_id_decisions)[pho];
        if(isPassLoose) setbit(tmpphoIDbit, 0);
        //cout<<"isPassLoose "<<isPassLoose<<endl;

        bool isPassMedium = (*medium_id_decisions)[pho];
        if(isPassMedium) setbit(tmpphoIDbit, 1);
        //cout<<"isPassMedium "<<isPassMedium<<endl;

        bool isPassTight  = (*tight_id_decisions)[pho];
        if(isPassTight) setbit(tmpphoIDbit, 2);
        //cout<<"isPassTight "<<isPassTight<<endl;

        phoIDMVA_.push_back((*mvaValues)[pho]);
      }

      phoIDbit_.push_back(tmpphoIDbit);      
      //cout<<"tmppho : phoIDbit: "<<tmpphoIDbit<<":"<<phoIDbit_[nPho_]<<endl;
    }

    nPho_++;
  }

  if (GEDIdTool) delete GEDIdTool;
}

void ggNtuplizer::cleanupPhotons()
{
  if (isAOD_ && runphoMVAID_) {
    delete tmvaReader_[0];
    delete tmvaReader_[1];
  }
}
