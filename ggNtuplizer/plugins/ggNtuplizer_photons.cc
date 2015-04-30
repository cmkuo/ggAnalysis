#include <TString.h>
#include <TMVA/Reader.h>
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "ggAnalysis/ggNtuplizer/interface/GEDPhoIDTools.h"
#include "ggAnalysis/ggNtuplizer/interface/ggNtuplizer.h"

using namespace std;

// (local) variables associated with tree branches
Int_t          nPho_;
vector<float>  phoE_;
vector<float>  phoEt_;
vector<float>  phoEta_;
vector<float>  phoPhi_;
vector<float>  phoSCE_;
vector<float>  phoSCRawE_;
vector<float>  phoESEn_;
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
vector<float>  phoSigmaIEtaIEta_2012_;
vector<float>  phoSigmaIEtaIPhi_2012_;
vector<float>  phoSigmaIPhiIPhi_2012_;
vector<float>  phoE1x3_2012_;
vector<float>  phoE2x2_2012_;
vector<float>  phoE2x5Max_2012_;
vector<float>  phoE5x5_2012_;
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
vector<float>  phoBC1E_;
vector<float>  phoBC1Eta_;
vector<float>  phoBC2E_;
vector<float>  phoBC2Eta_;

vector<float>  phoIDMVA_;

vector<float>  phoEcalRecHitSumEtConeDR03_;
vector<float>  phohcalDepth1TowerSumEtConeDR03_;
vector<float>  phohcalDepth2TowerSumEtConeDR03_;
vector<float>  phohcalTowerSumEtConeDR03_;
vector<float>  photrkSumPtHollowConeDR03_;

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

// TMVA Reader for applying MVA
TMVA::Reader *tmvaReader_[2];
TString methodName_[2];

void ggNtuplizer::branchesPhotons(TTree* tree)
{
  tree->Branch("nPho",                  &nPho_);
  tree->Branch("phoE",                  &phoE_);
  tree->Branch("phoEt",                 &phoEt_);
  tree->Branch("phoEta",                &phoEta_);
  tree->Branch("phoPhi",                &phoPhi_);
  tree->Branch("phoSCE",                &phoSCE_);
  tree->Branch("phoSCRawE",             &phoSCRawE_);
  tree->Branch("phoESEn",               &phoESEn_);
  tree->Branch("phoSCEta",              &phoSCEta_);
  tree->Branch("phoSCPhi",              &phoSCPhi_);
  tree->Branch("phoSCEtaWidth",         &phoSCEtaWidth_);
  tree->Branch("phoSCPhiWidth",         &phoSCPhiWidth_);
  tree->Branch("phoSCBrem",             &phoSCBrem_);
  tree->Branch("phohasPixelSeed",       &phohasPixelSeed_);
  tree->Branch("phoEleVeto",            &phoEleVeto_);
  tree->Branch("phoR9",                 &phoR9_);
  tree->Branch("phoHoverE",             &phoHoverE_);
  tree->Branch("phoSigmaIEtaIEta",      &phoSigmaIEtaIEta_);
  tree->Branch("phoSigmaIEtaIPhi",      &phoSigmaIEtaIPhi_);
  tree->Branch("phoSigmaIPhiIPhi",      &phoSigmaIPhiIPhi_);
  tree->Branch("phoE1x3",               &phoE1x3_);
  tree->Branch("phoE2x2",               &phoE2x2_);
  tree->Branch("phoE2x5Max",            &phoE2x5Max_);
  tree->Branch("phoE5x5",               &phoE5x5_);
  tree->Branch("phoESEffSigmaRR",       &phoESEffSigmaRR_);
  tree->Branch("phoSigmaIEtaIEta_2012", &phoSigmaIEtaIEta_2012_);
  tree->Branch("phoSigmaIEtaIPhi_2012", &phoSigmaIEtaIPhi_2012_);
  tree->Branch("phoSigmaIPhiIPhi_2012", &phoSigmaIPhiIPhi_2012_);
  tree->Branch("phoE1x3_2012",          &phoE1x3_2012_);
  tree->Branch("phoE2x2_2012",          &phoE2x2_2012_);
  tree->Branch("phoE2x5Max_2012",       &phoE2x5Max_2012_);
  tree->Branch("phoE5x5_2012",          &phoE5x5_2012_);
  tree->Branch("phoPFChIso",            &phoPFChIso_);
  tree->Branch("phoPFPhoIso",           &phoPFPhoIso_);
  tree->Branch("phoPFNeuIso",           &phoPFNeuIso_);
  tree->Branch("phoPFChWorstIso",       &phoPFChWorstIso_);
  tree->Branch("phoPFChIsoFrix1",       &phoPFChIsoFrix1_);
  tree->Branch("phoPFChIsoFrix2",       &phoPFChIsoFrix2_);
  tree->Branch("phoPFChIsoFrix3",       &phoPFChIsoFrix3_);
  tree->Branch("phoPFChIsoFrix4",       &phoPFChIsoFrix4_);
  tree->Branch("phoPFChIsoFrix5",       &phoPFChIsoFrix5_);
  tree->Branch("phoPFChIsoFrix6",       &phoPFChIsoFrix6_);
  tree->Branch("phoPFChIsoFrix7",       &phoPFChIsoFrix7_);
  tree->Branch("phoPFChIsoFrix8",       &phoPFChIsoFrix8_);
  tree->Branch("phoPFPhoIsoFrix1",      &phoPFPhoIsoFrix1_);
  tree->Branch("phoPFPhoIsoFrix2",      &phoPFPhoIsoFrix2_);
  tree->Branch("phoPFPhoIsoFrix3",      &phoPFPhoIsoFrix3_);
  tree->Branch("phoPFPhoIsoFrix4",      &phoPFPhoIsoFrix4_);
  tree->Branch("phoPFPhoIsoFrix5",      &phoPFPhoIsoFrix5_);
  tree->Branch("phoPFPhoIsoFrix6",      &phoPFPhoIsoFrix6_);
  tree->Branch("phoPFPhoIsoFrix7",      &phoPFPhoIsoFrix7_);
  tree->Branch("phoPFPhoIsoFrix8",      &phoPFPhoIsoFrix8_);
  tree->Branch("phoPFNeuIsoFrix1",      &phoPFNeuIsoFrix1_);
  tree->Branch("phoPFNeuIsoFrix2",      &phoPFNeuIsoFrix2_);
  tree->Branch("phoPFNeuIsoFrix3",      &phoPFNeuIsoFrix3_);
  tree->Branch("phoPFNeuIsoFrix4",      &phoPFNeuIsoFrix4_);
  tree->Branch("phoPFNeuIsoFrix5",      &phoPFNeuIsoFrix5_);
  tree->Branch("phoPFNeuIsoFrix6",      &phoPFNeuIsoFrix6_);
  tree->Branch("phoPFNeuIsoFrix7",      &phoPFNeuIsoFrix7_);
  tree->Branch("phoPFNeuIsoFrix8",      &phoPFNeuIsoFrix8_);
  tree->Branch("phoBC1E",               &phoBC1E_);
  tree->Branch("phoBC1Eta",             &phoBC1Eta_);
  tree->Branch("phoBC2E",               &phoBC2E_);
  tree->Branch("phoBC2Eta",             &phoBC2Eta_);

  tree->Branch("phoIDMVA",                        &phoIDMVA_);
  tree->Branch("phoEcalRecHitSumEtConeDR03",      &phoEcalRecHitSumEtConeDR03_);
  tree->Branch("phohcalDepth1TowerSumEtConeDR03", &phohcalDepth1TowerSumEtConeDR03_);
  tree->Branch("phohcalDepth2TowerSumEtConeDR03", &phohcalDepth2TowerSumEtConeDR03_);
  tree->Branch("phohcalTowerSumEtConeDR03",       &phohcalTowerSumEtConeDR03_);

  tree->Branch("photrkSumPtHollowConeDR03",       &photrkSumPtHollowConeDR03_);

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

void ggNtuplizer::fillPhotons(const edm::Event& e, const edm::EventSetup& es)
{

  // cleanup from previous execution
  phoE_                 .clear();
  phoEt_                .clear();
  phoEta_               .clear();
  phoPhi_               .clear();
  phoSCE_               .clear();
  phoSCRawE_            .clear();
  phoESEn_              .clear();
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
  phoSigmaIEtaIEta_2012_.clear();
  phoSigmaIEtaIPhi_2012_.clear();
  phoSigmaIPhiIPhi_2012_.clear();
  phoE1x3_2012_         .clear();
  phoE2x2_2012_         .clear();
  phoE2x5Max_2012_      .clear();
  phoE5x5_2012_         .clear();
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
  phoBC1E_              .clear();
  phoBC1Eta_            .clear();
  phoBC2E_              .clear();
  phoBC2Eta_            .clear();

  phoIDMVA_             .clear();

  phoEcalRecHitSumEtConeDR03_     .clear();
  phohcalDepth1TowerSumEtConeDR03_.clear();
  phohcalDepth2TowerSumEtConeDR03_.clear();
  phohcalTowerSumEtConeDR03_      .clear();
  photrkSumPtHollowConeDR03_      .clear();

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

  GEDPhoIDTools GEDIdTool(e);

  edm::Handle<reco::VertexCollection> recVtxs;
  e.getByToken(vtxLabel_, recVtxs);

  edm::Handle<reco::VertexCollection> recVtxsBS;
  e.getByToken(vtxBSLabel_, recVtxsBS);

  edm::Handle<reco::TrackCollection> tracksHandle;
  e.getByToken(tracklabel_, tracksHandle);

  edm::Handle<reco::GsfElectronCollection> gsfElectronHandle;
  e.getByToken(gsfElectronlabel_, gsfElectronHandle);

  edm::Handle<reco::PFCandidateCollection> pfAllCandidates;
  e.getByToken(pfAllParticles_, pfAllCandidates);

  edm::Handle<double> rhoHandle;
  e.getByToken(rhoLabel_, rhoHandle);
  double rho    = *(rhoHandle.product());

  cicPhotonId_->configure(recVtxsBS, tracksHandle, gsfElectronHandle, pfAllCandidates, rho);

  for (edm::View<pat::Photon>::const_iterator iPho = photonHandle->begin(); iPho != photonHandle->end(); ++iPho) {

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
    phoE1x3_          .push_back(lazyTool.e1x3(*((*iPho).superCluster()->seed())));
    phoE2x2_          .push_back(lazyTool.e2x2(*((*iPho).superCluster()->seed())));
    phoE2x5Max_       .push_back(lazyTool.e2x5Max(*((*iPho).superCluster()->seed())));
    phoE5x5_          .push_back(lazyTool.e5x5(*((*iPho).superCluster()->seed())));
    phoESEffSigmaRR_  .push_back(lazyTool.eseffsirir(*((*iPho).superCluster())));
    //phoPFChIso_       .push_back(iPho->chargedHadronIso());
    //phoPFPhoIso_      .push_back(iPho->photonIso());
    //phoPFNeuIso_      .push_back(iPho->neutralHadronIso());

    std::vector<float> vCov = lazyToolnoZS.localCovariances( *((*iPho).superCluster()->seed()) );
    const float see = (isnan(vCov[0]) ? 0. : sqrt(vCov[0]));
    const float spp = (isnan(vCov[2]) ? 0. : sqrt(vCov[2]));
    const float sep = vCov[1];

    phoSigmaIEtaIEta_2012_ .push_back(see);
    phoSigmaIEtaIPhi_2012_ .push_back(sep);
    phoSigmaIPhiIPhi_2012_ .push_back(spp);
    phoE1x3_2012_          .push_back(lazyToolnoZS.e1x3(*((*iPho).superCluster()->seed())));
    phoE2x2_2012_          .push_back(lazyToolnoZS.e2x2(*((*iPho).superCluster()->seed())));
    phoE2x5Max_2012_       .push_back(lazyToolnoZS.e2x5Max(*((*iPho).superCluster()->seed())));
    phoE5x5_2012_          .push_back(lazyToolnoZS.e5x5(*((*iPho).superCluster()->seed())));

    size_t rightRecoPho = -1;
    for (size_t iv = 0; iv < recoPhotonHandle->size(); ++iv) {
      reco::PhotonRef recophoRef2(recoPhotonHandle, iv);
      if (deltaR(iPho->eta(), iPho->phi(), recophoRef2->eta(), recophoRef2->phi()) < 0.01) rightRecoPho = iv;
    }
    reco::PhotonRef recophoRef(recoPhotonHandle, rightRecoPho);
    reco::Vertex pv = recVtxs->at(0);
    GEDIdTool.setPhotonP4(recophoRef, pv);

    phoPFChIso_       .push_back(GEDIdTool.SolidConeIso(0.3, reco::PFCandidate::h));
    phoPFPhoIso_      .push_back(GEDIdTool.SolidConeIso(0.3, reco::PFCandidate::gamma));
    phoPFNeuIso_      .push_back(GEDIdTool.SolidConeIso(0.3, reco::PFCandidate::h0));

    std::vector<double> IsoRings;
    GEDIdTool.FrixioneIso(0.1, 8, reco::PFCandidate::h, IsoRings);
    phoPFChIsoFrix1_.push_back(IsoRings[0]);
    phoPFChIsoFrix2_.push_back(IsoRings[1]);
    phoPFChIsoFrix3_.push_back(IsoRings[2]);
    phoPFChIsoFrix4_.push_back(IsoRings[3]);
    phoPFChIsoFrix5_.push_back(IsoRings[4]);
    phoPFChIsoFrix6_.push_back(IsoRings[5]);
    phoPFChIsoFrix7_.push_back(IsoRings[6]);
    phoPFChIsoFrix8_.push_back(IsoRings[7]);
    IsoRings.resize(0);

    GEDIdTool.FrixioneIso(0.1, 8, reco::PFCandidate::gamma, IsoRings);
    phoPFPhoIsoFrix1_.push_back(IsoRings[0]);
    phoPFPhoIsoFrix2_.push_back(IsoRings[1]);
    phoPFPhoIsoFrix3_.push_back(IsoRings[2]);
    phoPFPhoIsoFrix4_.push_back(IsoRings[3]);
    phoPFPhoIsoFrix5_.push_back(IsoRings[4]);
    phoPFPhoIsoFrix6_.push_back(IsoRings[5]);
    phoPFPhoIsoFrix7_.push_back(IsoRings[6]);
    phoPFPhoIsoFrix8_.push_back(IsoRings[7]);
    IsoRings.resize(0);

    GEDIdTool.FrixioneIso(0.1, 8, reco::PFCandidate::h0, IsoRings);
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
    varRho_          = rho;
    varPhoIsoRaw_    = phoPFPhoIso_[nPho_];
    varChIsoRaw_     = phoPFChIso_[nPho_];
    varWorstChRaw_   = phoPFNeuIso_[nPho_];
    // Declare spectator vars
    varPt_           = phoEt_[nPho_];
    varEta_          = phoEta_[nPho_];

    // 0=ECAL barrel or 1=ECAL endcaps
    int iBE = (fabs(phoSCEta_[nPho_]) < 1.479) ? 0 : 1;

    phoIDMVA_.push_back(tmvaReader_[iBE]->EvaluateMVA("BDT"));

    nPho_++;
  }

}

void ggNtuplizer::cleanupPhotons()
{
  delete tmvaReader_[0];
  delete tmvaReader_[1];
}
