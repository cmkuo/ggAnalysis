
#include "../interface/ggNtuplizer.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"



using namespace std;
using namespace edm;

Int_t          nPho_;
vector<float>  phoE_;
vector<float>  phoEt_;
vector<float>  phoEta_;
vector<float>  phoPhi_;
vector<int>  pho_genmatched_;
vector<float>  phoSigmaE_;
vector<float>  phoCalibE_;
vector<float>  phoCalibEt_;
vector<float>  phoSCE_;
vector<float>  phoSCRawE_;
vector<float>  phoESEn_;
vector<float>  phoESEnP1_;
vector<float>  phoESEnP2_;
vector<float>  phoESOverRawEn_;
vector<float>  phoSCEta_;
vector<float>  phoSCPhi_;
vector<float>  phoSCEtaWidth_;
vector<float>  phoSCPhiWidth_;
vector<float>  phoSCBrem_;
vector<int>    phohasPixelSeed_;
vector<int>    phoEleVeto_;
vector<float>  phoR9_;
vector<float>  phoHoverE_;
vector<float>  phoESEffSigmaRR_;
vector<float>  phoSigmaIEtaIEtaFull5x5_;
vector<float>  phoSigmaIEtaIPhiFull5x5_;
vector<float>  phoSigmaIPhiIPhiFull5x5_;
vector<float>  phoE2x2Full5x5_;
vector<float>  phoE5x5Full5x5_;
vector<float>  phoR9Full5x5_;
vector<float>  phoPFChIso_;
vector<float>  phoPFChPVIso_;
vector<float>  phoPFPhoIso_;
vector<float>  phoPFNeuIso_;
vector<float>  phoPFChWorstVetoIso_;
vector<float>  phoPFChWorstIso_;
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
vector<int>  phoSeediEta_;
vector<int>  phoSeediPhi_;
vector<vector<float>> phoEnergyMatrix_5x5_;
vector<vector<float>> phoEnergyMatrix_7x7_;
vector<vector<float>> phoEnergyMatrix_9x9_;
vector<vector<float>> phoEnergyMatrix_11x11_;
vector<vector<float>> phoEnergyMatrix_15x15_;
vector<vector<float>> phoEnergyMatrix_25x25_;

//Necessary for the Photon Footprint removal
template <class T, class U>
bool isInFootprint(const T& thefootprint, const U& theCandidate) {
  for ( auto itr = thefootprint.begin(); itr != thefootprint.end(); ++itr ) {

    if( itr.key() == theCandidate.key() ) return true;
    
  }
  return false;
}

void ggNtuplizer::branchesPhotons(TTree* tree) {

  //////////////////////////////Photon Branches/////////////////////////////////////////
  tree->Branch("nPho",                    &nPho_);
  tree->Branch("phoE",                    &phoE_);
  tree->Branch("phoEt",                   &phoEt_);
  tree->Branch("phoEta",                  &phoEta_);
  tree->Branch("phoPhi",                  &phoPhi_);
  tree->Branch("pho_genmatched",          &pho_genmatched_);
  //tree->Branch("phoSigmaE",               &phoSigmaE_);
  //tree->Branch("phoCalibE",               &phoCalibE_);
  //tree->Branch("phoCalibEt",              &phoCalibEt_);
  tree->Branch("phoSCE",                  &phoSCE_);
  tree->Branch("phoSCRawE",               &phoSCRawE_);
  tree->Branch("phoESEn",                 &phoESEn_);
  tree->Branch("phoESEnP1",               &phoESEnP1_);
  tree->Branch("phoESEnP2",               &phoESEnP2_);
  tree->Branch("phoESOverRawEn",          &phoESOverRawEn_);
  tree->Branch("phoSCEta",                &phoSCEta_);
  tree->Branch("phoSCPhi",                &phoSCPhi_);
  tree->Branch("phoSCEtaWidth",           &phoSCEtaWidth_);
  tree->Branch("phoSCPhiWidth",           &phoSCPhiWidth_);
  tree->Branch("phoSCBrem",               &phoSCBrem_);
  tree->Branch("phohasPixelSeed",         &phohasPixelSeed_);
  tree->Branch("phoEleVeto",              &phoEleVeto_);
  tree->Branch("phoR9",                   &phoR9_);
  tree->Branch("phoHoverE",               &phoHoverE_);
  tree->Branch("phoESEffSigmaRR",         &phoESEffSigmaRR_);
  tree->Branch("phoSigmaIEtaIEtaFull5x5", &phoSigmaIEtaIEtaFull5x5_);
  tree->Branch("phoSigmaIEtaIPhiFull5x5", &phoSigmaIEtaIPhiFull5x5_);
  tree->Branch("phoSigmaIPhiIPhiFull5x5", &phoSigmaIPhiIPhiFull5x5_);
  tree->Branch("phoE2x2Full5x5",          &phoE2x2Full5x5_);
  tree->Branch("phoE5x5Full5x5",          &phoE5x5Full5x5_);
  tree->Branch("phoR9Full5x5",            &phoR9Full5x5_);
  tree->Branch("phoPFChIso",              &phoPFChIso_);
  tree->Branch("phoPFChPVIso",            &phoPFChPVIso_);
  tree->Branch("phoPFPhoIso",             &phoPFPhoIso_);
  tree->Branch("phoPFNeuIso",             &phoPFNeuIso_);
  tree->Branch("phoPFChWorstIso",         &phoPFChWorstIso_);
  tree->Branch("phoPFChWorstVetoIso",     &phoPFChWorstVetoIso_);
  tree->Branch("phoEcalPFClusterIso",     &phoEcalPFClusterIso_);
  tree->Branch("phoHcalPFClusterIso",     &phoHcalPFClusterIso_);
  if(dumpCrystalinfo_){
  tree->Branch("phoSeedTime",             &phoSeedTime_);
  tree->Branch("phoSeedEnergy",           &phoSeedEnergy_);
  tree->Branch("phoSeediEta",             &phoSeediEta_);
  tree->Branch("phoSeediPhi",             &phoSeediPhi_);
  tree->Branch("phoEnergyMatrix_5x5",     &phoEnergyMatrix_5x5_);
  tree->Branch("phoEnergyMatrix_7x7",     &phoEnergyMatrix_7x7_);
  tree->Branch("phoEnergyMatrix_9x9",     &phoEnergyMatrix_9x9_);
  tree->Branch("phoEnergyMatrix_11x11",     &phoEnergyMatrix_11x11_);
  tree->Branch("phoEnergyMatrix_15x15",     &phoEnergyMatrix_15x15_);
  tree->Branch("phoEnergyMatrix_25x25",   &phoEnergyMatrix_25x25_);}

}

//
// member functions
//
// ------------ method called for each event  ------------
void ggNtuplizer::fillPhotons(const edm::Event& iEvent, const edm::EventSetup& iSetup) {


  //////////////////////////////Cleaning of previous event/////////////////////////////////////////
  phoE_                   .clear();
  phoSigmaE_              .clear();
  phoEt_                  .clear();
  phoEta_                 .clear();
  phoPhi_                 .clear();
  phoCalibE_              .clear();
  phoCalibEt_             .clear();
  phoSCE_                 .clear();
  phoSCRawE_              .clear();
  phoESEn_                .clear();
  phoESOverRawEn_         .clear();
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
  phoESEffSigmaRR_        .clear();
  phoSigmaIEtaIEtaFull5x5_.clear();
  phoSigmaIEtaIPhiFull5x5_.clear();
  phoSigmaIPhiIPhiFull5x5_.clear();
  phoE2x2Full5x5_         .clear();
  phoE5x5Full5x5_         .clear();
  phoR9Full5x5_           .clear();
  phoPFChIso_             .clear();
  phoPFChPVIso_           .clear();
  phoPFPhoIso_            .clear();
  phoPFNeuIso_            .clear();
  phoPFChWorstVetoIso_    .clear();
  phoPFChWorstIso_        .clear();
  phoEcalPFClusterIso_    .clear();
  phoHcalPFClusterIso_    .clear();
  phoSeedTime_            .clear();
  phoSeedEnergy_          .clear();
  phoSeediEta_            .clear();
  phoSeediPhi_            .clear();
  phoEnergyMatrix_5x5_    .clear();
  phoEnergyMatrix_7x7_    .clear();
  phoEnergyMatrix_9x9_    .clear();
  phoEnergyMatrix_11x11_  .clear();
  phoEnergyMatrix_15x15_  .clear();
  phoEnergyMatrix_25x25_ .clear();

  pho_genmatched_         .clear();
  nPho_ = 0;

  //////////////////////////////Filling of Variables///////////////////////////////////

  edm::Handle<edm::View<pat::Photon> > photonHandle;
  iEvent.getByToken(photonCollection_, photonHandle);

  // Get the RecHits from the event for barrel
  Handle<EcalRecHitCollection> pEBRecHits;
  iEvent.getByToken(ebReducedRecHitCollection_, pEBRecHits);

  // Get the RecHits from the event for endcap
  Handle<EcalRecHitCollection> pEERecHits;
  iEvent.getByToken(eeReducedRecHitCollection_, pEERecHits);


  if (!photonHandle.isValid()) {
    edm::LogWarning("ggNtuplizer") << "no pat::Photons in event";
    return;
  }

  edm::Handle<vector<reco::GenParticle> > genParticlesHandle;
  iEvent.getByToken(genParticlesCollection_, genParticlesHandle);


  EcalClusterLazyTools  lazyTool(iEvent, ecalClusterToolsESGetTokens_.get(iSetup), ebReducedRecHitCollection_, eeReducedRecHitCollection_, esReducedRecHitCollection_);
  
  noZS::EcalClusterLazyTools lazyToolnoZS(iEvent, ecalClusterToolsESGetTokens_.get(iSetup), ebReducedRecHitCollection_, eeReducedRecHitCollection_, esReducedRecHitCollection_);

  for (edm::View<pat::Photon>::const_iterator iPho = photonHandle->begin(); iPho != photonHandle->end(); ++iPho) {

    int genmatched = 0;
        
    if (genParticlesHandle.isValid()) {
      for (std::vector<reco::GenParticle>::const_iterator ip = genParticlesHandle->begin(); ip != genParticlesHandle->end(); ++ip) {
	const reco::Candidate *p = (const reco::Candidate*)&(*ip);
	if ( (p->status()==1) && abs(p->pdgId()) == 22 ) {
	  if ( ((reco::deltaR(*iPho,*p))<0.1) && ((iPho->pt()/p->pt())<0.2) ) genmatched=1;
	}  
      }
    }
    pho_genmatched_.push_back(genmatched);
    
    phoE_               .push_back(iPho->energy());
    phoEt_              .push_back(iPho->et());
    ////AR Requested UserFloat ecalEnergyPostCorr is not available! Possible UserFloats are: 
    ////AR PhotonMVAEstimatorRun2Spring16NonTrigV1Values PhotonMVAEstimatorRunIIFall17v1p1Values PhotonMVAEstimatorRunIIFall17v2Values 
    //phoCalibE_          .push_back(iPho->userFloat("ecalEnergyPostCorr"));
    //phoCalibEt_         .push_back(iPho->et()*iPho->userFloat("ecalEnergyPostCorr")/iPho->energy());
    //phoSigmaE_          .push_back(iPho->userFloat("ecalEnergyErrPostCorr"));
    phoEta_             .push_back(iPho->eta());
    phoPhi_             .push_back(iPho->phi());
    phoSCE_             .push_back((*iPho).superCluster()->energy());
    phoSCRawE_          .push_back((*iPho).superCluster()->rawEnergy());
    phoESEn_            .push_back((*iPho).superCluster()->preshowerEnergy());
    phoESEnP1_          .push_back((*iPho).superCluster()->preshowerEnergyPlane1());
    phoESEnP2_          .push_back((*iPho).superCluster()->preshowerEnergyPlane2());
    phoESOverRawEn_     .push_back((*iPho).superCluster()->preshowerEnergy()/(*iPho).superCluster()->rawEnergy());
    phoSCEta_           .push_back((*iPho).superCluster()->eta());
    phoSCPhi_           .push_back((*iPho).superCluster()->phi());
    phoSCEtaWidth_      .push_back((*iPho).superCluster()->etaWidth());
    phoSCPhiWidth_      .push_back((*iPho).superCluster()->phiWidth());
    phoSCBrem_          .push_back((*iPho).superCluster()->phiWidth()/(*iPho).superCluster()->etaWidth());
    phohasPixelSeed_    .push_back((Int_t)iPho->hasPixelSeed());
    phoEleVeto_         .push_back((Int_t)iPho->passElectronVeto());
    phoR9_              .push_back(iPho->r9());
    phoHoverE_          .push_back(iPho->hadTowOverEm());
    phoESEffSigmaRR_    .push_back(lazyTool.eseffsirir(*((*iPho).superCluster())));
    phoSigmaIEtaIEtaFull5x5_ .push_back(iPho->full5x5_sigmaIetaIeta());
    phoSigmaIEtaIPhiFull5x5_ .push_back(iPho->full5x5_showerShapeVariables().sigmaIetaIphi);
    phoSigmaIPhiIPhiFull5x5_ .push_back(iPho->full5x5_showerShapeVariables().sigmaIphiIphi);
    phoE2x2Full5x5_          .push_back(lazyToolnoZS.e2x2(*((*iPho).superCluster()->seed())));
    phoE5x5Full5x5_          .push_back(iPho->full5x5_e5x5());
    phoR9Full5x5_            .push_back(iPho->full5x5_r9());    
    phoPFChIso_         .push_back(iPho->chargedHadronIso()); //charged hadron isolation with dxy,dz match to pv
    phoPFChPVIso_       .push_back(iPho->chargedHadronPFPVIso()); //only considers particles assigned to the primary vertex (PV) by particle flow, corresponds to <10_6 chargedHadronIso
    phoPFPhoIso_        .push_back(iPho->photonIso());
    phoPFNeuIso_        .push_back(iPho->neutralHadronIso());
    phoPFChWorstIso_    .push_back(iPho->chargedHadronWorstVtxIso()); //max charged hadron isolation when dxy/dz matching to given vtx
    phoPFChWorstVetoIso_.push_back(iPho->chargedHadronWorstVtxGeomVetoIso()); //as chargedHadronWorstVtxIso but an additional geometry based veto cone
    phoEcalPFClusterIso_.push_back(iPho->ecalPFClusterIso());
    phoHcalPFClusterIso_.push_back(iPho->hcalPFClusterIso());

    ///////////////////////////////////////////////////////////////////////
    //const auto &seedSC = *(iPho->superCluster()->seed());
    ///////////////////////////////////////////////////////////////////////
    
    const reco::SuperCluster& superClus = *((*iPho).superCluster());
    const reco::CaloCluster & seedCluster = *superClus.seed();

    DetId seedDetId = (iPho->superCluster()->seed()->hitsAndFractions())[0].first;
    
    bool isBarrel = seedDetId.subdetId() == EcalBarrel;

    const EcalRecHitCollection& rechits = (isBarrel? *pEBRecHits : *pEERecHits); 
        
    EcalRecHitCollection::const_iterator theSeedHit = rechits.find(seedDetId);
    if (theSeedHit != rechits.end()) {
      //std::cout<<"(*theSeedHit).time()"<<(*theSeedHit).time()<<"seed energy: "<<(*theSeedHit).energy()<<std::endl;  
      if(dumpCrystalinfo_){ 
      phoSeedTime_  .push_back((*theSeedHit).time());
      phoSeedEnergy_.push_back((*theSeedHit).energy());}
      
      
      if (isBarrel) {
	EBDetId det = theSeedHit->id();
        if(dumpCrystalinfo_){
	phoSeediEta_  .push_back(det.ieta());
	phoSeediPhi_  .push_back(det.iphi());}
      }
      else {
	EEDetId det = theSeedHit->id();
        if(dumpCrystalinfo_){
	phoSeediEta_  .push_back(det.ix());
	phoSeediPhi_  .push_back(det.iy());}
      }

    } else{
      if(dumpCrystalinfo_){
      phoSeedTime_  .push_back(-99.);
      phoSeedEnergy_.push_back(-99.);
      phoSeediEta_  .push_back(-999);
      phoSeediPhi_  .push_back(-999);}
    }
    
    
    if(dumpCrystalinfo_){
    phoEnergyMatrix_5x5_.push_back(lazyToolnoZS.energyMatrix(seedCluster,2));
    phoEnergyMatrix_7x7_.push_back(lazyToolnoZS.energyMatrix(seedCluster,3));
    phoEnergyMatrix_9x9_.push_back(lazyToolnoZS.energyMatrix(seedCluster,4));
    phoEnergyMatrix_11x11_.push_back(lazyToolnoZS.energyMatrix(seedCluster,5));
    phoEnergyMatrix_15x15_.push_back(lazyToolnoZS.energyMatrix(seedCluster,7));
    phoEnergyMatrix_25x25_.push_back(lazyToolnoZS.energyMatrix(seedCluster, 12));}
    
    
    nPho_++;
  }

}
