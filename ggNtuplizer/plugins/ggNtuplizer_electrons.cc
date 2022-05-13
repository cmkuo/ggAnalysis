
#include "../interface/ggNtuplizer.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
using namespace std;
using namespace edm;

// (local) variables associated with tree branches
Int_t          nEle_;
vector<int>    eleCharge_;
vector<int>    eleChargeConsistent_;
vector<float>  eleEn_;
vector<float>  eleSCEn_;
vector<float>  eleEcalEn_;
vector<float>  eleESEnP1_;
vector<float>  eleESEnP2_;
vector<float>  eleESEnP1Raw_;
vector<float>  eleESEnP2Raw_;
vector<float>  eleD0_;
vector<float>  eleDz_;
vector<float>  eleSIP_;
vector<float>  elePt_;
vector<float>  elePtError_;
vector<float>  eleEta_;
vector<float>  elePhi_;
vector<float>  eleR9_;
vector<float>  eleCalibPt_;
vector<float>  eleCalibEn_;
vector<float>  eleSCEta_;
vector<float>  eleSCPhi_;
vector<float>  eleSCRawEn_;
vector<float>  eleSCEtaWidth_;
vector<float>  eleSCPhiWidth_;
vector<float>  eleHoverE_;
vector<float>  eleEoverP_;
vector<float>  eleEoverPout_;
vector<float>  eleEoverPInv_;
vector<float>  eleBrem_;
vector<float>  eledEtaAtVtx_;
vector<float>  eledPhiAtVtx_;
vector<float>  eleSigmaIEtaIEtaFull5x5_;
vector<float>  eleSigmaIPhiIPhiFull5x5_;
vector<int>    eleConvVeto_;
vector<int>    eleMissHits_;
vector<float>  eleESEffSigmaRR_;
vector<float>  elePFChIso_;
vector<float>  elePFPhoIso_;
vector<float>  elePFNeuIso_;
vector<float>  elePFPUIso_;
vector<float>  elePFClusEcalIso_;
vector<float>  elePFClusHcalIso_;
vector<float>  eleIDMVAIso_;
vector<float>  eleIDMVANoIso_;
vector<float>  eleR9Full5x5_;
vector<int>    eleEcalDrivenSeed_;
vector<int>  ele_genmatched_;
vector<vector<float>> eleEnergyMatrix_5x5_;
vector<vector<float>> eleEnergyMatrix_7x7_;
vector<vector<float>> eleEnergyMatrix_9x9_;
vector<vector<float>> eleEnergyMatrix_11x11_;
vector<vector<float>> eleEnergyMatrix_15x15_;
vector<vector<float>> eleEnergyMatrix_25x25_;


void ggNtuplizer::branchesElectrons(TTree* tree) {

  tree->Branch("nEle",                    &nEle_);
  tree->Branch("eleCharge",               &eleCharge_);
  tree->Branch("eleChargeConsistent",     &eleChargeConsistent_);
  tree->Branch("eleEn",                   &eleEn_);
  tree->Branch("eleSCEn",                 &eleSCEn_);
  tree->Branch("eleEcalEn",               &eleEcalEn_);
  tree->Branch("eleESEnP1",               &eleESEnP1_);
  tree->Branch("eleESEnP2",               &eleESEnP2_);
  tree->Branch("eleD0",                   &eleD0_);
  tree->Branch("eleDz",                   &eleDz_);
  tree->Branch("eleSIP",                  &eleSIP_);
  tree->Branch("elePt",                   &elePt_);
  tree->Branch("eleEta",                  &eleEta_);
  tree->Branch("elePhi",                  &elePhi_);
  tree->Branch("eleR9",                   &eleR9_);
  tree->Branch("ele_genmatched",         &ele_genmatched_);
  //tree->Branch("elePtError",              &elePtError_);
  //tree->Branch("eleCalibPt",              &eleCalibPt_);
  //tree->Branch("eleCalibEn",              &eleCalibEn_);
  tree->Branch("eleSCEta",                &eleSCEta_);
  tree->Branch("eleSCPhi",                &eleSCPhi_);
  tree->Branch("eleSCRawEn",              &eleSCRawEn_);
  tree->Branch("eleSCEtaWidth",           &eleSCEtaWidth_);
  tree->Branch("eleSCPhiWidth",           &eleSCPhiWidth_);
  tree->Branch("eleHoverE",               &eleHoverE_);
  tree->Branch("elePFClusEcalIso",        &elePFClusEcalIso_);
  tree->Branch("elePFClusHcalIso",        &elePFClusHcalIso_);

  tree->Branch("eleSigmaIEtaIEtaFull5x5", &eleSigmaIEtaIEtaFull5x5_);
  tree->Branch("eleSigmaIPhiIPhiFull5x5", &eleSigmaIPhiIPhiFull5x5_);
  tree->Branch("elePFChIso",              &elePFChIso_);
  tree->Branch("elePFPhoIso",             &elePFPhoIso_);
  tree->Branch("elePFNeuIso",             &elePFNeuIso_);
  tree->Branch("elePFPUIso",              &elePFPUIso_);
  tree->Branch("eleR9Full5x5",                &eleR9Full5x5_);
  tree->Branch("eleEcalDrivenSeed",           &eleEcalDrivenSeed_);
  if(dumpCrystalinfo_){
  tree->Branch("eleEnergyMatrix_5x5",     &eleEnergyMatrix_5x5_);
  tree->Branch("eleEnergyMatrix_7x7",     &eleEnergyMatrix_7x7_);
  tree->Branch("eleEnergyMatrix_9x9",     &eleEnergyMatrix_9x9_);
  tree->Branch("eleEnergyMatrix_11x11",   &eleEnergyMatrix_11x11_);
  tree->Branch("eleEnergyMatrix_15x15",   &eleEnergyMatrix_15x15_);
  tree->Branch("eleEnergyMatrix_25x25",   &eleEnergyMatrix_25x25_);}

}

//
// member functions
//

// ------------ method called for each event  ------------
void ggNtuplizer::fillElectrons(const edm::Event &iEvent, const edm::EventSetup &iSetup, math::XYZPoint &pv) {

  // cleanup of previous event
  eleCharge_                  .clear();
  eleChargeConsistent_        .clear();
  eleEn_                      .clear();
  eleSCEn_                    .clear();
  eleEcalEn_                  .clear();
  eleESEnP1_                  .clear();
  eleESEnP2_                  .clear();
  eleESEnP1Raw_               .clear();
  eleESEnP2Raw_               .clear();
  ele_genmatched_             .clear();
  /*  
  eleESEnEta_                 .clear();
  eleESEnPhi_                 .clear();
  eleESEnE_                   .clear();
  eleESEnZ_                   .clear();
  eleESEnP_                   .clear();
  eleESEnX_                   .clear();
  eleESEnY_                   .clear();
  eleESEnS_                   .clear();
  */
  eleD0_                      .clear();
  eleDz_                      .clear();
  eleSIP_                     .clear();
  elePt_                      .clear();
  elePtError_                 .clear();
  eleEta_                     .clear();
  elePhi_                     .clear();
  eleR9_                      .clear();
  eleCalibPt_                 .clear();
  eleCalibEn_                 .clear();
  eleSCEta_                   .clear();
  eleSCPhi_                   .clear();
  eleSCRawEn_                 .clear();
  eleSCEtaWidth_              .clear();
  eleSCPhiWidth_              .clear();
  eleHoverE_                  .clear();
  elePFClusEcalIso_           .clear();
  elePFClusHcalIso_           .clear();

  eleSigmaIEtaIEtaFull5x5_    .clear();
  eleSigmaIPhiIPhiFull5x5_    .clear();
  elePFChIso_                 .clear();
  elePFPhoIso_                .clear();
  elePFNeuIso_                .clear();
  elePFPUIso_                 .clear();
  eleR9Full5x5_               .clear();
  eleEcalDrivenSeed_          .clear();
  eleEnergyMatrix_5x5_        .clear();
  eleEnergyMatrix_7x7_        .clear();
  eleEnergyMatrix_9x9_        .clear();
  eleEnergyMatrix_11x11_      .clear(); 
  eleEnergyMatrix_15x15_      .clear();
  eleEnergyMatrix_25x25_       .clear();
  nEle_=0;

  //////////////////////////////Filling of Variables///////////////////////////////////
  edm::Handle<edm::View<pat::Electron> > electronHandle;
  iEvent.getByToken(electronCollection_, electronHandle);

  //edm::Handle<pat::PackedCandidateCollection> pfcands;
  //iEvent.getByToken(pckPFCandidateCollection_, pfcands);


  // Get the RecHits from the event for barrel                                                                                                
  //Handle<EcalRecHitCollection> pEBRecHits;
  //iEvent.getByToken(ebReducedRecHitCollection_, pEBRecHits);

  // Get the RecHits from the event for endcap                                                                                                
  //Handle<EcalRecHitCollection> pEERecHits;
  //iEvent.getByToken(eeReducedRecHitCollection_, pEERecHits);


  edm::Handle<vector<reco::GenParticle> > genParticlesHandle;
  iEvent.getByToken(genParticlesCollection_, genParticlesHandle);

  EcalClusterLazyTools lazyTool(iEvent, ecalClusterToolsESGetTokens_.get(iSetup), ebReducedRecHitCollection_, eeReducedRecHitCollection_, esReducedRecHitCollection_);
  noZS::EcalClusterLazyTools lazyToolnoZS(iEvent, ecalClusterToolsESGetTokens_.get(iSetup), ebReducedRecHitCollection_, eeReducedRecHitCollection_, esReducedRecHitCollection_);

  if (!electronHandle.isValid()) {
    edm::LogWarning("ggNtuplizer") << "no pat::Electrons in event";
    return;
  }

  for (edm::View<pat::Electron>::const_iterator iEle = electronHandle->begin(); iEle != electronHandle->end(); ++iEle) {


    int genmatched=0;
        
    if (genParticlesHandle.isValid()) {
      for (std::vector<reco::GenParticle>::const_iterator ip = genParticlesHandle->begin(); ip != genParticlesHandle->end(); ++ip) {
	const reco::Candidate *p = (const reco::Candidate*)&(*ip);
	if ( (p->status()==1) && abs(p->pdgId()) == 11 ) {
	  if ( ((reco::deltaR(*iEle,*p))<0.1) && ((iEle->pt()/p->pt())<0.2) ) genmatched=1;
	}  
      }
    }
    ele_genmatched_.push_back(genmatched);


    eleCharge_          .push_back(iEle->charge());
    eleChargeConsistent_.push_back((Int_t)iEle->isGsfCtfScPixChargeConsistent());
    eleEn_              .push_back(iEle->energy());
    eleD0_              .push_back(iEle->gsfTrack()->dxy(pv));
    eleDz_              .push_back(iEle->gsfTrack()->dz(pv));
    eleSIP_             .push_back(fabs(iEle->dB(pat::Electron::PV3D))/iEle->edB(pat::Electron::PV3D));
    elePt_              .push_back(iEle->pt());
    //elePtError_         .push_back(iEle->userFloat("ecalTrkEnergyErrPostCorr")*iEle->pt()/iEle->p());
    //eleCalibEn_         .push_back(iEle->userFloat("ecalEnergyPostCorr"));
    //eleCalibPt_         .push_back(iEle->userFloat("ecalTrkEnergyPostCorr")*iEle->pt()/iEle->p());
    eleEta_             .push_back(iEle->eta());
    elePhi_             .push_back(iEle->phi());
    eleR9_              .push_back(iEle->r9());
    eleSCEn_            .push_back(iEle->superCluster()->energy());
    eleEcalEn_          .push_back(iEle->ecalEnergy());
    eleESEnP1_          .push_back(iEle->superCluster()->preshowerEnergyPlane1());
    eleESEnP2_          .push_back(iEle->superCluster()->preshowerEnergyPlane2());
    eleSCEta_           .push_back(iEle->superCluster()->eta());
    eleSCPhi_           .push_back(iEle->superCluster()->phi());
    eleSCRawEn_         .push_back(iEle->superCluster()->rawEnergy());
    eleSCEtaWidth_      .push_back(iEle->superCluster()->etaWidth());
    eleSCPhiWidth_      .push_back(iEle->superCluster()->phiWidth());
    eleHoverE_          .push_back(iEle->hcalOverEcal());
    elePFClusEcalIso_   .push_back(iEle->ecalPFClusterIso());
    elePFClusHcalIso_   .push_back(iEle->hcalPFClusterIso());

    reco::GsfElectron::PflowIsolationVariables pfIso = iEle->pfIsolationVariables();
    elePFChIso_         .push_back(pfIso.sumChargedHadronPt);
    elePFPhoIso_        .push_back(pfIso.sumPhotonEt);
    elePFNeuIso_        .push_back(pfIso.sumNeutralHadronEt);
    elePFPUIso_         .push_back(pfIso.sumPUPt);

    eleSigmaIEtaIEtaFull5x5_.push_back(iEle->full5x5_sigmaIetaIeta());
    eleSigmaIPhiIPhiFull5x5_.push_back(iEle->full5x5_sigmaIphiIphi());
    eleR9Full5x5_           .push_back(iEle->full5x5_r9());
    eleEcalDrivenSeed_      .push_back(iEle->ecalDrivenSeed());


    ///////////////////////////////////////////////////////////
    //const auto &seedSC = *(iEle->superCluster()->seed());
    ///////////////////////////////////////////////////////////

    const reco::SuperCluster& superClus = *((*iEle).superCluster());
    const reco::CaloCluster &seedCluster = *superClus.seed();

    
    if(dumpCrystalinfo_){
    eleEnergyMatrix_5x5_.push_back(lazyToolnoZS.energyMatrix(seedCluster,2));
    eleEnergyMatrix_7x7_.push_back(lazyToolnoZS.energyMatrix(seedCluster,3));
    eleEnergyMatrix_9x9_.push_back(lazyToolnoZS.energyMatrix(seedCluster,4));
    eleEnergyMatrix_11x11_.push_back(lazyToolnoZS.energyMatrix(seedCluster,5));
    eleEnergyMatrix_15x15_.push_back(lazyToolnoZS.energyMatrix(seedCluster,7));
    eleEnergyMatrix_25x25_.push_back(lazyToolnoZS.energyMatrix(seedCluster,12));}

    nEle_++;
  }
}

