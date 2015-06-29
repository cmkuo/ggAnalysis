#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"

#include "ggAnalysis/ggNtuplizer/interface/ggNtuplizer.h"

using namespace std;

// (local) variables associated with tree branches
Int_t          nEle_;
vector<int>    eleCharge_;
vector<int>    eleChargeConsistent_;
vector<float>  eleEn_;
vector<float>  eleSCEn_;
vector<float>  eleESEn_;
vector<float>  eleD0_;
vector<float>  eleDz_;
vector<float>  elePt_;
vector<float>  eleEta_;
vector<float>  elePhi_;
vector<float>  eleR9_;
vector<float>  eleSCEta_;
vector<float>  eleSCPhi_;
vector<float>  eleSCRawEn_;
vector<float>  eleSCEtaWidth_;
vector<float>  eleSCPhiWidth_;
vector<float>  eleHoverE_;
vector<float>  eleEoverP_;
vector<float>  eleEoverPInv_;
vector<float>  eleBrem_;
vector<float>  eledEtaAtVtx_;
vector<float>  eledPhiAtVtx_;
vector<float>  eleSigmaIEtaIEta_;
vector<float>  eleSigmaIEtaIPhi_;
vector<float>  eleSigmaIPhiIPhi_;
vector<float>  eleSigmaIEtaIEtaFull5x5_;
vector<float>  eleSigmaIPhiIPhiFull5x5_;
vector<int>    eleConvVeto_;
vector<int>    eleMissHits_;
vector<float>  eleESEffSigmaRR_;
vector<float>  elePFChIso_;
vector<float>  elePFPhoIso_;
vector<float>  elePFNeuIso_;
vector<float>  elePFPUIso_;
vector<float>  eleIDMVANonTrg_;
vector<float>  eledEtaseedAtVtx_;
vector<float>  eleE1x5_;
vector<float>  eleE2x5_;
vector<float>  eleE5x5_;
vector<float>  eleE1x5Full5x5_;
vector<float>  eleE2x5Full5x5_;
vector<float>  eleE5x5Full5x5_;
vector<float>  eleR9Full5x5_;
vector<float>  eleRelIsoWithDBeta_;
vector<int>    eleEcalDrivenSeed_;
vector<float>  eleDr03EcalRecHitSumEt_;
vector<float>  eleDr03HcalDepth1TowerSumEt_;
vector<float>  eleDr03HcalDepth2TowerSumEt_;
vector<float>  eleDr03HcalTowerSumEt_;
vector<float>  eleDr03TkSumPt_;
vector<float>  elecaloEnergy_;
vector<float>  eleTrkdxy_;

//vector<ULong64_t> eleIDbit_;
vector<UShort_t> eleIDbit_;

void ggNtuplizer::branchesElectrons(TTree* tree) {


  tree->Branch("nEle",                    &nEle_);
  tree->Branch("eleCharge",               &eleCharge_);
  tree->Branch("eleChargeConsistent",     &eleChargeConsistent_);
  tree->Branch("eleEn",                   &eleEn_);
  tree->Branch("eleSCEn",                 &eleSCEn_);
  tree->Branch("eleESEn",                 &eleESEn_);
  tree->Branch("eleD0",                   &eleD0_);
  tree->Branch("eleDz",                   &eleDz_);
  tree->Branch("elePt",                   &elePt_);
  tree->Branch("eleEta",                  &eleEta_);
  tree->Branch("elePhi",                  &elePhi_);
  tree->Branch("eleR9",                   &eleR9_);
  tree->Branch("eleSCEta",                &eleSCEta_);
  tree->Branch("eleSCPhi",                &eleSCPhi_);
  tree->Branch("eleSCRawEn",              &eleSCRawEn_);
  tree->Branch("eleSCEtaWidth",           &eleSCEtaWidth_);
  tree->Branch("eleSCPhiWidth",           &eleSCPhiWidth_);
  tree->Branch("eleHoverE",               &eleHoverE_);
  tree->Branch("eleEoverP",               &eleEoverP_);
  tree->Branch("eleEoverPInv",            &eleEoverPInv_);
  tree->Branch("eleBrem",                 &eleBrem_);
  tree->Branch("eledEtaAtVtx",            &eledEtaAtVtx_);
  tree->Branch("eledPhiAtVtx",            &eledPhiAtVtx_);
  tree->Branch("eleSigmaIEtaIEta",        &eleSigmaIEtaIEta_);
  tree->Branch("eleSigmaIEtaIPhi",        &eleSigmaIEtaIPhi_);
  tree->Branch("eleSigmaIPhiIPhi",        &eleSigmaIPhiIPhi_);
  tree->Branch("eleSigmaIEtaIEtaFull5x5", &eleSigmaIEtaIEtaFull5x5_);
  tree->Branch("eleSigmaIPhiIPhiFull5x5", &eleSigmaIPhiIPhiFull5x5_);
  tree->Branch("eleConvVeto",             &eleConvVeto_);
  tree->Branch("eleMissHits",             &eleMissHits_);
  tree->Branch("eleESEffSigmaRR",         &eleESEffSigmaRR_);
  tree->Branch("elePFChIso",              &elePFChIso_);
  tree->Branch("elePFPhoIso",             &elePFPhoIso_);
  tree->Branch("elePFNeuIso",             &elePFNeuIso_);
  tree->Branch("elePFPUIso",              &elePFPUIso_);
  tree->Branch("eleIDMVANonTrg",          &eleIDMVANonTrg_);
  tree->Branch("eledEtaseedAtVtx",        &eledEtaseedAtVtx_);
  tree->Branch("eleE1x5",                 &eleE1x5_);
  tree->Branch("eleE2x5",                 &eleE2x5_);
  tree->Branch("eleE5x5",                 &eleE5x5_);
  tree->Branch("eleRelIsoWithDBeta",      &eleRelIsoWithDBeta_);
  tree->Branch("eleE1x5Full5x5",          &eleE1x5Full5x5_);
  tree->Branch("eleE2x5Full5x5",          &eleE2x5Full5x5_);
  tree->Branch("eleE5x5Full5x5",          &eleE5x5Full5x5_);
  tree->Branch("eleR9Full5x5",                &eleR9Full5x5_);
  tree->Branch("eleEcalDrivenSeed",           &eleEcalDrivenSeed_);
  tree->Branch("eleDr03EcalRecHitSumEt",      &eleDr03EcalRecHitSumEt_);
  tree->Branch("eleDr03HcalDepth1TowerSumEt", &eleDr03HcalDepth1TowerSumEt_);
  tree->Branch("eleDr03HcalDepth2TowerSumEt", &eleDr03HcalDepth2TowerSumEt_);
  tree->Branch("eleDr03HcalTowerSumEt",       &eleDr03HcalTowerSumEt_);
  tree->Branch("eleDr03TkSumPt",              &eleDr03TkSumPt_);
  tree->Branch("elecaloEnergy",               &elecaloEnergy_);
  tree->Branch("eleTrkdxy",                   &eleTrkdxy_);

  if(runeleIDVID_)
    tree->Branch("eleIDbit",                   &eleIDbit_);
  
}

void ggNtuplizer::fillElectrons(const edm::Event &e, const edm::EventSetup &es, math::XYZPoint &pv) {

    
  // cleanup from previous execution
  eleCharge_                  .clear();
  eleChargeConsistent_        .clear();
  eleEn_                      .clear();
  eleSCEn_                    .clear();
  eleESEn_                    .clear();
  eleD0_                      .clear();
  eleDz_                      .clear();
  elePt_                      .clear();
  eleEta_                     .clear();
  elePhi_                     .clear();
  eleR9_                      .clear();
  eleSCEta_                   .clear();
  eleSCPhi_                   .clear();
  eleSCRawEn_                 .clear();
  eleSCEtaWidth_              .clear();
  eleSCPhiWidth_              .clear();
  eleHoverE_                  .clear();
  eleEoverP_                  .clear();
  eleEoverPInv_               .clear();
  eleBrem_                    .clear();
  eledEtaAtVtx_               .clear();
  eledPhiAtVtx_               .clear();
  eleSigmaIEtaIEta_           .clear();
  eleSigmaIEtaIPhi_           .clear();
  eleSigmaIPhiIPhi_           .clear();
  eleSigmaIEtaIEtaFull5x5_    .clear();
  eleSigmaIPhiIPhiFull5x5_    .clear();
  eleConvVeto_                .clear();
  eleMissHits_                .clear();
  eleESEffSigmaRR_            .clear();
  elePFChIso_                 .clear();
  elePFPhoIso_                .clear();
  elePFNeuIso_                .clear();
  elePFPUIso_                 .clear();
  eleIDMVANonTrg_             .clear();
  eledEtaseedAtVtx_           .clear();
  eleE1x5_                    .clear();
  eleE2x5_                    .clear();
  eleE5x5_                    .clear();
  eleRelIsoWithDBeta_         .clear();
  eleEcalDrivenSeed_          .clear();
  eleDr03EcalRecHitSumEt_     .clear();
  eleDr03HcalDepth1TowerSumEt_.clear();
  eleDr03HcalDepth2TowerSumEt_.clear();
  eleDr03HcalTowerSumEt_      .clear();
  eleDr03TkSumPt_             .clear();
  eleE1x5Full5x5_             .clear();
  eleE2x5Full5x5_             .clear();
  eleE5x5Full5x5_             .clear();
  eleR9Full5x5_               .clear();
  elecaloEnergy_              .clear();
  eleTrkdxy_                  .clear();
  eleIDbit_                   .clear();

  nEle_ = 0;

  edm::Handle<edm::View<pat::Electron> > electronHandle;
  e.getByToken(electronCollection_, electronHandle);

  if (!electronHandle.isValid()) {
    edm::LogWarning("ggNtuplizer") << "no pat::Electrons in event";
    return;
  }

  edm::Handle<edm::ValueMap<bool> >  veto_id_decisions;
  edm::Handle<edm::ValueMap<bool> >  loose_id_decisions;
  edm::Handle<edm::ValueMap<bool> >  medium_id_decisions;
  edm::Handle<edm::ValueMap<bool> >  tight_id_decisions; 
  edm::Handle<edm::ValueMap<bool> >  heep_id_decisions;
  edm::Handle<edm::ValueMap<float> > eleMVAValues;

  if (runeleIDVID_) {
    e.getByToken(eleVetoIdMapToken_ ,   veto_id_decisions);
    e.getByToken(eleLooseIdMapToken_ ,  loose_id_decisions);
    e.getByToken(eleMediumIdMapToken_,  medium_id_decisions);
    e.getByToken(eleTightIdMapToken_,   tight_id_decisions);
    e.getByToken(eleHEEPIdMapToken_ ,   heep_id_decisions);
    e.getByToken(eleMVAValuesMapToken_, eleMVAValues);
  }

  edm::Handle<reco::VertexCollection> recVtxs;
  e.getByToken(vtxLabel_, recVtxs);
  
  EcalClusterLazyTools       lazyTool    (e, es, ebReducedRecHitCollection_, eeReducedRecHitCollection_, esReducedRecHitCollection_);
  noZS::EcalClusterLazyTools lazyToolnoZS(e, es, ebReducedRecHitCollection_, eeReducedRecHitCollection_, esReducedRecHitCollection_);
  
  for (edm::View<pat::Electron>::const_iterator iEle = electronHandle->begin(); iEle != electronHandle->end(); ++iEle) {

    eleCharge_          .push_back(iEle->charge());
    eleChargeConsistent_.push_back((Int_t)iEle->isGsfCtfScPixChargeConsistent());
    eleEn_              .push_back(iEle->energy());
    eleD0_              .push_back(iEle->gsfTrack()->dxy(pv));
    eleDz_              .push_back(iEle->gsfTrack()->dz(pv));
    elePt_              .push_back(iEle->pt());
    eleEta_             .push_back(iEle->eta());
    elePhi_             .push_back(iEle->phi());
    eleR9_              .push_back(iEle->r9());
    eleSCEn_            .push_back(iEle->superCluster()->energy());
    eleESEn_            .push_back(iEle->superCluster()->preshowerEnergy());
    eleSCEta_           .push_back(iEle->superCluster()->eta());
    eleSCPhi_           .push_back(iEle->superCluster()->phi());
    eleSCRawEn_         .push_back(iEle->superCluster()->rawEnergy());
    eleSCEtaWidth_      .push_back(iEle->superCluster()->etaWidth());
    eleSCPhiWidth_      .push_back(iEle->superCluster()->phiWidth());
    eleHoverE_          .push_back(iEle->hcalOverEcal());

    ///https://cmssdt.cern.ch/SDT/doxygen/CMSSW_7_2_2/doc/html/d8/dac/GsfElectron_8h_source.html
    eleEoverP_          .push_back(iEle->eSuperClusterOverP());
    eleBrem_            .push_back(iEle->fbrem());
    eledEtaAtVtx_       .push_back(iEle->deltaEtaSuperClusterTrackAtVtx());
    eledPhiAtVtx_       .push_back(iEle->deltaPhiSuperClusterTrackAtVtx());
    eleSigmaIEtaIEta_   .push_back(iEle->sigmaIetaIeta()); ///new sigmaietaieta
    eleSigmaIEtaIPhi_   .push_back(iEle->sigmaIetaIphi());
    eleSigmaIPhiIPhi_   .push_back(iEle->sigmaIphiIphi());
    eleConvVeto_        .push_back((Int_t)iEle->passConversionVeto()); // ConvVtxFit || missHit == 0
    eleMissHits_        .push_back(iEle->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS));
    eleESEffSigmaRR_    .push_back(lazyTool.eseffsirir(*((*iEle).superCluster())));

    if (iEle->ecalEnergy() == 0) {
      eleEoverPInv_.push_back(1e30);
    } else if (!std::isfinite(iEle->ecalEnergy())) {
      eleEoverPInv_.push_back(1e30);
    } else {
      eleEoverPInv_.push_back(1./iEle->ecalEnergy() - 1./iEle->p());
    }

    ///HEEP ID
    double eledEtaseedAtVtx = iEle->superCluster().isNonnull() && iEle->superCluster()->seed().isNonnull() ?
      iEle->deltaEtaSuperClusterTrackAtVtx() - iEle->superCluster()->eta() + iEle->superCluster()->seed()->eta() : std::numeric_limits<float>::max();

    eledEtaseedAtVtx_   .push_back(eledEtaseedAtVtx);

    eleE1x5_            .push_back(iEle->e1x5());
    eleE2x5_            .push_back(iEle->e2x5Max());
    eleE5x5_            .push_back(iEle->e5x5());

    reco::GsfElectron::PflowIsolationVariables pfIso = iEle->pfIsolationVariables();
    elePFChIso_         .push_back(pfIso.sumChargedHadronPt);
    elePFPhoIso_        .push_back(pfIso.sumPhotonEt);
    elePFNeuIso_        .push_back(pfIso.sumNeutralHadronEt);
    elePFPUIso_         .push_back(pfIso.sumPUPt);

    elecaloEnergy_      .push_back(iEle->caloEnergy());
    ///SJ - https://github.com/lgray/cmssw/blob/common_isolation_selection_70X/TestElectronID/ElectronIDAnalyzer/plugins/ElectronIDAnalyzer.cc

    float absiso = elePFChIso_[nEle_] + std::max(0.0 , elePFNeuIso_[nEle_] + elePFPhoIso_[nEle_] - 0.5 * elePFPUIso_[nEle_]);
    eleRelIsoWithDBeta_ .push_back(absiso/elePt_[nEle_]);

    /////quantities which were used for Run1 - these do not
    ///calculated through PF (meaning no energy is subtracted
    ///using PF)
    ///https://cmssdt.cern.ch/SDT/doxygen/CMSSW_7_2_2/doc/html/d9/d44/ElectronIDValueMapProducer_8cc_source.html
    ///line 120

    eleSigmaIEtaIEtaFull5x5_    .push_back(iEle->full5x5_sigmaIetaIeta());
    eleSigmaIPhiIPhiFull5x5_    .push_back(iEle->full5x5_sigmaIphiIphi());
    eleE1x5Full5x5_             .push_back(iEle->full5x5_e1x5());
    eleE2x5Full5x5_             .push_back(iEle->full5x5_e2x5Max());
    eleE5x5Full5x5_             .push_back(iEle->full5x5_e5x5());
    eleR9Full5x5_               .push_back(iEle->full5x5_r9());

    ///For HEEP ID
    eleEcalDrivenSeed_          .push_back(iEle->ecalDrivenSeed());
    eleDr03EcalRecHitSumEt_     .push_back(iEle->dr03EcalRecHitSumEt());
    eleDr03HcalDepth1TowerSumEt_.push_back(iEle->dr03HcalDepth1TowerSumEt());
    eleDr03HcalDepth2TowerSumEt_.push_back(iEle->dr03HcalDepth2TowerSumEt());
    eleDr03HcalTowerSumEt_      .push_back(iEle->dr03HcalTowerSumEt());
    eleDr03TkSumPt_             .push_back(iEle->dr03TkSumPt());

    reco::GsfTrackRef trackref = iEle->gsfTrack();
    if (iEle->gsfTrack().isNonnull()) {
      if (recVtxs->size() > 0)
        eleTrkdxy_        .push_back(trackref->dxy(recVtxs->front().position()));
    }
    
    //
    // Look up and save the ID decisions
    // 
    
    if (runeleIDVID_) {
      
      //edm::Ptr<reco::GsfElectron> recoEl(iEle);
      
      //const auto el = electrons->ptrAt(nEle_);
      const auto el = electronHandle->ptrAt(nEle_);
      
      UShort_t tmpeleIDbit = 0;

      if (isAOD_) {
	///el->electronID("cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-veto") also works
	bool isPassVeto  = (*veto_id_decisions)[el->originalObjectRef()];
	if(isPassVeto) setbit(tmpeleIDbit, 0);
	//cout<<"isVeto: "<<isPassVeto<<endl;
	
	bool isPassLoose  = (*loose_id_decisions)[el->originalObjectRef()];
	if(isPassLoose) setbit(tmpeleIDbit, 1);
	//cout<<"isLoose: "<<isPassLoose<<endl;
	
	bool isPassMedium = (*medium_id_decisions)[el->originalObjectRef()];
	if(isPassMedium) setbit(tmpeleIDbit, 2);
	//cout<<"isMedium: "<<isPassMedium<<endl;
	
	bool isPassTight  = (*tight_id_decisions)[el->originalObjectRef()];
	if(isPassTight) setbit(tmpeleIDbit, 3);
	//cout<<"isTight: "<<isPassTight<<endl;      

	bool isPassHEEP = (*heep_id_decisions)[el->originalObjectRef()];
	if(isPassHEEP) setbit(tmpeleIDbit, 4);
	//cout<<"isHeep: "<<isPassHEEP<<endl;

	eleIDMVANonTrg_.push_back((*eleMVAValues)[el->originalObjectRef()]);
      }

      if (!isAOD_) {
	///el->electronID("cutBasedElectronID-PHYS14-PU20bx25-V2-standalone-veto") also works
	bool isPassVeto  = (*veto_id_decisions)[el];
	if(isPassVeto) setbit(tmpeleIDbit, 0);
	//cout<<"isVeto: "<<isPassVeto<<endl;
	
	bool isPassLoose  = (*loose_id_decisions)[el];
	if(isPassLoose) setbit(tmpeleIDbit, 1);
	//cout<<"isLoose: "<<isPassLoose<<endl;
	
	bool isPassMedium = (*medium_id_decisions)[el];
	if(isPassMedium) setbit(tmpeleIDbit, 2);
	//cout<<"isMedium: "<<isPassMedium<<endl;
	
	bool isPassTight  = (*tight_id_decisions)[el];
	if(isPassTight) setbit(tmpeleIDbit, 3);
	//cout<<"isTight: "<<isPassTight<<endl;
      
	bool isPassHEEP = (*heep_id_decisions)[el];
	if(isPassHEEP) setbit(tmpeleIDbit, 4);
	//cout<<"isHeep: "<<isPassHEEP<<endl;

	eleIDMVANonTrg_.push_back((*eleMVAValues)[el]);
      }

      eleIDbit_.push_back(tmpeleIDbit);
      
      //cout<<"tmpele : eleIDbit: "<<tmpeleIDbit<<":"<<eleIDbit_[nEle_]<<endl;
    }//if(runeleIDVID_)

    nEle_++;
  }

}
