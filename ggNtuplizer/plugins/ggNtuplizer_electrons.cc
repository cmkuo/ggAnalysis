#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/EcalDetId/interface/ESDetId.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "ggAnalysis/ggNtuplizer/interface/ggNtuplizer.h"

using namespace std;
using namespace reco;

// (local) variables associated with tree branches
Int_t          nEle_;
vector<int>    eleCharge_;
vector<int>    eleChargeConsistent_;
vector<float>  eleEn_;
vector<float>  eleSCEn_;
vector<float>  eleESEn_;
vector<float>  eleESEnP1_;
vector<float>  eleESEnP2_;
vector<float>  eleESEnP1Raw_;
vector<float>  eleESEnP2Raw_;
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
vector<float>  eleEoverPout_;
vector<float>  eleEoverPInv_;
vector<float>  eleBrem_;
vector<float>  eledEtaAtVtx_;
vector<float>  eledPhiAtVtx_;
vector<float>  eledEtaAtCalo_;
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
vector<float>  eleIDMVATrg_;
vector<float>  eledEtaseedAtVtx_;
vector<float>  eleE1x5_;
vector<float>  eleE2x5_;
vector<float>  eleE5x5_;
vector<float>  eleE1x5Full5x5_;
vector<float>  eleE2x5Full5x5_;
vector<float>  eleE5x5Full5x5_;
vector<float>  eleR9Full5x5_;
vector<int>    eleEcalDrivenSeed_;
vector<float>  eleDr03EcalRecHitSumEt_;
vector<float>  eleDr03HcalDepth1TowerSumEt_;
vector<float>  eleDr03HcalDepth2TowerSumEt_;
vector<float>  eleDr03HcalTowerSumEt_;
vector<float>  eleDr03TkSumPt_;
vector<float>  elecaloEnergy_;
vector<float>  eleTrkdxy_;
vector<float>  eleKFHits_;
vector<float>  eleKFChi2_;
vector<float>  eleGSFChi2_;
vector<int>    eleFiredTrgs_;

vector<UShort_t> eleIDbit_;

vector<vector<float> > eleGSFPt_;
vector<vector<float> > eleGSFEta_;
vector<vector<float> > eleGSFPhi_;
vector<vector<float> > eleGSFCharge_;
vector<vector<int> >   eleGSFHits_;
vector<vector<int> >   eleGSFMissHits_;
vector<vector<int> >   eleGSFNHitsMax_;
vector<vector<float> > eleGSFVtxProb_;
vector<vector<float> > eleGSFlxyPV_;
vector<vector<float> > eleGSFlxyBS_;
vector<vector<float> > eleBCEn_;
vector<vector<float> > eleBCEta_;
vector<vector<float> > eleBCPhi_;
vector<vector<float> > eleBCS25_;
vector<vector<float> > eleBCS15_;
vector<vector<float> > eleBCSieie_;
vector<vector<float> > eleBCSieip_;
vector<vector<float> > eleBCSipip_;

Int_t nGSFTrk_;
vector<float> gsfPt_;
vector<float> gsfEta_;
vector<float> gsfPhi_;

void ggNtuplizer::branchesElectrons(TTree* tree) {

  tree->Branch("nEle",                    &nEle_);
  tree->Branch("eleCharge",               &eleCharge_);
  tree->Branch("eleChargeConsistent",     &eleChargeConsistent_);
  tree->Branch("eleEn",                   &eleEn_);
  tree->Branch("eleSCEn",                 &eleSCEn_);
  tree->Branch("eleESEn",                 &eleESEn_);
  tree->Branch("eleESEnP1",               &eleESEnP1_);
  tree->Branch("eleESEnP2",               &eleESEnP2_);
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
  tree->Branch("eleEoverPout",            &eleEoverPout_);
  tree->Branch("eleEoverPInv",            &eleEoverPInv_);
  tree->Branch("eleBrem",                 &eleBrem_);
  tree->Branch("eledEtaAtVtx",            &eledEtaAtVtx_);
  tree->Branch("eledPhiAtVtx",            &eledPhiAtVtx_);
  tree->Branch("eledEtaAtCalo",           &eledEtaAtCalo_);
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
  tree->Branch("eleIDMVATrg",             &eleIDMVATrg_);
  tree->Branch("eledEtaseedAtVtx",        &eledEtaseedAtVtx_);
  tree->Branch("eleE1x5",                 &eleE1x5_);
  tree->Branch("eleE2x5",                 &eleE2x5_);
  tree->Branch("eleE5x5",                 &eleE5x5_);
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
  tree->Branch("eleKFHits",                   &eleKFHits_);
  tree->Branch("eleKFChi2",                   &eleKFChi2_);
  tree->Branch("eleGSFPt",                    &eleGSFPt_);
  tree->Branch("eleGSFEta",                   &eleGSFEta_);
  tree->Branch("eleGSFPhi",                   &eleGSFPhi_);
  tree->Branch("eleGSFCharge",                &eleGSFCharge_);
  tree->Branch("eleGSFHits",                  &eleGSFHits_);
  tree->Branch("eleGSFMissHits",              &eleGSFMissHits_);
  tree->Branch("eleGSFNHitsMax",              &eleGSFNHitsMax_);
  tree->Branch("eleGSFVtxProb",               &eleGSFVtxProb_);
  tree->Branch("eleGSFlxyPV",                 &eleGSFlxyPV_);
  tree->Branch("eleGSFlxyBS",                 &eleGSFlxyBS_);
  tree->Branch("eleBCEn",                     &eleBCEn_);
  tree->Branch("eleBCEta",                    &eleBCEta_);
  tree->Branch("eleBCPhi",                    &eleBCPhi_);
  tree->Branch("eleBCS25",                    &eleBCS25_);
  tree->Branch("eleBCS15",                    &eleBCS15_);
  tree->Branch("eleBCSieie",                  &eleBCSieie_);
  tree->Branch("eleBCSieip",                  &eleBCSieip_);
  tree->Branch("eleBCSipip",                  &eleBCSipip_);
  tree->Branch("eleFiredTrgs",                &eleFiredTrgs_);

  if (runeleIDVID_) tree->Branch("eleIDbit",  &eleIDbit_);

  if (development_) {
    tree->Branch("eleESEnP1Raw",              &eleESEnP1Raw_);
    tree->Branch("eleESEnP2Raw",              &eleESEnP2Raw_);
    tree->Branch("nGSFTrk",                   &nGSFTrk_);
    tree->Branch("gsfPt",                     &gsfPt_);
    tree->Branch("gsfEta",                    &gsfEta_);
    tree->Branch("gsfPhi",                    &gsfPhi_);
  }
  
}

void ggNtuplizer::fillElectrons(const edm::Event &e, const edm::EventSetup &es, math::XYZPoint &pv) {
    
  // cleanup from previous execution
  eleCharge_                  .clear();
  eleChargeConsistent_        .clear();
  eleEn_                      .clear();
  eleSCEn_                    .clear();
  eleESEn_                    .clear();
  eleESEnP1_                  .clear();
  eleESEnP2_                  .clear();
  eleESEnP1Raw_               .clear();
  eleESEnP2Raw_               .clear();
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
  eleEoverPout_               .clear();
  eleEoverPInv_               .clear();
  eleBrem_                    .clear();
  eledEtaAtVtx_               .clear();
  eledPhiAtVtx_               .clear();
  eledEtaAtCalo_              .clear();
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
  eleIDMVATrg_                .clear();
  eledEtaseedAtVtx_           .clear();
  eleE1x5_                    .clear();
  eleE2x5_                    .clear();
  eleE5x5_                    .clear();
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
  eleKFHits_                  .clear();
  eleKFChi2_                  .clear();
  eleGSFChi2_                 .clear();
  eleGSFPt_                   .clear();
  eleGSFEta_                  .clear();
  eleGSFPhi_                  .clear();
  eleGSFCharge_               .clear();
  eleGSFHits_                 .clear();
  eleGSFMissHits_             .clear();
  eleGSFNHitsMax_             .clear();
  eleGSFVtxProb_              .clear();
  eleGSFlxyPV_                .clear();
  eleGSFlxyBS_                .clear();
  eleBCEn_                    .clear();
  eleBCEta_                   .clear();
  eleBCPhi_                   .clear();
  eleBCS25_                   .clear();
  eleBCS15_                   .clear();
  eleBCSieie_                 .clear();
  eleBCSieip_                 .clear();
  eleBCSipip_                 .clear();
  eleFiredTrgs_               .clear();
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
  edm::Handle<edm::ValueMap<float> > eleNonTrgMVAValues;
  edm::Handle<edm::ValueMap<float> > eleTrgMVAValues;

  if (runeleIDVID_) {
    e.getByToken(eleVetoIdMapToken_ ,         veto_id_decisions);
    e.getByToken(eleLooseIdMapToken_ ,        loose_id_decisions);
    e.getByToken(eleMediumIdMapToken_,        medium_id_decisions);
    e.getByToken(eleTightIdMapToken_,         tight_id_decisions);
    e.getByToken(eleHEEPIdMapToken_ ,         heep_id_decisions);
    e.getByToken(eleNonTrgMVAValuesMapToken_, eleNonTrgMVAValues);
    e.getByToken(eleTrgMVAValuesMapToken_,    eleTrgMVAValues);
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
    eleESEnP1_          .push_back(iEle->superCluster()->preshowerEnergyPlane1());
    eleESEnP2_          .push_back(iEle->superCluster()->preshowerEnergyPlane2());
    eleSCEta_           .push_back(iEle->superCluster()->eta());
    eleSCPhi_           .push_back(iEle->superCluster()->phi());
    eleSCRawEn_         .push_back(iEle->superCluster()->rawEnergy());
    eleSCEtaWidth_      .push_back(iEle->superCluster()->etaWidth());
    eleSCPhiWidth_      .push_back(iEle->superCluster()->phiWidth());
    eleHoverE_          .push_back(iEle->hcalOverEcal());

    eleFiredTrgs_       .push_back(matchElectronTriggerFilters(iEle->pt(), iEle->eta(), iEle->phi()));

    ///https://cmssdt.cern.ch/SDT/doxygen/CMSSW_7_2_2/doc/html/d8/dac/GsfElectron_8h_source.html
    eleEoverP_          .push_back(iEle->eSuperClusterOverP());
    eleEoverPout_       .push_back(iEle->eEleClusterOverPout());
    eleBrem_            .push_back(iEle->fbrem());
    eledEtaAtVtx_       .push_back(iEle->deltaEtaSuperClusterTrackAtVtx());
    eledPhiAtVtx_       .push_back(iEle->deltaPhiSuperClusterTrackAtVtx());
    eledEtaAtCalo_      .push_back(iEle->deltaEtaSeedClusterTrackAtCalo());
    eleSigmaIEtaIEta_   .push_back(iEle->sigmaIetaIeta()); ///new sigmaietaieta
    eleSigmaIEtaIPhi_   .push_back(iEle->sigmaIetaIphi());
    eleSigmaIPhiIPhi_   .push_back(iEle->sigmaIphiIphi());
    eleConvVeto_        .push_back((Int_t)iEle->passConversionVeto()); // ConvVtxFit || missHit == 0
    eleMissHits_        .push_back(iEle->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS));
    eleESEffSigmaRR_    .push_back(lazyTool.eseffsirir(*((*iEle).superCluster())));

    // VID calculation of (1/E - 1/p)
    if (iEle->ecalEnergy() == 0)   eleEoverPInv_.push_back(1e30);
    else if (!std::isfinite(iEle->ecalEnergy()))  eleEoverPInv_.push_back(1e30);
    else  eleEoverPInv_.push_back((1.0 - iEle->eSuperClusterOverP())/iEle->ecalEnergy());

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

    reco::GsfTrackRef gsfTrackRef = iEle->gsfTrack();
    if (iEle->gsfTrack().isNonnull()) {
      eleGSFChi2_.push_back(gsfTrackRef->normalizedChi2());
      if (recVtxs->size() > 0)
        eleTrkdxy_.push_back(gsfTrackRef->dxy(recVtxs->front().position()));
      else
	eleTrkdxy_.push_back(-999);
    } else {
      eleGSFChi2_.push_back(0.);
      eleTrkdxy_.push_back(-999);
    }
    
    reco::TrackRef kfTrackRef = iEle->closestCtfTrackRef();
    if (kfTrackRef.isAvailable() && kfTrackRef.isNonnull()) {
      eleKFHits_.push_back(kfTrackRef->hitPattern().trackerLayersWithMeasurement());
      eleKFChi2_.push_back(kfTrackRef->normalizedChi2());
    } else {
      eleKFHits_.push_back(-1.);
      eleKFChi2_.push_back(0.);
    }

    if (isAOD_) {

      edm::Handle<reco::BeamSpot> bsHandle;
      e.getByLabel("offlineBeamSpot", bsHandle);
      const reco::BeamSpot &beamspot = *bsHandle.product();

      edm::Handle<reco::ConversionCollection> hConversions;
      e.getByLabel("allConversions", hConversions);      

      // Close-by GSF tracks
      vector<float> eleGSFPt;
      vector<float> eleGSFEta;
      vector<float> eleGSFPhi;
      vector<float> eleGSFCharge;
      vector<int>   eleGSFHits;  
      vector<int>   eleGSFMissHits;
      vector<int>   eleGSFNHitsMax;
      vector<float> eleGSFVtxProb;
      vector<float> eleGSFlxyPV;
      vector<float> eleGSFlxyBS;
      eleGSFPt        .push_back(iEle->gsfTrack()->pt());
      eleGSFEta       .push_back(iEle->gsfTrack()->eta());
      eleGSFPhi       .push_back(iEle->gsfTrack()->phi());
      eleGSFCharge    .push_back(iEle->gsfTrack()->charge());
      eleGSFHits      .push_back(iEle->gsfTrack()->hitPattern().trackerLayersWithMeasurement());
      eleGSFMissHits  .push_back(iEle->gsfTrack()->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS));
      
      int   nHitsMax = -99;
      float lxyBS    = -99.;
      float lxyPV    = -99.;
      float vtxProb  = -99.;
      for (reco::ConversionCollection::const_iterator cnv = hConversions->begin(); cnv!= hConversions->end(); ++cnv) {
	reco::Vertex vtx1 = cnv->conversionVertex();
	if (vtx1.isValid()) {
	  if (ConversionTools::matchesConversion(iEle->gsfTrack(), *cnv)) {
	    vtxProb = TMath::Prob(vtx1.chi2(), vtx1.ndof());
	    
	    math::XYZVector mom1(cnv->refittedPairMomentum());
	    double dbsx1 = vtx1.x() - beamspot.position().x();   
	    double dbsy1 = vtx1.y() - beamspot.position().y();
	    lxyBS = (mom1.x()*dbsx1 + mom1.y()*dbsy1)/mom1.rho();
	    
	    double dpvx1 = vtx1.x() - pv.x();   
	    double dpvy1 = vtx1.y() - pv.y();
	    lxyPV = (mom1.x()*dpvx1 + mom1.y()*dpvy1)/mom1.rho();
	    
	    for (vector<uint8_t>::const_iterator it = cnv->nHitsBeforeVtx().begin(); it!=cnv->nHitsBeforeVtx().end(); ++it) {
	      if ((*it)>nHitsMax) nHitsMax = (*it);
	    }
	   
	    break; 
	  }
	}
      }
      eleGSFNHitsMax.push_back(nHitsMax);
      eleGSFlxyPV.push_back(lxyPV);
      eleGSFlxyBS.push_back(lxyBS);
      eleGSFVtxProb.push_back(vtxProb);
      
      for (GsfTrackRefVector::const_iterator igsf = iEle->ambiguousGsfTracksBegin(); igsf != iEle->ambiguousGsfTracksEnd(); ++igsf) {
	eleGSFPt        .push_back((*igsf)->pt());
	eleGSFEta       .push_back((*igsf)->eta());
	eleGSFPhi       .push_back((*igsf)->phi());
	eleGSFCharge    .push_back((*igsf)->charge());
	eleGSFHits      .push_back((*igsf)->hitPattern().trackerLayersWithMeasurement());
	eleGSFMissHits  .push_back((*igsf)->hitPattern().numberOfHits(reco::HitPattern::MISSING_INNER_HITS));

	nHitsMax = -99;
	lxyBS    = -99.;
	lxyPV    = -99.;
	vtxProb  = -99.;
	for (reco::ConversionCollection::const_iterator conv = hConversions->begin(); conv!= hConversions->end(); ++conv) {
	  reco::Vertex vtx = conv->conversionVertex();
	  if (vtx.isValid()) {
	    if (ConversionTools::matchesConversion(*igsf, *conv)) {
	      vtxProb = TMath::Prob( vtx.chi2(), vtx.ndof() );

	      math::XYZVector mom(conv->refittedPairMomentum());
	      double dbsx = vtx.x() - beamspot.position().x();   
	      double dbsy = vtx.y() - beamspot.position().y();
	      lxyBS = (mom.x()*dbsx + mom.y()*dbsy)/mom.rho();

	      double dpvx = vtx.x() - pv.x();   
	      double dpvy = vtx.y() - pv.y();
	      lxyPV = (mom.x()*dpvx + mom.y()*dpvy)/mom.rho();

	      for (std::vector<uint8_t>::const_iterator it = conv->nHitsBeforeVtx().begin(); it!=conv->nHitsBeforeVtx().end(); ++it) {
		if ((*it)>nHitsMax) nHitsMax = (*it);
	      }

	      break; 
	    }
	  }
	}
	eleGSFNHitsMax.push_back(nHitsMax);
	eleGSFlxyPV.push_back(lxyPV);
	eleGSFlxyBS.push_back(lxyBS);
	eleGSFVtxProb.push_back(vtxProb);
      }

      eleGSFPt_        .push_back(eleGSFPt);
      eleGSFEta_       .push_back(eleGSFEta);
      eleGSFPhi_       .push_back(eleGSFPhi);
      eleGSFCharge_    .push_back(eleGSFCharge);
      eleGSFHits_      .push_back(eleGSFHits);
      eleGSFMissHits_  .push_back(eleGSFMissHits);
      eleGSFNHitsMax_  .push_back(eleGSFNHitsMax);
      eleGSFlxyPV_     .push_back(eleGSFlxyPV);
      eleGSFlxyBS_     .push_back(eleGSFlxyBS);
      eleGSFVtxProb_   .push_back(eleGSFVtxProb);

      vector<float> eleBCEn;
      vector<float> eleBCEta;
      vector<float> eleBCPhi;
      vector<float> eleBCS25;
      vector<float> eleBCS15;
      vector<float> eleBCSieie;
      vector<float> eleBCSieip;
      vector<float> eleBCSipip;
      for (CaloCluster_iterator itbc = iEle->superCluster()->clustersBegin(); itbc != iEle->superCluster()->clustersEnd(); ++itbc) {
	eleBCEn   .push_back((*itbc)->energy());
	eleBCEta  .push_back((*itbc)->eta());
	eleBCPhi  .push_back((*itbc)->phi());
	eleBCS25  .push_back(lazyToolnoZS.e2x5Max(**itbc)/lazyToolnoZS.e5x5(**itbc));
	eleBCS15  .push_back(lazyToolnoZS.e1x5(**itbc)/lazyToolnoZS.e5x5(**itbc));
	vector<float> eleCov;
	eleCov.clear();
	eleCov = lazyToolnoZS.localCovariances(**itbc);
	eleBCSieie.push_back(sqrt(eleCov[0]));
	eleBCSieip.push_back(eleCov[1]);
	eleBCSipip.push_back(eleCov[2]);
      }
      eleBCEn_   .push_back(eleBCEn);
      eleBCEta_  .push_back(eleBCEta);
      eleBCPhi_  .push_back(eleBCPhi);
      eleBCS25_  .push_back(eleBCS25);
      eleBCS15_  .push_back(eleBCS15);
      eleBCSieie_.push_back(eleBCSieie);
      eleBCSieip_.push_back(eleBCSieip);
      eleBCSipip_.push_back(eleBCSipip);
    }

    if (development_) {

      Float_t ESp1 = 0;
      Float_t ESp2 = 0;
      for (CaloClusterPtrVector::const_iterator ips = iEle->superCluster()->preshowerClustersBegin(); ips != iEle->superCluster()->preshowerClustersEnd(); ++ips) {

	ESDetId esid = ESDetId((*ips)->seed());
	if (esid.plane() == 1) ESp1 += (*ips)->energy();
	if (esid.plane() == 2) ESp2 += (*ips)->energy();
      }

      eleESEnP1Raw_.push_back(ESp1);
      eleESEnP2Raw_.push_back(ESp2);
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

        eleIDMVANonTrg_.push_back((*eleNonTrgMVAValues)[el->originalObjectRef()]);
        eleIDMVATrg_   .push_back((*eleTrgMVAValues)[el->originalObjectRef()]);
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

        eleIDMVANonTrg_.push_back((*eleNonTrgMVAValues)[el]);
        eleIDMVATrg_   .push_back((*eleTrgMVAValues)[el]);
      }

      eleIDbit_.push_back(tmpeleIDbit);
      
      //cout<<"tmpele : eleIDbit: "<<tmpeleIDbit<<":"<<eleIDbit_[nEle_]<<endl;
    }//if(runeleIDVID_)

    nEle_++;
  }

  if (development_) {

    edm::Handle<edm::View<reco::GsfTrack> > GsfTrackHandle;
    e.getByToken(gsfTracks_, GsfTrackHandle);

    nGSFTrk_ = 0;
    gsfPt_ .clear();
    gsfEta_.clear();
    gsfPhi_.clear();

    for (edm::View<reco::GsfTrack>::const_iterator ig = GsfTrackHandle->begin(); ig != GsfTrackHandle->end(); ++ig) {
      gsfPt_ .push_back(ig->pt());
      gsfEta_.push_back(ig->eta());
      gsfPhi_.push_back(ig->phi());
      nGSFTrk_++;
    }



  }

}
