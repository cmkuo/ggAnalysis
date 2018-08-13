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
vector<float>  eleTrkdxy_;
vector<float>  eleKFHits_;
vector<float>  eleKFChi2_;
vector<float>  eleGSFChi2_;
vector<ULong64_t> eleFiredSingleTrgs_;
vector<ULong64_t> eleFiredDoubleTrgs_;
vector<ULong64_t> eleFiredL1Trgs_;
vector<UShort_t>  eleIDbit_;
vector<float>  eleScale_stat_up_;
vector<float>  eleScale_stat_dn_;
vector<float>  eleScale_syst_up_;
vector<float>  eleScale_syst_dn_;
vector<float>  eleScale_gain_up_;
vector<float>  eleScale_gain_dn_;
vector<float>  eleResol_rho_up_;
vector<float>  eleResol_rho_dn_;
vector<float>  eleResol_phi_up_;
vector<float>  eleResol_phi_dn_;

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

vector<vector<float> > eleESEnEta_;
vector<vector<float> > eleESEnPhi_;
vector<vector<int> >   eleESEnZ_;
vector<vector<int> >   eleESEnP_;
vector<vector<int> >   eleESEnX_;
vector<vector<int> >   eleESEnY_;
vector<vector<int> >   eleESEnS_;
vector<vector<float> > eleESEnE_;

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
  tree->Branch("eleEcalEn",               &eleEcalEn_);
  tree->Branch("eleESEnP1",               &eleESEnP1_);
  tree->Branch("eleESEnP2",               &eleESEnP2_);
  tree->Branch("eleD0",                   &eleD0_);
  tree->Branch("eleDz",                   &eleDz_);
  tree->Branch("eleSIP",                  &eleSIP_);
  tree->Branch("elePt",                   &elePt_);
  tree->Branch("elePtError",              &elePtError_);
  tree->Branch("eleEta",                  &eleEta_);
  tree->Branch("elePhi",                  &elePhi_);
  tree->Branch("eleR9",                   &eleR9_);
  tree->Branch("eleCalibPt",              &eleCalibPt_);
  tree->Branch("eleCalibEn",              &eleCalibEn_);
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
  tree->Branch("eleSigmaIEtaIEtaFull5x5", &eleSigmaIEtaIEtaFull5x5_);
  tree->Branch("eleSigmaIPhiIPhiFull5x5", &eleSigmaIPhiIPhiFull5x5_);
  tree->Branch("eleConvVeto",             &eleConvVeto_);
  tree->Branch("eleMissHits",             &eleMissHits_);
  tree->Branch("eleESEffSigmaRR",         &eleESEffSigmaRR_);
  tree->Branch("elePFChIso",              &elePFChIso_);
  tree->Branch("elePFPhoIso",             &elePFPhoIso_);
  tree->Branch("elePFNeuIso",             &elePFNeuIso_);
  tree->Branch("elePFPUIso",              &elePFPUIso_);
  tree->Branch("elePFClusEcalIso",        &elePFClusEcalIso_);
  tree->Branch("elePFClusHcalIso",        &elePFClusHcalIso_);
  tree->Branch("eleIDMVAIso",             &eleIDMVAIso_);
  tree->Branch("eleIDMVANoIso",           &eleIDMVANoIso_);
  tree->Branch("eleR9Full5x5",                &eleR9Full5x5_);
  tree->Branch("eleEcalDrivenSeed",           &eleEcalDrivenSeed_);
  tree->Branch("eleTrkdxy",                   &eleTrkdxy_);
  tree->Branch("eleKFHits",                   &eleKFHits_);
  tree->Branch("eleKFChi2",                   &eleKFChi2_);
  tree->Branch("eleGSFChi2",                  &eleGSFChi2_);
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
  tree->Branch("eleFiredSingleTrgs",          &eleFiredSingleTrgs_);
  tree->Branch("eleFiredDoubleTrgs",          &eleFiredDoubleTrgs_);
  tree->Branch("eleFiredL1Trgs",              &eleFiredL1Trgs_);
  tree->Branch("eleIDbit",                    &eleIDbit_);
  tree->Branch("eleScale_stat_up",            &eleScale_stat_up_);
  tree->Branch("eleScale_stat_dn",            &eleScale_stat_dn_);
  tree->Branch("eleScale_syst_up",            &eleScale_syst_up_);
  tree->Branch("eleScale_syst_dn",            &eleScale_syst_dn_);
  tree->Branch("eleScale_gain_up",            &eleScale_gain_up_);
  tree->Branch("eleScale_gain_dn",            &eleScale_gain_dn_);
  tree->Branch("eleResol_rho_up",             &eleResol_rho_up_);
  tree->Branch("eleResol_rho_dn",             &eleResol_rho_dn_);
  tree->Branch("eleResol_phi_up",             &eleResol_phi_up_);
  tree->Branch("eleResol_phi_dn",             &eleResol_phi_dn_);

  if (development_) {
    tree->Branch("eleESEnP1Raw",              &eleESEnP1Raw_);
    tree->Branch("eleESEnP2Raw",              &eleESEnP2Raw_);
    tree->Branch("eleESEnEta",                &eleESEnEta_);
    tree->Branch("eleESEnPhi",                &eleESEnPhi_);
    tree->Branch("eleESEnZ",                  &eleESEnZ_);
    tree->Branch("eleESEnP",                  &eleESEnP_);
    tree->Branch("eleESEnX",                  &eleESEnX_);
    tree->Branch("eleESEnY",                  &eleESEnY_);
    tree->Branch("eleESEnS",                  &eleESEnS_);
    tree->Branch("eleESEnE",                  &eleESEnE_);
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
  eleEcalEn_                  .clear();
  eleESEnP1_                  .clear();
  eleESEnP2_                  .clear();
  eleESEnP1Raw_               .clear();
  eleESEnP2Raw_               .clear();
  eleESEnEta_                 .clear();
  eleESEnPhi_                 .clear();
  eleESEnE_                   .clear();
  eleESEnZ_                   .clear();
  eleESEnP_                   .clear();
  eleESEnX_                   .clear();
  eleESEnY_                   .clear();
  eleESEnS_                   .clear();
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
  eleEoverP_                  .clear();
  eleEoverPout_               .clear();
  eleEoverPInv_               .clear();
  eleBrem_                    .clear();
  eledEtaAtVtx_               .clear();
  eledPhiAtVtx_               .clear();
  eleSigmaIEtaIEtaFull5x5_    .clear();
  eleSigmaIPhiIPhiFull5x5_    .clear();
  eleConvVeto_                .clear();
  eleMissHits_                .clear();
  eleESEffSigmaRR_            .clear();
  elePFChIso_                 .clear();
  elePFPhoIso_                .clear();
  elePFNeuIso_                .clear();
  elePFPUIso_                 .clear();
  elePFClusEcalIso_           .clear();
  elePFClusHcalIso_           .clear();
  eleIDMVAIso_                .clear();
  eleIDMVANoIso_              .clear();
  eleEcalDrivenSeed_          .clear();
  eleR9Full5x5_               .clear();
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
  eleFiredSingleTrgs_         .clear();
  eleFiredDoubleTrgs_         .clear();
  eleFiredL1Trgs_             .clear();
  eleIDbit_                   .clear();
  eleScale_stat_up_           .clear();
  eleScale_stat_dn_           .clear();
  eleScale_syst_up_           .clear();
  eleScale_syst_dn_           .clear();
  eleScale_gain_up_           .clear();
  eleScale_gain_dn_           .clear();
  eleResol_rho_up_            .clear();
  eleResol_rho_dn_            .clear();
  eleResol_phi_up_            .clear();
  eleResol_phi_dn_            .clear();

  nEle_ = 0;

  edm::Handle<edm::View<pat::Electron> > electronHandle;
  e.getByToken(electronCollection_, electronHandle);

  edm::Handle<pat::PackedCandidateCollection> pfcands;
  e.getByToken(pckPFCandidateCollection_, pfcands);

  if (!electronHandle.isValid()) {
    edm::LogWarning("ggNtuplizer") << "no pat::Electrons in event";
    return;
  }

  edm::Handle<reco::VertexCollection> recVtxs;
  e.getByToken(vtxLabel_, recVtxs);

  EcalClusterLazyTools       lazyTool    (e, es, ebReducedRecHitCollection_, eeReducedRecHitCollection_, esReducedRecHitCollection_);
  noZS::EcalClusterLazyTools lazyToolnoZS(e, es, ebReducedRecHitCollection_, eeReducedRecHitCollection_, esReducedRecHitCollection_);

  for (edm::View<pat::Electron>::const_iterator iEle = electronHandle->begin(); iEle != electronHandle->end(); ++iEle) {

    eleCharge_          .push_back(iEle->charge());
    eleChargeConsistent_.push_back((Int_t)iEle->isGsfCtfScPixChargeConsistent());
    eleEn_              .push_back(iEle->energy());
    eleCalibEn_         .push_back(iEle->userFloat("ecalEnergyPostCorr"));
    eleD0_              .push_back(iEle->gsfTrack()->dxy(pv));
    eleDz_              .push_back(iEle->gsfTrack()->dz(pv));
    eleSIP_             .push_back(fabs(iEle->dB(pat::Electron::PV3D))/iEle->edB(pat::Electron::PV3D));
    elePt_              .push_back(iEle->pt());
    eleCalibPt_         .push_back(iEle->userFloat("ecalTrkEnergyPostCorr")*iEle->pt()/iEle->p());
    elePtError_         .push_back(iEle->userFloat("ecalTrkEnergyErrPostCorr")*iEle->pt()/iEle->p());
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

    eleFiredSingleTrgs_ .push_back(matchSingleElectronTriggerFilters(iEle->pt(), iEle->eta(), iEle->phi()));
    eleFiredDoubleTrgs_ .push_back(matchDoubleElectronTriggerFilters(iEle->pt(), iEle->eta(), iEle->phi()));
    eleFiredL1Trgs_     .push_back(matchL1TriggerFilters(iEle->pt(), iEle->eta(), iEle->phi()));

    ///https://cmssdt.cern.ch/SDT/doxygen/CMSSW_7_2_2/doc/html/d8/dac/GsfElectron_8h_source.html
    eleEoverP_          .push_back(iEle->eSuperClusterOverP());
    eleEoverPout_       .push_back(iEle->eEleClusterOverPout());
    eleBrem_            .push_back(iEle->fbrem());
    eledEtaAtVtx_       .push_back(iEle->deltaEtaSuperClusterTrackAtVtx());
    eledPhiAtVtx_       .push_back(iEle->deltaPhiSuperClusterTrackAtVtx());
    eleConvVeto_        .push_back((Int_t)iEle->passConversionVeto()); // ConvVtxFit || missHit == 0
    eleMissHits_        .push_back(iEle->gsfTrack()->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS));
    eleESEffSigmaRR_    .push_back(lazyTool.eseffsirir(*((*iEle).superCluster())));

    // VID calculation of (1/E - 1/p)
    if (iEle->ecalEnergy() == 0)   eleEoverPInv_.push_back(1e30);
    else if (!std::isfinite(iEle->ecalEnergy()))  eleEoverPInv_.push_back(1e30);
    else  eleEoverPInv_.push_back((1.0 - iEle->eSuperClusterOverP())/iEle->ecalEnergy());

    reco::GsfElectron::PflowIsolationVariables pfIso = iEle->pfIsolationVariables();
    elePFChIso_         .push_back(pfIso.sumChargedHadronPt);
    elePFPhoIso_        .push_back(pfIso.sumPhotonEt);
    elePFNeuIso_        .push_back(pfIso.sumNeutralHadronEt);
    elePFPUIso_         .push_back(pfIso.sumPUPt);

    eleSigmaIEtaIEtaFull5x5_.push_back(iEle->full5x5_sigmaIetaIeta());
    eleSigmaIPhiIPhiFull5x5_.push_back(iEle->full5x5_sigmaIphiIphi());
    eleR9Full5x5_           .push_back(iEle->full5x5_r9());
    eleEcalDrivenSeed_      .push_back(iEle->ecalDrivenSeed());

    eleScale_stat_up_.push_back(iEle->userFloat("energyScaleStatUp"));
    eleScale_stat_dn_.push_back(iEle->userFloat("energyScaleStatDown"));
    eleScale_syst_up_.push_back(iEle->userFloat("energyScaleSystUp"));
    eleScale_syst_dn_.push_back(iEle->userFloat("energyScaleSystDown"));
    eleScale_gain_up_.push_back(iEle->userFloat("energyScaleGainUp"));
    eleScale_gain_dn_.push_back(iEle->userFloat("energyScaleGainDown"));
    eleResol_rho_up_ .push_back(iEle->userFloat("energySigmaRhoUp"));
    eleResol_rho_dn_ .push_back(iEle->userFloat("energySigmaRhoDown"));
    eleResol_phi_up_ .push_back(iEle->userFloat("energySigmaPhiUp"));
    eleResol_phi_dn_ .push_back(iEle->userFloat("energySigmaPhiDown"));

    reco::GsfTrackRef gsfTrackRef = iEle->gsfTrack();
    if (iEle->gsfTrack().isNonnull()) {
      eleGSFChi2_.push_back(gsfTrackRef->normalizedChi2());
      if (recVtxs->size() > 0)
        eleTrkdxy_.push_back(gsfTrackRef->dxy(recVtxs->front().position()));
      else
	eleTrkdxy_.push_back(-999);
    } else {
      eleGSFChi2_.push_back(999.);
      eleTrkdxy_.push_back(-999);
    }
    
    reco::TrackRef kfTrackRef = iEle->closestCtfTrackRef();
    if (kfTrackRef.isAvailable() && kfTrackRef.isNonnull()) {
      eleKFHits_.push_back(kfTrackRef->hitPattern().trackerLayersWithMeasurement());
      eleKFChi2_.push_back(kfTrackRef->normalizedChi2());
    } else {
      eleKFHits_.push_back(-1.);
      eleKFChi2_.push_back(999.);
    }

    if (development_) {

      Float_t ESp1 = 0;
      Float_t ESp2 = 0;
      vector<float> ESEta; 
      vector<float> ESPhi;
      vector<int>   ESZ;
      vector<int>   ESP;
      vector<int>   ESX;
      vector<int>   ESY;
      vector<int>   ESS;
      vector<float> ESE;
      for (CaloClusterPtrVector::const_iterator ips = iEle->superCluster()->preshowerClustersBegin(); ips != iEle->superCluster()->preshowerClustersEnd(); ++ips) {

	ESDetId esid = ESDetId((*ips)->seed());
	if (esid.plane() == 1) ESp1 += (*ips)->energy();
	if (esid.plane() == 2) ESp2 += (*ips)->energy();

	ESZ.push_back(esid.zside());
	ESP.push_back(esid.plane());
	ESX.push_back(esid.six());
	ESY.push_back(esid.siy());
	ESS.push_back(esid.strip());

	ESEta.push_back((*ips)->eta());
	ESPhi.push_back((*ips)->phi());
	ESE.push_back((*ips)->energy());
      }

      eleESEnP1Raw_.push_back(ESp1);
      eleESEnP2Raw_.push_back(ESp2);
      eleESEnEta_.push_back(ESEta);
      eleESEnPhi_.push_back(ESPhi);
      eleESEnZ_.push_back(ESZ);
      eleESEnP_.push_back(ESP);
      eleESEnX_.push_back(ESX);
      eleESEnY_.push_back(ESY);
      eleESEnS_.push_back(ESS);
      eleESEnE_.push_back(ESE);
    }

    // VID decisions 
    UShort_t tmpeleIDbit = 0;   
    bool isPassVeto   = iEle->electronID("cutBasedElectronID-Fall17-94X-V1-veto");
    if (isPassVeto)   setbit(tmpeleIDbit, 0);    
    bool isPassLoose  = iEle->electronID("cutBasedElectronID-Fall17-94X-V1-loose");
    if (isPassLoose)  setbit(tmpeleIDbit, 1);   
    bool isPassMedium = iEle->electronID("cutBasedElectronID-Fall17-94X-V1-medium");
    if (isPassMedium) setbit(tmpeleIDbit, 2);    
    bool isPassTight  = iEle->electronID("cutBasedElectronID-Fall17-94X-V1-tight");
    if (isPassTight)  setbit(tmpeleIDbit, 3);    
    bool isPassHEEP   = iEle->electronID("heepElectronID-HEEPV70");
    if (isPassHEEP)   setbit(tmpeleIDbit, 4);
    
    eleIDMVAIso_  .push_back(iEle->userFloat("ElectronMVAEstimatorRun2Fall17IsoV1Values"));
    eleIDMVANoIso_.push_back(iEle->userFloat("ElectronMVAEstimatorRun2Fall17NoIsoV1Values"));

    elePFClusEcalIso_.push_back(iEle->ecalPFClusterIso());
    elePFClusHcalIso_.push_back(iEle->hcalPFClusterIso());
    
    eleIDbit_.push_back(tmpeleIDbit);

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
