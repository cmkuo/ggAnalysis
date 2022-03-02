#include <TString.h>
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "ggAnalysis/ggNtuplizer/interface/ggNtuplizer.h"

using namespace std;

Int_t         nConv_;
vector<int>   convNTrks_;
vector<float> convVtxRadius_;
vector<float> convVtxX_;
vector<float> convVtxY_;
vector<float> convVtxZ_;
vector<float> convTrksPin0X_;
vector<float> convTrksPin0Y_;
vector<float> convTrksPin0Z_;
vector<float> convFitPairM_;
vector<float> convFitPairPX_;
vector<float> convFitPairPY_;
vector<float> convFitPairPZ_;
vector<float> convFitProb_;
vector<float> convD0_;
vector<float> convDz_;
vector<float> convL0_;
vector<float> convLz_;

void ggNtuplizer::branchesConversions(TTree* tree) {
  
  tree->Branch("nConv",                     &nConv_);
  tree->Branch("convNTrks",                 &convNTrks_);
  tree->Branch("convVtxRadius",             &convVtxRadius_);
  tree->Branch("convVtxX",                  &convVtxX_);
  tree->Branch("convVtxY",                  &convVtxY_);
  tree->Branch("convVtxZ",                  &convVtxZ_);
  tree->Branch("convTrksPin0X",             &convTrksPin0X_);
  tree->Branch("convTrksPin0Y",             &convTrksPin0Y_);
  tree->Branch("convTrksPin0Z",             &convTrksPin0Z_);
  tree->Branch("convFitPairPX",             &convFitPairPX_);
  tree->Branch("convFitPairPY",             &convFitPairPY_);
  tree->Branch("convFitPairPZ",             &convFitPairPZ_);
  tree->Branch("convFitPairM",              &convFitPairM_);
  tree->Branch("convFitProb",               &convFitProb_); 
  tree->Branch("convD0",                    &convD0_);             
  tree->Branch("convDz",                    &convDz_);             
  tree->Branch("convL0",                    &convL0_);             
  tree->Branch("convLz",                    &convLz_);             

}

void ggNtuplizer::fillConversions(const edm::Event& e, const edm::EventSetup& es) {
  
  // cleanup from previous execution
  convNTrks_               .clear();
  convVtxRadius_           .clear();
  convVtxX_                .clear();
  convVtxY_                .clear();
  convVtxZ_                .clear();
  convTrksPin0X_           .clear();
  convTrksPin0Y_           .clear();
  convTrksPin0Z_           .clear();
  convFitPairPX_           .clear();
  convFitPairPY_           .clear();
  convFitPairPZ_           .clear();
  convFitPairM_            .clear();
  convFitProb_             .clear();
  convD0_                  .clear();
  convDz_                  .clear();
  convL0_                  .clear();
  convLz_                  .clear();

  nConv_ = 0;

  edm::Handle<reco::ConversionCollection> ConversionHandle;
  e.getByToken(conversionsCollection_, ConversionHandle);

  edm::Handle<reco::ConversionCollection> SLConversionHandle;
  e.getByToken(conversionsCollectionSL_, SLConversionHandle);

  edm::Handle<reco::BeamSpot> beamSpotHandle;
  e.getByToken(beamSpot_, beamSpotHandle);
  const reco::BeamSpot &beamspot = *beamSpotHandle.product();

  if (!ConversionHandle.isValid()) {
    edm::LogWarning("ggNtuplizer") << "no pat::conversions in event";
    return;
  }

  if (!SLConversionHandle.isValid()) {
    edm::LogWarning("ggNtuplizer") << "no SL pat::conversions in event";
    return;
  }

  //for (edm::Ptr<reco::Conversion>::const_iterator iConv = ConversionHandle->begin(); iConv != ConversionHandle->end(); ++iConv) {
  for (reco::ConversionCollection::const_iterator iConv = ConversionHandle->begin(); iConv!= ConversionHandle->end(); ++iConv) {
  //for (const auto& iConv : *ConversionHandle) {

    convNTrks_                     .push_back(iConv->nTracks());
    convVtxRadius_                 .push_back(iConv->conversionVertex().position().perp2());
    convVtxX_                      .push_back(iConv->conversionVertex().x());
    convVtxY_                      .push_back(iConv->conversionVertex().y());
    convVtxZ_                      .push_back(iConv->conversionVertex().z());
    convTrksPin0X_                 .push_back(iConv->tracksPin()[0].x());
    convTrksPin0Y_                 .push_back(iConv->tracksPin()[0].y());
    convTrksPin0Z_                 .push_back(iConv->tracksPin()[0].z());
    convFitPairPX_                 .push_back(iConv->refittedPairMomentum().x());
    convFitPairPY_                 .push_back(iConv->refittedPairMomentum().y());
    convFitPairPZ_                 .push_back(iConv->refittedPairMomentum().z());
    convFitPairM_                  .push_back(iConv->refittedPair4Momentum().M());
    convFitProb_                   .push_back(TMath::Prob(iConv->conversionVertex().chi2(), iConv->conversionVertex().ndof()));
    convD0_                        .push_back(iConv->dxy(beamspot.position()));
    convDz_                        .push_back(iConv->dz(beamspot.position()));
    convL0_                        .push_back(iConv->lxy(beamspot.position()));
    convLz_                        .push_back(iConv->lz(beamspot.position()));

    nConv_++;
  }

  //for (edm::Ptr<reco::Conversion>::const_iterator iSLConv = SLConversionHandle->begin(); iSLConv != SLConversionHandle->end(); ++iSLConv) {
  for (reco::ConversionCollection::const_iterator iSLConv = SLConversionHandle->begin(); iSLConv!= SLConversionHandle->end(); ++iSLConv) {
    //for (const auto& iSLConv : *SLConversionHandle) {

    convNTrks_                     .push_back(iSLConv->nTracks());
    convVtxRadius_                 .push_back(iSLConv->conversionVertex().position().perp2());
    convVtxX_                      .push_back(iSLConv->conversionVertex().x());
    convVtxY_                      .push_back(iSLConv->conversionVertex().y());
    convVtxZ_                      .push_back(iSLConv->conversionVertex().z());
    convTrksPin0X_                 .push_back(iSLConv->tracksPin()[0].x());
    convTrksPin0Y_                 .push_back(iSLConv->tracksPin()[0].y());
    convTrksPin0Z_                 .push_back(iSLConv->tracksPin()[0].z());
    convFitPairPX_                 .push_back(iSLConv->refittedPairMomentum().x());
    convFitPairPY_                 .push_back(iSLConv->refittedPairMomentum().y());
    convFitPairPZ_                 .push_back(iSLConv->refittedPairMomentum().z());
    convFitPairM_                  .push_back(iSLConv->refittedPair4Momentum().M());
    convFitProb_                   .push_back(TMath::Prob(iSLConv->conversionVertex().chi2(), iSLConv->conversionVertex().ndof()));
    convD0_                        .push_back(iSLConv->dxy(beamspot.position()));
    convDz_                        .push_back(iSLConv->dz(beamspot.position()));
    convL0_                        .push_back(iSLConv->lxy(beamspot.position()));
    convLz_                        .push_back(iSLConv->lz(beamspot.position()));

    nConv_++;
  }

}

void ggNtuplizer::cleanupConversions() {

}
