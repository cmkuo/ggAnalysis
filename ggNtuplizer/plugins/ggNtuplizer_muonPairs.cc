#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "ggAnalysis/ggNtuplizer/interface/ggNtuplizer.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "Math/VectorUtil.h"

using namespace std;

// (local) variables associated with tree branches
Int_t            nMuPair_;
vector<int>      muPairIndex1_;
vector<int>      muPairIndex2_;
vector<int>      VtxIsValid_;
vector<float>   VtxPosX_;
vector<float>   VtxPosY_;
vector<float>   VtxPosZ_;
vector<float>   VtxProb_;
vector<float>   VtxDistXY_;
vector<float>   CosAlpha_; //the angle between the reconstructed momentum vector of the dimuon system and the vector from the PV to the dimuon vertex
vector<float>   Lxy_; // Transverse decay length (Lxy) between the dimuon vertex and the primary vertex
vector<float>   Rxy_;
vector<float>   eLxy_;
vector<float>   SLxy_; // Significance of Lxy ( Lxy divided by its uncertainty )
vector<float>   ctau_; // lifetime 


void ggNtuplizer::branchesMuonPairs(TTree* tree) {

  tree->Branch("nMuPair",       &nMuPair_);
  tree->Branch("muPairIndex1",   &muPairIndex1_);
  tree->Branch("muPairIndex2",   &muPairIndex2_);
  tree->Branch("VtxIsValid",   &VtxIsValid_);
  tree->Branch("VtxPosX",   &VtxPosX_);
  tree->Branch("VtxPosY",   &VtxPosY_);
  tree->Branch("VtxPosZ",   &VtxPosZ_);
  tree->Branch("VtxProb",   &VtxProb_);
  tree->Branch("CosAlpha",   &CosAlpha_);
  tree->Branch("Lxy",   &Lxy_);
  tree->Branch("Rxy",   &Rxy_);
  tree->Branch("eLxy",   &eLxy_);
  tree->Branch("SLxy",   &SLxy_);
  tree->Branch("ctau",   &ctau_);
}

void ggNtuplizer::fillMuonsPairs(const edm::Event& e, const edm::EventSetup& es, math::XYZPoint& pv, reco::Vertex vtx) {

  // cleanup from previous execution
  muPairIndex1_.clear();
  muPairIndex2_.clear();
  VtxIsValid_.clear();
  VtxPosX_.clear();
  VtxPosY_.clear();
  VtxPosZ_.clear();
  VtxProb_.clear();
  CosAlpha_.clear();
  Lxy_.clear();
  Rxy_.clear();
  eLxy_.clear();
  SLxy_.clear();
  ctau_.clear();
  nMuPair_ = 0;

  edm::Handle<edm::View<pat::Muon> > muonHandle;
  e.getByToken(muonCollection_, muonHandle);

  edm::Handle<pat::PackedCandidateCollection> pfcands;
  e.getByToken(pckPFCandidateCollection_, pfcands);

  edm::ESHandle<TransientTrackBuilder> transientTrackBuilder;
  es.get<TransientTrackRecord>().get("TransientTrackBuilder", transientTrackBuilder);

  if (!muonHandle.isValid()) {
    edm::LogWarning("ggNtuplizer") << "no pat::Muons in event";
    return;
  }

  int tmpMu1 = 0;
  int tmpMu2 = 0;

  for (edm::View<pat::Muon>::const_iterator iMu = muonHandle->begin(); iMu != muonHandle->end(); ++iMu) {

    // Build transientTrack
    const reco::TransientTrack &tt1 = transientTrackBuilder->build(iMu->bestTrack());

    for (edm::View<pat::Muon>::const_iterator jMu = iMu; jMu != muonHandle->end(); ++jMu) {

      if (iMu == jMu) {
        tmpMu2++;
        continue;
      }
      if (iMu->pt() < 3 || jMu->pt() < 3) continue;
      if (! (iMu->isPFMuon() || iMu->isGlobalMuon() || iMu->isTrackerMuon())) continue;
      if (! (jMu->isPFMuon() || jMu->isGlobalMuon() || jMu->isTrackerMuon())) continue;

      muPairIndex1_.push_back(tmpMu1);
      muPairIndex2_.push_back(tmpMu2);

      const reco::TransientTrack &tt2 = transientTrackBuilder->build(jMu->bestTrack());
      vector<reco::TransientTrack> t_tks = {tt1,tt2};

      KalmanVertexFitter fitter;
      TransientVertex tmpVertex = fitter.vertex(t_tks);

      if(tmpVertex.isValid()){
        VtxIsValid_.push_back(1);
        VtxPosX_.push_back(tmpVertex.position().x());
        VtxPosY_.push_back(tmpVertex.position().y());
        VtxPosZ_.push_back(tmpVertex.position().z());
        VtxProb_.push_back(ChiSquaredProbability(tmpVertex.totalChiSquared(), tmpVertex.degreesOfFreedom()));

        // Distance and significance
        TVector3 disp(tmpVertex.position().x() - vtx.x(), tmpVertex.position().y() - vtx.y(), 0);
        TLorentzVector dimuCand, mu1, mu2;
        mu1.SetPtEtaPhiM(iMu->pt(), iMu->eta(), iMu->phi(), 0.1057);
        mu2.SetPtEtaPhiM(jMu->pt(), jMu->eta(), jMu->phi(), 0.1057);
        dimuCand = mu1 + mu2;
        TVector3 cmom(dimuCand.Px(), dimuCand.Py(), 0 );
        double cosAlpha = disp.Dot( cmom ) / ( disp.Perp() * cmom.Perp() );
        double mass = dimuCand.M();
        AlgebraicVector3 vmom( dimuCand.Px(), dimuCand.Py(), 0 );
        VertexDistanceXY vdistXY;
        Measurement1D distXY = vdistXY.distance(tmpVertex,vtx);
        double ctauPV = (100./3.) * distXY.value() * cosAlpha * mass / cmom.Perp(); // unit : pico seconda

        CosAlpha_.push_back(cosAlpha);
        Lxy_.push_back(disp.Dot(cmom) / cmom.Perp());
        Rxy_.push_back(distXY.value());
        eLxy_.push_back(distXY.error());
        SLxy_.push_back(distXY.value() / distXY.error());
        ctau_.push_back(ctauPV);
      }
      else{
        VtxIsValid_.push_back(0);
        VtxPosX_.push_back(-999);
        VtxPosY_.push_back(-999);
        VtxPosZ_.push_back(-999);
        VtxProb_.push_back(-999);
        CosAlpha_.push_back(-999);
        Lxy_.push_back(-999);
        Rxy_.push_back(-999);
        eLxy_.push_back(-999);
        SLxy_.push_back(-999);
        ctau_.push_back(-999);
      }

      nMuPair_++;
      tmpMu2++;
    }
    tmpMu1++;
    tmpMu2 = tmpMu1;
  }
}
