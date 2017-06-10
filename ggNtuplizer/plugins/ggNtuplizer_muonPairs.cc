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
Int_t            ndiMu_;
vector<int>      diMuIndex1_;
vector<int>      diMuIndex2_;
vector<int>      diMuVtxIsValid_;
vector<float>    diMuVtx_;
vector<float>    diMuVty_;
vector<float>    diMuVtz_;
vector<float>    diMuVtxProb_;
vector<float>    diMu_CosAlpha_; //the angle between the reconstructed momentum vector of the dimuon system and the vector from the PV to the dimuon vertex
vector<float>    diMu_Lxy_; // Transverse decay length (Lxy) between the dimuon vertex and the primary vertex
vector<float>    diMu_Rxy_;
vector<float>    diMu_eLxy_;
vector<float>    diMu_SLxy_; // Significance of Lxy ( Lxy divided by its uncertainty )
vector<float>    diMu_ctau_; // lifetime 


void ggNtuplizer::branchesMuonPairs(TTree* tree) {

  tree->Branch("ndiMu",       &ndiMu_);
  tree->Branch("diMuIndex1",   &diMuIndex1_);
  tree->Branch("diMuIndex2",   &diMuIndex2_);
  tree->Branch("diMuVtxIsValid",   &diMuVtxIsValid_);
  tree->Branch("diMuVtx",   &diMuVtx_);
  tree->Branch("diMuVty",   &diMuVty_);
  tree->Branch("diMuVtz",   &diMuVtz_);
  tree->Branch("diMuVtxProb",   &diMuVtxProb_);
  tree->Branch("diMu_CosAlpha",   &diMu_CosAlpha_);
  tree->Branch("diMu_Lxy",   &diMu_Lxy_);
  tree->Branch("diMu_Rxy",   &diMu_Rxy_);
  tree->Branch("diMu_eLxy",   &diMu_eLxy_);
  tree->Branch("diMu_SLxy",   &diMu_SLxy_);
  tree->Branch("diMu_ctau",   &diMu_ctau_);
}

void ggNtuplizer::fillMuonsPairs(const edm::Event& e, const edm::EventSetup& es, math::XYZPoint& pv, reco::Vertex vtx) {

  // cleanup from previous execution
  diMuIndex1_.clear();
  diMuIndex2_.clear();
  diMuVtxIsValid_.clear();
  diMuVtx_.clear();
  diMuVty_.clear();
  diMuVtz_.clear();
  diMuVtxProb_.clear();
  diMu_CosAlpha_.clear();
  diMu_Lxy_.clear();
  diMu_Rxy_.clear();
  diMu_eLxy_.clear();
  diMu_SLxy_.clear();
  diMu_ctau_.clear();
  ndiMu_ = 0;

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

      diMuIndex1_.push_back(tmpMu1);
      diMuIndex2_.push_back(tmpMu2);

      const reco::TransientTrack &tt2 = transientTrackBuilder->build(jMu->bestTrack());
      vector<reco::TransientTrack> t_tks = {tt1,tt2};

      KalmanVertexFitter fitter;
      TransientVertex tmpVertex = fitter.vertex(t_tks);

      if(tmpVertex.isValid()){
        diMuVtxIsValid_.push_back(1);
        diMuVtx_.push_back(tmpVertex.position().x());
        diMuVty_.push_back(tmpVertex.position().y());
        diMuVtz_.push_back(tmpVertex.position().z());
        diMuVtxProb_.push_back(ChiSquaredProbability(tmpVertex.totalChiSquared(), tmpVertex.degreesOfFreedom()));

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

        diMu_CosAlpha_.push_back(cosAlpha);
        diMu_Lxy_.push_back(disp.Dot(cmom) / cmom.Perp());
        diMu_Rxy_.push_back(distXY.value());
        diMu_eLxy_.push_back(distXY.error());
        diMu_SLxy_.push_back(distXY.value() / distXY.error());
        diMu_ctau_.push_back(ctauPV);
      }
      else{
        diMuVtxIsValid_.push_back(0);
        diMuVtx_.push_back(-999);
        diMuVty_.push_back(-999);
        diMuVtz_.push_back(-999);
        diMuVtxProb_.push_back(-999);
        diMu_CosAlpha_.push_back(-999);
        diMu_Lxy_.push_back(-999);
        diMu_Rxy_.push_back(-999);
        diMu_eLxy_.push_back(-999);
        diMu_SLxy_.push_back(-999);
        diMu_ctau_.push_back(-999);
      }

      ndiMu_++;
      tmpMu2++;
    }
    tmpMu1++;
    tmpMu2 = tmpMu1;
  }
}
