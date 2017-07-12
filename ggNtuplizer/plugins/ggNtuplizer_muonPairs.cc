#include "ggAnalysis/ggNtuplizer/interface/ggNtuplizer.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "Math/VectorUtil.h"
// user include files
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
///For kinematic fit:
#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/TrackReco/interface/TrackFwd.h>
#include <DataFormats/TrackReco/interface/Track.h>
#include <DataFormats/VertexReco/interface/VertexFwd.h>
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/UserData.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "MagneticField/Engine/interface/MagneticField.h"            
#include "CommonTools/Statistics/interface/ChiSquaredProbability.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"
#include "RecoVertex/VertexTools/interface/VertexDistanceXY.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "TMath.h"
#include "Math/VectorUtil.h"
#include "TVector3.h"
#include <boost/foreach.hpp>
#include <string>

using namespace std;

// (local) variables associated with tree branches
Int_t            ndiMu_;
vector<int>      diMuIndex1_;
vector<int>      diMuIndex2_;
vector<int>      diMuVtxIsValid_;
vector<float>    diMuVtx_;
vector<float>    diMuVty_;
vector<float>    diMuVtz_;
vector<float>    diMuChi2_;
vector<int>      diMuNDF_;
vector<float>    diMuVtxProb_;
vector<float>    diMu_CosAlpha_; //the angle between the reconstructed momentum vector of the dimuon system and the vector from the PV to the dimuon vertex
vector<float>    diMu_Lxy_; // Transverse decay length (Lxy) between the dimuon vertex and the primary vertex
vector<float>    diMu_Rxy_;
vector<float>    diMu_eLxy_;
vector<float>    diMu_SLxy_; // Significance of Lxy ( Lxy divided by its uncertainty )
vector<float>    diMu_ctau_; // lifetime 
vector<float>    diMu_ctauErr_;
// Kinematic fit (dimuon fit with mass constraint)
vector<int>      diMuVtxIsValid_KinFit_;
vector<float>    diMuMass_KinFit_;
vector<float>    diMuPx_KinFit_;
vector<float>    diMuPy_KinFit_;
vector<float>    diMuPz_KinFit_;
vector<float>    diMuVtx_KinFit_;
vector<float>    diMuVty_KinFit_;
vector<float>    diMuVtz_KinFit_;
vector<float>    diMuChi2_KinFit_;
vector<int>      diMuNDF_KinFit_;
vector<float>    diMuVtxProb_KinFit_;
vector<float>    diMuCosAlpha_KinFit_; //the angle between the reconstructed momentum vector of the dimuon system and the vector from the PV to the dimuon vertex
vector<float>    diMu_Lxy_KinFit_; // Transverse decay length (Lxy) between the dimuon vertex and the primary vertex
vector<float>    diMu_Rxy_KinFit_;
vector<float>    diMu_eLxy_KinFit_;
vector<float>    diMu_SLxy_KinFit_; // Significance of Lxy ( Lxy divided by its uncertainty )
vector<float>    diMu_ctau_KinFit_; // lifetime 
vector<float>    diMu_ctauErr_KinFit_;
vector<float>    diMu_Mmu1_KinFit_;
vector<float>    diMu_Mmu2_KinFit_;
vector<float>    diMu_mu1Px_KinFit_;
vector<float>    diMu_mu1Py_KinFit_;
vector<float>    diMu_mu1Pz_KinFit_;
vector<float>    diMu_mu2Px_KinFit_;
vector<float>    diMu_mu2Py_KinFit_;
vector<float>    diMu_mu2Pz_KinFit_;
vector<float>    diMu_mu1En_KinFit_;
vector<float>    diMu_mu2En_KinFit_;
vector<int>      diMu_mu1Charge_KinFit_;
vector<int>      diMu_mu2Charge_KinFit_;

void ggNtuplizer::branchesMuonPairs(TTree* tree) {

  tree->Branch("ndiMu",       &ndiMu_);
  tree->Branch("diMuIndex1",   &diMuIndex1_);
  tree->Branch("diMuIndex2",   &diMuIndex2_);
  tree->Branch("diMuVtxIsValid",   &diMuVtxIsValid_);
  tree->Branch("diMuVtx",   &diMuVtx_);
  tree->Branch("diMuVty",   &diMuVty_);
  tree->Branch("diMuVtz",   &diMuVtz_);
  tree->Branch("diMuChi2",   &diMuChi2_);
  tree->Branch("diMuNDF",   &diMuNDF_);
  tree->Branch("diMuVtxProb",   &diMuVtxProb_);
  tree->Branch("diMu_CosAlpha",   &diMu_CosAlpha_);
  tree->Branch("diMu_Lxy",   &diMu_Lxy_);
  tree->Branch("diMu_Rxy",   &diMu_Rxy_);
  tree->Branch("diMu_eLxy",   &diMu_eLxy_);
  tree->Branch("diMu_SLxy",   &diMu_SLxy_);
  tree->Branch("diMu_ctau",   &diMu_ctau_);
  tree->Branch("diMu_ctauErr",   &diMu_ctauErr_);
  // Kinematic fit
  tree->Branch("diMuVtxIsValid_KinFit",   &diMuVtxIsValid_KinFit_);
  tree->Branch("diMuMass_KinFit",   &diMuMass_KinFit_);
  tree->Branch("diMuPx_KinFit",   &diMuPx_KinFit_);
  tree->Branch("diMuPy_KinFit",   &diMuPy_KinFit_);
  tree->Branch("diMuPz_KinFit",   &diMuPz_KinFit_);
  tree->Branch("diMuVtx_KinFit",   &diMuVtx_KinFit_);
  tree->Branch("diMuVty_KinFit",   &diMuVty_KinFit_);
  tree->Branch("diMuVtz_KinFit",   &diMuVtz_KinFit_);
  tree->Branch("diMuChi2_KinFit",   &diMuChi2_KinFit_);
  tree->Branch("diMuNDF_KinFit",   &diMuNDF_KinFit_);
  tree->Branch("diMuVtxProb_KinFit",   &diMuVtxProb_KinFit_);
  tree->Branch("diMuCosAlpha_KinFit",   &diMuCosAlpha_KinFit_);
  tree->Branch("diMu_Lxy_KinFit",   &diMu_Lxy_KinFit_);
  tree->Branch("diMu_Rxy_KinFit",   &diMu_Rxy_KinFit_);
  tree->Branch("diMu_eLxy_KinFit",   &diMu_eLxy_KinFit_);
  tree->Branch("diMu_SLxy_KinFit",   &diMu_SLxy_KinFit_);
  tree->Branch("diMu_ctau_KinFit",   &diMu_ctau_KinFit_);
  tree->Branch("diMu_ctauErr_KinFit",   &diMu_ctauErr_KinFit_);
  tree->Branch("diMu_Mmu1_KinFit",   &diMu_Mmu1_KinFit_);
  tree->Branch("diMu_Mmu2_KinFit",   &diMu_Mmu2_KinFit_);
  tree->Branch("diMu_mu1Px_KinFit",   &diMu_mu1Px_KinFit_);
  tree->Branch("diMu_mu1Py_KinFit",   &diMu_mu1Py_KinFit_);
  tree->Branch("diMu_mu1Pz_KinFit",   &diMu_mu1Pz_KinFit_);
  tree->Branch("diMu_mu2Px_KinFit",   &diMu_mu2Px_KinFit_);
  tree->Branch("diMu_mu2Py_KinFit",   &diMu_mu2Py_KinFit_);
  tree->Branch("diMu_mu2Pz_KinFit",   &diMu_mu2Pz_KinFit_);
  tree->Branch("diMu_mu1En_KinFit",   &diMu_mu1En_KinFit_);
  tree->Branch("diMu_mu2En_KinFit",   &diMu_mu2En_KinFit_);
  tree->Branch("diMu_mu1Charge_KinFit",   &diMu_mu1Charge_KinFit_);
  tree->Branch("diMu_mu2Charge_KinFit",   &diMu_mu2Charge_KinFit_);

}

void ggNtuplizer::fillMuonsPairs(const edm::Event& e, const edm::EventSetup& es, math::XYZPoint& pv, reco::Vertex vtx) {

  // cleanup from previous execution
  diMuIndex1_.clear();
  diMuIndex2_.clear();
  diMuVtxIsValid_.clear();
  diMuVtx_.clear();
  diMuVty_.clear();
  diMuVtz_.clear();
  diMuChi2_.clear();
  diMuNDF_.clear();
  diMuVtxProb_.clear();
  diMu_CosAlpha_.clear();
  diMu_Lxy_.clear();
  diMu_Rxy_.clear();
  diMu_eLxy_.clear();
  diMu_SLxy_.clear();
  diMu_ctau_.clear();
  diMu_ctauErr_.clear();
  // Kinematic fit 
  diMuVtxIsValid_KinFit_.clear();
  diMuMass_KinFit_.clear();
  diMuPx_KinFit_.clear();
  diMuPy_KinFit_.clear();
  diMuPz_KinFit_.clear();
  diMuVtx_KinFit_.clear();
  diMuVty_KinFit_.clear();
  diMuVtz_KinFit_.clear();
  diMuChi2_KinFit_.clear();
  diMuNDF_KinFit_.clear();
  diMuVtxProb_KinFit_.clear();
  diMuCosAlpha_KinFit_.clear();
  diMu_Lxy_KinFit_.clear();
  diMu_Rxy_KinFit_.clear();
  diMu_eLxy_KinFit_.clear();
  diMu_SLxy_KinFit_.clear();
  diMu_ctau_KinFit_.clear();
  diMu_ctauErr_KinFit_.clear();
  diMu_Mmu1_KinFit_.clear();
  diMu_Mmu2_KinFit_.clear();
  diMu_mu1Px_KinFit_.clear();
  diMu_mu2Px_KinFit_.clear();
  diMu_mu1Py_KinFit_.clear();
  diMu_mu2Py_KinFit_.clear();
  diMu_mu1Pz_KinFit_.clear();
  diMu_mu2Pz_KinFit_.clear();
  diMu_mu1En_KinFit_.clear();
  diMu_mu2En_KinFit_.clear();
  diMu_mu1Charge_KinFit_.clear();
  diMu_mu2Charge_KinFit_.clear();
  
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
        diMuChi2_.push_back(tmpVertex.totalChiSquared());
        diMuNDF_.push_back(tmpVertex.degreesOfFreedom());
        diMuVtxProb_.push_back(ChiSquaredProbability(tmpVertex.totalChiSquared(), tmpVertex.degreesOfFreedom()));

        // Distance and significance
        TVector3 disp(tmpVertex.position().x() - vtx.x(), tmpVertex.position().y() - vtx.y(), 0);
        TLorentzVector dimuCand, mu1, mu2;
        mu1.SetPtEtaPhiM(iMu->pt(), iMu->eta(), iMu->phi(), 0.1057);
        mu2.SetPtEtaPhiM(jMu->pt(), jMu->eta(), jMu->phi(), 0.1057);
        dimuCand = mu1 + mu2;
        TVector3 cmom(dimuCand.Px(), dimuCand.Py(), 0 );
        float cosAlpha = disp.Dot( cmom ) / ( disp.Perp() * cmom.Perp() );
        float mass = dimuCand.M();
        AlgebraicVector3 vmom( dimuCand.Px(), dimuCand.Py(), 0 );
        VertexDistanceXY vdistXY;
        Measurement1D distXY = vdistXY.distance(tmpVertex,vtx);
        float ctauPV = (100./3.) * distXY.value() * cosAlpha * mass / cmom.Perp(); // unit : pico seconda
        GlobalError sve = (reco::Vertex(tmpVertex)).error();
        GlobalError pve = vtx.error();
        AlgebraicSymMatrix33 vXYe = sve.matrix() + pve.matrix();
        float ctauErrPV = sqrt( ROOT::Math::Similarity( vmom, vXYe ) ) * mass / cmom.Perp2();

        diMu_CosAlpha_.push_back(cosAlpha);
        diMu_Lxy_.push_back(disp.Dot(cmom) / cmom.Perp());
        diMu_Rxy_.push_back(distXY.value());
        diMu_eLxy_.push_back(distXY.error());
        diMu_SLxy_.push_back(distXY.value() / distXY.error());
        diMu_ctau_.push_back(ctauPV);
        diMu_ctauErr_.push_back(ctauErrPV);
      }
      else{
        diMuVtxIsValid_.push_back(0);
        diMuVtx_.push_back(-999.);
        diMuVty_.push_back(-999.);
        diMuVtz_.push_back(-999.);
        diMuChi2_.push_back(-999.);
        diMuNDF_.push_back(-999.);
        diMuVtxProb_.push_back(-999.);
        diMu_CosAlpha_.push_back(-999.);
        diMu_Lxy_.push_back(-999.);
        diMu_Rxy_.push_back(-999.);
        diMu_eLxy_.push_back(-999.);
        diMu_SLxy_.push_back(-999.);
        diMu_ctau_.push_back(-999.);
        diMu_ctauErr_.push_back(-999.);
      }

      // Kinematic Fit
      const reco::TrackRef innmu1 = iMu->innerTrack();
      const reco::TrackRef innmu2 = jMu->innerTrack();
      if(!innmu1.isNull() && !innmu2.isNull()){
        const reco::TransientTrack &tt1_KinFit = transientTrackBuilder->build(iMu->innerTrack());
        const reco::TransientTrack &tt2_KinFit = transientTrackBuilder->build(jMu->innerTrack());
        vector<reco::TransientTrack> t_tks_KinFit = {tt1_KinFit,tt2_KinFit};

        KinematicParticleFactoryFromTransientTrack pFactory;
     
        const ParticleMass muonMass(0.1056583);
        float muonSigma = muonMass*1E-6;

        vector<RefCountedKinematicParticle> allDaughters;
        allDaughters.push_back(pFactory.particle(t_tks_KinFit[0], muonMass, float(0), float(0), muonSigma));
        allDaughters.push_back(pFactory.particle(t_tks_KinFit[1], muonMass, float(0), float(0), muonSigma));

        KinematicConstrainedVertexFitter constVertexFitter;
        double massConstraint = 3.0969;
        MultiTrackKinematicConstraint *onia_mtc = new TwoTrackMassKinematicConstraint(massConstraint); // Jpsi mass in GeV
        RefCountedKinematicTree diMuTree = constVertexFitter.fit(allDaughters,onia_mtc);

        if(!diMuTree->isEmpty()){
          diMuTree->movePointerToTheTop();
          RefCountedKinematicParticle fitDimu = diMuTree->currentParticle();
          RefCountedKinematicVertex diMuCompVertex = diMuTree->currentDecayVertex();
          // Get diMuon candidate reffited
          if (fitDimu->currentState().isValid()){
            diMuVtxIsValid_KinFit_.push_back(1);
            diMuMass_KinFit_.push_back(fitDimu->currentState().mass());
            diMuChi2_KinFit_.push_back(diMuCompVertex->chiSquared());
            diMuNDF_KinFit_.push_back(diMuCompVertex->degreesOfFreedom());
            diMuVtxProb_KinFit_.push_back(ChiSquaredProbability((double)diMuCompVertex->chiSquared(), (double)(diMuCompVertex->degreesOfFreedom())));
            diMuVtx_KinFit_.push_back(diMuCompVertex->position().x());
            diMuVty_KinFit_.push_back(diMuCompVertex->position().y());
            diMuVtz_KinFit_.push_back(diMuCompVertex->position().z());

            TVector3 disp_KinFit(diMuCompVertex->position().x() - vtx.x(), diMuCompVertex->position().y() - vtx.y(), 0);
            TLorentzVector dimuCand_KinFit;
            float diMu_mass_fit = fitDimu->currentState().mass();
            float diMu_px_fit = fitDimu->currentState().kinematicParameters().momentum().x();
            float diMu_py_fit = fitDimu->currentState().kinematicParameters().momentum().y();
            float diMu_pz_fit = fitDimu->currentState().kinematicParameters().momentum().z();
            float diMu_en_fit = sqrt(diMu_mass_fit*diMu_mass_fit+diMu_px_fit*diMu_px_fit+
                                      diMu_py_fit*diMu_py_fit+diMu_pz_fit*diMu_pz_fit);
            dimuCand_KinFit.SetPxPyPzE(diMu_px_fit, diMu_py_fit, diMu_pz_fit, diMu_en_fit);
            TVector3 cmom_KinFit(dimuCand_KinFit.Px(), dimuCand_KinFit.Py(), 0);
            float cosAlpha_KinFit = disp_KinFit.Dot( cmom_KinFit ) / ( disp_KinFit.Perp() * cmom_KinFit.Perp() );
            AlgebraicVector3 vmom_KinFit(dimuCand_KinFit.Px(), dimuCand_KinFit.Py(), 0);
            VertexDistanceXY vdistXY_KinFit;
            Measurement1D distXY_KinFit = vdistXY_KinFit.distance(reco::Vertex(*diMuCompVertex),vtx);
            float ctauPV_KinFit = (100./3.) * distXY_KinFit.value() * cosAlpha_KinFit * diMu_mass_fit / cmom_KinFit.Perp(); // unit : pico seconda
            GlobalError v1e = (reco::Vertex(*diMuCompVertex)).error();
            GlobalError v2e = vtx.error();
            AlgebraicSymMatrix33 vXYe = v1e.matrix()+ v2e.matrix();
            double ctauErrPV_KinFit = sqrt(ROOT::Math::Similarity(vmom_KinFit,vXYe))*diMu_mass_fit/(cmom_KinFit.Perp2());

            diMuPx_KinFit_.push_back(diMu_px_fit);
            diMuPy_KinFit_.push_back(diMu_py_fit);
            diMuPz_KinFit_.push_back(diMu_pz_fit);
            diMuCosAlpha_KinFit_.push_back(cosAlpha_KinFit);
            diMu_Lxy_KinFit_.push_back(disp_KinFit.Dot(cmom_KinFit) / cmom_KinFit.Perp());
            diMu_Rxy_KinFit_.push_back(distXY_KinFit.value());
            diMu_eLxy_KinFit_.push_back(distXY_KinFit.error());
            diMu_SLxy_KinFit_.push_back(distXY_KinFit.value() / distXY_KinFit.error());
            diMu_ctau_KinFit_.push_back(ctauPV_KinFit);
            diMu_ctauErr_KinFit_.push_back(ctauErrPV_KinFit);

            bool child = diMuTree->movePointerToTheFirstChild();
            RefCountedKinematicParticle fitMu1 = diMuTree->currentParticle();
            if(child){
              float m1_ma_fit = fitMu1->currentState().mass();
              int   m1_ch_fit = fitMu1->currentState().particleCharge();
              float m1_px_fit = fitMu1->currentState().kinematicParameters().momentum().x();
              float m1_py_fit = fitMu1->currentState().kinematicParameters().momentum().y();
              float m1_pz_fit = fitMu1->currentState().kinematicParameters().momentum().z();
              float m1_en_fit = sqrt(m1_ma_fit*m1_ma_fit+m1_px_fit*m1_px_fit+m1_py_fit*m1_py_fit+m1_pz_fit*m1_pz_fit);
              diMu_Mmu1_KinFit_.push_back(m1_ma_fit);
              diMu_mu1Px_KinFit_.push_back(m1_px_fit);
              diMu_mu1Py_KinFit_.push_back(m1_py_fit);
              diMu_mu1Pz_KinFit_.push_back(m1_pz_fit);
              diMu_mu1En_KinFit_.push_back(m1_en_fit);
              diMu_mu1Charge_KinFit_.push_back(m1_ch_fit);
            }
            else{
              diMu_Mmu1_KinFit_.push_back(-999.);
              diMu_mu1Px_KinFit_.push_back(-999.);
              diMu_mu1Py_KinFit_.push_back(-999.);
              diMu_mu1Pz_KinFit_.push_back(-999.);
              diMu_mu1En_KinFit_.push_back(-999.);
              diMu_mu1Charge_KinFit_.push_back(-999.);
            }

            child = diMuTree->movePointerToTheNextChild();
            RefCountedKinematicParticle fitMu2 = diMuTree->currentParticle();
            if(child){
              float m2_ma_fit = fitMu2->currentState().mass();
              int   m2_ch_fit = fitMu2->currentState().particleCharge();
              float m2_px_fit = fitMu2->currentState().kinematicParameters().momentum().x();
              float m2_py_fit = fitMu2->currentState().kinematicParameters().momentum().y();
              float m2_pz_fit = fitMu2->currentState().kinematicParameters().momentum().z();
              float m2_en_fit = sqrt(m2_ma_fit*m2_ma_fit+m2_px_fit*m2_px_fit+m2_py_fit*m2_py_fit+m2_pz_fit*m2_pz_fit);
              diMu_Mmu2_KinFit_.push_back(m2_ma_fit);
              diMu_mu2Px_KinFit_.push_back(m2_px_fit);
              diMu_mu2Py_KinFit_.push_back(m2_py_fit);
              diMu_mu2Pz_KinFit_.push_back(m2_pz_fit);
              diMu_mu2En_KinFit_.push_back(m2_en_fit);
              diMu_mu2Charge_KinFit_.push_back(m2_ch_fit);
            }
            else{
              diMu_Mmu2_KinFit_.push_back(-999.);
              diMu_mu2Px_KinFit_.push_back(-999.);
              diMu_mu2Py_KinFit_.push_back(-999.);
              diMu_mu2Pz_KinFit_.push_back(-999.);
              diMu_mu2En_KinFit_.push_back(-999.);
              diMu_mu2Charge_KinFit_.push_back(-999.);
            }
          }
          else{
            diMuVtxIsValid_KinFit_.push_back(0);
            diMuVtx_KinFit_.push_back(-999.);
            diMuVty_KinFit_.push_back(-999.);
            diMuVtz_KinFit_.push_back(-999.);
            diMuChi2_KinFit_.push_back(-999.);
            diMuNDF_KinFit_.push_back(-999.);
            diMuPx_KinFit_.push_back(-999);
            diMuPy_KinFit_.push_back(-999);
            diMuPz_KinFit_.push_back(-999);
            diMuVtxProb_KinFit_.push_back(-999.);
            diMuCosAlpha_KinFit_.push_back(-999.);
            diMu_Lxy_KinFit_.push_back(-999.);
            diMu_Rxy_KinFit_.push_back(-999.);
            diMu_eLxy_KinFit_.push_back(-999.);
            diMu_SLxy_KinFit_.push_back(-999.);
            diMu_ctau_KinFit_.push_back(-999.);
            diMu_ctauErr_KinFit_.push_back(-999.);
          }
        }
        else{
          diMu_Mmu1_KinFit_.push_back(-999.);
          diMu_mu1Px_KinFit_.push_back(-999.);
          diMu_mu1Py_KinFit_.push_back(-999.);
          diMu_mu1Pz_KinFit_.push_back(-999.);
          diMu_mu1En_KinFit_.push_back(-999.);
          diMu_mu1Charge_KinFit_.push_back(-999.);
          diMu_Mmu2_KinFit_.push_back(-999.);
          diMu_mu2Px_KinFit_.push_back(-999.);
          diMu_mu2Py_KinFit_.push_back(-999.);
          diMu_mu2Pz_KinFit_.push_back(-999.);
          diMu_mu2En_KinFit_.push_back(-999.);
          diMu_mu2Charge_KinFit_.push_back(-999.);
          diMuVtxIsValid_KinFit_.push_back(0);
          diMuVtx_KinFit_.push_back(-999.);
          diMuVty_KinFit_.push_back(-999.);
          diMuVtz_KinFit_.push_back(-999.);
          diMuChi2_KinFit_.push_back(-999.);
          diMuNDF_KinFit_.push_back(-999.);
          diMuPx_KinFit_.push_back(-999);
          diMuPy_KinFit_.push_back(-999);
          diMuPz_KinFit_.push_back(-999);
          diMuVtxProb_KinFit_.push_back(-999.);
          diMuCosAlpha_KinFit_.push_back(-999.);
          diMu_Lxy_KinFit_.push_back(-999.);
          diMu_Rxy_KinFit_.push_back(-999.);
          diMu_eLxy_KinFit_.push_back(-999.);
          diMu_SLxy_KinFit_.push_back(-999.);
          diMu_ctau_KinFit_.push_back(-999.);
          diMu_ctauErr_KinFit_.push_back(-999.);
        }
      }
      else{
        diMu_Mmu1_KinFit_.push_back(-999.);
        diMu_mu1Px_KinFit_.push_back(-999.);
        diMu_mu1Py_KinFit_.push_back(-999.);
        diMu_mu1Pz_KinFit_.push_back(-999.);
        diMu_mu1En_KinFit_.push_back(-999.);
        diMu_mu1Charge_KinFit_.push_back(-999.);
        diMu_Mmu2_KinFit_.push_back(-999.);
        diMu_mu2Px_KinFit_.push_back(-999.);
        diMu_mu2Py_KinFit_.push_back(-999.);
        diMu_mu2Pz_KinFit_.push_back(-999.);
        diMu_mu2En_KinFit_.push_back(-999.);
        diMu_mu2Charge_KinFit_.push_back(-999.);
        diMuVtxIsValid_KinFit_.push_back(0);
        diMuVtx_KinFit_.push_back(-999.);
        diMuVty_KinFit_.push_back(-999.);
        diMuVtz_KinFit_.push_back(-999.);
        diMuChi2_KinFit_.push_back(-999.);
        diMuNDF_KinFit_.push_back(-999.);
        diMuPx_KinFit_.push_back(-999);
        diMuPy_KinFit_.push_back(-999);
        diMuPz_KinFit_.push_back(-999);
        diMuVtxProb_KinFit_.push_back(-999.);
        diMuCosAlpha_KinFit_.push_back(-999.);
        diMu_Lxy_KinFit_.push_back(-999.);
        diMu_Rxy_KinFit_.push_back(-999.);
        diMu_eLxy_KinFit_.push_back(-999.);
        diMu_SLxy_KinFit_.push_back(-999.);
        diMu_ctau_KinFit_.push_back(-999.);
        diMu_ctauErr_KinFit_.push_back(-999.);
      }
      ndiMu_++;
      tmpMu2++;
    }
    tmpMu1++;
    tmpMu2 = tmpMu1;
  }
}
