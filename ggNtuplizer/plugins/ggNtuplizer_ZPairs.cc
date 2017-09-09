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

Int_t            nDiLep_;
vector<int>      diLepFlavor_;
vector<float>    diLepKFVtx_;
vector<float>    diLepKFVty_;
vector<float>    diLepKFVtz_;
vector<int>      diLepKFVtxIsValid_;
vector<int>      diLep1Index_;
vector<float>    diLep1KFPx_;
vector<float>    diLep1KFPy_;
vector<float>    diLep1KFPz_;
vector<int>      diLep2Index_;
vector<float>    diLep2KFPx_;
vector<float>    diLep2KFPy_;
vector<float>    diLep2KFPz_;

void ggNtuplizer::branchesZPairs(TTree* tree) {

  tree->Branch("nDiLep",            &nDiLep_);
  tree->Branch("diLepFlavor",       &diLepFlavor_);
  tree->Branch("diLepKFVtx",        &diLepKFVtx_);
  tree->Branch("diLepKFVty",        &diLepKFVty_);
  tree->Branch("diLepKFVtz",        &diLepKFVtz_);
  tree->Branch("diLepKFVtxIsValid", &diLepKFVtxIsValid_);
  tree->Branch("diLep1Index",       &diLep1Index_);
  tree->Branch("diLep1KFPx",        &diLep1KFPx_);
  tree->Branch("diLep1KFPy",        &diLep1KFPy_);
  tree->Branch("diLep1KFPz",        &diLep1KFPz_);
  tree->Branch("diLep2Index",       &diLep2Index_);
  tree->Branch("diLep2KFPx",        &diLep2KFPx_);
  tree->Branch("diLep2KFPy",        &diLep2KFPy_);
  tree->Branch("diLep2KFPz",        &diLep2KFPz_);

}

void ggNtuplizer::fillZPairs(const edm::Event& e, const edm::EventSetup& es, math::XYZPoint& pv, reco::Vertex vtx) {

  diLepFlavor_.clear();
  diLepKFVtx_.clear();
  diLepKFVty_.clear();
  diLepKFVtz_.clear();
  diLepKFVtxIsValid_.clear();
  diLep1Index_.clear();
  diLep1KFPx_.clear();
  diLep1KFPy_.clear();
  diLep1KFPz_.clear();
  diLep2Index_.clear();
  diLep2KFPx_.clear();
  diLep2KFPy_.clear();
  diLep2KFPz_.clear();
  
  nDiLep_ = 0;

  edm::Handle<edm::View<pat::Electron> > electronHandle;
  e.getByToken(electronCollection_, electronHandle);

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

  double massConstraint = 91.1876;
  int    tmpLep1 = 0;
  int    tmpLep2 = 0;

  for (edm::View<pat::Electron>::const_iterator iEle = electronHandle->begin(); iEle != electronHandle->end(); ++iEle) {

    tmpLep2 = tmpLep1;
    for (edm::View<pat::Electron>::const_iterator jEle = iEle; jEle != electronHandle->end(); ++jEle) {

      if (iEle == jEle) {
	tmpLep2++;
	continue;
      }

      diLepFlavor_.push_back(11);
      diLep1Index_.push_back(tmpLep1);
      diLep2Index_.push_back(tmpLep2);

      // Kinematic Fit
      const reco::TransientTrack &eltrk1_KinFit = transientTrackBuilder->build(iEle->bestTrack());
      const reco::TransientTrack &eltrk2_KinFit = transientTrackBuilder->build(jEle->bestTrack());
      vector<reco::TransientTrack> eltks_KinFit = {eltrk1_KinFit, eltrk2_KinFit};
      
      KinematicParticleFactoryFromTransientTrack elFactory;
      
      const ParticleMass electronMass(0.0005109989461);
      float electronSigma = electronMass*1E-6;
      
      vector<RefCountedKinematicParticle> allElDaughters;
      allElDaughters.push_back(elFactory.particle(eltks_KinFit[0], electronMass, float(0), float(0), electronSigma));
      allElDaughters.push_back(elFactory.particle(eltks_KinFit[1], electronMass, float(0), float(0), electronSigma));
      
      KinematicConstrainedVertexFitter constElVertexFitter;
      
      MultiTrackKinematicConstraint *Zee = new TwoTrackMassKinematicConstraint(massConstraint);
      RefCountedKinematicTree diElTree = constElVertexFitter.fit(allElDaughters, Zee);
      
      int   elKFVtxIsValid_ = 0;
      float elKFVtx_        = -999.;
      float elKFVty_        = -999.;
      float elKFVtz_        = -999.;
      float el1Px_          = -999.;
      float el1Py_          = -999.;
      float el1Pz_          = -999.;
      float el2Px_          = -999.;
      float el2Py_          = -999.;
      float el2Pz_          = -999.;

      if (!diElTree->isEmpty()){

	diElTree->movePointerToTheTop();
	RefCountedKinematicParticle fitDiel = diElTree->currentParticle();
	RefCountedKinematicVertex diElCompVertex = diElTree->currentDecayVertex();

	// Get refitted diMuon candidate
	if (fitDiel->currentState().isValid()){
	  
	  elKFVtxIsValid_ = 1;
	  elKFVtx_        = diElCompVertex->position().x();
	  elKFVty_        = diElCompVertex->position().y();
	  elKFVtz_        = diElCompVertex->position().z();	 

	  bool childEl = diElTree->movePointerToTheFirstChild();
	  RefCountedKinematicParticle fitEl1 = diElTree->currentParticle();
	  if (childEl) {
	    el1Px_ = fitEl1->currentState().kinematicParameters().momentum().x();
	    el1Py_ = fitEl1->currentState().kinematicParameters().momentum().y();
	    el1Pz_ = fitEl1->currentState().kinematicParameters().momentum().z();
	  } 
	  
	  childEl = diElTree->movePointerToTheNextChild();
	  RefCountedKinematicParticle fitEl2 = diElTree->currentParticle();
	  if (childEl) {
	    el2Px_ = fitEl2->currentState().kinematicParameters().momentum().x();
	    el2Py_ = fitEl2->currentState().kinematicParameters().momentum().y();
	    el2Pz_ = fitEl2->currentState().kinematicParameters().momentum().z();
	  }
	}
      }
      //}
      
      diLepKFVtxIsValid_.push_back(elKFVtxIsValid_);
      diLepKFVtx_.push_back(elKFVtx_);
      diLepKFVty_.push_back(elKFVty_);
      diLepKFVtz_.push_back(elKFVtz_);
      diLep1KFPx_.push_back(el1Px_);
      diLep1KFPy_.push_back(el1Py_);
      diLep1KFPz_.push_back(el1Pz_);
      diLep2KFPx_.push_back(el2Px_);
      diLep2KFPy_.push_back(el2Py_);
      diLep2KFPz_.push_back(el2Pz_);      

      nDiLep_++;
      tmpLep2++;
    } // jEle

    tmpLep1++;
  } // iEle

  tmpLep1 = 0;
  tmpLep2 = 0;

  for (edm::View<pat::Muon>::const_iterator iMu = muonHandle->begin(); iMu != muonHandle->end(); ++iMu) {

    if (iMu->pt() < 3) continue;
    if (! (iMu->isPFMuon() || iMu->isGlobalMuon() || iMu->isTrackerMuon())) continue;
    //cout<<"lep 1 : "<<tmpLep1<<endl;

    tmpLep2 = tmpLep1;
    for (edm::View<pat::Muon>::const_iterator jMu = iMu; jMu != muonHandle->end(); ++jMu) {

      if (iMu == jMu) {
	tmpLep2++;
	continue;
      }
      if (jMu->pt() < 3) continue;
      if (! (jMu->isPFMuon() || jMu->isGlobalMuon() || jMu->isTrackerMuon())) continue;

      //cout<<"lep pair : "<<tmpLep1<<" "<<tmpLep2<<endl;
      //cout<<"mu1      : "<<iMu->pt()<<" "<<iMu->eta()<<endl;
      //cout<<"mu2      : "<<jMu->pt()<<" "<<jMu->eta()<<endl;
      diLepFlavor_.push_back(13);
      diLep1Index_.push_back(tmpLep1);
      diLep2Index_.push_back(tmpLep2);

      // Kinematic Fit
      //const reco::TrackRef innmu1 = iMu->innerTrack();
      //const reco::TrackRef innmu2 = jMu->innerTrack();
      //if (!innmu1.isNull() && !innmu2.isNull()){
      const reco::TransientTrack &mutrk1_KinFit = transientTrackBuilder->build(iMu->bestTrack());
      const reco::TransientTrack &mutrk2_KinFit = transientTrackBuilder->build(jMu->bestTrack());
      vector<reco::TransientTrack> mutks_KinFit = {mutrk1_KinFit, mutrk2_KinFit};
      
      KinematicParticleFactoryFromTransientTrack muFactory;
      
      const ParticleMass muonMass(0.1056583);
      float muonSigma = muonMass*1E-6;
      
      vector<RefCountedKinematicParticle> allMuDaughters;
      allMuDaughters.push_back(muFactory.particle(mutks_KinFit[0], muonMass, float(0), float(0), muonSigma));
      allMuDaughters.push_back(muFactory.particle(mutks_KinFit[1], muonMass, float(0), float(0), muonSigma));
      
      KinematicConstrainedVertexFitter constMuVertexFitter;
      
      MultiTrackKinematicConstraint *Zmm = new TwoTrackMassKinematicConstraint(massConstraint);
      RefCountedKinematicTree diMuTree = constMuVertexFitter.fit(allMuDaughters, Zmm);
      
      int   muKFVtxIsValid_ = 0;
      float muKFVtx_        = -999.;
      float muKFVty_        = -999.;
      float muKFVtz_        = -999.;
      float mu1Px_          = -999.;
      float mu1Py_          = -999.;
      float mu1Pz_          = -999.;
      float mu2Px_          = -999.;
      float mu2Py_          = -999.;
      float mu2Pz_          = -999.;

      if (!diMuTree->isEmpty()){

	diMuTree->movePointerToTheTop();
	RefCountedKinematicParticle fitDimu = diMuTree->currentParticle();
	RefCountedKinematicVertex diMuCompVertex = diMuTree->currentDecayVertex();

	// Get refitted diMuon candidate
	if (fitDimu->currentState().isValid()){
	  
	  muKFVtxIsValid_ = 1;
	  muKFVtx_        = diMuCompVertex->position().x();
	  muKFVty_        = diMuCompVertex->position().y();
	  muKFVtz_        = diMuCompVertex->position().z();	 

	  bool childMu = diMuTree->movePointerToTheFirstChild();
	  RefCountedKinematicParticle fitMu1 = diMuTree->currentParticle();
	  if (childMu) {
	    mu1Px_ = fitMu1->currentState().kinematicParameters().momentum().x();
	    mu1Py_ = fitMu1->currentState().kinematicParameters().momentum().y();
	    mu1Pz_ = fitMu1->currentState().kinematicParameters().momentum().z();
	  } 
	  
	  childMu = diMuTree->movePointerToTheNextChild();
	  RefCountedKinematicParticle fitMu2 = diMuTree->currentParticle();
	  if (childMu) {
	    mu2Px_ = fitMu2->currentState().kinematicParameters().momentum().x();
	    mu2Py_ = fitMu2->currentState().kinematicParameters().momentum().y();
	    mu2Pz_ = fitMu2->currentState().kinematicParameters().momentum().z();
	  }
	}
      }
      //}
      
      diLepKFVtxIsValid_.push_back(muKFVtxIsValid_);
      diLepKFVtx_.push_back(muKFVtx_);
      diLepKFVty_.push_back(muKFVty_);
      diLepKFVtz_.push_back(muKFVtz_);
      diLep1KFPx_.push_back(mu1Px_);
      diLep1KFPy_.push_back(mu1Py_);
      diLep1KFPz_.push_back(mu1Pz_);
      diLep2KFPx_.push_back(mu2Px_);
      diLep2KFPy_.push_back(mu2Py_);
      diLep2KFPz_.push_back(mu2Pz_);
      
      //cout<<"vertex : "<<muKFVtxIsValid_<<" "<<sqrt(mu1Px_*mu1Px_+mu1Py_*mu1Py_)<<" "<<sqrt(mu2Px_*mu2Px_+mu2Py_*mu2Py_)<<endl;

      nDiLep_++;
      tmpLep2++;
    } // jMu

    tmpLep1++;
  } // iMu

  //cout<<"# of dilepton : "<<nDiLep_<<endl;

}
