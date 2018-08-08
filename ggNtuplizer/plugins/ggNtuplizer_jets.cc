#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "ggAnalysis/ggNtuplizer/interface/ggNtuplizer.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandGauss.h"

using namespace std;
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;

Int_t          nJet_;
vector<float>  jetPt_;
vector<float>  jetEn_;
vector<float>  jetEta_;
vector<float>  jetPhi_;
vector<float>  jetRawPt_;
vector<float>  jetRawEn_;
vector<float>  jetMt_;
vector<float>  jetArea_;
vector<float>  jetLeadTrackPt_;
vector<float>  jetLeadTrackEta_;
vector<float>  jetLeadTrackPhi_;
vector<int>    jetLepTrackPID_;
vector<float>  jetLepTrackPt_;
vector<float>  jetLepTrackEta_;
vector<float>  jetLepTrackPhi_;
vector<float>  jetCHF_;
vector<float>  jetNHF_;
vector<float>  jetCEF_;
vector<float>  jetNEF_;
vector<int>    jetNCH_;
vector<int>    jetNNP_;
vector<float>  jetMUF_;
vector<float>  jetHFHAE_;
vector<float>  jetHFEME_;
vector<int>    jetNConstituents_;
vector<float>  jetVtxPt_;
vector<float>  jetVtxMass_;
vector<float>  jetVtxNtrks_;
vector<float>  jetVtx3DVal_;
vector<float>  jetVtx3DSig_;
vector<float>  jetCSV2BJetTags_;
vector<float>  jetDeepCSVTags_b_;
vector<float>  jetDeepCSVTags_bb_;
vector<float>  jetDeepCSVTags_c_;
vector<float>  jetDeepCSVTags_udsg_;
vector<int>    jetPartonID_;
vector<int>    jetHadFlvr_;
vector<bool>   jetPFLooseId_;
vector<int>    jetID_; 
vector<float>  jetPUID_;
vector<int>    jetPUFullID_;
vector<float>  jetJECUnc_;
vector<float>  jetP4Smear_;
vector<float>  jetP4SmearUp_;
vector<float>  jetP4SmearDo_;
vector<ULong64_t> jetFiredTrgs_;
//gen-info for ak4
vector<float>  jetGenJetEn_;
vector<float>  jetGenJetPt_;
vector<float>  jetGenJetEta_;
vector<float>  jetGenJetPhi_;
vector<int>    jetGenPartonID_;
vector<float>  jetGenEn_;
vector<float>  jetGenPt_;
vector<float>  jetGenEta_;
vector<float>  jetGenPhi_;
vector<int>    jetGenPartonMomID_;

void ggNtuplizer::branchesJets(TTree* tree) {
  
  tree->Branch("nJet",                &nJet_);
  tree->Branch("jetPt",               &jetPt_);
  tree->Branch("jetEn",               &jetEn_);
  tree->Branch("jetEta",              &jetEta_);
  tree->Branch("jetPhi",              &jetPhi_);
  tree->Branch("jetRawPt",            &jetRawPt_);
  tree->Branch("jetRawEn",            &jetRawEn_);
  tree->Branch("jetMt",               &jetMt_);
  tree->Branch("jetArea",             &jetArea_);
  tree->Branch("jetLeadTrackPt",      &jetLeadTrackPt_);
  tree->Branch("jetLeadTrackEta",     &jetLeadTrackEta_);
  tree->Branch("jetLeadTrackPhi",     &jetLeadTrackPhi_);
  tree->Branch("jetLepTrackPID",      &jetLepTrackPID_);
  tree->Branch("jetLepTrackPt",       &jetLepTrackPt_);
  tree->Branch("jetLepTrackEta",      &jetLepTrackEta_);
  tree->Branch("jetLepTrackPhi",      &jetLepTrackPhi_);
  tree->Branch("jetCSV2BJetTags",     &jetCSV2BJetTags_);
  tree->Branch("jetDeepCSVTags_b",    &jetDeepCSVTags_b_);
  tree->Branch("jetDeepCSVTags_bb",   &jetDeepCSVTags_bb_);
  tree->Branch("jetDeepCSVTags_c",    &jetDeepCSVTags_c_);
  tree->Branch("jetDeepCSVTags_udsg", &jetDeepCSVTags_udsg_);
  if (doGenParticles_){
    tree->Branch("jetPartonID",       &jetPartonID_);
    tree->Branch("jetHadFlvr",        &jetHadFlvr_);
    tree->Branch("jetGenJetEn",       &jetGenJetEn_);
    tree->Branch("jetGenJetPt",       &jetGenJetPt_);
    tree->Branch("jetGenJetEta",      &jetGenJetEta_);
    tree->Branch("jetGenJetPhi",      &jetGenJetPhi_);
    tree->Branch("jetGenPartonID",    &jetGenPartonID_);
    tree->Branch("jetGenEn",          &jetGenEn_);
    tree->Branch("jetGenPt",          &jetGenPt_);
    tree->Branch("jetGenEta",         &jetGenEta_);
    tree->Branch("jetGenPhi",         &jetGenPhi_);
    tree->Branch("jetGenPartonMomID", &jetGenPartonMomID_);
    tree->Branch("jetP4Smear",        &jetP4Smear_);
    tree->Branch("jetP4SmearUp",      &jetP4SmearUp_);
    tree->Branch("jetP4SmearDo",      &jetP4SmearDo_);
  }  
  tree->Branch("jetPFLooseId", &jetPFLooseId_);
  tree->Branch("jetID",        &jetID_);
  tree->Branch("jetPUID",      &jetPUID_);
  tree->Branch("jetPUFullID",  &jetPUFullID_);
  tree->Branch("jetJECUnc",    &jetJECUnc_);
  tree->Branch("jetFiredTrgs", &jetFiredTrgs_);
  tree->Branch("jetCHF",       &jetCHF_);
  tree->Branch("jetNHF",       &jetNHF_);
  tree->Branch("jetCEF",       &jetCEF_);
  tree->Branch("jetNEF",       &jetNEF_);
  tree->Branch("jetNCH",       &jetNCH_);
  tree->Branch("jetNNP",       &jetNNP_);
  tree->Branch("jetMUF",       &jetMUF_);
  tree->Branch("jetVtxPt",     &jetVtxPt_);
  tree->Branch("jetVtxMass",   &jetVtxMass_);
  tree->Branch("jetVtxNtrks",  &jetVtxNtrks_);
  tree->Branch("jetVtx3DVal",  &jetVtx3DVal_);
  tree->Branch("jetVtx3DSig",  &jetVtx3DSig_);
  if (development_) {
    tree->Branch("jetHFHAE",         &jetHFHAE_);
    tree->Branch("jetHFEME",         &jetHFEME_);
    tree->Branch("jetNConstituents", &jetNConstituents_);
  }

}

void ggNtuplizer::fillJets(const edm::Event& e, const edm::EventSetup& es) {

  jetPt_                                  .clear();
  jetEn_                                  .clear();
  jetEta_                                 .clear();
  jetPhi_                                 .clear();
  jetRawPt_                               .clear();
  jetRawEn_                               .clear();
  jetMt_                                  .clear();
  jetArea_                                .clear();
  jetLeadTrackPt_                         .clear();
  jetLeadTrackEta_                        .clear();
  jetLeadTrackPhi_                        .clear();
  jetLepTrackPt_                          .clear();
  jetLepTrackPID_                         .clear();
  jetLepTrackEta_                         .clear();
  jetLepTrackPhi_                         .clear();
  jetCSV2BJetTags_                        .clear();
  jetDeepCSVTags_b_                       .clear();
  jetDeepCSVTags_bb_                      .clear();
  jetDeepCSVTags_c_                       .clear();
  jetDeepCSVTags_udsg_                    .clear();
  jetPartonID_                            .clear();
  jetHadFlvr_                             .clear();
  jetPFLooseId_                           .clear();
  jetID_                                  .clear();
  jetPUID_                                .clear();
  jetPUFullID_                            .clear();
  jetJECUnc_                              .clear();
  jetP4Smear_                             .clear();
  jetP4SmearUp_                           .clear();
  jetP4SmearDo_                           .clear();
  jetFiredTrgs_                           .clear();
  jetCHF_                                 .clear();
  jetNHF_                                 .clear();
  jetCEF_                                 .clear();
  jetNEF_                                 .clear();
  jetNCH_                                 .clear();
  jetNNP_                                 .clear();
  jetMUF_                                 .clear();
  jetVtxPt_                               .clear();
  jetVtxMass_                             .clear();
  jetVtxNtrks_                            .clear();
  jetVtx3DVal_                            .clear();
  jetVtx3DSig_                            .clear();
  if (development_) {
    jetHFHAE_                               .clear();
    jetHFEME_                               .clear();
    jetNConstituents_                       .clear();
  }
  jetGenJetEn_.clear();
  jetGenJetPt_.clear();
  jetGenJetEta_.clear();
  jetGenJetPhi_.clear();
  jetGenPartonID_.clear();
  jetGenEn_.clear();
  jetGenPt_.clear();
  jetGenEta_.clear();
  jetGenPhi_.clear();
  jetGenPartonMomID_.clear();

  nJet_ = 0;

  edm::Handle<edm::View<pat::Jet> > jetHandle;
  e.getByToken(jetsAK4Label_, jetHandle);

  if (!jetHandle.isValid()) {
    edm::LogWarning("ggNtuplizer") << "no pat::Jets (AK4) in event";
    return;
  }

  edm::Handle<vector<reco::GenParticle> > genParticlesHandle;
  if(doGenParticles_)e.getByToken(genParticlesCollection_, genParticlesHandle);
  
  //edm::Handle<double> rhoHandle;
  //e.getByToken(rhoLabel_, rhoHandle);
  //float rho = *(rhoHandle.product());
  
  edm::Handle<reco::VertexCollection> vtxHandle;
  e.getByToken(vtxLabel_, vtxHandle);
  if (!vtxHandle.isValid()) edm::LogWarning("ggNtuplizer") << "Primary vertices info not unavailable";
  
  // Accessing the JEC uncertainties 
  edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
  es.get<JetCorrectionsRecord>().get("AK4PFchs",JetCorParColl); 
  JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
  JetCorrectionUncertainty *jecUnc=0;
  jecUnc = new JetCorrectionUncertainty(JetCorPar);

  for (edm::View<pat::Jet>::const_iterator iJet = jetHandle->begin(); iJet != jetHandle->end(); ++iJet) {

    if (iJet->pt() < 20) continue;
    jetPt_.push_back(    iJet->pt());
    jetEn_.push_back(    iJet->energy());
    jetEta_.push_back(   iJet->eta());
    jetPhi_.push_back(   iJet->phi());
    jetRawPt_.push_back( (*iJet).correctedJet("Uncorrected").pt());
    jetRawEn_.push_back( (*iJet).correctedJet("Uncorrected").energy());
    jetMt_.push_back(    iJet->mt());
    jetArea_.push_back(  iJet->jetArea());
    jetCEF_.push_back(   iJet->chargedEmEnergyFraction());
    jetNEF_.push_back(   iJet->neutralEmEnergyFraction());
    jetCHF_.push_back(   iJet->chargedHadronEnergyFraction());
    jetNHF_.push_back(   iJet->neutralHadronEnergyFraction());
    jetNCH_.push_back(   iJet->chargedMultiplicity());
    jetNNP_.push_back(   iJet->neutralMultiplicity());
    jetMUF_.push_back(   iJet->muonEnergyFraction());
    if (development_) {
      jetHFHAE_.push_back( iJet->HFHadronEnergy());
      jetHFEME_.push_back( iJet->HFEMEnergy());
      jetNConstituents_.push_back(iJet->numberOfDaughters());
    }

    if (fabs(iJet->eta()) < 5.2) {
      jecUnc->setJetEta(iJet->eta());
      jecUnc->setJetPt(iJet->pt()); // here you must use the CORRECTED jet pt
      jetJECUnc_.push_back(jecUnc->getUncertainty(true));
    } else {
      jetJECUnc_.push_back(-1.);
    }

    jetFiredTrgs_.push_back(matchJetTriggerFilters(iJet->pt(), iJet->eta(), iJet->phi()));    

    //Searching for leading track and lepton
    float leadTrkPt  = -99;
    float leadTrkEta = -99;
    float leadTrkPhi = -99;
    int   lepTrkPID  = -99;
    float lepTrkPt   = -99;
    float lepTrkEta  = -99;
    float lepTrkPhi  = -99;

    for (unsigned id = 0; id < iJet->getJetConstituents().size(); id++) {

      const edm::Ptr<reco::Candidate> daughter = iJet->getJetConstituents().at(id);

      if (daughter.isNonnull() && daughter.isAvailable()) {
	if (daughter->charge() != 0 && daughter->pt() > leadTrkPt) {
	  leadTrkPt  = daughter->pt();
	  leadTrkEta = daughter->eta();
	  leadTrkPhi = daughter->phi();
	}

	if (abs(daughter->pdgId()) == 11 || abs(daughter->pdgId()) == 13) {
	  if (daughter->pt() > lepTrkPt) {
	    lepTrkPID = daughter->pdgId();
	    lepTrkPt  = daughter->pt();
	    lepTrkEta = daughter->eta();
	    lepTrkPhi = daughter->phi();
	  }
	}
      }
    }

    jetLeadTrackPt_ .push_back(leadTrkPt);
    jetLeadTrackEta_.push_back(leadTrkEta);
    jetLeadTrackPhi_.push_back(leadTrkPhi);
    jetLepTrackPID_ .push_back(lepTrkPID);
    jetLepTrackPt_  .push_back(lepTrkPt);
    jetLepTrackEta_ .push_back(lepTrkEta);
    jetLepTrackPhi_ .push_back(lepTrkPhi);    
    //jetVtxPt_       .push_back(sqrt(pow(iJet->userFloat("vtxPx"),2)+pow(iJet->userFloat("vtxPy"),2)));
    //jetVtxMass_     .push_back(iJet->userFloat("vtxMass"));
    //jetVtxNtrks_    .push_back(iJet->userFloat("vtxNtracks"));
    //jetVtx3DVal_    .push_back(iJet->userFloat("vtx3DVal"));
    //jetVtx3DSig_    .push_back(iJet->userFloat("vtx3DSig"));
    
    //b/c-tagging
    jetCSV2BJetTags_    .push_back(iJet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
    jetDeepCSVTags_b_   .push_back(iJet->bDiscriminator("pfDeepCSVJetTags:probb"));
    jetDeepCSVTags_bb_  .push_back(iJet->bDiscriminator("pfDeepCSVJetTags:probbb"));
    jetDeepCSVTags_c_   .push_back(iJet->bDiscriminator("pfDeepCSVJetTags:probc"));
    jetDeepCSVTags_udsg_.push_back(iJet->bDiscriminator("pfDeepCSVJetTags:probudsg"));
  
    //parton id
    jetPartonID_.push_back(iJet->partonFlavour());
    jetHadFlvr_.push_back(iJet->hadronFlavour());

    //jet PF Loose ID
    double NHF      = iJet->neutralHadronEnergyFraction();
    double NEMF     = iJet->neutralEmEnergyFraction();
    double NumConst = iJet->chargedMultiplicity()+iJet->neutralMultiplicity();
    double CHF      = iJet->chargedHadronEnergyFraction();
    double CHM      = iJet->chargedMultiplicity();
    double CEMF     = iJet->chargedEmEnergyFraction();
    double NNP      = iJet->neutralMultiplicity();

    bool looseJetID = false;    
    bool tightJetID = false;
    if (fabs(iJet->eta()) <= 2.7) {
      looseJetID = (NHF<0.99 && NEMF<0.99 && NumConst>1) && ((fabs(iJet->eta())<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || fabs(iJet->eta())>2.4);
      tightJetID = (NHF<0.90 && NEMF<0.90 && NumConst>1) && ((fabs(iJet->eta())<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || fabs(iJet->eta())>2.4);
    } else if (fabs(iJet->eta()) <= 3.0) {
      looseJetID = (NEMF>0.01 && NHF<0.98 && NNP>2);
      tightJetID = (NEMF>0.01 && NHF<0.98 && NNP>2);
    } else {
      looseJetID = (NEMF<0.90 && NNP>10); 
      tightJetID = (NEMF<0.90 && NNP>10);
    }
    jetPFLooseId_.push_back(looseJetID);
    Int_t jetIDdecision = 0;
    if (looseJetID) jetIDdecision += pow(2, 1);
    if (tightJetID) jetIDdecision += pow(2, 2);
    jetID_.push_back(jetIDdecision);    

    // PUJet ID from slimmedJets
    jetPUID_.push_back(iJet->userFloat("pileupJetId:fullDiscriminant"));
    jetPUFullID_.push_back(iJet->userInt("pileupJetId:fullId"));

    // gen jet and parton
    if (doGenParticles_ && genParticlesHandle.isValid()) {
      int jetGenPartonID    = -99;
      int jetGenPartonMomID = -99;
      float jetGenEn        = -999.;
      float jetGenPt        = -999.;
      float jetGenEta       = -999.;
      float jetGenPhi       = -999.;      
      if ((*iJet).genParton()) {
	jetGenPartonID = (*iJet).genParton()->pdgId();
	jetGenEn = (*iJet).genParton()->energy();
	jetGenPt = (*iJet).genParton()->pt();
	jetGenEta = (*iJet).genParton()->eta();
	jetGenPhi = (*iJet).genParton()->phi();
	if ((*iJet).genParton()->mother()) {
	  jetGenPartonMomID = (*iJet).genParton()->mother()->pdgId();
	}
      }
      
      jetGenPartonID_.push_back(jetGenPartonID);
      jetGenPartonMomID_.push_back(jetGenPartonMomID);
      jetGenEn_ .push_back(jetGenEn);
      jetGenPt_ .push_back(jetGenPt);
      jetGenEta_ .push_back(jetGenEta);
      jetGenPhi_ .push_back(jetGenPhi);
      
      float jetGenJetEn  = -999.;
      float jetGenJetPt  = -999.;
      float jetGenJetEta = -999.;
      float jetGenJetPhi = -999.;
      if ((*iJet).genJet()) {
	jetGenJetEn = (*iJet).genJet()->energy();
	jetGenJetPt = (*iJet).genJet()->pt();
	jetGenJetEta = (*iJet).genJet()->eta();
	jetGenJetPhi = (*iJet).genJet()->phi();
      }
      jetGenJetEn_.push_back(jetGenJetEn);
      jetGenJetPt_.push_back(jetGenJetPt);
      jetGenJetEta_.push_back(jetGenJetEta);
      jetGenJetPhi_.push_back(jetGenJetPhi);
      
      // access jet resolution       
      /*
      JME::JetParameters parameters;
      parameters.setJetPt(iJet->pt()).setJetEta(iJet->eta()).setRho(rho);
      float jetResolution = jetResolution_.getResolution(parameters);

      edm::Service<edm::RandomNumberGenerator> rng;
      if (!rng.isAvailable()) edm::LogError("JET : random number generator is missing !");
      CLHEP::HepRandomEngine & engine = rng->getEngine( e.streamID() );
      float rnd = CLHEP::RandGauss::shoot(&engine, 0., jetResolution);

      float jetResolutionSF   = jetResolutionSF_.getScaleFactor(parameters);
      float jetResolutionSFUp = jetResolutionSF_.getScaleFactor(parameters, Variation::UP);
      float jetResolutionSFDo = jetResolutionSF_.getScaleFactor(parameters, Variation::DOWN);

      float jetP4Smear   = -1.;
      float jetP4SmearUp = -1.;
      float jetP4SmearDo = -1.;
      if (jetGenJetPt > 0 && deltaR(iJet->eta(), iJet->phi(), jetGenJetEta, jetGenJetPhi) < 0.2 && fabs(iJet->pt()-jetGenJetPt) < 3*jetResolution*iJet->pt()) {
	jetP4Smear   = 1. + (jetResolutionSF   - 1.)*(iJet->pt() - jetGenJetPt)/iJet->pt();
	jetP4SmearUp = 1. + (jetResolutionSFUp - 1.)*(iJet->pt() - jetGenJetPt)/iJet->pt();
	jetP4SmearDo = 1. + (jetResolutionSFDo - 1.)*(iJet->pt() - jetGenJetPt)/iJet->pt();
      } else {
	jetP4Smear   = 1. + rnd*sqrt(max(pow(jetResolutionSF,   2)-1, 0.));
        jetP4SmearUp = 1. + rnd*sqrt(max(pow(jetResolutionSFUp, 2)-1, 0.));
	jetP4SmearDo = 1. + rnd*sqrt(max(pow(jetResolutionSFDo, 2)-1, 0.));
      }
      jetP4Smear_  .push_back(jetP4Smear);
      jetP4SmearUp_.push_back(jetP4SmearUp);
      jetP4SmearDo_.push_back(jetP4SmearDo);
      */
    }
    
    nJet_++;
  }
  
  delete jecUnc;
}
