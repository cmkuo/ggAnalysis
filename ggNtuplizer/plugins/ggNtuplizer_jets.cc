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

// ak4 jets
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
vector<float>  jetCSV2BJetTags_; // recommended
vector<float>  jetJetProbabilityBJetTags_;
vector<float>  jetpfCombinedMVAV2BJetTags_;
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
vector<UInt_t> jetFiredTrgs_;
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
//fat-jets (ak8)
Int_t         nAK8Jet_;
vector<float> AK8JetPt_;
vector<float> AK8JetEn_;
vector<float> AK8JetRawPt_;
vector<float> AK8JetRawEn_;
vector<float> AK8JetEta_;
vector<float> AK8JetPhi_;
vector<float> AK8JetMass_;
vector<float> AK8Jet_tau1_;
vector<float> AK8Jet_tau2_;
vector<float> AK8Jet_tau3_;
vector<float> AK8JetCHF_;
vector<float> AK8JetNHF_;
vector<float> AK8JetCEF_;
vector<float> AK8JetNEF_;
vector<int>   AK8JetNCH_;
vector<int>   AK8JetNNP_;
vector<float> AK8JetMUF_;
vector<int>   AK8Jetnconstituents_;
vector<bool>  AK8JetPFLooseId_;
vector<bool>  AK8JetPFTightLepVetoId_;
vector<float> AK8JetSoftDropMass_;
vector<float> AK8JetSoftDropMassCorr_;
vector<float> AK8JetPrunedMass_;
vector<float> AK8JetPrunedMassCorr_;
vector<float> AK8JetpfBoostedDSVBTag_;
vector<float> AK8JetDSVnewV4_;
vector<float> AK8JetCSV_;
vector<float> AK8JetJECUnc_;
vector<float> AK8JetL2L3corr_;

//gen-info for ak8
vector<int>   AK8JetPartonID_;
vector<int>   AK8JetHadFlvr_;
vector<int>   AK8JetGenJetIndex_;
vector<float> AK8JetGenJetEn_;
vector<float> AK8JetGenJetPt_;
vector<float> AK8JetGenJetEta_;
vector<float> AK8JetGenJetPhi_;
vector<int>   AK8JetGenPartonID_;
vector<float> AK8JetGenEn_;
vector<float> AK8JetGenPt_;
vector<float> AK8JetGenEta_;
vector<float> AK8JetGenPhi_;
vector<int>   AK8JetGenPartonMomID_;
vector<float> AK8JetP4Smear_;
vector<float> AK8JetP4SmearUp_;
vector<float> AK8JetP4SmearDo_;
//soft drop subjets
vector<int>             nAK8SDSJ_ ;
vector< vector<float> > AK8SDSJPt_ ;
vector< vector<float> > AK8SDSJEta_ ;
vector< vector<float> > AK8SDSJMass_ ;
vector< vector<float> > AK8SDSJPhi_ ;
vector< vector<float> > AK8SDSJE_ ;
vector< vector<int > >  AK8SDSJCharge_ ;
vector< vector<int > >  AK8SDSJFlavour_;
vector< vector<float> > AK8SDSJCSV_ ;
//puppi
vector<float> AK8puppiPt_;
vector<float> AK8puppiMass_;
vector<float> AK8puppiEta_;
vector<float> AK8puppiPhi_;
vector<float> AK8puppiTau1_;
vector<float> AK8puppiTau2_;
vector<float> AK8puppiTau3_;
vector<float> AK8puppiSDL2L3corr_;
vector<float> AK8puppiSDMass_;
vector<float> AK8puppiSDMassL2L3Corr_;

//puppi + soft drop subjets
vector<int>             nAK8puppiSDSJ_ ;
vector< vector<float> > AK8puppiSDSJPt_ ;
vector< vector<float> > AK8puppiSDSJEta_ ;
vector< vector<float> > AK8puppiSDSJMass_ ;
vector< vector<float> > AK8puppiSDSJPhi_ ;
vector< vector<float> > AK8puppiSDSJE_ ;
vector< vector<int > >  AK8puppiSDSJCharge_ ;
vector< vector<int > >  AK8puppiSDSJFlavour_;
vector< vector<float> > AK8puppiSDSJCSV_ ;

void ggNtuplizer::branchesJets(TTree* tree) {
  
  tree->Branch("nJet",            &nJet_);
  tree->Branch("jetPt",           &jetPt_);
  tree->Branch("jetEn",           &jetEn_);
  tree->Branch("jetEta",          &jetEta_);
  tree->Branch("jetPhi",          &jetPhi_);
  tree->Branch("jetRawPt",        &jetRawPt_);
  tree->Branch("jetRawEn",        &jetRawEn_);
  tree->Branch("jetMt",           &jetMt_);
  tree->Branch("jetArea",         &jetArea_);
  tree->Branch("jetLeadTrackPt",  &jetLeadTrackPt_);
  tree->Branch("jetLeadTrackEta", &jetLeadTrackEta_);
  tree->Branch("jetLeadTrackPhi", &jetLeadTrackPhi_);
  tree->Branch("jetLepTrackPID",  &jetLepTrackPID_);
  tree->Branch("jetLepTrackPt",   &jetLepTrackPt_);
  tree->Branch("jetLepTrackEta",  &jetLepTrackEta_);
  tree->Branch("jetLepTrackPhi",  &jetLepTrackPhi_);
  tree->Branch("jetCSV2BJetTags", &jetCSV2BJetTags_);
  tree->Branch("jetJetProbabilityBJetTags", &jetJetProbabilityBJetTags_);
  tree->Branch("jetpfCombinedMVAV2BJetTags", &jetpfCombinedMVAV2BJetTags_);
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

  // SubJet
  if (dumpSubJets_) {
    tree->Branch("nAK8Jet",                  &nAK8Jet_);
    tree->Branch("AK8JetPt",                 &AK8JetPt_);
    tree->Branch("AK8JetEn",                 &AK8JetEn_);
    tree->Branch("AK8JetRawPt",              &AK8JetRawPt_);
    tree->Branch("AK8JetRawEn",              &AK8JetRawEn_);
    tree->Branch("AK8JetEta",                &AK8JetEta_);
    tree->Branch("AK8JetPhi",                &AK8JetPhi_);
    tree->Branch("AK8JetMass",               &AK8JetMass_);
    tree->Branch("AK8Jet_tau1",              &AK8Jet_tau1_);
    tree->Branch("AK8Jet_tau2",              &AK8Jet_tau2_);
    tree->Branch("AK8Jet_tau3",              &AK8Jet_tau3_);
    tree->Branch("AK8JetCHF",                &AK8JetCHF_);
    tree->Branch("AK8JetNHF",                &AK8JetNHF_);
    tree->Branch("AK8JetCEF",                &AK8JetCEF_);
    tree->Branch("AK8JetNEF",                &AK8JetNEF_);
    tree->Branch("AK8JetNCH",                &AK8JetNCH_);
    tree->Branch("AK8JetNNP",                &AK8JetNNP_);
    tree->Branch("AK8JetMUF",                &AK8JetMUF_);
    tree->Branch("AK8Jetnconstituents",      &AK8Jetnconstituents_);
    tree->Branch("AK8JetPFLooseId",          &AK8JetPFLooseId_);
    tree->Branch("AK8JetPFTightLepVetoId",   &AK8JetPFTightLepVetoId_);
    tree->Branch("AK8JetSoftDropMass",    &AK8JetSoftDropMass_);
    tree->Branch("AK8JetSoftDropMassCorr",&AK8JetSoftDropMassCorr_);
    tree->Branch("AK8JetPrunedMass",      &AK8JetPrunedMass_);
    tree->Branch("AK8JetPrunedMassCorr",      &AK8JetPrunedMassCorr_);
    tree->Branch("AK8JetpfBoostedDSVBTag",   &AK8JetpfBoostedDSVBTag_);
    tree->Branch("AK8JetDSVnewV4",   &AK8JetDSVnewV4_);
    tree->Branch("AK8JetCSV",        &AK8JetCSV_);
    tree->Branch("AK8JetJECUnc",             &AK8JetJECUnc_);
    tree->Branch("AK8JetL2L3corr",             &AK8JetL2L3corr_);
    tree->Branch("AK8puppiPt", &AK8puppiPt_);
    tree->Branch("AK8puppiMass", &AK8puppiMass_);
    tree->Branch("AK8puppiEta", &AK8puppiEta_);
    tree->Branch("AK8puppiPhi", &AK8puppiPhi_);
    tree->Branch("AK8puppiTau1", &AK8puppiTau1_);
    tree->Branch("AK8puppiTau2", &AK8puppiTau2_);
    tree->Branch("AK8puppiTau3", &AK8puppiTau3_);
    tree->Branch("AK8puppiSDL2L3corr", &AK8puppiSDL2L3corr_);
    tree->Branch("AK8puppiSDMass", &AK8puppiSDMass_);
    tree->Branch("AK8puppiSDMassL2L3Corr", &AK8puppiSDMassL2L3Corr_);

    if (doGenParticles_){
      tree->Branch("AK8JetPartonID", &AK8JetPartonID_);
      tree->Branch("AK8JetHadFlvr", &AK8JetHadFlvr_);
      tree->Branch("AK8JetGenJetIndex", &AK8JetGenJetIndex_);
      tree->Branch("AK8JetGenJetEn", &AK8JetGenJetEn_);
      tree->Branch("AK8JetGenJetPt", &AK8JetGenJetPt_);
      tree->Branch("AK8JetGenJetEta", &AK8JetGenJetEta_);
      tree->Branch("AK8JetGenJetPhi", &AK8JetGenJetPhi_);
      tree->Branch("AK8JetGenPartonID", &AK8JetGenPartonID_);
      tree->Branch("AK8JetGenEn", &AK8JetGenEn_);
      tree->Branch("AK8JetGenPt", &AK8JetGenPt_);
      tree->Branch("AK8JetGenEta", &AK8JetGenEta_);
      tree->Branch("AK8JetGenPhi", &AK8JetGenPhi_);
      tree->Branch("AK8JetGenPartonMomID", &AK8JetGenPartonMomID_);
      tree->Branch("AK8JetP4Smear",        &AK8JetP4Smear_);
      tree->Branch("AK8JetP4SmearUp",      &AK8JetP4SmearUp_);
      tree->Branch("AK8JetP4SmearDo",      &AK8JetP4SmearDo_);
    }
    tree->Branch("nAK8SDSJ",       &nAK8SDSJ_);
    tree->Branch("AK8SDSJPt",      &AK8SDSJPt_);
    tree->Branch("AK8SDSJEta",     &AK8SDSJEta_);
    tree->Branch("AK8SDSJPhi",     &AK8SDSJPhi_);
    tree->Branch("AK8SDSJMass",    &AK8SDSJMass_);
    tree->Branch("AK8SDSJE",       &AK8SDSJE_);
    tree->Branch("AK8SDSJCharge",  &AK8SDSJCharge_);
    tree->Branch("AK8SDSJFlavour", &AK8SDSJFlavour_);
    tree->Branch("AK8SDSJCSV",     &AK8SDSJCSV_);
    tree->Branch("nAK8puppiSDSJ",       &nAK8puppiSDSJ_);
    tree->Branch("AK8puppiSDSJPt",      &AK8puppiSDSJPt_);
    tree->Branch("AK8puppiSDSJEta",     &AK8puppiSDSJEta_);
    tree->Branch("AK8puppiSDSJPhi",     &AK8puppiSDSJPhi_);
    tree->Branch("AK8puppiSDSJMass",    &AK8puppiSDSJMass_);
    tree->Branch("AK8puppiSDSJE",       &AK8puppiSDSJE_);
    tree->Branch("AK8puppiSDSJCharge",  &AK8puppiSDSJCharge_);
    tree->Branch("AK8puppiSDSJFlavour", &AK8puppiSDSJFlavour_);
    tree->Branch("AK8puppiSDSJCSV",     &AK8puppiSDSJCSV_);
  }
}


void ggNtuplizer::fillJets(const edm::Event& e, const edm::EventSetup& es) {

  // cleanup from previous execution
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
  jetJetProbabilityBJetTags_              .clear();
  jetpfCombinedMVAV2BJetTags_             .clear();
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
  
  // SubJet
  AK8JetPt_           .clear();
  AK8JetEn_           .clear();
  AK8JetRawPt_        .clear();
  AK8JetRawEn_        .clear();
  AK8JetEta_          .clear();
  AK8JetPhi_          .clear();
  AK8JetMass_         .clear();
  AK8Jet_tau1_        .clear();
  AK8Jet_tau2_        .clear();
  AK8Jet_tau3_        .clear();
  AK8JetCHF_          .clear();
  AK8JetNHF_          .clear();
  AK8JetCEF_          .clear();
  AK8JetNEF_          .clear();
  AK8JetNCH_            .clear();
  AK8JetNNP_            .clear();
  AK8Jetnconstituents_     .clear();
  AK8JetMUF_          .clear();
  AK8JetPFLooseId_         .clear();
  AK8JetPFTightLepVetoId_  .clear();
  AK8JetSoftDropMass_   .clear();
  AK8JetSoftDropMassCorr_   .clear();
  AK8JetPrunedMass_   .clear();
  AK8JetPrunedMassCorr_   .clear();
  AK8JetpfBoostedDSVBTag_  .clear();
  AK8JetDSVnewV4_  .clear();
  AK8JetCSV_   .clear();
  AK8JetJECUnc_            .clear();
  AK8JetL2L3corr_            .clear();

  AK8puppiPt_                .clear();
  AK8puppiMass_                .clear();
  AK8puppiEta_                .clear();
  AK8puppiPhi_                .clear();
  AK8puppiTau1_                .clear();
  AK8puppiTau2_                .clear();
  AK8puppiTau3_                .clear();
  AK8puppiSDL2L3corr_                .clear();
  AK8puppiSDMass_                .clear();
  AK8puppiSDMassL2L3Corr_                .clear();

  AK8JetPartonID_ .clear();
  AK8JetHadFlvr_ .clear();
  AK8JetGenJetIndex_.clear();
  AK8JetGenJetEn_.clear();
  AK8JetGenJetPt_.clear();
  AK8JetGenJetEta_.clear();
  AK8JetGenJetPhi_.clear();
  AK8JetGenPartonID_.clear();
  AK8JetGenEn_.clear();
  AK8JetGenPt_.clear();
  AK8JetGenEta_.clear();
  AK8JetGenPhi_.clear();
  AK8JetGenPartonMomID_.clear();
  AK8JetP4Smear_.clear();
  AK8JetP4SmearUp_.clear();
  AK8JetP4SmearDo_.clear();

  nAK8SDSJ_ .clear();
  AK8SDSJPt_ .clear();
  AK8SDSJEta_ .clear();
  AK8SDSJPhi_ .clear();
  AK8SDSJMass_ .clear();
  AK8SDSJCharge_ .clear();
  AK8SDSJE_ .clear();
  AK8SDSJFlavour_ .clear();
  AK8SDSJCSV_ .clear();

  nAK8puppiSDSJ_ .clear();
  AK8puppiSDSJPt_ .clear();
  AK8puppiSDSJEta_ .clear();
  AK8puppiSDSJPhi_ .clear();
  AK8puppiSDSJMass_ .clear();
  AK8puppiSDSJCharge_ .clear();
  AK8puppiSDSJE_ .clear();
  AK8puppiSDSJFlavour_ .clear();
  AK8puppiSDSJCSV_ .clear();

  nJet_ = 0;

  edm::Handle<edm::View<pat::Jet> > jetHandle;
  e.getByToken(jetsAK4Label_, jetHandle);

  if (!jetHandle.isValid()) {
    edm::LogWarning("ggNtuplizer") << "no pat::Jets (AK4) in event";
    return;
  }

  edm::Handle<vector<reco::GenParticle> > genParticlesHandle;
  if(doGenParticles_)e.getByToken(genParticlesCollection_, genParticlesHandle);
  
  edm::Handle<double> rhoHandle;
  e.getByToken(rhoLabel_, rhoHandle);
  float rho = *(rhoHandle.product());
  
  edm::Handle<reco::VertexCollection> vtxHandle;
  e.getByToken(vtxLabel_, vtxHandle);
  if (!vtxHandle.isValid()) edm::LogWarning("ggNtuplizer") << "Primary vertices info not unavailable";
  
  // Accessing the JEC uncertainties 
  //ak4  
  edm::ESHandle<JetCorrectorParametersCollection> JetCorParColl;
  es.get<JetCorrectionsRecord>().get("AK4PFchs",JetCorParColl); 
  JetCorrectorParameters const & JetCorPar = (*JetCorParColl)["Uncertainty"];
  JetCorrectionUncertainty *jecUnc=0;
  jecUnc = new JetCorrectionUncertainty(JetCorPar);
  //ak8
  edm::ESHandle<JetCorrectorParametersCollection> AK8JetCorParColl;
  es.get<JetCorrectionsRecord>().get("AK8PFchs",AK8JetCorParColl);
  JetCorrectorParameters const & AK8JetCorPar = (*AK8JetCorParColl)["Uncertainty"];
  JetCorrectionUncertainty *AK8jecUnc=0;
  AK8jecUnc = new JetCorrectionUncertainty(AK8JetCorPar);
  
  //start jets Lvdp
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
    
    //b-tagging
    jetCSV2BJetTags_           .push_back(iJet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
    jetJetProbabilityBJetTags_ .push_back(iJet->bDiscriminator("pfJetProbabilityBJetTags"));
    jetpfCombinedMVAV2BJetTags_.push_back(iJet->bDiscriminator("pfCombinedMVAV2BJetTags"));
  
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
  
  if (dumpSubJets_) {
    edm::Handle<edm::View<pat::Jet> > jetsAK8;
    e.getByToken(jetsAK8Label_, jetsAK8);
    
    if (!jetsAK8.isValid()) {
      edm::LogWarning("ggNtuplizer") << "no pat::Jets (AK8AK8) in event";
      return;
    }
    
    //Double B-tagger 
//    edm::Handle<reco::JetTagCollection> pfBoostedDoubleSecondaryVertex; 
//    e.getByToken(boostedDoubleSVLabel_, pfBoostedDoubleSecondaryVertex); 
    
    //https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetWtagging#Recipes_to_apply_JEC_on_the_prun
    
    std::vector<JetCorrectorParameters> vPar;
    for ( std::vector<std::string>::const_iterator payloadBegin = jecAK8PayloadNames_.begin(), payloadEnd = jecAK8PayloadNames_.end(), ipayload = payloadBegin; ipayload != payloadEnd; ++ipayload ) {
      JetCorrectorParameters pars(*ipayload);
    vPar.push_back(pars);
    }
    
    // Make the FactorizedJetCorrector
    jecAK8_ = boost::shared_ptr<FactorizedJetCorrector> ( new FactorizedJetCorrector(vPar) );
    jecAK8pSD_ = boost::shared_ptr<FactorizedJetCorrector> ( new FactorizedJetCorrector(vPar) );
       
    nAK8Jet_ = 0;
    //jet substructure
    int nsubjets = 0;
    std::vector<float> vecSDSJcsv ;
    std::vector<float> vecSDSJpt ;
    std::vector<float> vecSDSJeta ;
    std::vector<float> vecSDSJmass ;
    std::vector<float> vecSDSJphi ;
    std::vector<float> vecSDSJe ;
    std::vector<int > vecSDSJcharge ;
    std::vector<int > vecSDSJflavour;

    int nPuppiSJs = 0;
    std::vector<float> vecPuppiSDSJcsv ;
    std::vector<float> vecPuppiSDSJpt ;
    std::vector<float> vecPuppiSDSJeta ;
    std::vector<float> vecPuppiSDSJmass ;
    std::vector<float> vecPuppiSDSJphi ;
    std::vector<float> vecPuppiSDSJe ;
    std::vector<int > vecPuppiSDSJcharge ;
    std::vector<int > vecPuppiSDSJflavour;

    edm::View<pat::Jet>::const_iterator beginAK8 = jetsAK8->begin();
    edm::View<pat::Jet>::const_iterator endAK8 = jetsAK8->end();
    edm::View<pat::Jet>::const_iterator ijetAK8 = beginAK8;
    int ijetRef = -1;
    // Loop over the "hard" jets
    for(ijetAK8 = beginAK8; ijetAK8 != endAK8; ++ijetAK8 ) {
      ijetRef++;
      if( ijetAK8->pt() < 30.0 ) continue;
      nAK8Jet_++;
      AK8JetPt_.push_back( ijetAK8->pt() );
      AK8JetEn_.push_back( ijetAK8->energy() );
      AK8JetMass_.push_back( ijetAK8->mass() );
      AK8JetRawPt_.push_back( (*ijetAK8).correctedJet("Uncorrected").pt() );
      AK8JetRawEn_.push_back( (*ijetAK8).correctedJet("Uncorrected").energy() );
      AK8JetEta_.push_back( ijetAK8->eta() );
      AK8JetPhi_.push_back( ijetAK8->phi() );
      AK8Jet_tau1_.push_back( ijetAK8->userFloat("NjettinessAK8:tau1") );
      AK8Jet_tau2_.push_back( ijetAK8->userFloat("NjettinessAK8:tau2") );
      AK8Jet_tau3_.push_back( ijetAK8->userFloat("NjettinessAK8:tau3") );
      AK8JetCHF_.push_back( ijetAK8->chargedHadronEnergyFraction()); // 0.0
      AK8JetNHF_.push_back( ijetAK8->neutralHadronEnergyFraction()); //0.99
      AK8JetCEF_.push_back( ijetAK8->chargedEmEnergyFraction()); //0.99
      AK8JetNEF_.push_back( ijetAK8->neutralEmEnergyFraction()); //0.99
      AK8JetNCH_.push_back( ijetAK8->chargedMultiplicity()); //0
      AK8JetNNP_.push_back( ijetAK8->neutralMultiplicity()); //0
      AK8Jetnconstituents_.push_back( ijetAK8->chargedMultiplicity() + ijetAK8->neutralMultiplicity()); //1  
      AK8JetMUF_.push_back(ijetAK8->muonEnergyFraction());
      //https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID#Recommendations_for_13_TeV_data
      //LooseID
      bool AK8jetID = true;
      if (fabs(ijetAK8->eta()) <= 2.7) {
        if (!(ijetAK8->neutralHadronEnergyFraction() < 0.99))                       AK8jetID = false;
        if (!(ijetAK8->neutralEmEnergyFraction() < 0.99))                           AK8jetID = false;
        if (!((ijetAK8->chargedMultiplicity() + ijetAK8->neutralMultiplicity()) > 1))  AK8jetID = false;
        if (fabs(ijetAK8->eta()) <= 2.4) {
          if (!(ijetAK8->chargedHadronEnergyFraction() > 0))  AK8jetID = false;
          if (!(ijetAK8->chargedMultiplicity() > 0))          AK8jetID = false;
          if (!(ijetAK8->chargedEmEnergyFraction() < 0.99))   AK8jetID = false;
        }
      }
      else if (fabs(ijetAK8->eta()) > 2.7 && fabs(ijetAK8->eta()) <= 3.0) {
        if (!(ijetAK8->neutralEmEnergyFraction() > 0.01))     AK8jetID = false;
	if (!(ijetAK8->neutralHadronEnergyFraction() < 0.98)) AK8jetID = false;
        if (!(ijetAK8->neutralMultiplicity() > 2))            AK8jetID = false;
      }
      else if (fabs(ijetAK8->eta()) > 3.0) {
        if (!(ijetAK8->neutralEmEnergyFraction() < 0.90))  AK8jetID = false;
        if (!(ijetAK8->neutralMultiplicity() > 10))        AK8jetID = false;
      }
      AK8JetPFLooseId_.push_back(AK8jetID);
      //TightIDMuon Fraction
      bool AK8jetIDTightLepVeto = true;
      if (fabs(ijetAK8->eta()) <= 3.0) {
        if (!(ijetAK8->neutralHadronEnergyFraction() < 0.90))                       AK8jetIDTightLepVeto = false;
        if (!(ijetAK8->neutralEmEnergyFraction() < 0.90))                           AK8jetIDTightLepVeto = false;
        if (!((ijetAK8->chargedMultiplicity() + ijetAK8->neutralMultiplicity()) > 1))  AK8jetIDTightLepVeto = false;
        if (!(ijetAK8->muonEnergyFraction()<0.8))                                   AK8jetIDTightLepVeto = false;
        if (fabs(ijetAK8->eta()) <= 2.4) {
          if (!(ijetAK8->chargedHadronEnergyFraction() > 0))  AK8jetIDTightLepVeto = false;
          if (!(ijetAK8->chargedMultiplicity() > 0))          AK8jetIDTightLepVeto = false;
          if (!(ijetAK8->chargedEmEnergyFraction() < 0.90))   AK8jetIDTightLepVeto = false;
        }
      }
      AK8JetPFTightLepVetoId_.push_back(AK8jetIDTightLepVeto);

      AK8JetpfBoostedDSVBTag_.push_back(ijetAK8->bDiscriminator("pfBoostedDoubleSecondaryVertexAK8BJetTags"));
      AK8JetDSVnewV4_.push_back(ijetAK8->bDiscriminator("newV4pfBoostedDoubleSecondaryVertexAK8BJetTags"));
      AK8JetCSV_.push_back(ijetAK8->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
      //AK8JetpfBoostedDSVBTag_.push_back(     (*pfBoostedDoubleSecondaryVertex).value(ijetRef));

      const LorentzVector uncorrJet = (*ijetAK8).correctedP4(0);
      jecAK8_->setJetEta(uncorrJet.eta());
      jecAK8_->setJetPt ( uncorrJet.pt() );
      jecAK8_->setJetE  ( uncorrJet.energy() );
      jecAK8_->setJetA  ( (*ijetAK8).jetArea() );
      jecAK8_->setRho   ( rho );
      jecAK8_->setNPV   ( vtxHandle->size() );

      float corr = jecAK8_->getCorrection();
      AK8JetL2L3corr_.push_back(corr);
	    
      AK8JetSoftDropMass_.push_back(ijetAK8->userFloat("ak8PFJetsCHSSoftDropMass"));
      AK8JetPrunedMass_.push_back(ijetAK8->userFloat("ak8PFJetsCHSPrunedMass"));
      AK8JetSoftDropMassCorr_.push_back(corr*(ijetAK8->userFloat("ak8PFJetsCHSSoftDropMass")));
      AK8JetPrunedMassCorr_.push_back(corr*(ijetAK8->userFloat("ak8PFJetsCHSPrunedMass")));

      //JEC uncertainty
      if (fabs(ijetAK8->eta()) < 5.2) {
        AK8jecUnc->setJetEta(ijetAK8->eta());
        AK8jecUnc->setJetPt(ijetAK8->pt()); // here you must use the CORRECTED jet pt
        AK8JetJECUnc_.push_back(AK8jecUnc->getUncertainty(true));
      } else {
        AK8JetJECUnc_.push_back(-1.);
      }

      //save gen-info for ak8 jets
      //parton id                                                                                                                                                           
      AK8JetPartonID_.push_back(ijetAK8->partonFlavour());
      AK8JetHadFlvr_.push_back(ijetAK8->hadronFlavour());
      int AK8JetGenPartonID = -99;
      int AK8JetGenPartonMomID = -99;
      float AK8JetGenEn = -999.;
      float AK8JetGenPt = -999.;
      float AK8JetGenEta = -999.;
      float AK8JetGenPhi = -999.;
      if (doGenParticles_ && genParticlesHandle.isValid() ) {
        if ((*ijetAK8).genParton()) {
          AK8JetGenPartonID = (*ijetAK8).genParton()->pdgId();
          AK8JetGenEn = (*ijetAK8).genParton()->energy();
          AK8JetGenPt = (*ijetAK8).genParton()->pt();
          AK8JetGenEta = (*ijetAK8).genParton()->eta();
          AK8JetGenPhi = (*ijetAK8).genParton()->phi();
          if ((*ijetAK8).genParton()->mother()) {
            AK8JetGenPartonMomID = (*ijetAK8).genParton()->mother()->pdgId();
          }
        }
      }
      AK8JetGenPartonID_.push_back(AK8JetGenPartonID);
      AK8JetGenPartonMomID_.push_back(AK8JetGenPartonMomID);
      AK8JetGenEn_ .push_back(AK8JetGenEn);
      AK8JetGenPt_ .push_back(AK8JetGenPt);
      AK8JetGenEta_ .push_back(AK8JetGenEta);
      AK8JetGenPhi_ .push_back(AK8JetGenPhi);
      int AK8JetGenJetIndex = -1;
      float AK8JetGenJetEn = -999.;
      float AK8JetGenJetPt = -999.;
      float AK8JetGenJetEta = -999.;
      float AK8JetGenJetPhi = -999.;
      if (doGenParticles_ && genParticlesHandle.isValid() ) {
        if ((*ijetAK8).genJet()) {
          AK8JetGenJetIndex = 1;
          AK8JetGenJetEn = (*ijetAK8).genJet()->energy();
          AK8JetGenJetPt = (*ijetAK8).genJet()->pt();
          AK8JetGenJetEta = (*ijetAK8).genJet()->eta();
          AK8JetGenJetPhi = (*ijetAK8).genJet()->phi();
        }
        // access AK8jet resolution       
        JME::JetParameters AK8parameters;
        AK8parameters.setJetPt(ijetAK8->pt()).setJetEta(ijetAK8->eta()).setRho(rho);
        float AK8jetResolution = AK8jetResolution_.getResolution(AK8parameters);

        edm::Service<edm::RandomNumberGenerator> rng;
        if (!rng.isAvailable()) edm::LogError("JET : random number generator is missing !");
        CLHEP::HepRandomEngine & engine = rng->getEngine( e.streamID() );
        float rnd = CLHEP::RandGauss::shoot(&engine, 0., AK8jetResolution);

        float AK8jetResolutionSF   = AK8jetResolutionSF_.getScaleFactor(AK8parameters);
        float AK8jetResolutionSFUp = AK8jetResolutionSF_.getScaleFactor(AK8parameters, Variation::UP);
        float AK8jetResolutionSFDo = AK8jetResolutionSF_.getScaleFactor(AK8parameters, Variation::DOWN);

        float AK8JetP4Smear   = -1.;
        float AK8JetP4SmearUp = -1.;
        float AK8JetP4SmearDo = -1.;
        if (AK8JetGenJetPt > 0 && deltaR(ijetAK8->eta(), ijetAK8->phi(), AK8JetGenJetEta, AK8JetGenJetPhi) < 0.4 && fabs(ijetAK8->pt()-AK8JetGenJetPt) < 3*AK8jetResolution*ijetAK8->pt()) {
          AK8JetP4Smear   = 1. + (AK8jetResolutionSF   - 1.)*(ijetAK8->pt() - AK8JetGenJetPt)/ijetAK8->pt();
          AK8JetP4SmearUp = 1. + (AK8jetResolutionSFUp - 1.)*(ijetAK8->pt() - AK8JetGenJetPt)/ijetAK8->pt();
          AK8JetP4SmearDo = 1. + (AK8jetResolutionSFDo - 1.)*(ijetAK8->pt() - AK8JetGenJetPt)/ijetAK8->pt();
        } else {
          AK8JetP4Smear   = 1. + rnd*sqrt(max(pow(AK8jetResolutionSF,   2)-1, 0.));
          AK8JetP4SmearUp = 1. + rnd*sqrt(max(pow(AK8jetResolutionSFUp, 2)-1, 0.));
          AK8JetP4SmearDo = 1. + rnd*sqrt(max(pow(AK8jetResolutionSFDo, 2)-1, 0.));
        }
        AK8JetP4Smear_  .push_back(AK8JetP4Smear);
        AK8JetP4SmearUp_.push_back(AK8JetP4SmearUp);
        AK8JetP4SmearDo_.push_back(AK8JetP4SmearDo);
      }
      AK8JetGenJetIndex_.push_back(AK8JetGenJetIndex);
      AK8JetGenJetEn_.push_back(AK8JetGenJetEn);
      AK8JetGenJetPt_.push_back(AK8JetGenJetPt);
      AK8JetGenJetEta_.push_back(AK8JetGenJetEta);
      AK8JetGenJetPhi_.push_back(AK8JetGenJetPhi);
      
      //save Softdrop subjet info Lvdp
      vecSDSJcsv.clear();
      vecSDSJpt.clear();
      vecSDSJeta.clear();
      vecSDSJmass.clear();
      vecSDSJphi.clear();
      vecSDSJe.clear();
      vecSDSJcharge.clear();
      vecSDSJflavour.clear();
      nsubjets = 0;
      if(dumpSoftDrop_) {
        auto const & sdSubjets = ijetAK8->subjets("SoftDrop");
        for ( auto const & SDSJ : sdSubjets ) {
          nsubjets++;
          vecSDSJpt.push_back(SDSJ->pt());
          vecSDSJeta.push_back(SDSJ->eta());
          vecSDSJmass.push_back(SDSJ->mass());
          vecSDSJphi.push_back(SDSJ->phi());
          vecSDSJe.push_back(SDSJ->energy());
          vecSDSJflavour.push_back(abs(SDSJ->partonFlavour()));
          vecSDSJcharge.push_back(SDSJ->charge());
          vecSDSJcsv.push_back(SDSJ->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") );
        }
      }
      nAK8SDSJ_.push_back(nsubjets);
      AK8SDSJPt_.push_back(vecSDSJpt);
      AK8SDSJEta_.push_back(vecSDSJeta);
      AK8SDSJPhi_.push_back(vecSDSJphi);
      AK8SDSJMass_.push_back(vecSDSJmass);
      AK8SDSJE_.push_back(vecSDSJe);
      AK8SDSJCharge_.push_back(vecSDSJcharge);
      AK8SDSJFlavour_.push_back(vecSDSJflavour);
      AK8SDSJCSV_.push_back(vecSDSJcsv);


//for some r&d on puppi + softdrop

      AK8puppiPt_.push_back( ijetAK8->userFloat("ak8PFJetsPuppiValueMap:pt"));
      AK8puppiMass_.push_back( ijetAK8->userFloat("ak8PFJetsPuppiValueMap:mass"));
      AK8puppiEta_.push_back( ijetAK8->userFloat("ak8PFJetsPuppiValueMap:eta"));
      AK8puppiPhi_.push_back( ijetAK8->userFloat("ak8PFJetsPuppiValueMap:phi"));
      AK8puppiTau1_.push_back( ijetAK8->userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau1"));
      AK8puppiTau2_.push_back( ijetAK8->userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau2"));
      AK8puppiTau3_.push_back( ijetAK8->userFloat("ak8PFJetsPuppiValueMap:NjettinessAK8PuppiTau3"));

      //save puppi-Softdrop subjet info Lvdp
      vecPuppiSDSJcsv.clear();
      vecPuppiSDSJpt.clear();
      vecPuppiSDSJeta.clear();
      vecPuppiSDSJmass.clear();
      vecPuppiSDSJphi.clear();
      vecPuppiSDSJe.clear();
      vecPuppiSDSJcharge.clear();
      vecPuppiSDSJflavour.clear();
      nPuppiSJs = 0;

      TLorentzVector puppi_softdrop, puppi_softdrop_subjet;
      auto const & sdSubjetsPuppi = ijetAK8->subjets("SoftDropPuppi");
      for ( auto const & puppiSDSJ : sdSubjetsPuppi ) {
        nPuppiSJs++;
        vecPuppiSDSJpt.push_back(puppiSDSJ->pt());
        vecPuppiSDSJeta.push_back(puppiSDSJ->eta());
        vecPuppiSDSJmass.push_back(puppiSDSJ->mass());
        vecPuppiSDSJphi.push_back(puppiSDSJ->phi());
        vecPuppiSDSJe.push_back(puppiSDSJ->energy());
        vecPuppiSDSJflavour.push_back(abs(puppiSDSJ->partonFlavour()));
        vecPuppiSDSJcharge.push_back(puppiSDSJ->charge());
        vecPuppiSDSJcsv.push_back(puppiSDSJ->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") );

        puppi_softdrop_subjet.SetPtEtaPhiM(puppiSDSJ->correctedP4(0).pt(),puppiSDSJ->correctedP4(0).eta(),puppiSDSJ->correctedP4(0).phi(),puppiSDSJ->correctedP4(0).mass());
        puppi_softdrop+=puppi_softdrop_subjet;
      }
//fir L2L3 corrections
      jecAK8pSD_->setJetEta( puppi_softdrop.Eta() );
      jecAK8pSD_->setJetPt ( puppi_softdrop.Pt() );
      jecAK8pSD_->setJetE  ( puppi_softdrop.E() );
      jecAK8pSD_->setJetA  ( (*ijetAK8).jetArea() );
      jecAK8pSD_->setRho   ( rho );
      jecAK8pSD_->setNPV   ( vtxHandle->size() );

      float corr_puppiSD = jecAK8pSD_->getCorrection();
      AK8puppiSDL2L3corr_.push_back(corr_puppiSD);
      AK8puppiSDMass_.push_back(puppi_softdrop.M());
      AK8puppiSDMassL2L3Corr_.push_back(corr_puppiSD*puppi_softdrop.M());

      nAK8puppiSDSJ_.push_back(nPuppiSJs);
      AK8puppiSDSJPt_.push_back(vecPuppiSDSJpt);
      AK8puppiSDSJEta_.push_back(vecPuppiSDSJeta);
      AK8puppiSDSJPhi_.push_back(vecPuppiSDSJphi);
      AK8puppiSDSJMass_.push_back(vecPuppiSDSJmass);
      AK8puppiSDSJE_.push_back(vecPuppiSDSJe);
      AK8puppiSDSJCharge_.push_back(vecPuppiSDSJcharge);
      AK8puppiSDSJFlavour_.push_back(vecPuppiSDSJflavour);
      AK8puppiSDSJCSV_.push_back(vecPuppiSDSJcsv);
    }
  }
  delete jecUnc;
  delete AK8jecUnc;
}
