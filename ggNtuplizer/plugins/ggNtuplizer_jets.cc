#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "ggAnalysis/ggNtuplizer/interface/ggNtuplizer.h"

using namespace std;

// (local) variables associated with tree branches
Int_t          nJet_;
vector<float>  jetPt_;
vector<float>  jetEn_;
vector<float>  jetEta_;
vector<float>  jetPhi_;
vector<float>  jetCHF_;
vector<float>  jetNHF_;
vector<float>  jetCEF_;
vector<float>  jetNEF_;
vector<int>    jetNCH_;
vector<float>  jetHFHAE_;
vector<float>  jetHFEME_;
vector<int>    jetNConstituents_;
vector<float>  jetCombinedSecondaryVtxBJetTags_; // recommended
vector<float>  jetJetProbabilityBJetTags_;
vector<float>  jetJetBProbabilityBJetTags_;
vector<float>  jetTrackCountingHighPurBJetTags_;
vector<float>  jetTrackCountingHighEffBJetTags_;
vector<float>  jetSimpleSecondaryVertexHighEffBJetTags_;
vector<float>  jetSimpleSecondaryVertexHighPurBJetTags_;
vector<int>    jetPartonID_;
vector<bool>   jetPFLooseId_;

//SubJet
Int_t          nAK8Jet_;
vector<float>  AK8JetPt_;
vector<float>  AK8JetEn_;
vector<float>  AK8JetEta_;
vector<float>  AK8JetPhi_;
vector<float>  AK8JetMass_;
vector<float>  AK8Jet_tau1_;
vector<float>  AK8Jet_tau2_;
vector<float>  AK8Jet_tau3_;
vector<float>  AK8JetCHF_;
vector<float>  AK8JetNHF_;
vector<float>  AK8JetCEF_;
vector<float>  AK8JetNEF_;
vector<int>    AK8JetNCH_;
vector<int>    AK8Jetnconstituents_;
vector<float>  AK8CHSSoftDropJetMass_;

void ggNtuplizer::branchesJets(TTree* tree)
{
  tree->Branch("nJet",   &nJet_);
  tree->Branch("jetPt",  &jetPt_);
  tree->Branch("jetEn",  &jetEn_);
  tree->Branch("jetEta", &jetEta_);
  tree->Branch("jetPhi", &jetPhi_);
  tree->Branch("jetCHF", &jetCHF_);
  tree->Branch("jetNHF", &jetNHF_);
  tree->Branch("jetCEF", &jetCEF_);
  tree->Branch("jetNEF", &jetNEF_);
  tree->Branch("jetNCH", &jetNCH_);
  tree->Branch("jetHFHAE", &jetHFHAE_);
  tree->Branch("jetHFEME", &jetHFEME_);
  tree->Branch("jetNConstituents", &jetNConstituents_);
  tree->Branch("jetCombinedSecondaryVtxBJetTags", &jetCombinedSecondaryVtxBJetTags_);
  tree->Branch("jetJetProbabilityBJetTags", &jetJetProbabilityBJetTags_);
  tree->Branch("jetJetBProbabilityBJetTags", &jetJetBProbabilityBJetTags_);
  tree->Branch("jetTrackCountingHighPurBJetTags", &jetTrackCountingHighPurBJetTags_);
  tree->Branch("jetTrackCountingHighEffBJetTags", &jetTrackCountingHighEffBJetTags_);
  tree->Branch("jetSimpleSecondaryVertexHighEffBJetTags", &jetSimpleSecondaryVertexHighEffBJetTags_);
  tree->Branch("jetSimpleSecondaryVertexHighPurBJetTags", &jetSimpleSecondaryVertexHighPurBJetTags_);
  if (doGenParticles_) tree->Branch("jetPartonID", &jetPartonID_);
  tree->Branch("jetPFLooseId", &jetPFLooseId_);

  // SubJet
  if (dumpSubJets_) {
    tree->Branch("nAK8Jet",             &nAK8Jet_);
    tree->Branch("AK8JetPt",            &AK8JetPt_);
    tree->Branch("AK8JetEn",            &AK8JetEn_);
    tree->Branch("AK8JetEta",           &AK8JetEta_);
    tree->Branch("AK8JetPhi",           &AK8JetPhi_);
    tree->Branch("AK8JetMass",          &AK8JetMass_);
    tree->Branch("AK8Jet_tau1",         &AK8Jet_tau1_);
    tree->Branch("AK8Jet_tau2",         &AK8Jet_tau2_);
    tree->Branch("AK8Jet_tau3",         &AK8Jet_tau3_);
    tree->Branch("AK8JetCHF",           &AK8JetCHF_);
    tree->Branch("AK8JetNHF",           &AK8JetNHF_);
    tree->Branch("AK8JetCEF",           &AK8JetCEF_);
    tree->Branch("AK8JetNEF",           &AK8JetNEF_);
    tree->Branch("AK8JetNCH",           &AK8JetNCH_);
    tree->Branch("AK8Jetnconstituents", &AK8Jetnconstituents_);
    tree->Branch("AK8CHSSoftDropJetMass",    &AK8CHSSoftDropJetMass_);
  }

}

void ggNtuplizer::fillJets(const edm::Event& e)
{

  // cleanup from previous execution
  jetPt_                                  .clear();
  jetEn_                                  .clear();
  jetEta_                                 .clear();
  jetPhi_                                 .clear();
  jetCHF_                                 .clear();
  jetNHF_                                 .clear();
  jetCEF_                                 .clear();
  jetNEF_                                 .clear();
  jetNCH_                                 .clear();
  jetHFHAE_                               .clear();
  jetHFEME_                               .clear();
  jetNConstituents_                       .clear();
  jetCombinedSecondaryVtxBJetTags_        .clear();
  jetJetProbabilityBJetTags_              .clear();
  jetJetBProbabilityBJetTags_             .clear();
  jetTrackCountingHighPurBJetTags_        .clear();
  jetTrackCountingHighEffBJetTags_        .clear();
  jetSimpleSecondaryVertexHighEffBJetTags_.clear();
  jetSimpleSecondaryVertexHighPurBJetTags_.clear();
  jetPartonID_                            .clear();
  jetPFLooseId_                           .clear();

  // SubJet
  AK8JetPt_           .clear();
  AK8JetEn_           .clear();
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
  AK8JetNCH_          .clear();
  AK8Jetnconstituents_.clear();
  AK8CHSSoftDropJetMass_   .clear();

  nJet_ = 0;

  edm::Handle<edm::View<pat::Jet> > jetHandle;
  e.getByToken(jetsAK4Label_, jetHandle);

  if (!jetHandle.isValid()) {
    edm::LogWarning("ggNtuplizer") << "no pat::Jets (AK4) in event";
    return;
  }

  //start jets Lvdp
  for (edm::View<pat::Jet>::const_iterator iJet = jetHandle->begin(); iJet != jetHandle->end(); ++iJet) {
    //    cout<<iJet->pt() <<endl;
    jetPt_.push_back(    iJet->pt());
    jetEn_.push_back(    iJet->energy());
    jetEta_.push_back(   iJet->eta());
    jetPhi_.push_back(   iJet->phi());
    jetCEF_.push_back(   iJet->chargedEmEnergyFraction());
    jetNEF_.push_back(   iJet->neutralEmEnergyFraction());
    jetCHF_.push_back(   iJet->chargedHadronEnergyFraction());
    jetNHF_.push_back(   iJet->neutralHadronEnergyFraction());
    jetHFHAE_.push_back( iJet->HFHadronEnergy());
    jetHFEME_.push_back( iJet->HFEMEnergy());
    jetNCH_.push_back(   iJet->chargedMultiplicity());
    jetNConstituents_.push_back(iJet->numberOfDaughters());
    //b-tagging
    jetCombinedSecondaryVtxBJetTags_.push_back(iJet->bDiscriminator("combinedSecondaryVertexBJetTags"));
    jetJetProbabilityBJetTags_.push_back(iJet->bDiscriminator("pfJetProbabilityBJetTags"));
    jetJetBProbabilityBJetTags_.push_back(iJet->bDiscriminator("pfJetBProbabilityBJetTags"));
    jetTrackCountingHighPurBJetTags_.push_back(iJet->bDiscriminator("pfTrackCountingHighPurBJetTags"));
    jetTrackCountingHighEffBJetTags_.push_back(iJet->bDiscriminator("pfTrackCountingHighEffBJetTags"));
    jetSimpleSecondaryVertexHighEffBJetTags_.push_back(iJet->bDiscriminator("pfSimpleSecondaryVertexHighEffBJetTags"));
    jetSimpleSecondaryVertexHighPurBJetTags_.push_back(iJet->bDiscriminator("pfSimpleSecondaryVertexHighPurBJetTags"));
  
    //parton id
    jetPartonID_.push_back(iJet->partonFlavour());
    //jet PF Loose ID
    pat::strbitset retjet = pfLooseId_.getBitTemplate();
    jetPFLooseId_.push_back(pfLooseId_(*iJet, retjet));
    nJet_++;
  }

  if(dumpSubJets_) {
    edm::Handle<edm::View<pat::Jet> > jetsAK8;
    e.getByToken(jetsAK8Label_, jetsAK8);

    if (!jetsAK8.isValid()) {
      edm::LogWarning("ggNtuplizer") << "no pat::Jets (AK8AK8) in event";
      return;
    }

    nAK8Jet_ = 0;
    //jet substructure

    edm::View<pat::Jet>::const_iterator beginAK8 = jetsAK8->begin();
    edm::View<pat::Jet>::const_iterator endAK8 = jetsAK8->end();
    edm::View<pat::Jet>::const_iterator ijetAK8 = beginAK8;

    // Loop over the "hard" jets
    for(ijetAK8 = beginAK8; ijetAK8 != endAK8; ++ijetAK8 ) {
      if( ijetAK8->pt() < 30.0 ) continue;

      nAK8Jet_++;
      AK8JetPt_.push_back( ijetAK8->pt() );
      AK8JetEn_.push_back( ijetAK8->energy() );
      AK8JetEta_.push_back( ijetAK8->eta() );
      AK8JetPhi_.push_back( ijetAK8->phi() );
      AK8JetMass_.push_back( ijetAK8->mass() );
      if (isAOD_ )
       {    AK8Jet_tau1_.push_back( ijetAK8->userFloat("NjettinessAK8CHS:tau1") );
            AK8Jet_tau2_.push_back( ijetAK8->userFloat("NjettinessAK8CHS:tau2") );
            AK8Jet_tau3_.push_back( ijetAK8->userFloat("NjettinessAK8CHS:tau3") );
       }
      else 
       {    AK8Jet_tau1_.push_back( ijetAK8->userFloat("NjettinessAK8:tau1") );
            AK8Jet_tau2_.push_back( ijetAK8->userFloat("NjettinessAK8:tau2") );
            AK8Jet_tau3_.push_back( ijetAK8->userFloat("NjettinessAK8:tau3") );
       }
      AK8JetCHF_.push_back( ijetAK8->chargedHadronEnergyFraction()); // 0.0
      AK8JetNHF_.push_back( ( ijetAK8->neutralHadronEnergy() + ijetAK8->HFHadronEnergy() ) / ijetAK8->energy()); //0.99
      AK8JetCEF_.push_back( ijetAK8->chargedEmEnergyFraction()); //0.99
      AK8JetNEF_.push_back( ijetAK8->neutralEmEnergyFraction()); //0.99
      AK8JetNCH_.push_back( ijetAK8->chargedMultiplicity()); //0
      AK8Jetnconstituents_.push_back( ijetAK8->numberOfDaughters()); //1  
      AK8CHSSoftDropJetMass_.push_back(ijetAK8->userFloat("ak8PFJetsCHSSoftDropMass")); //new miniAOD
//      AK8CHSSoftDropJetMass_.push_back(ijetAK8->userFloat("ak8PFJetsCHSPrunedLinks")); //phys14

    }
  }
}
