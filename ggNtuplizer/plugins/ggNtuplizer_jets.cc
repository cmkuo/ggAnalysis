#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "ggAnalysis/ggNtuplizer/interface/ggNtuplizer.h"

using namespace std;

// (local) variables associated with tree branches
Int_t          nJet_;
vector<float>  jetPt_;
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
Int_t          nCA8Jet_;
vector<float>  CA8JetPt_;
vector<float>  CA8JetEta_;
vector<float>  CA8JetPhi_;
vector<float>  CA8JetMass_;
vector<float>  CA8JetArea_;
vector<float>  CA8Jet_tau1_;
vector<float>  CA8Jet_tau2_;
vector<float>  CA8Jet_tau3_;
vector<float>  CA8JetCHF_;
vector<float>  CA8JetNHF_;
vector<float>  CA8JetCEF_;
vector<float>  CA8JetNEF_;
vector<int>    CA8JetNCH_;
vector<int>    CA8Jetnconstituents_;
vector<float>  CA8prunedJetMass_;

void ggNtuplizer::branchesJets(TTree* tree)
{
  tree->Branch("nJet",   &nJet_);
  tree->Branch("jetPt",  &jetPt_);
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
    tree->Branch("nCA8Jet",             &nCA8Jet_);
    tree->Branch("CA8JetPt",            &CA8JetPt_);
    tree->Branch("CA8JetEta",           &CA8JetEta_);
    tree->Branch("CA8JetPhi",           &CA8JetPhi_);
    tree->Branch("CA8JetMass",          &CA8JetMass_);
    tree->Branch("CA8JetArea",          &CA8JetArea_);
    tree->Branch("CA8Jet_tau1",         &CA8Jet_tau1_);
    tree->Branch("CA8Jet_tau2",         &CA8Jet_tau2_);
    tree->Branch("CA8Jet_tau3",         &CA8Jet_tau3_);
    tree->Branch("CA8JetCHF",           &CA8JetCHF_);
    tree->Branch("CA8JetNHF",           &CA8JetNHF_);
    tree->Branch("CA8JetCEF",           &CA8JetCEF_);
    tree->Branch("CA8JetNEF",           &CA8JetNEF_);
    tree->Branch("CA8JetNCH",           &CA8JetNCH_);
    tree->Branch("CA8Jetnconstituents", &CA8Jetnconstituents_);
    tree->Branch("CA8prunedJetMass",    &CA8prunedJetMass_);
  }

}

void ggNtuplizer::fillJets(const edm::Event& e)
{

  // cleanup from previous execution
  jetPt_                                  .clear();
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
  CA8JetPt_           .clear();
  CA8JetEta_          .clear();
  CA8JetPhi_          .clear();
  CA8JetMass_         .clear();
  CA8JetArea_         .clear();
  CA8Jet_tau1_        .clear();
  CA8Jet_tau2_        .clear();
  CA8Jet_tau3_        .clear();
  CA8JetCHF_          .clear();
  CA8JetNHF_          .clear();
  CA8JetCEF_          .clear();
  CA8JetNEF_          .clear();
  CA8JetNCH_          .clear();
  CA8Jetnconstituents_.clear();
  CA8prunedJetMass_   .clear();

  nJet_ = 0;

  edm::Handle<edm::View<pat::Jet> > jetHandle;
  e.getByToken(jetCollection_, jetHandle);

  if (!jetHandle.isValid()) {
    edm::LogWarning("ggNtuplizer") << "no pat::Jets in event";
    return;
  }

  //start jets Lvdp
  for (edm::View<pat::Jet>::const_iterator iJet = jetHandle->begin(); iJet != jetHandle->end(); ++iJet) {
    jetPt_.push_back(    iJet->pt());
    jetEta_.push_back(   iJet->eta());
    jetPhi_.push_back(   iJet->phi());
    jetCEF_.push_back(   iJet->chargedEmEnergyFraction());
    jetNEF_.push_back(   iJet->neutralEmEnergyFraction());
    jetCHF_.push_back(   iJet->chargedHadronEnergyFraction());
    jetNHF_.push_back(   iJet->neutralHadronEnergyFraction());
    jetHFHAE_.push_back( iJet->HFHadronEnergy());
    jetHFEME_.push_back( iJet->HFEMEnergy());
    jetNCH_.push_back(   iJet->chargedMultiplicity());
    jetNConstituents_.push_back(iJet->getPFConstituents().size());
    //b-tagging
    jetCombinedSecondaryVtxBJetTags_.push_back(iJet->bDiscriminator("combinedSecondaryVertexBJetTags"));
    jetJetProbabilityBJetTags_.push_back(iJet->bDiscriminator("jetProbabilityBJetTags"));
    jetJetBProbabilityBJetTags_.push_back(iJet->bDiscriminator("jetBProbabilityBJetTags"));
    jetTrackCountingHighPurBJetTags_.push_back(iJet->bDiscriminator("trackCountingHighPurBJetTags"));
    jetTrackCountingHighEffBJetTags_.push_back(iJet->bDiscriminator("trackCountingHighEffBJetTags"));
    jetSimpleSecondaryVertexHighEffBJetTags_.push_back(iJet->bDiscriminator("simpleSecondaryVertexHighEffBJetTags"));
    jetSimpleSecondaryVertexHighPurBJetTags_.push_back(iJet->bDiscriminator("simpleSecondaryVertexHighPurBJetTags"));
    //parton id
    jetPartonID_.push_back(iJet->partonFlavour());
    //jet PF Loose ID
    pat::strbitset retjet = pfLooseId_.getBitTemplate();
    jetPFLooseId_.push_back(pfLooseId_(*iJet, retjet));
    nJet_++;
  }

  if(dumpSubJets_) {
    edm::Handle<edm::View<pat::Jet> > jetsCHS;
    e.getByToken(jetsCHSLabel_, jetsCHS);

    if (!jetsCHS.isValid()) {
      edm::LogWarning("ggNtuplizer") << "no pat::Jets (CHS) in event";
      return;
    }

    nCA8Jet_ = 0;
    //jet substructure

    edm::View<pat::Jet>::const_iterator beginCHS = jetsCHS->begin();
    edm::View<pat::Jet>::const_iterator endCHS = jetsCHS->end();
    edm::View<pat::Jet>::const_iterator ijetCHS = beginCHS;

    // Loop over the "hard" jets
    for(ijetCHS = beginCHS; ijetCHS != endCHS; ++ijetCHS ) {
      if( ijetCHS->pt() < 30.0 ) continue;
      nCA8Jet_++;
      CA8JetPt_.push_back( ijetCHS->pt() );
      CA8JetEta_.push_back( ijetCHS->eta() );
      CA8JetPhi_.push_back( ijetCHS->phi() );
      CA8JetMass_.push_back( ijetCHS->mass() );
      CA8Jet_tau1_.push_back( ijetCHS->userFloat("NjettinessCA8:tau1") );
      CA8Jet_tau2_.push_back( ijetCHS->userFloat("NjettinessCA8:tau2") );
      CA8Jet_tau3_.push_back( ijetCHS->userFloat("NjettinessCA8:tau3") );

      CA8JetCHF_.push_back( ijetCHS->chargedHadronEnergyFraction()); // 0.0
      CA8JetNHF_.push_back( ( ijetCHS->neutralHadronEnergy() + ijetCHS->HFHadronEnergy() ) / ijetCHS->energy()); //0.99
      CA8JetCEF_.push_back( ijetCHS->chargedEmEnergyFraction()); //0.99
      CA8JetNEF_.push_back( ijetCHS->neutralEmEnergyFraction()); //0.99
      CA8JetNCH_.push_back( ijetCHS->chargedMultiplicity()); //0
      CA8Jetnconstituents_.push_back( ijetCHS->numberOfDaughters()); //1
      CA8prunedJetMass_.push_back(ijetCHS->userFloat("ca8PFJetsCHSPrunedLinks"));
    }
  }
}
