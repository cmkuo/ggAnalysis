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
vector<float>  jetpfCombinedInclusiveSecondaryVertexV2BJetTags_; // recommended
vector<float>  jetJetProbabilityBJetTags_;
vector<float>  jetpfCombinedMVABJetTags_;
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

vector<int> nAK8softdropSubjet_ ;
vector< vector<float> > AK8softdropSubjetPt_ ;
vector< vector<float> > AK8softdropSubjetEta_ ;
vector< vector<float> > AK8softdropSubjetMass_ ;
vector< vector<float> > AK8softdropSubjetPhi_ ;
vector< vector<float> > AK8softdropSubjetE_ ;
vector< vector<int > > AK8softdropSubjetCharge_ ;
vector< vector<int > > AK8softdropSubjetFlavour_;
vector< vector<float> > AK8softdropSubjetCSV_ ;

//
void ggNtuplizer::branchesJets(TTree* tree)
{
  tree->Branch("nJet",   &nJet_);
  tree->Branch("jetPt",  &jetPt_);
  //tree->Branch("jetEn",  &jetEn_);
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
  tree->Branch("jetpfCombinedInclusiveSecondaryVertexV2BJetTags", &jetpfCombinedInclusiveSecondaryVertexV2BJetTags_);
  tree->Branch("jetJetProbabilityBJetTags", &jetJetProbabilityBJetTags_);
  tree->Branch("jetpfCombinedMVABJetTags", &jetpfCombinedMVABJetTags_);
  if (doGenParticles_) tree->Branch("jetPartonID", &jetPartonID_);
  tree->Branch("jetPFLooseId", &jetPFLooseId_);

  // SubJet
  if (dumpSubJets_) {
    tree->Branch("nAK8Jet",             &nAK8Jet_);
    tree->Branch("AK8JetPt",            &AK8JetPt_);
    //tree->Branch("AK8JetEn",            &AK8JetEn_);
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
    tree->Branch("nAK8softdropSubjet",    &nAK8softdropSubjet_);
    tree->Branch("AK8softdropSubjetPt",    &AK8softdropSubjetPt_);
    tree->Branch("AK8softdropSubjetEta",    &AK8softdropSubjetEta_);
    tree->Branch("AK8softdropSubjetPhi",    &AK8softdropSubjetPhi_);
    tree->Branch("AK8softdropSubjetMass",    &AK8softdropSubjetMass_);
    tree->Branch("AK8softdropSubjetE",    &AK8softdropSubjetE_);
    tree->Branch("AK8softdropSubjetCharge",    &AK8softdropSubjetCharge_);
    tree->Branch("AK8softdropSubjetFlavour",    &AK8softdropSubjetFlavour_);
    tree->Branch("AK8softdropSubjetCSV",    &AK8softdropSubjetCSV_);


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
  jetpfCombinedInclusiveSecondaryVertexV2BJetTags_.clear();
  jetJetProbabilityBJetTags_              .clear();
  jetpfCombinedMVABJetTags_             .clear();
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
  nAK8softdropSubjet_ .clear();
  AK8softdropSubjetPt_ .clear();
  AK8softdropSubjetEta_ .clear();
  AK8softdropSubjetPhi_ .clear();
  AK8softdropSubjetMass_ .clear();
  AK8softdropSubjetCharge_ .clear();
  AK8softdropSubjetE_ .clear();
  AK8softdropSubjetFlavour_ .clear();
  AK8softdropSubjetCSV_ .clear();

  nJet_ = 0;

  edm::Handle<edm::View<pat::Jet> > jetHandle;
  e.getByToken(jetsAK4Label_, jetHandle);

  if (!jetHandle.isValid()) {
    edm::LogWarning("ggNtuplizer") << "no pat::Jets (AK4) in event";
    return;
  }

  //start jets Lvdp
  for (edm::View<pat::Jet>::const_iterator iJet = jetHandle->begin(); iJet != jetHandle->end(); ++iJet) {
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
    jetpfCombinedInclusiveSecondaryVertexV2BJetTags_.push_back(iJet->bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
    jetJetProbabilityBJetTags_.push_back(iJet->bDiscriminator("pfJetProbabilityBJetTags"));
    jetpfCombinedMVABJetTags_.push_back(iJet->bDiscriminator("pfCombinedMVABJetTags"));
  
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
    int nsubjets = 0;
    std::vector<float> vecSoftdropSubjetcsv ;
    std::vector<float> vecSoftdropSubjetpt ;
    std::vector<float> vecSoftdropSubjeteta ;
    std::vector<float> vecSoftdropSubjetmass ;
    std::vector<float> vecSoftdropSubjetphi ;
    std::vector<float> vecSoftdropSubjete ;
    std::vector<int > vecSoftdropSubjetcharge ;
    std::vector<int > vecSoftdropSubjetflavour;
    
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

//save Softdrop subjet info Lvdp
      vecSoftdropSubjetcsv.clear();
      vecSoftdropSubjetpt.clear();
      vecSoftdropSubjeteta.clear();
      vecSoftdropSubjetmass.clear();
      vecSoftdropSubjetphi.clear();
      vecSoftdropSubjete.clear();
      vecSoftdropSubjetcharge.clear();
      vecSoftdropSubjetflavour.clear();
      nsubjets = 0;
      const std::vector<edm::Ptr<pat::Jet> > &wSubjets = ijetAK8->subjets("SoftDrop");
      if(ijetAK8->subjets("SoftDrop").size()>0){
	for ( const pat::Jet & softdropsubjet : wSubjets ) {
	  if( softdropsubjet.pt() < 0.01 ) continue;
	  nsubjets++;
	  vecSoftdropSubjetpt.push_back(softdropsubjet.pt());
	  vecSoftdropSubjeteta.push_back(softdropsubjet.eta());
	  vecSoftdropSubjetmass.push_back(softdropsubjet.mass());
	  vecSoftdropSubjetphi.push_back(softdropsubjet.phi());
	  vecSoftdropSubjete.push_back(softdropsubjet.energy());
	  vecSoftdropSubjetflavour.push_back(abs(softdropsubjet.partonFlavour()));
	  vecSoftdropSubjetcharge.push_back(softdropsubjet.charge());
	  vecSoftdropSubjetcsv.push_back(softdropsubjet.bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags") );
	}
      }
      nAK8softdropSubjet_.push_back(nsubjets);
      AK8softdropSubjetPt_.push_back(vecSoftdropSubjetpt);
      AK8softdropSubjetEta_.push_back(vecSoftdropSubjeteta);
      AK8softdropSubjetPhi_.push_back(vecSoftdropSubjetphi);
      AK8softdropSubjetMass_.push_back(vecSoftdropSubjetmass);
      AK8softdropSubjetE_.push_back(vecSoftdropSubjete);
      AK8softdropSubjetCharge_.push_back(vecSoftdropSubjetcharge);
      AK8softdropSubjetFlavour_.push_back(vecSoftdropSubjetflavour);
      AK8softdropSubjetCSV_.push_back(vecSoftdropSubjetcsv);
    }
  }
}
