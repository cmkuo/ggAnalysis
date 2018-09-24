#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "ATGC-Analysis/Analyzer/interface/Analyzer.h"
#include "ATGC-Analysis/Analyzer/interface/GenParticleParentage.h"
using namespace std;

// (local) variables associated with tree branches
Int_t            nGenJets_;
vector<float>  genJetEn_;
vector<float>  genJetPt_;
vector<float>  genJetEta_;
vector<float>  genJetPhi_;
//vector<int>  genJetPartonID_;


void Analyzer::branchesGenJetPart(TTree* tree) {

  tree->Branch("nGenJets",     &nGenJets_);
  tree->Branch("genJetEn",       &genJetEn_);
  tree->Branch("genJetPt",       &genJetPt_);
  tree->Branch("genJetEta",      &genJetEta_);
  tree->Branch("genJetPhi",      &genJetPhi_);
  //tree->Branch("genJetPartonID",      &genJetPartonID_);
  
}

void Analyzer::fillGenJetInfo(const edm::Event& e) {

  nGenJets_ = -99;
  genJetEn_.clear();
  genJetPt_.clear();
  genJetEta_.clear();
  genJetPhi_.clear();
  //genJetPartonID_.clear();
  
  if(doGenJets_){

    // Get GenJets
    edm::Handle<std::vector<reco::GenJet>> genjets;
    e.getByToken(GenJetLabel_,genjets);
    if(!genjets.isValid()) {
      edm::LogWarning("Analyzer") << "Could not find GenJet vector named ";


      return;
      }

    
    nGenJets_= genjets->size();

    for (vector<reco::GenJet>::const_iterator ip = genjets->begin(); ip != genjets->end(); ++ip) {
      genJetEn_.push_back(ip->energy());
      genJetPt_.push_back(ip->pt());
      genJetEta_.push_back(ip->eta());
      genJetPhi_.push_back(ip->phi());
      //genJetPartonID_.push_back(ip->genParton()->pdgId());
    }
      
  }//if(doGenJets_)





  //////////////////////////////////////////////////////
  // cleanup from previous execution

}


