#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/EgammaReco/interface/HFEMClusterShapeAssociation.h"
#include "DataFormats/EgammaReco/interface/HFEMClusterShape.h"

#include "ggAnalysis/ggNtuplizer/interface/ggNtuplizer.h"

using namespace std;
using namespace reco;

// (local) variables associated with tree branches
Int_t          nHFEle_;
vector<float>  hfeleEn_;
vector<float>  hfelePt_;
vector<float>  hfeleEta_;
vector<float>  hfelePhi_;
vector<float>  hfeleE1x1_;
vector<float>  hfeleE3x3_;
vector<float>  hfeleE5x5_;
vector<float>  hfeleELong1x1_;
vector<float>  hfeleELong3x3_;
vector<float>  hfeleELong5x5_;
vector<float>  hfeleEShort1x1_;
vector<float>  hfeleEShort3x3_;
vector<float>  hfeleEShort5x5_;
vector<float>  hfeleE9E25_;
vector<float>  hfeleECoreE9_;
vector<float>  hfeleECore_;
vector<float>  hfeleESeL_;

void ggNtuplizer::branchesHFElectrons(TTree* tree) {

  tree->Branch("nHFEle",         &nHFEle_);
  tree->Branch("hfeleEn",        &hfeleEn_);
  tree->Branch("hfelePt",        &hfelePt_);
  tree->Branch("hfeleEta",       &hfeleEta_);
  tree->Branch("hfelePhi",       &hfelePhi_);
  tree->Branch("hfeleE1x1",      &hfeleE1x1_);
  tree->Branch("hfeleE3x3",      &hfeleE3x3_);
  tree->Branch("hfeleE5x5",      &hfeleE5x5_);
  tree->Branch("hfeleELong1x1",  &hfeleELong1x1_);
  tree->Branch("hfeleELong3x3",  &hfeleELong3x3_);
  tree->Branch("hfeleELong5x5",  &hfeleELong5x5_);
  tree->Branch("hfeleEShort1x1", &hfeleEShort1x1_);
  tree->Branch("hfeleEShort3x3", &hfeleEShort3x3_);
  tree->Branch("hfeleEShort5x5", &hfeleEShort5x5_);
  tree->Branch("hfeleE9E25",     &hfeleE9E25_);
  tree->Branch("hfeleECoreE9",   &hfeleECoreE9_);
  tree->Branch("hfeleECore",     &hfeleECore_);
  tree->Branch("hfeleESeL",      &hfeleESeL_);
  
}

void ggNtuplizer::fillHFElectrons(const edm::Event &e) {
    
  // cleanup from previous execution
  hfeleEn_       .clear();
  hfelePt_       .clear();
  hfeleEta_      .clear();
  hfelePhi_      .clear();

  hfeleE1x1_     .clear();
  hfeleE3x3_     .clear();
  hfeleE5x5_     .clear();
  hfeleELong1x1_ .clear();
  hfeleELong3x3_ .clear();
  hfeleELong5x5_ .clear();
  hfeleEShort1x1_.clear();
  hfeleEShort3x3_.clear();
  hfeleEShort5x5_.clear();
  hfeleE9E25_    .clear();
  hfeleECoreE9_  .clear();
  hfeleECore_    .clear();
  hfeleESeL_     .clear();

  nHFEle_ = 0;

  edm::Handle<std::vector<reco::RecoEcalCandidate> > hfElectronHandle;
  e.getByToken(hfElectronCollection_, hfElectronHandle);

  if (!hfElectronHandle.isValid()) {
    edm::LogWarning("ggNtuplizer") << "no reco::RecoEcalCandidate in event";
    return;
  }

  edm::Handle<reco::HFEMClusterShapeAssociationCollection> hfClusterMapHandle;
  e.getByToken(hfClusterMapCollection_, hfClusterMapHandle);

  if (!hfClusterMapHandle.isValid()) {
    edm::LogWarning("ggNtuplizer") << "no reco::HFEMClusterShapeAssociationCollection in event";
    return;
  }
  
  for(unsigned int i=0; i!= hfElectronHandle->size(); ++i) {
    reco::RecoEcalCandidate hfele = hfElectronHandle->at(i);
    hfeleEn_ .push_back(hfele.energy());
    hfelePt_ .push_back(hfele.pt());
    hfeleEta_.push_back(hfele.eta());
    hfelePhi_.push_back(hfele.phi());
    
    edm::Ref<std::vector<reco::HFEMClusterShape> > hfCluster;
    hfCluster = hfClusterMapHandle->find(hfele.superCluster())->val;
    hfeleE1x1_.push_back(hfCluster->e1x1());
    hfeleE3x3_.push_back(hfCluster->e3x3());
    hfeleE5x5_.push_back(hfCluster->e5x5());
    hfeleELong1x1_.push_back(hfCluster->eLong1x1());
    hfeleELong3x3_.push_back(hfCluster->eLong3x3());
    hfeleELong5x5_.push_back(hfCluster->eLong5x5());
    hfeleEShort1x1_.push_back(hfCluster->eShort1x1());
    hfeleEShort3x3_.push_back(hfCluster->eShort3x3());
    hfeleEShort5x5_.push_back(hfCluster->eShort5x5());
    hfeleE9E25_    .push_back(hfCluster->e9e25());
    hfeleECoreE9_  .push_back(hfCluster->eCOREe9());
    hfeleECore_    .push_back(hfCluster->eCore());
    hfeleESeL_     .push_back(hfCluster->eSeL());

    ++nHFEle_;
  }

}
