#ifndef TRACKMET_H
#define TRACKMET_H

#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TLorentzVector.h"

class trackMET {

 public:
  
  explicit trackMET(const edm::ParameterSet&);
  ~trackMET() {};

  void configure(edm::Handle<reco::VertexCollection>, edm::Handle<reco::PFCandidateCollection>);
  std::vector<float> trackMETWithPV(float ptMin, float dzMax, float dxyMax, int ivtx = 0);
  std::vector<float> trackMETWithVertex(float ptMin, float dzMax, float dxyMax);
  
 private:

  edm::Handle<reco::VertexCollection> vtxHandle;
  edm::Handle<reco::PFCandidateCollection> pfHandle;

};

#endif
