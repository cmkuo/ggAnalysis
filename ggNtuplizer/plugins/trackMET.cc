#define MSDEBUG 0
#include "ggAnalysis/ggNtuplizer/interface/trackMET.h"
#include "DataFormats/Math/interface/deltaR.h"

trackMET::trackMET(const edm::ParameterSet& config) {
}

void trackMET::configure(edm::Handle<reco::VertexCollection> vtx, edm::Handle<reco::PFCandidateCollection> pf) {

  vtxHandle = vtx;
  pfHandle  = pf;

}

std::vector<float> trackMET::trackMETWithPV(float ptMin, float dzMax, float dxyMax, int ivtx) {
  
  std::vector<float> result;
  const reco::PFCandidateCollection* chargedPFParticlesForMET = pfHandle.product();

  // Calculate track MET for PV
  reco::VertexRef vtx(vtxHandle, ivtx);
  
  float px     = 0;
  float py     = 0;
  float metx   = 0;
  float mety   = 0;
  float metphi = 0;
  float met    = 0;
  
  // Loop over the PFCandidates
  for (unsigned i=0; i<chargedPFParticlesForMET->size(); ++i) {
    
    const reco::PFCandidate& pfcand = (*chargedPFParticlesForMET)[i];
    
    // require that PFCandidate is a charged particle
    if (pfcand.trackRef().isNonnull()) {
      
      if (pfcand.pt() < ptMin) continue;
      
      float dz = fabs(pfcand.trackRef()->dz(vtx->position()));
      if (dz > dzMax) continue;
      
      float dxy = fabs(pfcand.trackRef()->dxy(vtx->position()));
      if (dxy > dxyMax) continue;
      
      px += pfcand.px();
      py += pfcand.py();
    }
  }
  
  metx   = -px;
  mety   = -py;
  metphi = atan2(-py, -px);
  met    = sqrt(px*px + py*py);
  
  result.push_back(metx);
  result.push_back(mety);
  result.push_back(metphi);
  result.push_back(met);

  return result;
}

std::vector<float> trackMET::trackMETWithVertex(float ptMin, float dzMax, float dxyMax) {
  
  std::vector<float> result;
  const reco::PFCandidateCollection* chargedPFParticlesForMET = pfHandle.product();

  // Calculate track MET separately for each vertex
  for (unsigned int ivtx=0; ivtx<vtxHandle->size(); ++ivtx) {
    
    reco::VertexRef vtx(vtxHandle, ivtx);
    
    float px     = 0;
    float py     = 0;
    float metx   = 0;
    float mety   = 0;
    float metphi = 0;
    float met    = 0;

    // Loop over the PFCandidates
    for (unsigned i=0; i<chargedPFParticlesForMET->size(); ++i) {
     
      const reco::PFCandidate& pfcand = (*chargedPFParticlesForMET)[i];

      // require that PFCandidate is a charged particle
      if (pfcand.trackRef().isNonnull()) {
	
	if (pfcand.pt() < ptMin) continue;
	
	float dz = fabs(pfcand.trackRef()->dz(vtx->position()));
	if (dz > dzMax) continue;
	
	float dxy = fabs(pfcand.trackRef()->dxy(vtx->position()));
	if (dxy > dxyMax) continue;
	
	px += pfcand.px();
	py += pfcand.py();
      }
    }

    metx   = -px;
    mety   = -py;
    metphi = atan2(-py, -px);
    met    = sqrt(px*px + py*py);

    result.push_back(metx);
    result.push_back(mety);
    result.push_back(metphi);
    result.push_back(met);
  }
  
  return result;
}

