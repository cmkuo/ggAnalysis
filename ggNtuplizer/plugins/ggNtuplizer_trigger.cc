#include "map"
#include "ggAnalysis/ggNtuplizer/interface/ggNtuplizer.h"

using namespace std;

// local variables: per-filter per-electron/muon/photon arrays of matched trigger objects
// NOTE: number of elements in the arrays equals sizeof(Int_t)
vector<double> trgElePt[32], trgEleEta[32];
vector<double> trgPhoPt[32], trgPhoEta[32];
vector<double> trgMuPt[32],  trgMuEta[32];

void ggNtuplizer::initTriggerFilters(const edm::Event &e) {
  // Fills the arrays above.

  // cleanup from previous execution
  for (size_t i = 0; i < 32; i++) {
    trgElePt [i].clear();
    trgEleEta[i].clear();
    trgPhoPt [i].clear();
    trgPhoEta[i].clear();
    trgMuPt  [i].clear();
    trgMuEta [i].clear();
  }

  // AOD vs miniAOD
  edm::Handle<trigger::TriggerEvent> triggerHandle;
  edm::Handle<pat::TriggerObjectStandAloneCollection> triggerHandleMiniAOD;
  if (isAOD_)
    e.getByToken(trgEventLabel_, triggerHandle);
  else
    e.getByToken(triggerObjectsLabel_, triggerHandleMiniAOD);

  // filter => index (in trg*[] arrays) mappings
  std::map<string,size_t> eleFilters;
  std::map<string,size_t> phoFilters;
  std::map<string,size_t> muFilters;

  // FIXME: define actual filters
  // FIXME: using only the latest filter in a HLT is not sufficient
  eleFilters["hltSingleEle22WPLooseGsfTrackIsoFilter"] = 0;
  eleFilters["hltL1sL1SingleEG20ORL1SingleEG15"] = 1;
  eleFilters["hltEle25WP60SC4HcalIsoFilter"] = 2;

  if (isAOD_) {
    const trigger::TriggerObjectCollection& trgObjects = triggerHandle->getObjects();

    // loop over particular filters (and not over full HLTs)
    for (trigger::size_type iF = 0; iF != triggerHandle->sizeFilters(); ++iF) {
      // full filter name and its keys each corresponding to a matched (pt, eta, phi, ...) object
      string const&        label = triggerHandle->filterTag(iF).label();
      const trigger::Keys& keys  = triggerHandle->filterKeys(iF);

      std::map<string,size_t>::iterator idxEle = eleFilters.find(label);
      std::map<string,size_t>::iterator idxPho = phoFilters.find(label);
      std::map<string,size_t>::iterator idxMu  = muFilters.find(label);

      // electron filters
      if (idxEle != eleFilters.end()) {
        size_t idx = idxEle->second;

        for (size_t iK = 0; iK < keys.size(); ++iK) {
          const trigger::TriggerObject& trgV = trgObjects.at(keys[iK]);
          trgElePt [idx].push_back(trgV.pt());
          trgEleEta[idx].push_back(trgV.eta());
        }
      }

      // photon filters
      if (idxPho != phoFilters.end()) {
        size_t idx = idxPho->second;

        for (size_t iK = 0; iK < keys.size(); ++iK) {
          const trigger::TriggerObject& trgV = trgObjects.at(keys[iK]);
          trgPhoPt [idx].push_back(trgV.pt());
          trgPhoEta[idx].push_back(trgV.eta());
        }
      }

      // muon filters
      if (idxMu != muFilters.end()) {
        size_t idx = idxMu->second;

        for (size_t iK = 0; iK < keys.size(); ++iK) {
          const trigger::TriggerObject& trgV = trgObjects.at(keys[iK]);
          trgMuPt [idx].push_back(trgV.pt());
          trgMuEta[idx].push_back(trgV.eta());
        }
      }
    } // HLT filter loop

    return;
  } // if AOD

  // miniAOD
  // TODO
}

Int_t ggNtuplizer::matchElectronTriggerFilters(double pt, double eta) {

  // bits in the return value correspond to decisions from filters defined above
  Int_t result = 0;

  for (size_t f = 0; f < 32; ++f)
    for (size_t v = 0; v < trgElePt[f].size(); ++v)
      if (fabs(pt - trgElePt[f][v])/trgElePt[f][v] < trgFilterDeltaPtCut_ &&
          fabs(trgEleEta[f][v] - eta) < trgFilterDeltaEtaCut_) {
        result |= (1<<f);
        break;
      }

  return result;
}

Int_t ggNtuplizer::matchPhotonTriggerFilters(double pt, double eta) {

  // bits in the return value correspond to decisions from filters defined above
  Int_t result = 0;

  for (size_t f = 0; f < 32; ++f)
    for (size_t v = 0; v < trgPhoPt[f].size(); ++v)
      if (fabs(pt - trgPhoPt[f][v])/trgPhoPt[f][v] < trgFilterDeltaPtCut_ &&
          fabs(trgPhoEta[f][v] - eta) < trgFilterDeltaEtaCut_) {
        result |= (1<<f);
        break;
      }

  return result;
}

Int_t ggNtuplizer::matchMuonTriggerFilters(double pt, double eta) {

  // bits in the return value correspond to decisions from filters defined above
  Int_t result = 0;

  for (size_t f = 0; f < 32; ++f)
    for (size_t v = 0; v < trgMuPt[f].size(); ++v)
      if (fabs(pt - trgMuPt[f][v])/trgMuPt[f][v] < trgFilterDeltaPtCut_ &&
          fabs(trgMuEta[f][v] - eta) < trgFilterDeltaEtaCut_) {
        result |= (1<<f);
        break;
      }

  return result;
}
