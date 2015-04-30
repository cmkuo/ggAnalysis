#ifndef ggNtuplizer_h
#define ggNtuplizer_h

#include "TTree.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/CiCPhotonID.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"

using namespace std;

class ggNtuplizer : public edm::EDAnalyzer {
 public:

  explicit ggNtuplizer(const edm::ParameterSet&);
  ~ggNtuplizer();
  
//   static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
 private:
  
//   virtual void beginJob() {};
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
//   virtual void endJob() {};

  void branchesGlobalEvent(TTree*);
  void branchesGenInfo    (TTree*, edm::Service<TFileService>&);
  void branchesGenPart    (TTree*);
  void branchesPhotons    (TTree*);
  void branchesElectrons  (TTree*);
  void branchesMuons      (TTree*);
  void branchesTaus       (TTree*);
  void branchesJets       (TTree*);

  void fillGlobalEvent(const edm::Event&);
  void fillGenInfo    (const edm::Event&);
  void fillGenPart    (const edm::Event&);
  void fillPhotons    (const edm::Event&, const edm::EventSetup&);
  void fillElectrons  (const edm::Event&, const edm::EventSetup&, math::XYZPoint&);
  void fillMuons      (const edm::Event&, math::XYZPoint&);
  void fillTaus       (const edm::Event&);
  void fillJets       (const edm::Event&);

  void cleanupPhotons();
  
  Bool_t doGenParticles_;
  Bool_t runOnParticleGun_;
  Bool_t dumpTaus_;
  Bool_t dumpJets_;
  Bool_t dumpSubJets_;

  vector<int> newparticles_;

  edm::EDGetTokenT<reco::VertexCollection>      vtxLabel_;
  edm::EDGetTokenT<reco::VertexCollection>      vtxBSLabel_;
  edm::EDGetTokenT<double>                      rhoLabel_;
  edm::EDGetTokenT<GenEventInfoProduct>         generatorLabel_;
  edm::EDGetTokenT<vector<PileupSummaryInfo> >  puCollection_;
  edm::EDGetTokenT<vector<reco::GenParticle> >  genParticlesCollection_;
  edm::EDGetTokenT<edm::View<pat::MET> >        pfMETlabel_;
  edm::EDGetTokenT<edm::View<pat::Electron> >   electronCollection_;
  edm::EDGetTokenT<edm::View<pat::Photon> >     photonCollection_;
  edm::EDGetTokenT<edm::View<pat::Muon> >       muonCollection_;
  edm::EDGetTokenT<vector<pat::Tau> >           tauCollection_;
  edm::EDGetTokenT<EcalRecHitCollection>        ebReducedRecHitCollection_;
  edm::EDGetTokenT<EcalRecHitCollection>        eeReducedRecHitCollection_;
  //edm::EDGetTokenT<EcalRecHitCollection> esReducedRecHitCollection_; // FIXME: not used anymore
  edm::EDGetTokenT<reco::PhotonCollection>      recophotonCollection_;
  edm::EDGetTokenT<reco::TrackCollection>       tracklabel_;
  edm::EDGetTokenT<reco::GsfElectronCollection> gsfElectronlabel_;
  edm::EDGetTokenT<reco::PFCandidateCollection> pfAllParticles_;
  edm::EDGetTokenT<edm::View<pat::Jet> >        jetCollection_;
  edm::EDGetTokenT<edm::View<pat::Jet> >        jetsCHSLabel_;

  TTree   *tree_;
  TH1F    *hEvents_;

  CiCPhotonID    *cicPhotonId_;

  PFJetIDSelectionFunctor pfLooseId_;

};

#endif
