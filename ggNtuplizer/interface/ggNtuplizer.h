#ifndef ggNtuplizer_h
#define ggNtuplizer_h

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Particle.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/CiCPhotonID.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "RecoEgamma/PhotonIdentification/interface/GEDPhoIDTools.h"

#include "TTree.h"
#include "TH1F.h"

#include <memory>
#include <fstream>
#include <map>

using namespace std;
using namespace edm;
using namespace reco;
using namespace pat;

class ggNtuplizer : public edm::EDAnalyzer {

   public:

  explicit ggNtuplizer(const edm::ParameterSet&);
  ~ggNtuplizer();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
 private:
  
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  void clearVectors();
  
  void getHandles(const edm::Event & event,
		  edm::Handle<std::vector<reco::GenParticle> > & genParticles,
		  edm::Handle<VertexCollection>                & recVtxs,
		  edm::Handle<VertexCollection>                & recVtxsBS,
		  edm::Handle<double>                          & rhoHandle,
		  edm::Handle<View<pat::MET> >                 & pfMETHandle_,
		  edm::Handle<edm::View<pat::Electron> >       & electronHandle,
		  edm::Handle<edm::View<pat::Photon> >         & photonHandle,
		  edm::Handle<edm::View<pat::Muon> >           & muonHandle,
		  edm::Handle<EcalRecHitCollection>            & EBReducedRecHits,
		  edm::Handle<EcalRecHitCollection>            & EEReducedRecHits,
		  edm::Handle<EcalRecHitCollection>            & ESRecHits,
		  edm::Handle<reco::PhotonCollection>          & recoPhotonHandle,
		  edm::Handle<TrackCollection>                 & tracksHandle,
		  edm::Handle<GsfElectronCollection>           & gsfElectronHandle,
		  edm::Handle<PFCandidateCollection>           & pfAllCandidates
		  );

  float getGenCalIso(edm::Handle<reco::GenParticleCollection> handle, reco::GenParticleCollection::const_iterator thisPho, 
		     const Float_t dRMax=0.4, Bool_t removeMu=true, Bool_t removeNu=false);
  float getGenTrkIso(edm::Handle<reco::GenParticleCollection> handle, reco::GenParticleCollection::const_iterator thisPho, const Float_t dRMax=0.4);
  
  Bool_t   doGenParticles_;
  Bool_t   runOnParticleGun_;  

  InputTag vtxLabel_;  
  InputTag vtxBSLabel_;
  InputTag rhoLabel_;
  InputTag generatorLabel_;
  InputTag puCollection_;
  InputTag genParticlesCollection_;
  InputTag pfMETlabel_;
  InputTag electronCollection_;
  InputTag photonCollection_;
  InputTag muonCollection_;
  InputTag ebReducedRecHitCollection_;
  InputTag eeReducedRecHitCollection_;
  InputTag esReducedRecHitCollection_;
  InputTag recophotonCollection_;
  InputTag tracklabel_;
  InputTag gsfElectronlabel_;
  InputTag pfAllParticles_;

  TTree   *tree_;
  TH1F    *hEvents_; 
  TH1F    *hPU_;
  TH1F    *hPUTrue_;

  EcalClusterLazyTools       *lazyTool;
  noZS::EcalClusterLazyTools *lazyToolnoZS;
  CiCPhotonID          *cicPhotonId_;

  Int_t          run_;
  Long64_t       event_;
  Int_t          lumis_;
  Bool_t         isData_;
  vector<float>  pdf_;
  Float_t        pthat_;
  Float_t        processID_;
  Int_t          nVtx_;
  Int_t          nTrks_;
  Float_t        rho_;
  // genParticle
  Int_t          nMC_;
  vector<int>    mcPID;
  vector<float>  mcVtx_x;
  vector<float>  mcVtx_y;
  vector<float>  mcVtx_z;
  vector<float>  mcPt;
  vector<float>  mcMass;
  vector<float>  mcEta;
  vector<float>  mcPhi;
  vector<float>  mcE;
  vector<float>  mcEt;
  vector<int>    mcGMomPID;
  vector<int>    mcMomPID;
  vector<float>  mcMomPt;
  vector<float>  mcMomMass;
  vector<float>  mcMomEta;
  vector<float>  mcMomPhi;
  vector<int>    mcIndex;
  vector<int>    mcDecayType;
  vector<int>    mcParentage;
  vector<int>    mcStatus;
  vector<float>  mcCalIsoDR03;
  vector<float>  mcTrkIsoDR03;
  vector<float>  mcCalIsoDR04;
  vector<float>  mcTrkIsoDR04;
  // PU
  Int_t          nPUInfo_;  
  vector<int>    nPU_;
  vector<int>    puBX_;
  vector<float>  puTrue_;
  // pfMET Type1 
  Float_t        genMET_;
  Float_t        genMETPhi_;
  Float_t        pfMET_;
  Float_t        pfMETPhi_;
  Float_t        pfMETsumEt_;
  Float_t        pfMETmEtSig_;
  Float_t        pfMETSig_;
  // Electron
  Int_t          nEle_;
  vector<int>    eleCharge_;
  vector<int>    eleChargeConsistent_;
  vector<float>  eleEn_;
  vector<float>  eleSCEn_;
  vector<float>  eleESEn_;
  vector<float>  eleD0_;
  vector<float>  eleDz_;
  vector<float>  elePt_;
  vector<float>  eleEta_;
  vector<float>  elePhi_;
  vector<float>  eleSCEta_;
  vector<float>  eleSCPhi_;
  vector<float>  eleSCRawEn_;
  vector<float>  eleSCEtaWidth_;
  vector<float>  eleSCPhiWidth_;
  vector<float>  eleHoverE_;
  vector<float>  eleEoverP_;
  vector<float>  eleEoverPInv_;
  vector<float>  eleBrem_;
  vector<float>  eledEtaAtVtx_;
  vector<float>  eledPhiAtVtx_;
  vector<float>  eleSigmaIEtaIEta_;
  vector<float>  eleSigmaIEtaIEta_2012_;
  vector<float>  eleSigmaIEtaIPhi_;
  vector<float>  eleSigmaIPhiIPhi_;
  vector<int>    eleConvVeto_;
  vector<int>    eleMissHits_;
  vector<float>  eleESEffSigmaRR_;
  vector<float>  elePFChIso_;
  vector<float>  elePFPhoIso_;
  vector<float>  elePFNeuIso_;
  vector<float>  elePFPUIso_;
  vector<float>  eleBC1E_;
  vector<float>  eleBC1Eta_;
  vector<float>  eleBC2E_;
  vector<float>  eleBC2Eta_;
  // Photon
  Int_t          nPho_;      
  vector<float>  phoE_;
  vector<float>  phoEt_;
  vector<float>  phoEta_;
  vector<float>  phoPhi_;
  vector<float>  phoSCE_;
  vector<float>  phoSCRawE_;
  vector<float>  phoESEn_;
  vector<float>  phoSCEta_;
  vector<float>  phoSCPhi_;
  vector<float>  phoSCEtaWidth_;
  vector<float>  phoSCPhiWidth_;
  vector<float>  phoSCBrem_;
  vector<int>    phohasPixelSeed_;
  vector<int>    phoEleVeto_;
  vector<float>  phoR9_;
  vector<float>  phoHoverE_;
  vector<float>  phoSigmaIEtaIEta_;
  vector<float>  phoSigmaIEtaIPhi_;
  vector<float>  phoSigmaIPhiIPhi_;
  vector<float>  phoE1x3_;
  vector<float>  phoE2x2_;
  vector<float>  phoE2x5Max_;
  vector<float>  phoE5x5_;
  vector<float>  phoESEffSigmaRR_;
  vector<float>  phoSigmaIEtaIEta_2012_;
  vector<float>  phoSigmaIEtaIPhi_2012_;
  vector<float>  phoSigmaIPhiIPhi_2012_;
  vector<float>  phoE1x3_2012_;
  vector<float>  phoE2x2_2012_;
  vector<float>  phoE2x5Max_2012_;
  vector<float>  phoE5x5_2012_;
  vector<float>  phoPFChIso_;
  vector<float>  phoPFPhoIso_;
  vector<float>  phoPFNeuIso_;
  vector<float>  phoPFChWorstIso_;
  vector<float>  phoPFChIsoFrix1_;
  vector<float>  phoPFChIsoFrix2_;
  vector<float>  phoPFChIsoFrix3_;
  vector<float>  phoPFChIsoFrix4_;
  vector<float>  phoPFChIsoFrix5_;
  vector<float>  phoPFChIsoFrix6_;
  vector<float>  phoPFChIsoFrix7_;
  vector<float>  phoPFChIsoFrix8_;
  vector<float>  phoPFPhoIsoFrix1_;
  vector<float>  phoPFPhoIsoFrix2_;
  vector<float>  phoPFPhoIsoFrix3_;
  vector<float>  phoPFPhoIsoFrix4_;
  vector<float>  phoPFPhoIsoFrix5_;
  vector<float>  phoPFPhoIsoFrix6_;
  vector<float>  phoPFPhoIsoFrix7_;
  vector<float>  phoPFPhoIsoFrix8_;
  vector<float>  phoPFNeuIsoFrix1_;
  vector<float>  phoPFNeuIsoFrix2_;
  vector<float>  phoPFNeuIsoFrix3_;
  vector<float>  phoPFNeuIsoFrix4_;
  vector<float>  phoPFNeuIsoFrix5_;
  vector<float>  phoPFNeuIsoFrix6_;
  vector<float>  phoPFNeuIsoFrix7_;
  vector<float>  phoPFNeuIsoFrix8_;
  vector<float>  phoBC1E_;
  vector<float>  phoBC1Eta_;
  vector<float>  phoBC2E_;
  vector<float>  phoBC2Eta_;
  // muon
  Int_t          nMu_;
  vector<float>  muPt_;
  vector<float>  muEta_;
  vector<float>  muPhi_;
  vector<int>    muCharge_;
  vector<int>    muType_;
  vector<int>    muIsGood_;
  //vector<int>    muID_;
  vector<float>  muD0_;
  vector<float>  muDz_;
  vector<float>  muChi2NDF_;
  vector<float>  muInnerD0_;
  vector<float>  muInnerDz_;
  vector<int>    muTrkLayers_; 
  vector<int>    muPixelLayers_;
  vector<int>    muPixelHits_;
  vector<int>    muMuonHits_;
  vector<int>    muStations_;
  vector<int>    muTrkQuality_;
  vector<float>  muIsoTrk_;
  vector<float>  muPFChIso_;
  vector<float>  muPFPhoIso_;
  vector<float>  muPFNeuIso_;
  vector<float>  muPFPUIso_;

  // Physics objects handles
  Handle<std::vector<reco::GenParticle> > genParticlesHandle_;
  Handle<VertexCollection>                recVtxs_;
  Handle<VertexCollection>                recVtxsBS_;
  Handle<double>                          rhoHandle_;
  Handle<View<pat::MET> >                 pfMETHandle_;
  Handle<View<pat::Electron> >            electronHandle_;
  Handle<View<pat::Photon> >              photonHandle_;
  Handle<View<pat::Muon> >                muonHandle_;
  Handle<EcalRecHitCollection>            EBReducedRecHits_;
  Handle<EcalRecHitCollection>            EEReducedRecHits_;
  Handle<EcalRecHitCollection>            ESReducedRecHits_;
  Handle<reco::PhotonCollection>          recoPhotonHandle_;
  Handle<TrackCollection>                 tracksHandle_;
  Handle<GsfElectronCollection>           gsfElectronHandle_;
  Handle<PFCandidateCollection>           pfAllCandidates_; 

};

#endif
