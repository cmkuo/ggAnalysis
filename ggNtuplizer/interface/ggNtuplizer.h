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
#include "ggAnalysis/ggNtuplizer/interface/GEDPhoIDTools.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "TTree.h"
#include "TH1F.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"
#include <memory>
#include <fstream>
#include <map>

#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

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
  
  void makeBranchesGenParticles(TTree* tree);
  void fillGenParticles(const edm::Event& e);

  void clearVectors();
  
  void getHandles(const edm::Event & event,
		  edm::Handle<VertexCollection>                & recVtxs,
		  edm::Handle<VertexCollection>                & recVtxsBS,
		  edm::Handle<double>                          & rhoHandle,
		  edm::Handle<View<pat::MET> >                 & pfMETHandle_,
		  edm::Handle<edm::View<pat::Electron> >       & electronHandle,
		  edm::Handle<edm::View<pat::Photon> >         & photonHandle,
		  edm::Handle<edm::View<pat::Muon> >           & muonHandle,
		  edm::Handle<vector<pat::Tau> >               & tauHandle,  //lvdp
		  edm::Handle<EcalRecHitCollection>            & EBReducedRecHits,
		  edm::Handle<EcalRecHitCollection>            & EEReducedRecHits,
		  edm::Handle<EcalRecHitCollection>            & ESRecHits,
		  edm::Handle<reco::PhotonCollection>          & recoPhotonHandle,
		  edm::Handle<TrackCollection>                 & tracksHandle,
		  edm::Handle<GsfElectronCollection>           & gsfElectronHandle,
		  edm::Handle<PFCandidateCollection>           & pfAllCandidates,
                  edm::Handle<edm::View<pat::Jet> >            & jetHandle   //lvdp		  
		  );

  Bool_t   doGenParticles_;
  Bool_t   runOnParticleGun_;  
  vector<int> newparticles_;

  InputTag vtxLabel_;  
  InputTag vtxBSLabel_;
  InputTag rhoLabel_;
  InputTag generatorLabel_;
  InputTag puCollection_;
  edm::EDGetTokenT<vector<reco::GenParticle> > genParticlesCollection_;
  InputTag pfMETlabel_;
  InputTag electronCollection_;
  InputTag photonCollection_;
  InputTag muonCollection_;
  InputTag tauCollection_;
  edm::EDGetTokenT<EcalRecHitCollection> ebReducedRecHitCollection_;
  edm::EDGetTokenT<EcalRecHitCollection> eeReducedRecHitCollection_;
  edm::EDGetTokenT<EcalRecHitCollection> esReducedRecHitCollection_;
  InputTag recophotonCollection_;
  InputTag tracklabel_;
  InputTag gsfElectronlabel_;
  InputTag pfAllParticles_;
  InputTag jetsCHSLabel_;
  InputTag jetCollection_;


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
  vector<float>  eleTheta_;
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
  vector<float>  eleIDMVA_;

  vector<float>    eledEtaseedAtVtx_;
  vector<float>    eleE1x5_;
  vector<float>    eleE2x5_;
  vector<float>    eleE5x5_;

  vector<float>    eleE1x5_2012_;
  vector<float>    eleE2x5_2012_;
  vector<float>    eleE5x5_2012_;
  vector<float>    eleRelIsoWithDBeta_;

  vector<int>      eleEcalDrivenSeed_;
  vector<float>    eleDr03EcalRecHitSumEt_;
  vector<float>    eleDr03HcalDepth1TowerSumEt_;
  vector<float>    eleDr03HcalDepth2TowerSumEt_;
  vector<float>    eleDr03HcalTowerSumEt_;
  vector<float>    eleDr03TkSumPt_;
 
  vector<float>    elecaloEnergy_;
  vector<float>    eleTrkdxy_;
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

  vector<float>  phoIDMVA_;

  vector<float>  phoEcalRecHitSumEtConeDR03_;
  vector<float>  phohcalDepth1TowerSumEtConeDR03_;
  vector<float>  phohcalDepth2TowerSumEtConeDR03_;
  vector<float>  phohcalTowerSumEtConeDR03_;
  vector<float>  photrkSumPtHollowConeDR03_;
  
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

  ///SJ
  vector<float>  muInnervalidFraction_;
  vector<float>  musegmentCompatibility_;
  vector<float>  muchi2LocalPosition_;
  vector<float>  mutrkKink_;
  vector<float>  muBestTrkPtError_;
  vector<float>  muBestTrkPt_;

  //SubJet
  Int_t nCA8Jet_;
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
  vector<int>  CA8JetNCH_;
  vector<int>  CA8Jetnconstituents_;
  vector<float>  CA8prunedJetMass_;

  //jets
  Int_t    nJet_;
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
  vector<int> jetPartonID_;
  vector<bool> jetPFLooseId_;


//Taus:  Lvdp
 
  Int_t nTau_;
  // decay mode discriminators
  std::vector<bool>   tauByLooseElectronRejection_;
  std::vector<bool>   tauByMediumElectronRejection_;
  std::vector<bool>   tauByTightElectronRejection_;
  std::vector<bool>   tauByMVA5LooseElectronRejection_;
  std::vector<bool>   tauByMVA5MediumElectronRejection_;
  std::vector<bool>   tauByMVA5TightElectronRejection_;
  std::vector<bool>   tauByMVA5VTightElectronRejection_;
  std::vector<bool>   tauByLooseMuonRejection_;
  std::vector<bool>   tauByMediumMuonRejection_;
  std::vector<bool>   tauByTightMuonRejection_;
  std::vector<bool>   tauByLooseMuonRejection3_;
  std::vector<bool>   tauByTightMuonRejection3_;
  std::vector<bool>   tauByMVALooseMuonRejection_;
  std::vector<bool>   tauByMVAMediumMuonRejection_;
  std::vector<bool>   tauByMVATightMuonRejection_;
  std::vector<bool>   tauByMVArawMuonRejection_;
  std::vector<bool>   pfTausDiscriminationByDecayModeFinding_;
  std::vector<bool>   tauByVLooseIsolation_;
  std::vector<bool>   tauByVLooseCombinedIsolationDBSumPtCorr_;
  std::vector<bool>   tauByLooseCombinedIsolationDBSumPtCorr_;
  std::vector<bool>   tauByMediumCombinedIsolationDBSumPtCorr_;
  std::vector<bool>   tauByTightCombinedIsolationDBSumPtCorr_;
  std::vector<bool>   tauByLooseCombinedIsolationDBSumPtCorr3Hits_;
  std::vector<bool>   tauByMediumCombinedIsolationDBSumPtCorr3Hits_;
  std::vector<bool>   tauByTightCombinedIsolationDBSumPtCorr3Hits_;
  std::vector<bool>   tauByVLooseIsolationMVA3newDMwoLT_;
  std::vector<bool>   tauByLooseIsolationMVA3newDMwoLT_;
  std::vector<bool>   tauByMediumIsolationMVA3newDMwoLT_;
  std::vector<bool>   tauByTightIsolationMVA3newDMwoLT_;
  std::vector<bool>   tauByVTightIsolationMVA3newDMwoLT_;
  std::vector<bool>   tauByVVTightIsolationMVA3newDMwoLT_;
  std::vector<bool>   tauByIsolationMVA3newDMwoLTraw_;
  std::vector<bool>   tauByVLooseIsolationMVA3oldDMwLT_;
  std::vector<bool>   tauByLooseIsolationMVA3oldDMwLT_;
  std::vector<bool>   tauByMediumIsolationMVA3oldDMwLT_;
  std::vector<bool>   tauByTightIsolationMVA3oldDMwLT_;
  std::vector<bool>   tauByVTightIsolationMVA3oldDMwLT_;
  std::vector<bool>   tauByVVTightIsolationMVA3oldDMwLT_;
  std::vector<bool>   tauByIsolationMVA3oldDMwLTraw_;
  std::vector<bool>   tauByVLooseIsolationMVA3oldDMwoLT_;
  std::vector<bool>   tauByLooseIsolationMVA3oldDMwoLT_;
  std::vector<bool>   tauByTightIsolationMVA3oldDMwoLT_;
  std::vector<bool>   tauByVTightIsolationMVA3oldDMwoLT_;
  std::vector<bool>   tauByVVTightIsolationMVA3oldDMwoLT_;
  std::vector<bool>   tauByIsolationMVA3oldDMwoLTraw_;
  std::vector<bool>   tauByLooseIsolationMVA3newDMwLT_;
  std::vector<bool>   tauByVLooseIsolationMVA3newDMwLT_;
  std::vector<bool>   tauByMediumIsolationMVA3newDMwLT_;
  std::vector<bool>   tauByTightIsolationMVA3newDMwLT_;
  std::vector<bool>   tauByVTightIsolationMVA3newDMwLT_;
  std::vector<bool>   tauByVVTightIsolationMVA3newDMwLT_;
  std::vector<bool>   tauByIsolationMVA3newDMwLTraw_;



  // isolation
  vector<bool> tauCombinedIsolationDeltaBetaCorrRaw3Hits_;
  vector<bool> tauLooseCombinedIsolationDeltaBetaCorr3Hits_;
  vector<bool> tauMediumCombinedIsolationDeltaBetaCorr3Hits_;
  vector<bool> tauTightCombinedIsolationDeltaBetaCorr3Hits_;

  // kinematics
  vector<float> tauEta_;
  vector<float> tauPhi_;
  vector<float> tauPt_;
  vector<float> tauEt_;
  vector<float> tauCharge_;
  vector<int>   tauDecayMode_;
  vector<float> tauEMFraction_;
  vector<float> tauHCAL3x3OverPLead_;
  vector<float> tauHCALMaxOverPLead_;
  vector<float> tauHCALTotOverPLead_;
  vector<float> tauIsolationPFChargedHadrCandsPtSum_;
  vector<float> tauIsolationPFGammaCandsEtSum_;
  vector<float> tauLeadPFChargedHadrCandsignedSipt_;
  vector<bool>  tauLeadChargedHadronExists_;
  vector<float> tauLeadChargedHadronEta_;
  vector<float> tauLeadChargedHadronPhi_;
  vector<float> tauLeadChargedHadronPt_;


  // Physics objects handles
  Handle<VertexCollection>                recVtxs_;
  Handle<VertexCollection>                recVtxsBS_;
  Handle<double>                          rhoHandle_;
  Handle<View<pat::MET> >                 pfMETHandle_;
  Handle<View<pat::Electron> >            electronHandle_;
  Handle<View<pat::Photon> >              photonHandle_;
  Handle<View<pat::Muon> >                muonHandle_;
  Handle<vector<pat::Tau> >               tauHandle_;  //lvdp
  Handle<EcalRecHitCollection>            EBReducedRecHits_;
  Handle<EcalRecHitCollection>            EEReducedRecHits_;
  Handle<EcalRecHitCollection>            ESReducedRecHits_;
  Handle<reco::PhotonCollection>          recoPhotonHandle_;
  Handle<TrackCollection>                 tracksHandle_;
  Handle<GsfElectronCollection>           gsfElectronHandle_;
  Handle<PFCandidateCollection>           pfAllCandidates_; 
  Handle<View<pat::Jet> >       jetHandle_;
  Bool_t dumpSubJets_;
  Bool_t dumpJets_;
  Bool_t dumpTaus_;
  PFJetIDSelectionFunctor pfLooseId_;


  // Variables that will be containers on which TMVA Reader works
  // The variables
  float varPhi_;
  float varR9_; 
  float varSieie_;
  float varSieip_; 
  float varE1x3overE5x5_; 
  float varE2x2overE5x5_; 
  float varE2x5overE5x5_; 
  float varSCEta_; 
  float varRawE_; 
  float varSCEtaWidth_; 
  float varSCPhiWidth_; 
  float varRho_;
  float varPhoIsoRaw_;
  float varChIsoRaw_; 
  float varWorstChRaw_;
  float varESEnOverRawE_; // for endcap MVA only
  float varESEffSigmaRR_; // for endcap MVA only
  // The spectators
  float varPt_; 
  float varEta_;

  // TMVA Reader for applying MVA
  TMVA::Reader *tmvaReader_[2];
  TString methodName_[2];


};

#endif
