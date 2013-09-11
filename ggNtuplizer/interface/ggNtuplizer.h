#ifndef ggNtuplizer_h
#define ggNtuplizer_h

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CommonTools/Utils/interface/TFileDirectory.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
#include "DataFormats/HeavyIonEvent/interface/CentralityProvider.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementSuperCluster.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/Particle.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "CMGTools/External/interface/PileupJetIdAlgo.h"
#include "CMGTools/External/interface/PileupJetIdentifier.h"
#include "EGamma/EGammaAnalysisTools/interface/EGammaMvaEleEstimator.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloTopology/interface/EcalPreshowerTopology.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "RecoEgamma/EgammaTools/interface/EGEnergyCorrector.h"
#include "RecoEgamma/EgammaElectronAlgos/interface/ElectronHcalHelper.h"
#include "HiggsAnalysis/HiggsTo2photons/interface/CiCPhotonID.h"
#include "PhysicsTools/PatUtils/interface/TriggerHelper.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "ggAnalysis/ggNtuplizer/interface/ggPFIsolation.h"
#include "ggAnalysis/ggNtuplizer/interface/trackMET.h"

#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TString.h"

#include <memory>
#include <fstream>
#include <map>

using namespace edm;
using namespace std;
using namespace reco;
using namespace pat;
using namespace pat::helper;

const Int_t maxP = 600;

class ggNtuplizer : public EDProducer {

public:

  explicit ggNtuplizer(const ParameterSet&);
  virtual ~ggNtuplizer();

protected:

  virtual void beginJob();
  virtual void produce( Event &, const EventSetup & );
  virtual void endJob();

  // Get handles
  void getHandles(edm::Event  & event,
		  edm::Handle<std::vector<reco::GenParticle> > & genParticles,
		  edm::Handle<VertexCollection> &            recVtxs,
		  edm::Handle<VertexCollection> &            recVtxsBS,
		  edm::Handle<TriggerResults> &              trgResultsHandle,
		  edm::Handle<TriggerEvent> &                triggerEvent,
		  edm::Handle<TrackCollection> &             tracksHandle,
		  edm::Handle<GsfElectronCollection> &       gsfElectronHandle,
		  edm::Handle<edm::View<pat::MET> > &        METHandle,
		  edm::Handle<View<pat::MET> > &             pfMETHandle_,
		  edm::Handle<View<pat::MET> > &             pfType01METHandle_,
		  edm::Handle<reco::PFMETCollection> &       recoPfMETHandle_,
		  edm::Handle<edm::View<pat::Electron> > &   electronHandle,
		  edm::Handle<edm::View<pat::Photon> > &     photonHandle,
		  edm::Handle<reco::PhotonCollection> &      recoPhotonHandle,
		  edm::Handle<edm::View<pat::Muon> > &       muonHandle,
		  edm::Handle<edm::View<pat::Jet> > &        jetHandle,
		  edm::Handle<BeamSpot> &                    beamSpotHandle,
		  edm::Handle<EcalRecHitCollection> &        EBReducedRecHits,
		  edm::Handle<EcalRecHitCollection> &        EEReducedRecHits,
		  edm::Handle<EcalRecHitCollection> &        ESRecHits,
		  edm::Handle<CaloTowerCollection> &         towers, 
                  edm::Handle<PFCandidateCollection>&        pfAllCandidates,
                  edm::Handle<PFCandidateCollection>&        pfCandidates,
	          edm::Handle<PFCandidateCollection>&        pfCandidatesPhotons,
                  edm::Handle<reco::ConversionCollection>&   convH,
		  edm::Handle<reco::PhotonCollection>& pfPhoTranslator
		  //  edm::Handle<reco::SuperClusterCollection>& pfPhoSClusters,
		  // edm::Handle<reco::SuperClusterCollection>& pfEleSClusters,
		  //edm::Handle<reco::PhotonCoreCollection>& pfPhotonCore,
		  //edm::Handle<reco::GsfElectronCoreCollection>& pfElectronCore
		  );

  Double_t eT(Double_t pt1, Double_t pt2) const;
  Double_t massT(Double_t pt1, Double_t pt2, Double_t wpx, Double_t wpy) const;
  Double_t acop(Double_t phi1, Double_t phi2) const;
  
  float getPhotonTrkIso(double photonVz, double vbsx, double vbsy, double vbsz, double photonEta, double photonPhi, edm::Handle<reco::TrackCollection> tracksHandle_, double outercone, double innercone, double etastrip, double dZCut, double dxy, int Option);
  float getPhotondRtoTrk(edm::Handle<reco::GsfElectronCollection> gsfHandle, const pat::Photon & pho, float minPt, int maxMissingHits);

  float getGenCalIso(edm::Handle<reco::GenParticleCollection> handle, reco::GenParticleCollection::const_iterator thisPho, 
		     const Float_t dRMax=0.4, Bool_t removeMu=true, Bool_t removeNu=false);
  float getGenTrkIso(edm::Handle<reco::GenParticleCollection> handle, reco::GenParticleCollection::const_iterator thisPho, const Float_t dRMax=0.4);
  void fillBlockInfo(reco::PFCandidate PFCand, bool doPhotons, bool reclassify, const CaloSubdetectorGeometry* geomBar, const CaloSubdetectorGeometry* geomEnd, const CaloSubdetectorGeometry* geometry_p); 
  void mergeBlockInfo(reco::PFCandidate PFCandPhoton, reco::PFCandidate PFCandElectron, const CaloSubdetectorGeometry* geomBar, const CaloSubdetectorGeometry* geomEnd, const CaloSubdetectorGeometry* geometry_p); 
  bool SLPoint(const reco::Photon phot, reco::ConversionRefVector  &SLconversions, std::vector<float> &Zint);
  double fsrPhotonIso03(const reco::PFCandidateCollection* pfAllCandidates, const reco::PFCandidateCollection* pfCandidates, const reco::PFCandidate pfPhoton);

  //std::vector<float> getESHits(const reco::CaloClusterPtr bcPtr, map<DetId, EcalRecHit> rechits_map, const CaloSubdetectorGeometry*& geometry_p, CaloSubdetectorTopology *topology_p, int row=0);
  std::vector<float> getESHits(double X, double Y, double Z, map<DetId, EcalRecHit> rechits_map, const CaloSubdetectorGeometry*& geometry_p, CaloSubdetectorTopology *topology_p, int row=0);
  vector<float> getESEffSigmaRR(vector<float> ESHits0);
  vector<float> getESEn(vector<float> ESHits0);

  map<DetId, EcalRecHit> rechits_map_;

  typedef std::vector< edm::Handle< edm::ValueMap<reco::IsoDeposit> > > IsoDepositMaps;
  typedef std::vector< edm::Handle< edm::ValueMap<double> > > IsoDepositVals;

  std::vector<edm::InputTag> inputTagPhotonIsoDeposits_;
  std::vector<edm::InputTag> inputTagIsoDepElectrons_;
  std::vector<edm::InputTag> inputTagIsoDepPhotons_;
  std::vector<edm::InputTag> inputTagIsoValElectronsPFId_;
  std::vector<edm::InputTag> inputTagIsoValPhotonsPFId_;

  // SC footprint remover parameters
  edm::ParameterSet scRemover03Pset_;
  edm::ParameterSet scRemover04Pset_;
  InputTag genParticlesCollection_;
  InputTag vtxlabel_;
  InputTag tracklabel_;
  InputTag gsfElectronlabel_;
  InputTag trgResults_;
  InputTag trgEvent_;
  InputTag METCollection_;
  InputTag pfMETlabel_;
  InputTag pfType01METlabel_;
  InputTag recoPfMETlabel_;
  InputTag electronCollection_;
  InputTag photonCollection_;
  InputTag recophotonCollection_;
  InputTag muonCollection_;
  InputTag jetCollection_; 
  InputTag ebReducedRecHitCollection_;
  InputTag eeReducedRecHitCollection_;
  InputTag towerCollection_;
  InputTag esRecHitCollection_;
  InputTag beamSpotCollection_;
  InputTag puCollection_;
 
  InputTag rho2011Label_;
  InputTag rhoLepPFisoCollection_;
  InputTag rhoCollection25_;
  InputTag rhoCollection25_neu_;
  InputTag rhoCollection44_;
  InputTag pfAllParticles_;
  InputTag pfParticles_;
  InputTag allConversionsColl_;
  InputTag pfPhotonParticles_;
  InputTag pfPhotonCollection_;
  std::vector<edm::ParameterSet> jetMVAAlgos_; 
  PFJetIDSelectionFunctor pfLooseId_;
  std::vector<edm::InputTag> inputTagMuonTags_; // muon cocktail
  Bool_t doGenParticles_;
  Bool_t doSkim_;
  Bool_t getBlocks_;
  Bool_t runOnParticleGun_;
  Bool_t dumpESHits_;
  Bool_t dumpESClusterInfo_;
  Bool_t dumpTrks_;
  Bool_t dumpJets_;
  Bool_t useAllPF_;
  Bool_t develop_;
  Bool_t doCentrality_;
  Int_t  verbosity_;

  HLTConfigProvider hltConfigProvider_;
  EcalClusterLazyTools *lazyTool;
  EGEnergyCorrector egCorrEle_;
  EGEnergyCorrector egCorrPho_;
  CiCPhotonID *cicPhotonId_;
  trackMET *trackMET_;
  ElectronHcalHelper *hcalHelper, *hcalHelperPflow;
  ElectronHcalHelper::Configuration hcalCfg, hcalCfgPflow;
  CentralityProvider *centProvider_;

  TTree *tree_;
  TH1F  *hEvents_; 
  TH1F  *hPU_;
  TH1F  *hPUTrue_;

  Int_t    run_;
  Long64_t event_;
  Int_t    lumis_;
  Bool_t   isData_;
  Float_t  pdf_[7];
  Float_t  pthat_;
  Float_t  processID_;
  Int_t    nHLT_;
  Int_t    HLT_[maxP];
  Int_t    HLTIndex_[70];
  Float_t  bspotPos_[3];
  Int_t    nVtx_;
  Float_t  vtx_[150][3];
  Int_t    vtxNTrk_[150];
  Float_t  vtxNDF_[150];
  Float_t  vtxD0_[150];
  Int_t    IsVtxGood_;
  Int_t    nGoodVtx_;
  Int_t    nVtxBS_;
  Float_t  vtxbs_[150][3];
  Float_t  vtxbsPtMod_[150];
  Float_t  vtxbsSumPt2_[150];
  std::vector<std::vector<Int_t> >*   vtxbsTkIndex_;
  std::vector<std::vector<Float_t> >* vtxbsTkWeight_;
  Int_t    nTrk_;
  Float_t  trkP_[30000][3];
  Float_t  trkVtx_[30000][3];
  Float_t  trkd0_[30000];
  Float_t  trkd0Err_[30000];
  Float_t  trkdz_[30000];
  Float_t  trkdzErr_[30000];
  Float_t  trkPtErr_[30000];
  Int_t    trkQuality_[30000];
  Int_t    nGoodTrk_;
  Int_t    IsTracksGood_;
  Float_t  centrality_[5];
  // genParticle
  Int_t    nMC_;
  Int_t    mcPID[maxP];
  Float_t  mcVtx[maxP][3];
  Float_t  mcPt[maxP];
  Float_t  mcMass[maxP];
  Float_t  mcEta[maxP];
  Float_t  mcPhi[maxP];
  Float_t  mcE[maxP];
  Float_t  mcEt[maxP];
  Int_t    mcGMomPID[maxP];
  Int_t    mcMomPID[maxP];
  Float_t  mcMomPt[maxP];
  Float_t  mcMomMass[maxP];
  Float_t  mcMomEta[maxP];
  Float_t  mcMomPhi[maxP];
  Int_t    mcIndex[maxP];
  Int_t    mcDecayType[maxP];
  Int_t    mcParentage[maxP];
  Int_t    mcStatus[maxP];
  // PU
  Int_t    nPUInfo_;  
  Int_t    nPU_[maxP];  
  Int_t    puBX_[maxP];  
  Float_t  puTrue_[maxP];
  // Gen & Reco MET
  Float_t  genMET_;
  Float_t  genMETPhi_;
  Float_t  MET_;
  Float_t  METx_;
  Float_t  METy_;
  Float_t  METPhi_;
  Float_t  METsumEt_;
  // pfMET Type1
  Float_t  pfMET_;
  Float_t  pfMETPhi_;
  Float_t  pfMETsumEt_;
  Float_t  pfMETmEtSig_;
  Float_t  pfMETSig_;
  // pfMET (Type 0+1)
  Float_t  pfType01MET_;
  Float_t  pfType01METPhi_;
  Float_t  pfType01METsumEt_;
  Float_t  pfType01METmEtSig_;
  Float_t  pfType01METSig_;
  // reco pf met
  Float_t  recoPfMET_;
  Float_t  recoPfMETPhi_;
  Float_t  recoPfMETsumEt_;
  Float_t  recoPfMETmEtSig_;
  Float_t  recoPfMETSig_;
  // track MET
  Float_t  trkMETxPV_;
  Float_t  trkMETyPV_;
  Float_t  trkMETPhiPV_;
  Float_t  trkMETPV_;
  Float_t  trkMETx_[maxP];
  Float_t  trkMETy_[maxP];
  Float_t  trkMETPhi_[maxP];
  Float_t  trkMET_[maxP];
  // MET filters
  Int_t    metFilters_[10];
  // Electron
  Int_t    nEle_;
  Int_t    eleTrg_[maxP][16];
  Int_t    eleClass_[maxP];
  Int_t    eleIsEcalDriven_[maxP];
  Int_t    eleCharge_[maxP];
  Int_t    eleChargeConsistent_[maxP];
  Float_t  eleEn_[maxP];
  Float_t  eleEcalEn_[maxP];
  Float_t  eleSCEn_[maxP];
  Float_t  eleESEn_[maxP];
  Float_t  eleVtx_[maxP][3];
  Float_t  eleD0_[maxP];
  Float_t  eleDz_[maxP];
  Float_t  eleD0GV_[maxP];
  Float_t  eleDzGV_[maxP];
  Float_t  eleD0Vtx_[maxP][100];
  Float_t  eleDzVtx_[maxP][100];
  Float_t  elePt_[maxP];
  Float_t  eleEta_[maxP];
  Float_t  elePhi_[maxP];
  Float_t  eleSCEta_[maxP];
  Float_t  eleSCPhi_[maxP];
  Float_t  eleEtVtx_[maxP][100];
  Float_t  eleEtaVtx_[maxP][100];
  Float_t  elePhiVtx_[maxP][100];
  Float_t  eleSCRawEn_[maxP];
  Float_t  eleSCEtaWidth_[maxP];
  Float_t  eleSCPhiWidth_[maxP];
  Float_t  eleHoverE_[maxP];
  Float_t  eleHoverE12_[maxP]; 
  Float_t  eleEoverP_[maxP];
  Float_t  elePin_[maxP];
  Float_t  elePout_[maxP];
  Float_t  eleTrkMomErr_[maxP];
  Float_t  eleBrem_[maxP];
  Float_t  eledEtaAtVtx_[maxP];
  Float_t  eledPhiAtVtx_[maxP];
  Float_t  eleSigmaIEtaIEta_[maxP];
  Float_t  eleSigmaIEtaIPhi_[maxP];
  Float_t  eleSigmaIPhiIPhi_[maxP];
  Float_t  eleEmax_[maxP];
  Float_t  eleE1x5_[maxP];
  Float_t  eleE3x3_[maxP];
  Float_t  eleE5x5_[maxP];
  Float_t  eleE2x5Max_[maxP];
  Float_t  eleRegrE_[maxP];
  Float_t  eleRegrEerr_[maxP];
  Float_t  elePhoRegrE_[maxP];
  Float_t  elePhoRegrEerr_[maxP];
  Float_t  eleSeedTime_[maxP];
  Int_t    eleRecoFlag_[maxP];
  Int_t    elePos_[maxP];
  Int_t    eleGenIndex_[maxP];
  Int_t    eleGenGMomPID_[maxP];
  Int_t    eleGenMomPID_[maxP];
  Float_t  eleGenMomPt_[maxP];
  Float_t  eleIsoTrkDR03_[maxP];
  Float_t  eleIsoEcalDR03_[maxP];
  Float_t  eleIsoHcalDR03_[maxP];
  Float_t  eleIsoHcalDR0312_[maxP];
  Float_t  eleIsoTrkDR04_[maxP];
  Float_t  eleIsoEcalDR04_[maxP];
  Float_t  eleIsoHcalDR04_[maxP];
  Float_t  eleIsoHcalDR0412_[maxP];
  Float_t  eleModIsoTrk_[maxP];
  Float_t  eleModIsoEcal_[maxP];
  Float_t  eleModIsoHcal_[maxP];
  Float_t  eleChi2NDF_[maxP];
  Int_t    eleMissHits_[maxP];
  Float_t  eleConvDist_[maxP];
  Float_t  eleConvDcot_[maxP];
  Int_t    eleConvVtxFit_[maxP];
  Float_t  eleIP3D_[maxP];
  Float_t  eleIP3DErr_[maxP];
  Float_t  eleIDMVANonTrig_[maxP];
  Float_t  eleIDMVATrig_[maxP];
  Int_t    eleID2012_[maxP][4];
  Int_t    eleESDetId_[maxP][2];
  Float_t  eleESHits_[maxP][3][62];
  Float_t  eleESEffSigmaRR_[maxP][3];
  Float_t  eleESE1_[maxP][2];
  Float_t  eleESE3_[maxP][2];
  Float_t  eleESE5_[maxP][2];
  Float_t  eleESE7_[maxP][2];
  Float_t  eleESE11_[maxP][2];
  Float_t  eleESE21_[maxP][2];
  Int_t	   eleNBC_[maxP];
  Float_t  eleBrLinear_[maxP];
  Float_t  eleCetaCorrE_[maxP];
  Float_t  eleCetaCorrEt_[maxP];
  Float_t  eleBremCorrE_[maxP];
  Float_t  eleBremCorrEt_[maxP];
  Float_t  eleFullCorrE_[maxP];
  Float_t  eleFullCorrEt_[maxP];
  Float_t  elePFChIso03_[maxP];
  Float_t  elePFPhoIso03_[maxP];
  Float_t  elePFNeuIso03_[maxP];
  Float_t  elePFChIso04_[maxP];
  Float_t  elePFPhoIso04_[maxP];
  Float_t  elePFNeuIso04_[maxP];
  // Photon
  Int_t    nPho_;
  Int_t    phoTrg_[maxP][8];
  Int_t    phoTrgFilter_[maxP][50];
  Bool_t   phoIsPhoton_[maxP];
  Float_t  phoE_[maxP];
  Float_t  phoEt_[maxP];
  Float_t  phoEta_[maxP];
  Float_t  phoPhi_[maxP];
  Float_t  phoEtVtx_[maxP][100];
  Float_t  phoEtaVtx_[maxP][100];
  Float_t  phoPhiVtx_[maxP][100];
  Float_t  phoVtx_[maxP][3];
  Float_t  phoSCPos_[maxP][3];
  Float_t  phoCaloPos_[maxP][3];
  Float_t  phoR9_[maxP];
  Float_t  phoCetaCorrE_[maxP];
  Float_t  phoCetaCorrEt_[maxP];
  Float_t  phoBremCorrE_[maxP];
  Float_t  phoBremCorrEt_[maxP];
  Float_t  phoFullCorrE_[maxP];
  Float_t  phoFullCorrEt_[maxP];
  Float_t  phoTrkIsoSolidDR03_[maxP];
  Float_t  phoTrkIsoHollowDR03_[maxP];
  Float_t  phoEcalIsoDR03_[maxP];
  Float_t  phoHcalIsoSolidDR03_[maxP];
  Float_t  phoHcalIsoDR03_[maxP];
  Float_t  phoHcalIsoDR0312_[maxP];
  Float_t  phoTrkIsoSolidDR04_[maxP];
  Float_t  phoTrkIsoHollowDR04_[maxP];
  Float_t  phoEcalIsoDR04_[maxP];
  Float_t  phoHcalIsoDR04_[maxP];
  Float_t  phoHcalIsoDR0412_[maxP];
  Float_t  phoHcalIsoSolidDR04_[maxP];
  Float_t  phoCiCTrkIsoDR03_[maxP][100];
  Float_t  phoCiCTrkIsoDR04_[maxP][100];
  Float_t  phoCiCdRtoTrk_[maxP];
  Float_t  phoHoverE_[maxP];
  Float_t  phoHoverEBCdepth1_[maxP];
  Float_t  phoHoverEBCdepth2_[maxP];
  Float_t  phoHoverE12_[maxP];
  Float_t  phoSigmaIEtaIEta_[maxP];
  Float_t  phoSigmaIEtaIPhi_[maxP];
  Float_t  phoSigmaIPhiIPhi_[maxP];
  Float_t  phoEmax_[maxP];
  Float_t  phoE3x3_[maxP];
  Float_t  phoE5x5_[maxP];
  Float_t  phoE2x5Max_[maxP];
  Float_t  phoE5x1_[maxP];
  Float_t  phoE1x5_[maxP];
  Float_t  phoE3x1_[maxP];
  Float_t  phoE1x3_[maxP];
  Float_t  phoE2x2_[maxP];
  Float_t  phoRegrE_[maxP];
  Float_t  phoRegrEerr_[maxP];
  Float_t  phoPFChIso_[maxP];
  Float_t  phoPFPhoIso_[maxP];
  Float_t  phoPFNeuIso_[maxP];
  Float_t  phoSCRChIso_[maxP];
  Float_t  phoSCRPhoIso_[maxP];
  Float_t  phoSCRNeuIso_[maxP];
  Float_t  phoSCRChIso04_[maxP];
  Float_t  phoSCRPhoIso04_[maxP];
  Float_t  phoSCRNeuIso04_[maxP];
  Float_t  phoRandConeChIso_[maxP];
  Float_t  phoRandConePhoIso_[maxP];
  Float_t  phoRandConeNeuIso_[maxP];
  Float_t  phoRandConeChIso04_[maxP];
  Float_t  phoRandConePhoIso04_[maxP];
  Float_t  phoRandConeNeuIso04_[maxP];
  Float_t  phoSeedTime_[maxP];
  Float_t  phoLICTD_[maxP];
  Float_t  phoCiCPF4phopfIso005_[maxP];
  Float_t  phoCiCPF4phopfIso01_[maxP];
  Float_t  phoCiCPF4phopfIso02_[maxP];
  Float_t  phoCiCPF4phopfIso03_[maxP];
  Float_t  phoCiCPF4phopfIso04_[maxP];
  Float_t  phoCiCPF4phopfIso05_[maxP];
  Float_t  phoCiCPF4phopfIso06_[maxP];
  Float_t  phoCiCPF4phopfIso07_[maxP];
  Float_t  phoCiCPF4phopfIso08_[maxP];
  Float_t  phoCiCPF4chgpfIso005_[maxP][100];
  Float_t  phoCiCPF4chgpfIso01_[maxP][100];
  Float_t  phoCiCPF4chgpfIso02_[maxP][100];
  Float_t  phoCiCPF4chgpfIso03_[maxP][100];
  Float_t  phoCiCPF4chgpfIso04_[maxP][100];
  Float_t  phoCiCPF4chgpfIso05_[maxP][100];
  Float_t  phoCiCPF4chgpfIso06_[maxP][100];
  Float_t  phoCiCPF4chgpfIso07_[maxP][100];
  Float_t  phoCiCPF4chgpfIso08_[maxP][100];
  
  Float_t  phoCiCPF4phopfIsoNoVETO005_[maxP];
  Float_t  phoCiCPF4phopfIsoNoVETO01_[maxP];
  Float_t  phoCiCPF4phopfIsoNoVETO02_[maxP];
  Float_t  phoCiCPF4phopfIsoNoVETO03_[maxP];
  Float_t  phoCiCPF4phopfIsoNoVETO04_[maxP];
  Float_t  phoCiCPF4phopfIsoNoVETO05_[maxP];
  Float_t  phoCiCPF4phopfIsoNoVETO06_[maxP];
  Float_t  phoCiCPF4phopfIsoNoVETO07_[maxP];
  Float_t  phoCiCPF4phopfIsoNoVETO08_[maxP];
  Float_t  phoCiCPF4chgpfIsoNoVETO005_[maxP][100];
  Float_t  phoCiCPF4chgpfIsoNoVETO01_[maxP][100];
  Float_t  phoCiCPF4chgpfIsoNoVETO02_[maxP][100];
  Float_t  phoCiCPF4chgpfIsoNoVETO03_[maxP][100];
  Float_t  phoCiCPF4chgpfIsoNoVETO04_[maxP][100];
  Float_t  phoCiCPF4chgpfIsoNoVETO05_[maxP][100];
  Float_t  phoCiCPF4chgpfIsoNoVETO06_[maxP][100];
  Float_t  phoCiCPF4chgpfIsoNoVETO07_[maxP][100];
  Float_t  phoCiCPF4chgpfIsoNoVETO08_[maxP][100];
  
  Float_t  phoCiCPF4chgpfIso02AppVeto_[maxP][100];
  Float_t  phoCiCPF4chgpfIso03AppVeto_[maxP][100];
  Float_t  phoCiCPF4chgpfIso04AppVeto_[maxP][100];
  Float_t  phoCiCPF4phopfIso03AppVeto_[maxP];
  Float_t  phoCiCPF4phopfIso04AppVeto_[maxP];
  Float_t  phoCiCPF4phopfIso03AppVetoClean_[maxP];
  Float_t  phoCiCPF4phopfIso04AppVetoClean_[maxP];
  Float_t  phoCiCPF4chgpfIso03Mod_[maxP][100];
  Float_t  phoCiCPF4chgpfIso04Mod_[maxP][100];
  
  
  Int_t    phoSeedDetId1_[maxP];
  Int_t    phoSeedDetId2_[maxP];
  Int_t    phoRecoFlag_[maxP];
  Int_t    phoPos_[maxP];
  Int_t    phoGenIndex_[maxP];
  Int_t    phoGenGMomPID[maxP];
  Int_t    phoGenMomPID[maxP];
  Float_t  phoGenMomPt[maxP];
  Float_t  phoSCE_[maxP];
  Float_t  phoSCRawE_[maxP];
  Float_t  phoESEn_[maxP];
  Float_t  phoSCEt_[maxP];
  Float_t  phoSCEta_[maxP];
  Float_t  phoSCPhi_[maxP];
  Float_t  phoSCEtaWidth_[maxP];
  Float_t  phoSCPhiWidth_[maxP];
  Float_t  phoSCBrem_[maxP];
  Int_t    phoOverlap_[maxP];
  Int_t    phohasPixelSeed_[maxP];
  Int_t    phoIsConv_[maxP];
  Int_t    phoEleVeto_[maxP];
  Int_t    phoESDetId_[maxP][2];
  Float_t  phoESHits_[maxP][3][62];
  Float_t  phoESEffSigmaRR_[maxP][3];
  Float_t  phoESE1_[maxP][2];
  Float_t  phoESE3_[maxP][2];
  Float_t  phoESE5_[maxP][2];
  Float_t  phoESE7_[maxP][2];
  Float_t  phoESE11_[maxP][2];
  Float_t  phoESE21_[maxP][2];
  Int_t    phoNConv_[maxP];
  Int_t    phoConvNTrks_[maxP];
  Float_t  phoConvInvMass_[maxP];
  Float_t  phoConvCotTheta_[maxP];
  Float_t  phoConvEoverP_[maxP];
  Float_t  phoConvZofPVfromTrks_[maxP];
  Float_t  phoConvMinDist_[maxP];
  Float_t  phoConvdPhiAtVtx_[maxP];
  Float_t  phoConvdPhiAtCalo_[maxP];
  Float_t  phoConvdEtaAtCalo_[maxP];
  Float_t  phoConvTrkd0_[maxP][2];
  Float_t  phoConvTrkPin_[maxP][2];
  Float_t  phoConvTrkPout_[maxP][2];
  Float_t  phoConvTrkdz_[maxP][2];
  Float_t  phoConvTrkdzErr_[maxP][2];
  Float_t  phoConvChi2_[maxP];
  Float_t  phoConvChi2Prob_[maxP];
  Float_t  phoConvCharge_[maxP][2];
  Float_t  phoConvValidVtx_[maxP];
  Float_t  phoConvLikeLihood_[maxP];
  Float_t  phoConvP4_[maxP][4];
  Float_t  phoConvVtx_[maxP][3];
  Float_t  phoConvVtxErr_[maxP][3];
  Float_t  phoConvPairMomentum_[maxP][3];
  Float_t  phoConvRefittedMomentum_[maxP][3];
  Int_t    SingleLegConv_[maxP];
  Float_t  phoPFConvRefitMom_[maxP][3];
  Float_t  phoPFConvMom_[maxP][3];
  Float_t  phoPFConvVtx_[maxP][3];
  
  //new PF Variables stored when Matched to Reco:
  Int_t    PFRecoMatch_[maxP];
  Int_t    PFEleMatch_[maxP];
  Int_t    PFEleVeto_[maxP];    
  Float_t  PFPreShowerE1_[maxP];
  Float_t  PFPreShowerE2_[maxP];
  Float_t  MustacheEin_[maxP];
  Float_t  MustacheEOut_[maxP];
  Float_t  MustacheEtOut_[maxP];
  Float_t  PFLowestClustE_[maxP];
  Float_t  PFClustdEta_[maxP];
  Float_t  PFClustdPhi_[maxP];
  Float_t  PFClustRMSPhi_[maxP];
  Float_t  PFClustRMSPhiMust_[maxP];
  Float_t  PFClustEneCorr_[maxP];
  Int_t    pho_hasSLConvPf_[maxP];
  Int_t    pho_hasConvPf_[maxP];
  Float_t  pho_pfconvVtxZ_[maxP];
  Float_t  pho_pfconvVtxZErr_[maxP];
  Int_t    pho_nSLConv_[maxP];
  Float_t  pho_pfSLConvPos_[maxP][20][3];
  Float_t  pho_pfSLConvVtxZ_[maxP][20];
  // Muon
  Int_t    nMu_;
  Int_t    muTrg_[maxP][10];
  Float_t  muEta_[maxP];
  Float_t  muPhi_[maxP];
  Int_t    muCharge_[maxP];
  Float_t  muPt_[maxP];
  Float_t  muPz_[maxP];
  Float_t  muVtx_[maxP][3];
  Float_t  muVtxGlb_[maxP][3];
  Int_t    muGenIndex_[maxP];
  Float_t  mucktPt_[maxP];
  Float_t  mucktPtErr_[maxP];
  Float_t  mucktEta_[maxP];
  Float_t  mucktPhi_[maxP];
  Float_t  mucktdxy_[maxP];
  Float_t  mucktdz_[maxP];
  Float_t  muIsoTrk_[maxP];
  Float_t  muIsoCalo_[maxP];
  Float_t  muIsoEcal_[maxP];
  Float_t  muIsoHcal_[maxP];
  Float_t  muChi2NDF_[maxP];
  Float_t  muInnerChi2NDF_[maxP];
  Float_t  muPFIsoR04_CH_[maxP];
  Float_t  muPFIsoR04_NH_[maxP];
  Float_t  muPFIsoR04_Pho_[maxP];
  Float_t  muPFIsoR04_PU_[maxP];
  Float_t  muPFIsoR04_CPart_[maxP];
  Float_t  muPFIsoR04_NHHT_[maxP];
  Float_t  muPFIsoR04_PhoHT_[maxP];
  Float_t  muPFIsoR03_CH_[maxP];
  Float_t  muPFIsoR03_NH_[maxP];
  Float_t  muPFIsoR03_Pho_[maxP];
  Float_t  muPFIsoR03_PU_[maxP];
  Float_t  muPFIsoR03_CPart_[maxP];
  Float_t  muPFIsoR03_NHHT_[maxP];
  Float_t  muPFIsoR03_PhoHT_[maxP];
  Int_t    muType_[maxP];
  Float_t  muD0_[maxP];
  Float_t  muDz_[maxP];
  Float_t  muD0GV_[maxP];
  Float_t  muDzGV_[maxP];
  Float_t  muD0Vtx_[maxP][100];
  Float_t  muDzVtx_[maxP][100];
  Float_t  muInnerD0_[maxP];
  Float_t  muInnerDz_[maxP];
  Float_t  muInnerD0GV_[maxP];
  Float_t  muInnerDzGV_[maxP];
  Float_t  muInnerPt_[maxP];
  Float_t  muInnerPtErr_[maxP];
  Int_t    muNumberOfValidTrkLayers_[maxP]; 
  Int_t    muNumberOfValidTrkHits_[maxP];
  Int_t    muNumberOfValidPixelLayers_[maxP];
  Int_t    muNumberOfValidPixelHits_[maxP];
  Int_t    muNumberOfValidMuonHits_[maxP];
  Int_t    muStations_[maxP];
  Int_t    muChambers_[maxP];
  Float_t  muIP3D_[maxP];
  Float_t  muIP3DErr_[maxP];
  //PFPhotons  
  Int_t    nPFPho_;  
  Float_t  PFPhoE_[maxP];  
  Float_t  PFPhoEt_[maxP];  
  Float_t  PFPhoEta_[maxP];  
  Float_t  PFPhoPhi_[maxP];  
  Int_t    PFPhoType_[maxP];
  Float_t  PFPhoIso_[maxP];
  Int_t    nPFPhoClust_[maxP];
  Float_t  PFPhoMaxEtaWidth_[maxP];
  Float_t  PFPhoMaxPhiWidth_[maxP];  
  Int_t    nPFPho_clust; 
  Int_t    nPFPhoCrys_[maxP][20]; 
  Float_t  PFPho_clusteta_[maxP][20];  
  Float_t  PFPho_clustphi_[maxP][20];  
  Float_t  PFPho_clustE_[maxP][20];
  Float_t  PFPho_clustSCEfrac_[maxP][20];
  Float_t  PFPho_clustH_[maxP][20];
  Float_t  PFPho_clustEseed_[maxP][20];
  Float_t  PFPho_clustEtop_[maxP][20];
  Float_t  PFPho_clustEbottom_[maxP][20];
  Float_t  PFPho_clustEleft_[maxP][20];
  Float_t  PFPho_clustEright_[maxP][20];
  Float_t  PFPho_clustE3x3_[maxP][20];
  Float_t  PFPho_clustE3x1_[maxP][20];
  Float_t  PFPho_clustE1x3_[maxP][20];
  Float_t  PFPho_clustE5x5_[maxP][20];
  Float_t  PFPho_clustE5x1_[maxP][20];
  Float_t  PFPho_clustE1x5_[maxP][20];
  Float_t  PFPho_clustE2x5Max_[maxP][20];
  Float_t  PFPho_clustE2x5Top_[maxP][20];
  Float_t  PFPho_clustE2x5Bottom_[maxP][20];
  Float_t  PFPho_clustE2x5Left_[maxP][20];
  Float_t  PFPho_clustE2x5Right_[maxP][20];
  Float_t  PFPho_clustER9_[maxP][20];
  Float_t  PFPho_clustEt_[maxP][20];
  Float_t  PFPho_sigIetaIeta_[maxP][20];
  Float_t  PFPho_sigIphiIphi_[maxP][20];
  Float_t  PFPho_crysphi_[maxP][20];
  Float_t  PFPho_cryseta_[maxP][20];
  Float_t  PFPho_crysphifix_[maxP][20]; 	 
  Float_t  PFPho_crysetafix_[maxP][20]; 	 
  Float_t  PFPho_modphifix_[maxP][20]; 	 
  Float_t  PFPho_modetafix_[maxP][20]; 	 
  Float_t  PFPho_Smodphifix_[maxP][20]; 	 
  Float_t  PFPho_Smodetafix_[maxP][20];
  Float_t  PFPho_crysphiTilt_[maxP][20]; 	 
  Float_t  PFPho_crysthetaTilt_[maxP][20];

  Float_t  PFPho_crysxfix_[maxP][20]; 	 
  Float_t  PFPho_crysyfix_[maxP][20]; 	 
  Float_t  PFPho_modxfix_[maxP][20]; 	 
  Float_t  PFPho_modyfix_[maxP][20]; 	 
  Float_t  PFPho_Smodxfix_[maxP][20]; 	 
  Float_t  PFPho_Smodyfix_[maxP][20];
  Float_t  PFPho_ES1Energy_[maxP][20];
  Float_t  PFPho_ES2Energy_[maxP][20];
  Int_t    nPFPho_ES1clust;
  Int_t    nPFPhoES1Clust_[maxP];
  Float_t  PFPho_ES1clustE_[maxP][30];
  Float_t  PFPho_ES1Clust_perEEclust_[maxP][30];
  Int_t    nPFPho_ES1Clust_perEEclust_[maxP][20];
  Int_t    PFPho_ES1size_perEEclust_[maxP][20];
  Float_t  PFPho_ES1weightedX_perEEclust_[maxP][20]; 
  Float_t  PFPho_ES1linkD_perEEclust_[maxP][20];
  Float_t  PFPho_ES1clustx_[maxP][30];
  Float_t  PFPho_ES1clusty_[maxP][30];
  Float_t  PFPho_ES1clustLinkDist_[maxP][30];
  Float_t  PFPho_ES1clusteta_[maxP][30];
  Float_t  PFPho_ES1clustphi_[maxP][30];
  Float_t  PFPho_ES1clustz_[maxP][30];
  Int_t    PFPho_ES1_stripsDetId_[maxP][30][30];
  Float_t  PFPho_ES1_stripsX_[maxP][30][30];
  Float_t  PFPho_ES1_stripsY_[maxP][30][30];
  Float_t  PFPho_ES1_stripsZ_[maxP][30][30];
  Float_t  PFPho_ES1_stripsEta_[maxP][30][30];
  Float_t  PFPho_ES1_stripsPhi_[maxP][30][30];
  Float_t  PFPho_ES1_stripsFrac_[maxP][30][30];
  Float_t  PFPho_ES1_stripsE_[maxP][30][30];
  Int_t    PFPho_ES1size_[maxP][30];
  Float_t  PFPho_weightedXst_perES1clust_[maxP][30];

  Int_t    nPFPho_ES2clust;
  Int_t    nPFPhoES2Clust_[maxP];
  Float_t  PFPho_ES2clustE_[maxP][30];
  Float_t  PFPho_ES2Clust_perEEclust_[maxP][30];
  Int_t    nPFPho_ES2Clust_perEEclust_[maxP][20];
  Int_t    PFPho_ES2size_perEEclust_[maxP][20];
  Float_t  PFPho_ES2weightedY_perEEclust_[maxP][20];
  Float_t  PFPho_ES_siphiiphi_[maxP][20];
  Float_t  PFPho_ES_sietaieta_[maxP][20];
  Float_t  PFPho_ES2linkD_perEEclust_[maxP][20];
  Float_t  PFPho_ES2clustx_[maxP][30];
  Float_t  PFPho_ES2clusty_[maxP][30];
  Float_t  PFPho_ES2clustLinkDist_[maxP][30];
  Float_t  PFPho_ES2clusteta_[maxP][30];
  Float_t  PFPho_ES2clustphi_[maxP][30];
  Float_t  PFPho_ES2clustz_[maxP][30];
  Int_t    PFPho_ES2_stripsDetId_[maxP][30][30];
  Float_t  PFPho_ES2_stripsX_[maxP][30][30];
  Float_t  PFPho_ES2_stripsY_[maxP][30][30];
  Float_t  PFPho_ES2_stripsZ_[maxP][30][30];
  Float_t  PFPho_ES2_stripsEta_[maxP][30][30];
  Float_t  PFPho_ES2_stripsPhi_[maxP][30][30];
  Float_t  PFPho_ES2_stripsFrac_[maxP][30][30];
  Float_t  PFPho_ES2_stripsE_[maxP][30][30];
  Int_t    PFPho_ES2size_[maxP][30];
  Float_t  PFPho_weightedYst_perES2clust_[maxP][30];

  Int_t    PFPho_crysIphi_[maxP][20];
  Int_t    PFPho_crysIeta_[maxP][20];
  Int_t    PFPho_crysIX_[maxP][20];
  Int_t    PFPho_crysIY_[maxP][20];
  Float_t  PFPhoMustEout_[maxP];
  Float_t  PFPhoMustEin_[maxP];  
  Float_t  PFPhoMustEtOut_[maxP];
  Float_t  PFPhoMustExcl_[maxP];
  Int_t    hasGSF_[maxP];
  Float_t  PFPhoGSFPin_[maxP][3];  
  Float_t  PFPhoGSFPout_[maxP];  
  Float_t  PFPhoGSFChi2NDF_[maxP];  
  Int_t    PFPhoCharge_[maxP];  
  Float_t  PFPhoGsf_In_[maxP][3];
  Float_t  PFPhoGsf_Theta_[maxP];
  Float_t  PFPhoGsf_ThetaErr_[maxP];  
  Float_t  PFPhoGsfeta_In_[maxP];  
  Float_t  PFPhoGsfeta_Out_[maxP];  
  Float_t  PFPhoGsfphi_In_[maxP];  
  Float_t  PFPhoGsfphi_Out_[maxP];
  Int_t    nPFPho_tks;  
  Int_t    nPFPhoTks_[maxP];
  Int_t    PFPho_tkq_[maxP][20];  
  Float_t  PFPho_tketa_[maxP][20];  
  Float_t  PFPho_tkphi_[maxP][20];  
  Float_t  PFPho_tkpt_[maxP][20];  
  Float_t  PFPho_tkR_[maxP][20];  
  Float_t  PFPho_tkZ_[maxP][20];
  Float_t  PFPho_tkTheta_[maxP][20];
  Float_t  PFPho_tkerreta_[maxP][20];
  Float_t  PFPho_tkerrphi_[maxP][20];
  Float_t  PFPho_tkerrthet_[maxP][20];
  Float_t  PFPho_tkPos_[maxP][20][3];
  Int_t    PFPhoIsConv_[maxP];  
  Float_t  PFPhoTrkIsoHollowDR04_[maxP];  
  Float_t  PFPhoEcalIsoDR04_[maxP];  
  Float_t  PFPhoHcalIsoDR04_[maxP];  
  Float_t  PFPhoHoverE_[maxP];  
  Int_t    isIsoPFPho_[maxP];  
  Float_t  PFPhoR9_[maxP];
  Float_t  PFPhoIetaIeta_[maxP]; 
  Float_t  PFPhoE5x5_[maxP];
  Float_t  PFPhoE3x3_[maxP];
  Int_t    nPFPhoIso_[maxP];
  Float_t  PFPho_isoeta_[maxP][50];  
  Float_t  PFPho_isophi_[maxP][50];  
  Float_t  PFPho_isodeltaR_[maxP][50];  
  Float_t  PFPho_isovalue_[maxP][50];    
  Float_t  PFPho_isotype_[maxP][50]; 

  //PFElectrons:  
  Int_t    nPFEle_;
  Int_t    nPFEleClust_[maxP];
  Float_t  PFElePt_[maxP];  
  Float_t  PFEleEta_[maxP];  
  Float_t  PFElePhi_[maxP];
  Float_t  PFEleMaxEtaWidth_[maxP];
  Float_t  PFEleMaxPhiWidth_[maxP];
  Int_t    nPFEle_clust;  
  Int_t    nPFEleCrys_[maxP][20];
  Float_t  PFEle_clusteta_[maxP][20];  
  Float_t  PFEle_clustphi_[maxP][20];  
  Float_t  PFEle_clustE_[maxP][20];
  Float_t  PFEle_clustSCEfrac_[maxP][20];
  Float_t  PFEle_clustEt_[maxP][20];
  Float_t  PFEle_clustH_[maxP][20];
  Float_t  PFEle_clustEseed_[maxP][20];
  Float_t  PFEle_clustEtop_[maxP][20];
  Float_t  PFEle_clustEbottom_[maxP][20];
  Float_t  PFEle_clustEleft_[maxP][20];
  Float_t  PFEle_clustEright_[maxP][20];
  Float_t  PFEle_clustE3x3_[maxP][20];
  Float_t  PFEle_clustE3x1_[maxP][20];
  Float_t  PFEle_clustE1x3_[maxP][20];
  Float_t  PFEle_clustE5x5_[maxP][20];
  Float_t  PFEle_clustE5x1_[maxP][20];
  Float_t  PFEle_clustE1x5_[maxP][20];
  Float_t  PFEle_clustE2x5Max_[maxP][20];
  Float_t  PFEle_clustE2x5Top_[maxP][20];
  Float_t  PFEle_clustE2x5Bottom_[maxP][20];
  Float_t  PFEle_clustE2x5Left_[maxP][20];
  Float_t  PFEle_clustE2x5Right_[maxP][20];
  Float_t  PFEle_clustER9_[maxP][20];
  Int_t    PFEle_crysIphi_[maxP][20];
  Int_t    PFEle_crysIeta_[maxP][20];
  Int_t    PFEle_crysIX_[maxP][20];
  Int_t    PFEle_crysIY_[maxP][20];
  Float_t  PFEle_crysphiTilt_[maxP][20]; 	 
  Float_t  PFEle_crysthetaTilt_[maxP][20];
  Float_t  PFEle_sigIetaIeta_[maxP][20];
  Float_t  PFEle_sigIphiIphi_[maxP][20];
  Float_t  PFEle_crysx_[maxP][20];
  Float_t  PFEle_crysy_[maxP][20];
  Float_t  PFEle_crysxfix_[maxP][20]; 	 
  Float_t  PFEle_crysyfix_[maxP][20]; 	 
  Float_t  PFEle_modxfix_[maxP][20]; 	 
  Float_t  PFEle_modyfix_[maxP][20]; 	 
  Float_t  PFEle_Smodxfix_[maxP][20]; 	 
  Float_t  PFEle_Smodyfix_[maxP][20];
  Float_t  PFEle_cryseta_[maxP][20];
  Float_t  PFEle_crysphi_[maxP][20];
  Float_t  PFEle_crysetafix_[maxP][20]; 	 
  Float_t  PFEle_crysphifix_[maxP][20]; 	 
  Float_t  PFEle_modetafix_[maxP][20]; 	 
  Float_t  PFEle_modphifix_[maxP][20]; 	 
  Float_t  PFEle_Smodetafix_[maxP][20]; 	 
  Float_t  PFEle_Smodphifix_[maxP][20];
  Float_t  PFEle_ES1Energy_[maxP][20];
  Float_t  PFEle_ES2Energy_[maxP][20];
  Int_t    nPFEle_ES1clust;
  Int_t    nPFEleES1Clust_[maxP];
  Float_t  PFEle_ES1clustE_[maxP][30];
  Float_t  PFEle_ES1clusteta_[maxP][30];
  Float_t  PFEle_ES1clustphi_[maxP][30];
  Float_t  PFEle_ES1clustx_[maxP][30];
  Float_t  PFEle_ES1clusty_[maxP][30];
  Float_t  PFEle_ES1clustz_[maxP][30];
  Int_t    nPFEle_ES2clust;
  Int_t    nPFEleES2Clust_[maxP];
  Float_t  PFEle_ES2clustE_[maxP][30];
  Float_t  PFEle_ES2clusteta_[maxP][30];
  Float_t  PFEle_ES2clustphi_[maxP][30];
  Float_t  PFEle_ES2clustx_[maxP][30];
  Float_t  PFEle_ES2clusty_[maxP][30];
  Float_t  PFEle_ES2clustz_[maxP][30];
  Float_t  PFEleMustEtOut_[maxP];
  Float_t  PFEleMustExcl_[maxP];
  Float_t  PFEleMustEout_[maxP];
  Float_t  PFEleMustEin_[maxP];  
  Float_t  PFEleHoverE_[maxP];  
  Float_t  PFEleEoverP_[maxP];  
  Float_t  PFElePin_[maxP][3];  
  Float_t  PFElePout_[maxP];  
  Float_t  PFEleChi2NDF_[maxP];  
  Int_t    PFEleCharge_[maxP];  
  Float_t  PFEleEn_[maxP];  
  Float_t  PFEleSCEta_[maxP];  
  Float_t  PFEleSCPhi_[maxP];  
  Float_t  PFEleSCRawEn_[maxP];
  Float_t  PFGsf_In_[maxP][3];
  Float_t  PFGsf_Theta_[maxP];
  Float_t  PFGsf_ThetaErr_[maxP];  
  Float_t  PFGsfeta_In_[maxP];  
  Float_t  PFGsfeta_Out_[maxP];  
  Float_t  PFGsfphi_In_[maxP];  
  Float_t  PFGsfphi_Out_[maxP];  
  Int_t    nPFEleIso_[maxP];
  Float_t  PFEle_isoeta_[maxP][50];  
  Float_t  PFEle_isophi_[maxP][50];  
  Float_t  PFEle_isodeltaR_[maxP][50];  
  Float_t  PFEle_isovalue_[maxP][50];    
  Float_t  PFEle_isotype_[maxP][50]; 
  //PFHadrons:  //charged:  
  Int_t    nPFchad_;  
  Int_t    PFhad_charge_[maxP];  
  Float_t  PFchad_Eta_[maxP];  
  Float_t  PFchad_Phi_[maxP];  
  Float_t  PFchad_Pt_[maxP];  
  Float_t  PFchad_E_[maxP];  
  Float_t  PFchad_P_[maxP];  
  Float_t  PFchad_Ecal_[maxP];  
  Float_t  PFchad_Hcal_[maxP];  
  //neutral:  
  Int_t    nPFnhad_;  
  Float_t  PFnhad_Eta_[maxP];  
  Float_t  PFnhad_Phi_[maxP]; 
  Float_t  PFnhad_Pt_[maxP];  
  Float_t  PFnhad_E_[maxP];  
  Float_t  PFnhad_P_[maxP];  
  Float_t  PFnhad_Hcal_[maxP];    
  //PFClusters  //ECAL:  
  Int_t    nPFEcal_;  
  Float_t  PFEcal_eta_[maxP];  
  Float_t  PFEcal_phi_[maxP];  
  Float_t  PFEcal_energy_[maxP];    
  // rho correction
  Float_t  rho25_;
  Float_t  rho25_neu_;
  Float_t  rho25_muPFiso_;
  Float_t  rho25_elePFiso_;
  Float_t  rho2011_;
  Float_t  rho2012_;
  // Jet
  Int_t    nJet_;
  Int_t    jetTrg_[maxP][14];
  Int_t    jetAlgo_[maxP];
  Float_t  jetEn_[maxP];
  Float_t  jetPt_[maxP];
  Float_t  jetEta_[maxP];
  Float_t  jetPhi_[maxP];
  Float_t  jetEt_[maxP];
  Float_t  jetRawPt_[maxP];
  Float_t  jetRawEn_[maxP];
  Float_t  jetCharge_[maxP];
  Float_t  jetArea_[maxP];
  Float_t  jetCHF_[maxP];
  Float_t  jetNHF_[maxP];
  Float_t  jetCEF_[maxP];
  Float_t  jetNEF_[maxP];
  Int_t    jetNCH_[maxP];
  Float_t  jetHFHAE_[maxP];
  Float_t  jetHFEME_[maxP];
  Int_t    jetPartonID_[maxP];
  Int_t    jetNConstituents_[maxP];
  Float_t  jetCombinedSecondaryVtxBJetTags_[maxP]; ////  rob      ------------> recommended
  Float_t  jetCombinedSecondaryVtxMVABJetTags_[maxP]; ////  rob
  Float_t  jetBetaStar_[maxP][100];
  Int_t    jetGenJetIndex_[maxP];
  Float_t  jetGenJetEn_[maxP];
  Float_t  jetGenJetPt_[maxP];
  Float_t  jetGenJetEta_[maxP];
  Float_t  jetGenJetPhi_[maxP];
  Int_t    jetGenPartonID_[maxP];
  Float_t  jetGenEn_[maxP];
  Float_t  jetGenPt_[maxP];
  Float_t  jetGenEta_[maxP];
  Float_t  jetGenPhi_[maxP];
  Int_t    jetGenPartonMomID_[maxP];
  // Jet ID MVA variables
  Float_t  jetMVAs_[maxP][4];
  Int_t    jetWPLevels_[maxP][4];
  Float_t  jetMVAsExt_[maxP][4][100];
  Int_t    jetWPLevelsExt_[maxP][4][100];
  // std::vector<float * > jetMVAs_;
  // std::vector<int * >   jetWPLevels_;
  // std::vector<std::vector<std::vector<float> > * > jetMVAsExt_;
  // std::vector<std::vector<std::vector<int> > * > jetWPLevelsExt_;
  std::vector<PileupJetIdAlgo* > pujetIDalgos_;
  Float_t jetDRMean_[maxP];
  Float_t jetDR2Mean_[maxP];
  Float_t jetDZ_[maxP];
  Float_t jetFrac01_[maxP];
  Float_t jetFrac02_[maxP];
  Float_t jetFrac03_[maxP];
  Float_t jetFrac04_[maxP];
  Float_t jetFrac05_[maxP];
  Float_t jetFrac06_[maxP];
  Float_t jetFrac07_[maxP];
  Float_t jetBeta_[maxP];
  Float_t jetBetaStarCMG_[maxP];
  Float_t jetBetaStarClassic_[maxP];
  Float_t jetBetaExt_[maxP][100];
  Float_t jetBetaStarCMGExt_[maxP][100];
  Float_t jetBetaStarClassicExt_[maxP][100];
  Float_t jetNNeutrals_[maxP];
  Float_t jetNCharged_[maxP];
  Bool_t  jetPFLooseId_[maxP];
  // b-jet regression variables
  Float_t jetMt_[maxP];
  Float_t jetJECUnc_[maxP];
  Float_t jetLeadTrackPt_[maxP];
  Float_t jetVtxPt_[maxP];
  Float_t jetVtxMass_[maxP];
  Float_t jetVtx3dL_[maxP];
  Float_t jetVtx3deL_[maxP];
  Float_t jetSoftLeptPt_[maxP];
  Float_t jetSoftLeptPtRel_[maxP];
  Float_t jetSoftLeptdR_[maxP];
  Float_t jetSoftLeptIdlooseMu_[maxP];
  Float_t jetSoftLeptIdEle95_[maxP];
  Float_t jetDPhiMETJet_[maxP];
  Float_t jetPuJetIdL_[maxP];
  Float_t jetPuJetIdM_[maxP];
  Float_t jetPuJetIdT_[maxP];
  // Low Pt Jets
  Int_t   nLowPtJet_;
  Float_t jetLowPtEn_[maxP];
  Float_t jetLowPtPt_[maxP];
  Float_t jetLowPtEta_[maxP];
  Float_t jetLowPtPhi_[maxP];
  Float_t jetLowPtCharge_[maxP];
  Float_t jetLowPtEt_[maxP];
  Float_t jetLowPtRawPt_[maxP];
  Float_t jetLowPtRawEn_[maxP];
  Float_t jetLowPtArea_[maxP];
  Float_t jetLowPtGenJetEn_[maxP];
  Float_t jetLowPtGenJetPt_[maxP];
  Float_t jetLowPtGenJetEta_[maxP];
  Float_t jetLowPtGenJetPhi_[maxP];
  Int_t   jetLowPtPartonID_[maxP];
  Int_t   jetLowPtGenPartonID_[maxP];
  Float_t jetLowPtGenEn_[maxP];
  Float_t jetLowPtGenPt_[maxP];
  Float_t jetLowPtGenEta_[maxP];
  Float_t jetLowPtGenPhi_[maxP];
  // Converted Photon Collection
  Int_t    nConv_;
  Int_t    convBarrel_[maxP];
  Int_t    convScInd_[maxP];
  Float_t  convP4_[maxP][4];
  Float_t  convVtx_[maxP][3];
  Float_t  convVtxErr_[maxP][3];
  Float_t  convPairMomentum_[maxP][3];
  Float_t  convRefittedMomentum_[maxP][3];
  Int_t    convNTracks_[maxP];
  Float_t  convPairInvMass_[maxP];
  Float_t  convPairCotThetaSep_[maxP];
  Float_t  convEoverP_[maxP];
  Float_t  convDistOfMinApproach_[maxP];
  Float_t  convDPhiTrksAtVtx_[maxP];
  Float_t  convDPhiTrksAtEcal_[maxP];
  Float_t  convDEtaTrksAtEcal_[maxP];
  Float_t  convDxy_[maxP]; 
  Float_t  convDz_[maxP];    
  Float_t  convLxy_[maxP]; 
  Float_t  convLz_[maxP];    
  Float_t  convZofPrimVtxFromTrks_[maxP];
  Int_t    convNHitsBeforeVtx_[maxP][2];
  Int_t    convNSharedHits_[maxP];
  Int_t    convValidVtx_[maxP];
  Float_t  convMVALikelihood_[maxP];
  // vertex quantities 
  Float_t  convChi2_[maxP];
  Float_t  convChi2Probability_[maxP];
  // per track quantities
  Float_t  convTk1Dz_[maxP];
  Float_t  convTk2Dz_[maxP];
  Float_t  convTk1DzErr_[maxP];
  Float_t  convTk2DzErr_[maxP];
  //Int_t    convTk1Nh_[maxP];
  //Int_t    convTk2Nh_[maxP];
  Int_t    convCh1Ch2_[maxP];
  Float_t  convTk1D0_[maxP];
  Float_t  convTk1Pout_[maxP];
  Float_t  convTk1Pin_[maxP];
  Float_t  convTk2D0_[maxP];
  Float_t  convTk2Pout_[maxP];
  Float_t  convTk2Pin_[maxP];

  // Physics objects handles
  Handle<std::vector<reco::GenParticle> > genParticlesHandle_;
  Handle<VertexCollection>      recVtxs_;
  Handle<VertexCollection>      recVtxsBS_;
  Handle<TriggerResults>        trgResultsHandle_;
  Handle<TriggerEvent>          triggerEvent_;
  Handle<TrackCollection>       tracksHandle_;
  Handle<GsfElectronCollection> gsfElectronHandle_;
  Handle<View<pat::MET> >       METHandle_;
  Handle<View<pat::MET> >       tcMETHandle_;
  Handle<View<pat::MET> >       pfMETHandle_;
  Handle<View<pat::MET> >       pfType01METHandle_;
  Handle<reco::PFMETCollection> recoPfMETHandle_;
  Handle<View<pat::Electron> >  electronHandle_;
  Handle<View<pat::Photon> >    photonHandle_;
  Handle<View<pat::Muon> >      muonHandle_;
  Handle<View<pat::Jet> >       jetHandle_;
  Handle<BeamSpot>              beamSpotHandle_;
  Handle<EcalRecHitCollection>  EBReducedRecHits_;
  Handle<EcalRecHitCollection>  EEReducedRecHits_;
  Handle<EcalRecHitCollection>  ESRecHits_;
  Handle<CaloTowerCollection>   towers_;
  Handle<PFCandidateCollection> pfAllCandidates; 
  Handle<PFCandidateCollection> pfCandidates; 
  Handle<PFCandidateCollection> pfCandidatePhotons;   
  Handle<PFClusterCollection>   pfClust_ecal;  
  Handle<PFClusterCollection>   pfClust_hcal;
  Handle<reco::ConversionCollection> convH_;
  Handle<reco::PhotonCollection> pfPhoTranslator_;
  Handle<reco::PhotonCollection> recoPhotonHandle_;
  //Handle<reco::SuperClusterCollection> pfPhoSClusters_;
  // Handle<reco::SuperClusterCollection> pfEleSClusters_;
  // Handle<reco::PhotonCoreCollection>pfPhotonCore_;
  // Handle<reco::GsfElectronCoreCollection>pfElectronCore_;
  std::vector<std::vector<Float_t> > phoPFIsoChargedVec;
  std::vector<Float_t> phoPFIsoNeutralVec;
  std::vector<Float_t> phoPFIsoPhotonVec;
  
  ggPFIsolation isolation;

};

#endif

