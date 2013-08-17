#ifndef ggPFIsolation_h
#define ggPFIsolation_h

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PileUpPFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PileUpPFCandidateFwd.h"

#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"

#include "DataFormats/ParticleFlowReco/interface/PFBlock.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFBlockElementSuperCluster.h"

#include "DataFormats/RecoCandidate/interface/IsoDeposit.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/ValueMap.h"


#include "TMVA/Tools.h"
#include "TMVA/Reader.h"



#include <memory>
#include <fstream>
#include <map>

using namespace edm;
using namespace std;
using namespace reco;



const  float p00 = -0.107537;//-0.569806;
const  float p01 = 0.590969;//1.45503;

const float p02 = -0.076494;//-0.0898262;
const float p10 = -0.0268843;//-0.142451;
const float p11 = 0.147742;//0.363757;
const float p12 = -0.0191235;//-0.0224565;
  
const float w00 = -0.00571429;
const float w01 = -0.002;
const float w10 = 0.0135714;
const float w11 = 0.001;
const Int_t maxPho = 1000;

class ggPFIsolation  {

public:

  explicit ggPFIsolation();
  virtual ~ggPFIsolation();


  typedef std::vector< edm::Handle< edm::ValueMap<reco::IsoDeposit> > > IsoDepositMaps;

  virtual float mvaID(const reco::PFCandidateCollection* pfParticlesColl,const reco::Photon *, edm::Handle< reco::VertexCollection >);
  virtual float mvaIDWithIsodeps(const reco::PFCandidateCollection* pfParticlesColl,const reco::Photon *,reco::PFCandidatePtr PFCandPtr, IsoDepositMaps photonIsoDep, float );
  virtual float getMustacheValue(reco::PFCandidate );
  virtual vector<float >getPFIsolationNeutral(const reco::PFCandidateCollection* pfParticlesColl,int iPgdID,math::XYZPoint photon,bool doVeto);
  virtual vector<vector<float > > getPFIsolationCharged(const reco::PFCandidateCollection* pfParticlesColl,int iPgdID, math::XYZPoint photon,bool doVeto,edm::Handle< reco::VertexCollection > recVtxs_);
  virtual bool isType(int,int);
  virtual vector<float >  getIsolationFromIsodeps(reco::PFCandidatePtr PFCand,IsoDepositMaps photonIsoDep); 
  virtual float CiCPhotonIDPF4pfEcalIso(edm::Handle<reco::PFCandidateCollection> pfHandle, reco::PhotonRef localPho, float dRmax, float dRVetoBarrel, float dRVetoEndcap, float etaStrip, float phiStrip, float energyBarrel, float energyEndcap, std::vector<reco::PFCandidate::ParticleType> pVetoes, bool cleanfoot);
  
 virtual std::vector<float> CiCPhotonPF4pfTkIsoWithVertex(
							  edm::Handle<reco::PFCandidateCollection> pfHandle, edm::Handle<reco::VertexCollection> vtxHandle,
							  reco::PhotonRef localPho, float dRmax, float dRvetoBarrel, float dRvetoEndcap, 
							  float ptMin, float dzMax, float dxyMax,
							  std::vector<reco::PFCandidate::ParticleType> pVetoes);
    //==== Parameters for Mustache ID  =====================
 
  int imax; //=-1  //imax for each photon
  float etmax;// = 0; //etmax for each photon
  float phot_phi_maxcl; //phi for clust w/in a photon with max Et
  float phot_eta_maxcl; //eta for clust w/in a photon with  max Et
  
  Float_t  PFPho_clusteta_[20];  
  Float_t  PFPho_clustphi_[20];  
  Float_t  PFPho_clustEt_[20];


   

   //---------------- Variables for isolation cone --------------------------------
   TMVA::Reader *tmvaReaderID;
   Float_t isoBin01;
   Float_t isoBin02;
   Float_t isoBin03;
   Float_t isoBin04;
   Float_t isoBin05;
   Float_t isoBin06;
   Float_t isoBin07;
   Float_t isoBin08;

   Float_t SLisoBin01;
   Float_t SLisoBin02;
   Float_t SLisoBin03;
   Float_t SLisoBin04;
   Float_t SLisoBin05;
   Float_t SLisoBin06;
   Float_t SLisoBin07;
   Float_t SLisoBin08;

   Float_t rho25;

   Float_t phoPFIsoBDT,phoPtBDT,phoEtaBDT,phoR9BDT;
   Float_t fMustache;
};

#endif
