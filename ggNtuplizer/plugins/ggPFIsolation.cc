#include "ggAnalysis/ggNtuplizer/interface/ggPFIsolation.h"
#include "DataFormats/Math/interface/deltaR.h"


#include"TMath.h"
#include <iostream>

ggPFIsolation::ggPFIsolation(){

  
  //initialize TMVA for leading photon ID
  tmvaReaderID = new TMVA::Reader("!Color:Silent");
  tmvaReaderID->AddVariable( "isoBin01",&isoBin01);
  tmvaReaderID->AddVariable( "isoBin02",&isoBin02);
  tmvaReaderID->AddVariable( "isoBin03",&isoBin03);
  tmvaReaderID->AddVariable( "isoBin04",&isoBin04);
  tmvaReaderID->AddVariable( "isoBin05",&isoBin05);
  tmvaReaderID->AddVariable( "isoBin06",&isoBin06);
  tmvaReaderID->AddVariable( "isoBin07",&isoBin07);
  tmvaReaderID->AddVariable( "isoBin08",&isoBin08);

  tmvaReaderID->AddVariable( "rho25",&rho25);
  tmvaReaderID->AddVariable( "PFPho_excl_Et",&fMustache );

  tmvaReaderID->AddSpectator("phoPFIso := isoBin02+isoBin03+isoBin04",&phoPFIsoBDT);
  tmvaReaderID->AddSpectator( "phoPt",&phoPtBDT );
  tmvaReaderID->AddSpectator( "phoEta",&phoEtaBDT );
  tmvaReaderID->AddSpectator( "phoR9",&phoR9BDT );

 

  //tmvaReaderID->BookMVA("BDT","/afs/cern.ch/user/v/vchetlur/public/TMVAClassification_BDT.weights.xml");
 
}

ggPFIsolation::~ggPFIsolation() {
}


float ggPFIsolation::mvaID(const reco::PFCandidateCollection* pfParticlesColl,const reco::Photon *recoPhoton ,edm::Handle< reco::VertexCollection > recoVtx){

  //Position vector of the photon SC with respect to 0,0,0
  math::XYZVector phtotonSC(recoPhoton->caloPosition().x(),recoPhoton->caloPosition().y(),recoPhoton->caloPosition().z()); 

  //Find a PF particle(photon or a electron) match to the reco::Photon
  Int_t iMatch=-1;
  float_t fDeltaRClosest = 999;
  for(unsigned iPF=0; iPF<pfParticlesColl->size(); iPF++) {
    const reco::PFCandidate& pfParticle= (*pfParticlesColl)[iPF]; 
    if( !((pfParticle.pdgId()==22||abs(pfParticle.pdgId())==11)))
	continue;
    math::XYZVector photonWRTPFParticleVtx(recoPhoton->caloPosition().x()-pfParticle.vx(),recoPhoton->caloPosition().y()-pfParticle.vy(),recoPhoton->caloPosition().z()-pfParticle.vz());
    
    Float_t fDeltaR =  deltaR(pfParticle.eta(),pfParticle.phi(),photonWRTPFParticleVtx.Eta(),photonWRTPFParticleVtx.Phi());
    if(fDeltaR<fDeltaRClosest){
      iMatch = iPF;
      fDeltaRClosest = fDeltaR;
    }
   
  } 
  
  fMustache = -1;

  if(!(iMatch>=0&&fDeltaRClosest<0.1))
    return -10;
  
  fMustache = getMustacheValue( (*pfParticlesColl)[iMatch]);
  
  math::XYZPoint photon(recoPhoton->caloPosition().x(),recoPhoton->caloPosition().y(),recoPhoton->caloPosition().z());
  
  
  vector<float > fIsolationGamma = getPFIsolationNeutral(pfParticlesColl,0,photon,true);
  vector<float > fIsolationNeutral = getPFIsolationNeutral(pfParticlesColl,1,photon,false);
  vector<vector<float > > fIsolationCharged = getPFIsolationCharged(pfParticlesColl,2,photon,true,recoVtx);

  
  //code to initialize TMVA and get the ID
 
  return 1.0;
}

float ggPFIsolation::mvaIDWithIsodeps(const reco::PFCandidateCollection* pfParticlesColl,const reco::Photon *recoPhoton, reco::PFCandidatePtr PFCandPtr, IsoDepositMaps photonIsoDep, float fastJetCorr){
  //Disabling the 44X code 

/*
  math::XYZVector phtotonSC(recoPhoton->caloPosition().x(),recoPhoton->caloPosition().y(),recoPhoton->caloPosition().z()); 
  
  //Find a PF particle(photon or a electron) match to the reco::Photon
  Int_t iMatch=-1;
  float_t fDeltaRClosest = 999;
  for(unsigned iPF=0; iPF<pfParticlesColl->size(); iPF++) {
    const reco::PFCandidate& pfParticle= (*pfParticlesColl)[iPF]; 
    if( !((pfParticle.pdgId()==22||abs(pfParticle.pdgId())==11)))
	continue;
    math::XYZVector photonWRTPFParticleVtx(recoPhoton->caloPosition().x()-pfParticle.vx(),recoPhoton->caloPosition().y()-pfParticle.vy(),recoPhoton->caloPosition().z()-pfParticle.vz());
    
    Float_t fDeltaR =  deltaR(pfParticle.eta(),pfParticle.phi(),photonWRTPFParticleVtx.Eta(),photonWRTPFParticleVtx.Phi());
    if(fDeltaR<fDeltaRClosest){
      iMatch = iPF;
      fDeltaRClosest = fDeltaR;
    }
   
  } 
  
  fMustache = -1;

  if(!(iMatch>=0&&fDeltaRClosest<0.1))
    return -10;
  
  fMustache = getMustacheValue( (*pfParticlesColl)[iMatch]);


  vector<float > fIsolationIsodeps = getIsolationFromIsodeps(PFCandPtr,photonIsoDep);
  
 
  isoBin02 = fIsolationIsodeps [1];
  isoBin03 = fIsolationIsodeps [2];
  isoBin04 = fIsolationIsodeps [3];
  isoBin05 = fIsolationIsodeps [4];
  isoBin06 = fIsolationIsodeps [5];
  isoBin07 = fIsolationIsodeps [6];
  isoBin08 = fIsolationIsodeps [7];

  phoPtBDT = recoPhoton->et();
  phoEtaBDT = recoPhoton->eta();
  phoR9BDT = recoPhoton->r9();

  rho25 = fastJetCorr;
 
  return  tmvaReaderID->EvaluateMVA("BDT");
*/

  return -1;
}


float ggPFIsolation::getMustacheValue(reco::PFCandidate PFCand){
  
 
  etmax = 0;
  imax = -1;
  phot_eta_maxcl = 0.0;
  phot_phi_maxcl = 0.0;

  Float_t fExcluded = 0.0;

  float deta, dphi;
  float upper_cut, lower_cut;
  float b_upper, b_lower;
  float a_upper, a_lower;
  float curv_low, curv_up;
  float midpoint;
 
  PFCandidate::ElementsInBlocks eleInBlocks = PFCand.elementsInBlocks();
   
  Int_t nClust =0;
   
  for(unsigned i=0; i<eleInBlocks.size(); i++){
    
    PFBlockRef blockRef = eleInBlocks[i].first;
    unsigned indexInBlock = eleInBlocks[i].second;
    PFBlock::LinkData linkData    = blockRef->linkData();	      
    const edm::OwnVector< reco::PFBlockElement >&  elements=eleInBlocks[i].first->elements();  
    const reco::PFBlockElement& element = elements[indexInBlock];
    if(element.type()==reco::PFBlockElement::ECAL){
      PFPho_clusteta_[nClust]=element.clusterRef()->position().eta();
      PFPho_clustphi_[nClust]=element.clusterRef()->position().phi();
      PFPho_clustEt_[nClust]=element.clusterRef()->pt();
      nClust++;
    }
    
  }
  
  
  //loop over clusters
  for(int k=0; k<nClust; k++){
    
    //search for highest Et cluster, set phi and eta
    if(etmax < PFPho_clustEt_[k]){
      imax = k;
      etmax = PFPho_clustEt_[k];
      phot_eta_maxcl = PFPho_clusteta_[k];
      phot_phi_maxcl = PFPho_clustphi_[k];
    }//end search for highest Et cluster
    
  }//end loop over clusters
  
  for(int k=0; k<nClust; k++){
    deta = 0.0; //= sin(photeta[0])*(clust_eta-photeta[0])
    dphi = 0.0; //= (clust_phi-photphi[0])
    upper_cut = 0.0;
    lower_cut = 0.0;
    b_upper = 0.0;
    b_lower = 0.0;
    a_upper = 0.0;
    a_lower = 0.0;
    curv_low = 0.0;
    curv_up = 0.0;
    midpoint = 0.0;  

    	
    deta = sin(phot_eta_maxcl)*(PFPho_clusteta_[k]-phot_eta_maxcl);	
    dphi = PFPho_clustphi_[k]-phot_phi_maxcl;
    
    //2 parabolas (upper and lower) 
    //of the form: y = a*x*x + b      
    
    //b comes from a fit to the width
    //and has a slight dependence on Et on the upper edge
    b_lower = w00*sin(phot_eta_maxcl)*phot_eta_maxcl + w01 / sqrt(log10(PFPho_clustEt_[k])+1.1);
    b_upper = w10*sin(phot_eta_maxcl)*phot_eta_maxcl + w11  / sqrt(log10(PFPho_clustEt_[k])+1.1);
    
    //here make an adjustment to the width for the offset from 0.
    midpoint = b_upper - (b_upper-b_lower)/2.;
    b_lower = b_lower - midpoint;
    b_upper = b_upper - midpoint;
    
    //the curvature comes from a parabolic 
    //fit for many slices in eta given a 
    //slice -0.1 < log10(Et) < 0.1
    curv_up = p00*pow(phot_eta_maxcl*sin(phot_eta_maxcl),2)+p01*phot_eta_maxcl*sin(phot_eta_maxcl)+p02;
    curv_low = p10*pow(phot_eta_maxcl*sin(phot_eta_maxcl),2)+p11*phot_eta_maxcl*sin(phot_eta_maxcl)+p12;
    
    //solving for the curviness given the width of this particular point
    a_lower = (1/(4*curv_low))-fabs(b_lower);
    a_upper = (1/(4*curv_up))-fabs(b_upper);
    
    upper_cut =(1./(4.*a_upper))*pow(dphi,2)+b_upper;
    lower_cut =(1./(4.*a_lower))*pow(dphi,2)+b_lower;
    
  
    if (!(deta < upper_cut && deta > lower_cut)){
      fExcluded += PFPho_clustEt_[k];
    }
    
  }
  return fExcluded;
}

//Function to calculate the Isolation for neutral particles. It returns the isolation in 8 rings around the photon candidate.
vector<float > ggPFIsolation::getPFIsolationNeutral(const reco::PFCandidateCollection* pfParticlesColl,int iType, math::XYZPoint photon,bool doVeto){

   
  vector<float> fIsolation;

  for(int isoBin =0;isoBin<8;isoBin++){
    float fTemp = 0.0;
    fIsolation.push_back(fTemp);
  }
  
  for(unsigned iPF=0; iPF<pfParticlesColl->size(); iPF++) {

    const reco::PFCandidate& pfParticle= (*pfParticlesColl)[iPF]; 

    if(!isType(iType,pfParticle.pdgId()))
   	continue;

    math::XYZVector photonWRTPFParticleVtx(photon.x()-pfParticle.vx(),photon.y()-pfParticle.vy(),photon.z()-pfParticle.vz());
        
    float fDeltaR = deltaR(pfParticle.eta(),pfParticle.phi(),photonWRTPFParticleVtx.Eta(),photonWRTPFParticleVtx.Phi()); 
    
    if(fDeltaR>0.8)
      continue;

    if(doVeto){
      if(abs(photonWRTPFParticleVtx.Eta()-pfParticle.eta())<0.025&&deltaPhi(photonWRTPFParticleVtx.Phi(),pfParticle.phi())<0.4)
	continue;
    }
    
    int iBin = (int)(fDeltaR * 10.);
    
    if(iBin<8)
      fIsolation[iBin]+=pfParticle.pt();
    
  }
  
  return fIsolation;
}


//Function to calculate the Isolation for charged particles. It returns the isolation in 8 rings around the photon candidate for a collection of input vertecies.
vector<vector<float > > ggPFIsolation::getPFIsolationCharged(const reco::PFCandidateCollection* pfParticlesColl,int iType, math::XYZPoint photon,bool doVeto,edm::Handle< reco::VertexCollection > recVtxs_){
 
  
  vector<vector<float > > fIsolationWithRespectToVertex;
  
  //chack if the input vertex collection looks reasonable
  if(!recVtxs_.isValid()){
    cout<<"Invalid Vertex collection"<<endl;
  }
  
  //Read the vertex collection
  const reco::VertexCollection& vertices = *(recVtxs_.product());


  //Loop through different verticies
  for(reco::VertexCollection::const_iterator vertex=vertices.begin(); vertex!=vertices.end(); ++vertex) {
  
    vector<float> fIsolation;


    ///initialize the isolation vector 
    for(int isoBin =0;isoBin<8;isoBin++){
      float fTemp = 0.0;
      fIsolation.push_back(fTemp);
    }
    
 
    //Associate the photon to the given vertex
    math::XYZVector photonWrtVertex(photon.x() - vertex->x(), 
				    photon.y() - vertex->y(), 
				    photon.z() - vertex->z());
    
    //Loop through the photons
    for(unsigned iPF=0; iPF<pfParticlesColl->size(); iPF++) {
      
      const reco::PFCandidate& pfParticle= (*pfParticlesColl)[iPF]; 
    
      if(!isType(iType,pfParticle.pdgId()))
	continue;
     
      math::XYZPoint pfParticleVtx(pfParticle.momentum());
      
    
      float dz = fabs(pfParticle.vz() - vertex->z());
      if (dz > 1.)
	continue;
      
      double dxy = ( -(pfParticle.vx() - vertex->x())*pfParticle.py() + (pfParticle.vy() - vertex->y())*pfParticle.px()) / pfParticle.pt();
      if(fabs(dxy) > 0.1)
	continue;

  
      float fDeltaR = deltaR(pfParticle.eta(),pfParticle.phi(),photonWrtVertex.eta(),photonWrtVertex.phi()); 
      
      if(fDeltaR>0.8)
	continue;

      if(doVeto){
	if(abs(photonWrtVertex.eta()-pfParticle.eta())<0.025&&deltaPhi(photonWrtVertex.phi(),pfParticle.phi())<0.05)
	  continue;
      }
    
      
      int iBin = (int)(fDeltaR * 10.);

      if(iBin<8)
	fIsolation[iBin]+=pfParticle.pt();
    }

    
    fIsolationWithRespectToVertex.push_back(fIsolation);
  }

   
  return fIsolationWithRespectToVertex;
}   
   	 



vector<float > ggPFIsolation::getIsolationFromIsodeps(reco::PFCandidatePtr PFCand,IsoDepositMaps photonIsoDep){

  vector<float >  fIsolationWithIsoDeps;

  for(int isoBin =0;isoBin<8;isoBin++){
    float fTemp = 0.0;
    fIsolationWithIsoDeps.push_back(fTemp);
  }
  
  
  edm::Handle<edm::ValueMap<reco::PFCandidatePtr> > photonValMapH;
  const edm::ValueMap<reco::PFCandidatePtr> & myPhotonValMap(*photonValMapH); 
 
  unsigned nIsoDepTypes=photonIsoDep.size(); 

  for(unsigned ideptype=0; ideptype<nIsoDepTypes;++ideptype) {
    const reco::IsoDeposit & isoDep((*photonIsoDep[ideptype])[PFCand]);
    typedef reco::IsoDeposit::const_iterator IM;
    for(IM im=isoDep.begin(); im != isoDep.end(); ++im) {
      int iBin = (int)(im->dR()* 10.);
      if(iBin<8)
	fIsolationWithIsoDeps[iBin]+=im->value();
      
    }
  }
  return fIsolationWithIsoDeps;
}

float ggPFIsolation::CiCPhotonIDPF4pfEcalIso(edm::Handle<reco::PFCandidateCollection> pfHandle, reco::PhotonRef localPho, float dRmax, float dRVetoBarrel, float dRVetoEndcap, float etaStrip, float phiStrip, float energyBarrel, float energyEndcap, std::vector<reco::PFCandidate::ParticleType> pVetoes, bool cleanfoot){
  float sum = 0;
  float dRVeto;
  const reco::PFCandidateCollection* forIsolation = pfHandle.product();
  for(unsigned i=0; i<forIsolation->size(); i++) {
    
    const reco::PFCandidate& pfc = (*forIsolation)[i];
    
    bool process = false;
    for (std::vector<reco::PFCandidate::ParticleType>::const_iterator it = pVetoes.begin();
	 it != pVetoes.end(); ++it) {
      if (pfc.particleId() == *it) {
	process = true;
	break;
      }
    }
    if(cleanfoot){
      if(pfc.mva_nothing_gamma()>0.01)process=false;
    }
    if (process) {
      if (fabs(pfc.momentum().Eta()) < 1.479) {
	dRVeto = dRVetoBarrel;
	//if (fabs(pfc.pt()) < energyBarrel)
	//continue;
      } else {
	dRVeto = dRVetoEndcap;
	etaStrip=0; //no eta strip applied in the endcap
	//if (fabs(pfc.energy()) < energyEndcap)
	//continue;
      }
      math::XYZVector pfvtx(pfc.vx(), pfc.vy(), pfc.vz()); 
      math::XYZVector vCand;
      vCand = math::XYZVector(localPho->superCluster()->x()-pfvtx.X(), localPho->superCluster()->y()-pfvtx.Y(), localPho->superCluster()->z()-pfvtx.Z());
      
      //float r = vCand.R();
      //math::XYZVector pvm((pfc.momentum()*r/pfc.momentum().R()) + pfvtx);
      
      float dR = deltaR(vCand.Eta(), vCand.Phi(), pfc.momentum().Eta(), pfc.momentum().Phi());
      float dEta = fabs(vCand.Eta() - pfc.momentum().Eta());
      //double dPhi = fabs(vCand.Phi() - pvm.Phi());
      //if(dPhi >TMath::Pi())
      //dPhi = TMath::TwoPi() - dPhi;
      
      if (dEta < etaStrip)
	continue;
      
      //if (dPhi < phiStrip)
      //continue;
      
      if(dR > dRmax || dR < dRVeto)
	continue;
      sum += pfc.pt();
    }
  }
      
  return sum;
}  

std::vector<float> ggPFIsolation::CiCPhotonPF4pfTkIsoWithVertex(edm::Handle<reco::PFCandidateCollection> pfHandle,
						    edm::Handle<reco::VertexCollection> vtxHandle,
						    reco::PhotonRef localPho, float dRmax, float dRvetoBarrel, float dRvetoEndcap, 
						  float ptMin, float dzMax, float dxyMax,
						  std::vector<reco::PFCandidate::ParticleType> pVetoes) {
  
  float dRveto;
  if (localPho->isEB())
    dRveto = dRvetoBarrel;
  else
    dRveto = dRvetoEndcap;

  std::vector<float> result;
  const reco::PFCandidateCollection* forIsolation = pfHandle.product();

  for(unsigned int ivtx=0; ivtx<vtxHandle->size(); ++ivtx) {
    
    reco::VertexRef vtx(vtxHandle, ivtx);
    math::XYZVector vCand(localPho->superCluster()->x() - vtx->x(), 
			  localPho->superCluster()->y() - vtx->y(), 
			  localPho->superCluster()->z() - vtx->z());
    
    float sum = 0;
    for(unsigned i=0; i<forIsolation->size(); i++) {
    
      const reco::PFCandidate& pfc = (*forIsolation)[i];
      
      bool process = false;
      for (std::vector<reco::PFCandidate::ParticleType>::const_iterator it = pVetoes.begin();
	   it != pVetoes.end(); ++it) {
	if (pfc.particleId() == *it) {
	  process = true;
	  break;
	}
      }
      
      if (process) {
	if (pfc.pt() < ptMin)
	  continue;
	
	//float dz = fabs(pfc.vz() - vtx->z());
	float dz=pfc.trackRef()->dz(vtx->position());
	if (dz > dzMax)
	  continue;
	
	//	double dxy = ( -(pfc.vx() - vtx->x())*pfc.py() + (pfc.vy() - vtx->y())*pfc.px()) / pfc.pt();
	float dxy = fabs(pfc.trackRef()->dxy(vtx->position()));
	if(fabs(dxy) > dxyMax)
	  continue;
	
	math::XYZVector pvi(pfc.momentum());
	float dR = deltaR(vCand.Eta(), vCand.Phi(), pvi.Eta(), pvi.Phi());
	
	if(dR > dRmax || dR < dRveto)
	  continue;

	sum += pfc.pt();
      }
    }

    result.push_back(sum);
  }
  
  return result;
}

bool ggPFIsolation::isType(int iType,int iPdgID){

  
  if(iType==0){
    if(iPdgID==22)
      return true;
    else
      return false;
  }


  if(iType==1){
    if((abs(iPdgID))==130)
      return true;
    else
      return false;
  }        


  if(iType==2){
    if((abs(iPdgID)==211)||(abs(iPdgID)==11))
      return true;
    else
      return false;
  }
	   
  
  return false;

}


