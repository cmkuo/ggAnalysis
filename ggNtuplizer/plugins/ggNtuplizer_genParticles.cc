#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "ggAnalysis/ggNtuplizer/interface/ggNtuplizer.h"
#include "ggAnalysis/ggNtuplizer/interface/GenParticleParentage.h"
using namespace std;

// (local) variables associated with tree branches
vector<float>    pdf_;
float            pthat_;
float            processID_;
float            genWeight_;
float            genHT_;
float            genPho1_;
float            genPho2_;
TString          EventTag_;
float            pdfWeight_;     
vector<float>    pdfSystWeight_;
Int_t            nLHE_;
vector<int>      lhePID;
vector<float>    lhePx;
vector<float>    lhePy;
vector<float>    lhePz; 
vector<float>    lheE;

Int_t            nPUInfo_;
vector<int>      nPU_;
vector<int>      puBX_;
vector<float>    puTrue_;

Int_t            nMC_;
vector<int>      mcPID;
vector<float>    mcVtx;
vector<float>    mcVty;
vector<float>    mcVtz;
vector<float>    mcPt;
vector<float>    mcMass;
vector<float>    mcEta;
vector<float>    mcPhi;
vector<float>    mcE;
vector<float>    mcEt;
vector<int>      mcGMomPID;
vector<int>      mcMomPID;
vector<float>    mcMomPt;
vector<float>    mcMomMass;
vector<float>    mcMomEta;
vector<float>    mcMomPhi;
//vector<int>      mcIndex;
vector<UShort_t> mcStatusFlag;
vector<int>      mcParentage;
vector<int>      mcStatus;
vector<float>    mcCalIsoDR03;
vector<float>    mcTrkIsoDR03;
vector<float>    mcCalIsoDR04;
vector<float>    mcTrkIsoDR04;

float getGenCalIso(edm::Handle<reco::GenParticleCollection> handle,
                   reco::GenParticleCollection::const_iterator thisPart,
                   float dRMax, bool removeMu, bool removeNu) {

  // Returns Et sum.
  float etSum = 0;

  for (reco::GenParticleCollection::const_iterator p = handle->begin(); p != handle->end(); ++p) {
    if (p == thisPart) continue;
    if (p->status() != 1) continue;

    // has to come from the same collision
    if (thisPart->collisionId() != p->collisionId())
        continue;

    int pdgCode = abs(p->pdgId());

    // skip muons/neutrinos, if requested
    if (removeMu && pdgCode == 13) continue;
    if (removeNu && (pdgCode == 12 || pdgCode == 14 || pdgCode == 16)) continue;
    
    // pass a minimum Et threshold
    // if (p->et() < 0) continue;

    // must be within deltaR cone
    float dR = reco::deltaR(thisPart->momentum(), p->momentum());
    if (dR > dRMax) continue;

    etSum += p->et();
  }

  return etSum;
}

float getGenTrkIso(edm::Handle<reco::GenParticleCollection> handle,
                   reco::GenParticleCollection::const_iterator thisPart, float dRMax) {

   // Returns pT sum without counting neutral particles.
   float ptSum = 0;

   for (reco::GenParticleCollection::const_iterator p = handle->begin(); p != handle->end(); ++p) {
      if (p == thisPart) continue;
      if (p->status() != 1) continue;
      if (p->charge() == 0) continue;  // do not count neutral particles

      // has to come from the same collision
      if (thisPart->collisionId() != p->collisionId())
         continue;

      // pass a minimum pt threshold
      // if (p->pt() < 0) continue;

      // must be within deltaR cone
      float dR = reco::deltaR(thisPart->momentum(), p->momentum());
      if (dR > dRMax) continue;

      ptSum += p->pt();
   }

   return ptSum;
}

void ggNtuplizer::branchesGenInfo(TTree* tree, edm::Service<TFileService> &fs) {

  tree->Branch("pdf",           &pdf_);
  tree->Branch("pthat",         &pthat_);
  tree->Branch("processID",     &processID_);
  tree->Branch("genWeight",     &genWeight_);
  tree->Branch("genHT",         &genHT_);
  tree->Branch("genPho1",       &genPho1_);
  tree->Branch("genPho2",       &genPho2_);
  if (dumpPDFSystWeight_) {
    tree->Branch("pdfWeight",     &pdfWeight_);
    tree->Branch("pdfSystWeight", &pdfSystWeight_);
  }
  tree->Branch("EventTag",      &EventTag_);

  tree->Branch("nPUInfo",       &nPUInfo_);
  tree->Branch("nPU",           &nPU_);
  tree->Branch("puBX",          &puBX_);
  tree->Branch("puTrue",        &puTrue_);
  tree->Branch("nLHE",          &nLHE_);
  tree->Branch("lhePID",        &lhePID);
  tree->Branch("lhePx",         &lhePx);
  tree->Branch("lhePy",         &lhePy); 
  tree->Branch("lhePz",         &lhePz); 
  tree->Branch("lheE",          &lheE); 

  hPU_        = fs->make<TH1F>("hPU",        "number of pileup",      200,  0, 200);
  hPUTrue_    = fs->make<TH1F>("hPUTrue",    "number of true pilepu", 1000, 0, 200);
  hGenWeight_ = fs->make<TH1F>("hGenWeight", "Gen weights",           2,    0, 2);
  hSumGenWeight_ = fs->make<TH1F>("hSumGenWeight", "Sum of Gen weights",1,  0, 1);
}

void ggNtuplizer::branchesGenPart(TTree* tree) {

  tree->Branch("nMC",          &nMC_);
  tree->Branch("mcPID",        &mcPID);
  tree->Branch("mcVtx",        &mcVtx);
  tree->Branch("mcVty",        &mcVty);
  tree->Branch("mcVtz",        &mcVtz);
  tree->Branch("mcPt",         &mcPt);
  tree->Branch("mcMass",       &mcMass);
  tree->Branch("mcEta",        &mcEta);
  tree->Branch("mcPhi",        &mcPhi);
  tree->Branch("mcE",          &mcE);
  tree->Branch("mcEt",         &mcEt);
  tree->Branch("mcGMomPID",    &mcGMomPID);
  tree->Branch("mcMomPID",     &mcMomPID);
  tree->Branch("mcMomPt",      &mcMomPt);
  tree->Branch("mcMomMass",    &mcMomMass);
  tree->Branch("mcMomEta",     &mcMomEta);
  tree->Branch("mcMomPhi",     &mcMomPhi);
  //tree->Branch("mcIndex",      &mcIndex);
  tree->Branch("mcStatusFlag", &mcStatusFlag); //-999:non W or Z, 1:hardronic, 2:e, 3:mu, 4:tau
  tree->Branch("mcParentage",  &mcParentage);  // 16*lepton + 8*boson + 4*non-prompt + 2*qcd + exotics
  tree->Branch("mcStatus",     &mcStatus);     // status of the particle
  tree->Branch("mcCalIsoDR03", &mcCalIsoDR03);
  tree->Branch("mcTrkIsoDR03", &mcTrkIsoDR03);
  tree->Branch("mcCalIsoDR04", &mcCalIsoDR04);
  tree->Branch("mcTrkIsoDR04", &mcTrkIsoDR04);
}

void ggNtuplizer::fillGenInfo(const edm::Event& e) {

  // cleanup from previous execution
  pthat_     = -99;
  processID_ = -99;
  genWeight_ = -99;
  genHT_     = -99;
  genPho1_   = -99;
  genPho2_   = -99;
  nPUInfo_   = 0;
  pdfWeight_ = -99;
  EventTag_  = "";
  pdf_          .clear();
  pdfSystWeight_.clear();
  nPU_          .clear();
  puBX_         .clear();
  puTrue_       .clear();
  lhePID        .clear();
  lhePx         .clear();
  lhePy         .clear();
  lhePz         .clear();
  lheE          .clear();

  edm::Handle<GenEventInfoProduct> genEventInfoHandle;
  e.getByToken(generatorLabel_, genEventInfoHandle);

  if (genEventInfoHandle.isValid()) {

    if (genEventInfoHandle->pdf()) {
      pdf_.push_back(genEventInfoHandle->pdf()->id.first);    // PDG ID of incoming parton #1
      pdf_.push_back(genEventInfoHandle->pdf()->id.second);   // PDG ID of incoming parton #2
      pdf_.push_back(genEventInfoHandle->pdf()->x.first);     // x value of parton #1
      pdf_.push_back(genEventInfoHandle->pdf()->x.second);    // x value of parton #2
      pdf_.push_back(genEventInfoHandle->pdf()->xPDF.first);  // PDF weight for parton #1
      pdf_.push_back(genEventInfoHandle->pdf()->xPDF.second); // PDF weight for parton #2
      pdf_.push_back(genEventInfoHandle->pdf()->scalePDF);    // scale of the hard interaction
    }

    if (genEventInfoHandle->hasBinningValues())
      pthat_ = genEventInfoHandle->binningValues()[0];
    processID_ = genEventInfoHandle->signalProcessID();
    genWeight_ = genEventInfoHandle->weight();
    if (genWeight_ >= 0) hGenWeight_->Fill(0.5);    
    else hGenWeight_->Fill(1.5);
    if (abs(genWeight_)>1) hSumGenWeight_->Fill(0.5,genWeight_/abs(genWeight_));
    else hSumGenWeight_->Fill(0.5,genWeight_);
  } else
    edm::LogWarning("ggNtuplizer") << "no GenEventInfoProduct in event";
  
  // access generator level HT  
  edm::Handle<LHEEventProduct> lheEventProduct;
  e.getByToken(lheEventLabel_, lheEventProduct);
  
  nLHE_ = 0;  
  double lheHt   = 0.;
  double lhePho1 = 0.;
  double lhePho2 = 0.;
  if (lheEventProduct.isValid()){
    const lhef::HEPEUP& lheEvent = lheEventProduct->hepeup();
    std::vector<lhef::HEPEUP::FiveVector> lheParticles = lheEvent.PUP;
    size_t numParticles = lheParticles.size();
    int nMCPho = 0;
    for ( size_t idxParticle = 0; idxParticle < numParticles; ++idxParticle ) {
      int absPdgId = TMath::Abs(lheEvent.IDUP[idxParticle]);
      int status = lheEvent.ISTUP[idxParticle];

      if (status == 1) {
	lhePID.push_back(lheEvent.IDUP[idxParticle]);
	lhePx .push_back(lheParticles[idxParticle][0]);
	lhePy .push_back(lheParticles[idxParticle][1]);
	lhePz .push_back(lheParticles[idxParticle][2]);
	lheE  .push_back(lheParticles[idxParticle][3]);
	nLHE_++;
      } 

      if (status == 1 && ((absPdgId >= 1 && absPdgId <= 6) || absPdgId == 21) ) { // quarks and gluons
	lheHt += TMath::Sqrt(TMath::Power(lheParticles[idxParticle][0], 2.) + TMath::Power(lheParticles[idxParticle][1], 2.)); // first entry is px, second py
      } 
      if (status == 1 && absPdgId == 22 && nMCPho == 0) { // first photon
	lhePho1 = TMath::Sqrt(TMath::Power(lheParticles[idxParticle][0], 2.) + TMath::Power(lheParticles[idxParticle][1], 2.));
	nMCPho++;
      }
      if (status == 1 && absPdgId == 22 && nMCPho == 1) { // first photon
	lhePho2 = TMath::Sqrt(TMath::Power(lheParticles[idxParticle][0], 2.) + TMath::Power(lheParticles[idxParticle][1], 2.));
	nMCPho++;
      }

      typedef std::vector<std::string>::const_iterator comments_const_iterator;

      comments_const_iterator c_begin = lheEventProduct->comments_begin();
      comments_const_iterator c_end = lheEventProduct->comments_end();

      TString model_params;
      for(comments_const_iterator cit=c_begin; cit!=c_end; ++cit) {
	size_t found = (*cit).find("model");
	if(found != std::string::npos)   { 
	  model_params = *cit;
	}
      }
      EventTag_ = model_params;
    }

    if (dumpPDFSystWeight_) {
      pdfWeight_ = lheEventProduct->originalXWGTUP(); // PDF weight of this event !
      for (unsigned i = 0; i < lheEventProduct->weights().size(); ++i) {
	pdfSystWeight_.push_back(lheEventProduct->weights()[i].wgt);
      }
    }
  }
  genHT_   = lheHt;  
  genPho1_ = lhePho1;  
  genPho2_ = lhePho2;  

  edm::Handle<vector<PileupSummaryInfo> > genPileupHandle;
  e.getByToken(puCollection_, genPileupHandle);
  
  if (genPileupHandle.isValid()) {
    for (vector<PileupSummaryInfo>::const_iterator pu = genPileupHandle->begin(); pu != genPileupHandle->end(); ++pu) {
      if (pu->getBunchCrossing() == 0) {
        hPU_->Fill(pu->getPU_NumInteractions());
        hPUTrue_->Fill(pu->getTrueNumInteractions());
      }

      nPU_   .push_back(pu->getPU_NumInteractions());
      puTrue_.push_back(pu->getTrueNumInteractions());
      puBX_  .push_back(pu->getBunchCrossing());

      nPUInfo_++;
    }
  }
  else
    edm::LogWarning("ggNtuplizer") << "no PileupSummaryInfo in event";

}

void ggNtuplizer::fillGenPart(const edm::Event& e) {

  mcPID       .clear();
  mcVtx       .clear();
  mcVty       .clear();
  mcVtz       .clear();
  mcPt        .clear();
  mcMass      .clear();
  mcEta       .clear();
  mcPhi       .clear();
  mcE         .clear();
  mcEt        .clear();
  mcGMomPID   .clear();
  mcMomPID    .clear();
  mcMomPt     .clear();
  mcMomMass   .clear();
  mcMomEta    .clear();
  mcMomPhi    .clear();
  //mcIndex     .clear();
  mcStatusFlag.clear();
  mcParentage .clear();
  mcStatus    .clear();
  mcCalIsoDR03.clear();
  mcTrkIsoDR03.clear();
  mcCalIsoDR04.clear();
  mcTrkIsoDR04.clear();

  nMC_ = 0;

  edm::Handle<vector<reco::GenParticle> > genParticlesHandle;
  e.getByToken(genParticlesCollection_, genParticlesHandle);

  if (!genParticlesHandle.isValid()) {
    edm::LogWarning("ggNtuplizer") << "no reco::GenParticles in event";
    return;
  }

  int genIndex = 0;

  for (vector<reco::GenParticle>::const_iterator ip = genParticlesHandle->begin(); ip != genParticlesHandle->end(); ++ip) {
    genIndex++;

    int status = ip->status();
    //bool stableFinalStateParticle = status == 1 && ip->pt() > 5.0;
    
    //bool quarks = abs(ip->pdgId())<7;

    // keep non-FSR photons with pT > 5.0 and all leptons with pT > 3.0;
    bool photonOrLepton =
      (ip->pdgId() == 22 && (ip->isPromptFinalState() || ip->isLastCopy())) ||
      (status == 1 && abs(ip->pdgId()) == 11 && (ip->isPromptFinalState() || ip->isLastCopy())) || 
      (status == 1 && abs(ip->pdgId()) == 13 && (ip->isPromptFinalState() || ip->isLastCopy())) ||
      (status == 1 && (abs(ip->pdgId()) == 12 || abs(ip->pdgId()) == 14 || abs(ip->pdgId()) == 16)) ||
      (status == 1 && ( abs(ip->pdgId()) >= 11 && abs(ip->pdgId()) <= 16 ) && ip->pt() > 3.0)  ||
      (status < 10 && abs(ip->pdgId()) == 15 && ip->pt() > 3.0);
      
    // select also Z, W, H, top and b 
    bool heavyParticle =
      ((    ip->pdgId()  == 23 && ip->isLastCopy()) || 
       (abs(ip->pdgId()) == 24 && ip->isLastCopy()) || 
       (    ip->pdgId()  == 25 && ip->isLastCopy()) ||
       (abs(ip->pdgId()) ==  6 && ip->isLastCopy()) || 
       (abs(ip->pdgId()) ==  5 && ip->isLastCopy()));
    
    bool newParticle = false;
    for (size_t inp = 0; inp < newparticles_.size(); ++inp) {
      if (abs(ip->pdgId()) == newparticles_[inp]) newParticle = true;
    }
    
    //if ( heavyParticle || photonOrLepton || quarks || newParticle ) {
    if ( heavyParticle || photonOrLepton || newParticle ) {
      
      const reco::Candidate *p = (const reco::Candidate*)&(*ip);
      if (!runOnParticleGun_ && !p->mother()) continue;

      mcPID    .push_back(p->pdgId());
      mcVtx    .push_back(p->vx());
      mcVty    .push_back(p->vy());
      mcVtz    .push_back(p->vz());
      mcPt     .push_back(p->pt());
      mcMass   .push_back(p->mass());
      mcEta    .push_back(p->eta());
      mcPhi    .push_back(p->phi());
      mcE      .push_back(p->energy());
      mcEt     .push_back(p->et());
      mcStatus .push_back(p->status());
	
      UShort_t tmpStatusFlag = 0;
      if (ip->fromHardProcessFinalState()) setbit(tmpStatusFlag, 0);
      if (ip->isPromptFinalState())        setbit(tmpStatusFlag, 1);
      if (ip->isHardProcess())             setbit(tmpStatusFlag, 2);

      // if genParticle is W or Z, check its decay type
      if ( ip->pdgId() == 23 || abs(ip->pdgId()) == 24 ) {
        for (size_t k=0; k < p->numberOfDaughters(); ++k) {
          const reco::Candidate *dp = p->daughter(k);
          if (abs(dp->pdgId())<=6)                               setbit(tmpStatusFlag, 4);
          else if (abs(dp->pdgId())==11 || abs(dp->pdgId())==12) setbit(tmpStatusFlag, 5);
          else if (abs(dp->pdgId())==13 || abs(dp->pdgId())==14) setbit(tmpStatusFlag, 6);
          else if (abs(dp->pdgId())==15 || abs(dp->pdgId())==16) setbit(tmpStatusFlag, 7);
        }
      }
      //if (ip->isLastCopy())  setbit(tmpStatusFlag, 8);
      mcStatusFlag.push_back(tmpStatusFlag);

      int mcGMomPID_ = -999;
      int mcMomPID_  = -999;
      float mcMomPt_    = -999.;
      float mcMomMass_  = -999.;
      float mcMomEta_   = -999.;
      float mcMomPhi_   = -999.;
      if (!runOnSherpa_) {
	
	reco::GenParticleRef partRef = reco::GenParticleRef(genParticlesHandle,
							    ip-genParticlesHandle->begin());
	genpartparentage::GenParticleParentage particleHistory(partRef);
	
	mcParentage.push_back(particleHistory.hasLeptonParent()*16   +
			      particleHistory.hasBosonParent()*8     +
			      particleHistory.hasNonPromptParent()*4 +
			      particleHistory.hasQCDParent()*2       +
			      particleHistory.hasExoticParent());      
	
	if ( particleHistory.hasRealParent() ) {
	  reco::GenParticleRef momRef = particleHistory.parent();
	  if ( momRef.isNonnull() && momRef.isAvailable() ) {
	    mcMomPID_  = momRef->pdgId();
	    mcMomPt_   = momRef->pt();
	    mcMomMass_ = momRef->mass();
	    mcMomEta_  = momRef->eta();
	    mcMomPhi_  = momRef->phi();
	    
	    // get Granny
	    genpartparentage::GenParticleParentage motherParticle(momRef);
	    if ( motherParticle.hasRealParent() ) {
	      reco::GenParticleRef granny = motherParticle.parent();
	      mcGMomPID_ = granny->pdgId();
	    }
	  }
	}
	mcGMomPID.push_back(mcGMomPID_);
	mcMomPID.push_back(mcMomPID_);
	mcMomPt.push_back(mcMomPt_);
	mcMomMass.push_back(mcMomMass_);
	mcMomEta.push_back(mcMomEta_);
	mcMomPhi.push_back(mcMomPhi_);
      }

      //mcIndex.push_back(genIndex-1);

      if (photonOrLepton) {
	mcCalIsoDR03.push_back( getGenCalIso(genParticlesHandle, ip, 0.3, false, false) );
	mcTrkIsoDR03.push_back( getGenTrkIso(genParticlesHandle, ip, 0.3) );
	mcCalIsoDR04.push_back( getGenCalIso(genParticlesHandle, ip, 0.4, false, false) );
	mcTrkIsoDR04.push_back( getGenTrkIso(genParticlesHandle, ip, 0.4) );
      } else {
	mcCalIsoDR03.push_back( -999. );
	mcTrkIsoDR03.push_back( -999. );
	mcCalIsoDR04.push_back( -999. );
	mcTrkIsoDR04.push_back( -999. );
      }

      nMC_++;
    } // save info on particles of interest
  } // loop over gen-level particles

}
