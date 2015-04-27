#include "ggAnalysis/ggNtuplizer/interface/ggNtuplizer.h"
#include "ggAnalysis/ggNtuplizer/interface/GenParticleParentage.h"

using namespace std;

// (local) variables associated with tree branches
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

float getGenCalIso(edm::Handle<reco::GenParticleCollection> handle, reco::GenParticleCollection::const_iterator thisPho, const Float_t dRMax, bool removeMu, bool removeNu) {

  const Float_t etMin = 0.0;
  Float_t genCalIsoSum = 0.0;
  if (!handle.isValid()) return genCalIsoSum;

  for (reco::GenParticleCollection::const_iterator it_gen=handle->begin(); it_gen!=handle->end(); ++it_gen) {

    if (it_gen == thisPho) continue;        // can't be the original photon
    if (it_gen->status() != 1) continue;    // need to be a stable particle
    if (thisPho->collisionId() != it_gen->collisionId()) continue; // has to come from the same collision

    Int_t pdgCode = abs(it_gen->pdgId());
    // we should not count neutrinos, muons
    if (removeMu && pdgCode == 13 ) continue;
    if (removeNu && (pdgCode == 12 || pdgCode == 14 || pdgCode == 16)) continue;

    Float_t et = it_gen->et();
    if (et < etMin) continue; // pass a minimum et threshold, default 0

    Float_t dR = reco::deltaR(thisPho->momentum(), it_gen->momentum());
    if (dR > dRMax) continue; // within deltaR cone

    genCalIsoSum += et;

  }

  return genCalIsoSum;
}

float getGenTrkIso(edm::Handle<reco::GenParticleCollection> handle, reco::GenParticleCollection::const_iterator thisPho, const Float_t dRMax) {

  const Float_t ptMin = 0.0;
  Float_t genTrkIsoSum = 0.0;
  if (!handle.isValid()) return genTrkIsoSum;

  for (reco::GenParticleCollection::const_iterator it_gen=handle->begin(); it_gen!=handle->end(); ++it_gen){

    if (it_gen == thisPho) continue;        // can't be the original photon
    if (it_gen->status() != 1) continue;    // need to be a stable particle
    if (thisPho->collisionId() != it_gen->collisionId()) continue; // has to come from the same collision

    if (it_gen->charge() == 0) continue;    // we should not count neutral particles

    Float_t pt = it_gen->pt();
    if (pt < ptMin) continue; // pass a minimum pt threshold, default 0

    Float_t dR = reco::deltaR(thisPho->momentum(), it_gen->momentum());
    if (dR > dRMax) continue; // within deltaR cone
    genTrkIsoSum += pt;

  }// end of loop over gen particles

  return genTrkIsoSum;
}

void ggNtuplizer::makeBranchesGenParticles(TTree* tree_)
{
  tree_->Branch("nMC", &nMC_, "nMC/I");
  tree_->Branch("mcPID", &mcPID);
  tree_->Branch("mcVtx_x", &mcVtx_x);
  tree_->Branch("mcVtx_y", &mcVtx_y);
  tree_->Branch("mcVtx_z", &mcVtx_z);
  tree_->Branch("mcPt", &mcPt);
  tree_->Branch("mcMass", &mcMass);
  tree_->Branch("mcEta", &mcEta);
  tree_->Branch("mcPhi", &mcPhi);
  tree_->Branch("mcE", &mcE);
  tree_->Branch("mcEt", &mcEt);
  tree_->Branch("mcGMomPID", &mcGMomPID);
  tree_->Branch("mcMomPID", &mcMomPID);
  tree_->Branch("mcMomPt", &mcMomPt);
  tree_->Branch("mcMomMass", &mcMomMass);
  tree_->Branch("mcMomEta", &mcMomEta);
  tree_->Branch("mcMomPhi", &mcMomPhi);
  tree_->Branch("mcIndex", &mcIndex);
  tree_->Branch("mcDecayType", &mcDecayType); //-999:non W or Z, 1:hardronic, 2:e, 3:mu, 4:tau
  tree_->Branch("mcParentage", &mcParentage); // 16*lepton + 8*boson + 4*non-prompt + 2*qcd + exotics
  tree_->Branch("mcStatus", &mcStatus); // status of the particle
  tree_->Branch("mcCalIsoDR03", &mcCalIsoDR03);
  tree_->Branch("mcTrkIsoDR03", &mcTrkIsoDR03);
  tree_->Branch("mcCalIsoDR04", &mcCalIsoDR04);
  tree_->Branch("mcTrkIsoDR04", &mcTrkIsoDR04);
}

void ggNtuplizer::fillGenParticles(const edm::Event& e)
{
  // Fills tree branches with generated particle info.

  // cleanup from previous execution
  nMC_ = 0;
  int genIndex = 0;
  mcPID.clear();
  mcVtx_x.clear();
  mcVtx_y.clear();
  mcVtx_z.clear();
  mcPt.clear();
  mcMass.clear();
  mcEta.clear();
  mcPhi.clear();
  mcE.clear();
  mcEt.clear();
  mcGMomPID.clear();
  mcMomPID.clear();
  mcMomPt.clear();
  mcMomMass.clear();
  mcMomEta.clear();
  mcMomPhi.clear();
  mcIndex.clear();
  mcDecayType.clear();
  mcParentage.clear();
  mcStatus.clear();
  mcCalIsoDR03.clear();
  mcTrkIsoDR03.clear();
  mcCalIsoDR04.clear();
  mcTrkIsoDR04.clear();

  edm::Handle<vector<reco::GenParticle> > genParticlesHandle;
  e.getByToken(genParticlesCollection_, genParticlesHandle);

  if (!genParticlesHandle.isValid()) {
    // FIXME: print some warning message here?
    return;
  }

  for (vector<GenParticle>::const_iterator ip = genParticlesHandle->begin(); ip != genParticlesHandle->end(); ++ip) {
    genIndex++;

    int status = ip->status() - 10*(ip->status()/10);
    bool stableFinalStateParticle = status == 1 && ip->pt() > 5.0;

    // keep all the photons with pT > 5.0 and all leptons;
    bool photonOrLepton =
      (status == 1 && ip->pdgId() == 22 && ip->pt() > 5.0 ) ||
      (status == 1 && ( abs(ip->pdgId()) >= 11 && abs(ip->pdgId()) <= 16 ))  ||
      (status < 10 && abs(ip->pdgId()) == 15 );
    // select also Z, W, H, and top
    bool heavyParticle =
      (ip->pdgId() == 23 || abs(ip->pdgId()) == 24 || ip->pdgId() == 25 ||
        abs(ip->pdgId()) == 6 || abs(ip->pdgId()) == 5);

    if ( stableFinalStateParticle || heavyParticle || photonOrLepton ) {
      const Candidate *p = (const Candidate*)&(*ip);
      if (!runOnParticleGun_ && !p->mother()) continue;

      reco::GenParticleRef partRef = reco::GenParticleRef(genParticlesHandle,
                      ip-genParticlesHandle->begin());
      genpartparentage::GenParticleParentage particleHistory(partRef);

      mcPID    .push_back(p->pdgId());
      mcVtx_x  .push_back(p->vx());
      mcVtx_y  .push_back(p->vy());
      mcVtx_z  .push_back(p->vz());
      mcPt     .push_back(p->pt());
      mcMass   .push_back(p->mass());
      mcEta    .push_back(p->eta());
      mcPhi    .push_back(p->phi());
      mcE      .push_back(p->energy());
      mcEt     .push_back(p->et());
      mcParentage.push_back(particleHistory.hasLeptonParent()*16   +
      particleHistory.hasBosonParent()*8     +
      particleHistory.hasNonPromptParent()*4 +
      particleHistory.hasQCDParent()*2       +
      particleHistory.hasExoticParent());
      mcStatus.push_back(p->status());

      int mcDecayType_ = -999;
      // if genParticle is W or Z, check its decay type
      if ( ip->pdgId() == 23 || abs(ip->pdgId()) == 24 ) {
        for (size_t k=0; k < p->numberOfDaughters(); ++k) {
          const Candidate *dp = p->daughter(k);
          if (abs(dp->pdgId())<=6)
            mcDecayType_ = 1;
          else if (abs(dp->pdgId())==11 || abs(dp->pdgId())==12)
            mcDecayType_ = 2;
          else if (abs(dp->pdgId())==13 || abs(dp->pdgId())==14)
            mcDecayType_ = 3;
          else if (abs(dp->pdgId())==15 || abs(dp->pdgId())==16)
            mcDecayType_ = 4;
        }
      }
      mcDecayType.push_back(mcDecayType_);
      int mcGMomPID_ = -999;
      int mcMomPID_  = -999;
      float mcMomPt_    = -999.;
      float mcMomMass_  = -999.;
      float mcMomEta_   = -999.;
      float mcMomPhi_   = -999.;
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

      mcIndex.push_back(genIndex-1);

      mcCalIsoDR03.push_back( getGenCalIso(genParticlesHandle, ip, 0.3, false, false) );
            mcTrkIsoDR03.push_back( getGenTrkIso(genParticlesHandle, ip, 0.3) );
            mcCalIsoDR04.push_back( getGenCalIso(genParticlesHandle, ip, 0.4, false, false) );
      mcTrkIsoDR04.push_back( getGenTrkIso(genParticlesHandle, ip, 0.4) );

      nMC_++;
    } // save info on particles of interest
  } // loop over gen-level particles

}
