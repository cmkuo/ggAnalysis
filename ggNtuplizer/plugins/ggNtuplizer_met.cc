#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "ggAnalysis/ggNtuplizer/interface/ggNtuplizer.h"

using namespace std;

Int_t metFilters_;
float genMET_;
float genMETPhi_;
float pfMET_;
float pfMETPhi_;
float pfMETsumEt_;
float pfMETmEtSig_;
float pfMETSig_;
float pfMET_T1JERUp_;
float pfMET_T1JERDo_;
float pfMET_T1JESUp_;
float pfMET_T1JESDo_;
float pfMET_T1MESUp_;
float pfMET_T1MESDo_;
float pfMET_T1EESUp_;
float pfMET_T1EESDo_;
float pfMET_T1PESUp_;
float pfMET_T1PESDo_;
float pfMET_T1TESUp_;
float pfMET_T1TESDo_;
float pfMET_T1UESUp_;
float pfMET_T1UESDo_;
float pfMET_T1TxyPhi_;
float pfMET_T1TxyPt_;

/*
float     noHFMET_;
float     noHFMETPhi_;
float     noHFMETType1_;
float     noHF_T1JERUp_;
float     noHF_T1JERDo_;
float     noHF_T1JESUp_;
float     noHF_T1JESDo_;
float     noHF_T1MESUp_;
float     noHF_T1MESDo_;
float     noHF_T1EESUp_;
float     noHF_T1EESDo_;
float     noHF_T1TESUp_;
float     noHF_T1TESDo_;
float     noHF_T1UESUp_;
float     noHF_T1UESDo_;
float     noHF_T1Phi_;
float     noHF_T1Px_;
float     noHF_T1Py_;
float     noHF_T1SumEt_;
float     noHF_T1TxyPhi_;
float     noHF_T1TxyPt_;
*/

void ggNtuplizer::branchesMET(TTree* tree) {

  if (addFilterInfoAOD_ || addFilterInfoMINIAOD_) {
    tree->Branch("metFilters", &metFilters_);
  }
  if (doGenParticles_) {
    tree->Branch("genMET",      &genMET_);
    tree->Branch("genMETPhi",   &genMETPhi_);
  }
  tree->Branch("pfMET",         &pfMET_);
  tree->Branch("pfMETPhi",      &pfMETPhi_);
  tree->Branch("pfMETsumEt",    &pfMETsumEt_);
  tree->Branch("pfMETmEtSig",   &pfMETmEtSig_);
  tree->Branch("pfMETSig",      &pfMETSig_);
  tree->Branch("pfMET_T1JERUp", &pfMET_T1JERUp_);
  tree->Branch("pfMET_T1JERDo", &pfMET_T1JERDo_);
  tree->Branch("pfMET_T1JESUp", &pfMET_T1JESUp_);
  tree->Branch("pfMET_T1JESDo", &pfMET_T1JESDo_);
  tree->Branch("pfMET_T1MESUp", &pfMET_T1MESUp_);
  tree->Branch("pfMET_T1MESDo", &pfMET_T1MESDo_);
  tree->Branch("pfMET_T1EESUp", &pfMET_T1EESUp_);
  tree->Branch("pfMET_T1EESDo", &pfMET_T1EESDo_);
  tree->Branch("pfMET_T1PESUp", &pfMET_T1PESUp_);
  tree->Branch("pfMET_T1PESDo", &pfMET_T1PESDo_);
  tree->Branch("pfMET_T1TESUp", &pfMET_T1TESUp_);
  tree->Branch("pfMET_T1TESDo", &pfMET_T1TESDo_);
  tree->Branch("pfMET_T1UESUp", &pfMET_T1UESUp_);
  tree->Branch("pfMET_T1UESDo", &pfMET_T1UESDo_);

  /*
  if (doNoHFMET_){
    tree->Branch("noHFMET",      &noHFMET_);
    tree->Branch("noHFMETPhi",   &noHFMETPhi_);
    tree->Branch("noHFMETType1", &noHFMETType1_);     
    tree->Branch("noHF_T1JERUp", &noHF_T1JERUp_);
    tree->Branch("noHF_T1JERDo", &noHF_T1JERDo_);
    tree->Branch("noHF_T1JESUp", &noHF_T1JESUp_);
    tree->Branch("noHF_T1JESDo", &noHF_T1JESDo_);
    tree->Branch("noHF_T1MESUp", &noHF_T1MESUp_);
    tree->Branch("noHF_T1MESDo", &noHF_T1MESDo_);
    tree->Branch("noHF_T1EESUp", &noHF_T1EESUp_);
    tree->Branch("noHF_T1EESDo", &noHF_T1EESDo_);
    tree->Branch("noHF_T1TESUp", &noHF_T1TESUp_);
    tree->Branch("noHF_T1TESDo", &noHF_T1TESDo_);
    tree->Branch("noHF_T1UESUp", &noHF_T1UESUp_);
    tree->Branch("noHF_T1UESDo", &noHF_T1UESDo_);
    tree->Branch("noHF_T1Phi",   &noHF_T1Phi_);
    tree->Branch("noHF_T1Px",    &noHF_T1Px_);
    tree->Branch("noHF_T1Py",    &noHF_T1Py_);
    tree->Branch("noHF_T1SumEt", &noHF_T1SumEt_);
    tree->Branch("noHF_T1TxyPhi", &noHF_T1TxyPhi_);
    tree->Branch("noHF_T1TxyPt", &noHF_T1TxyPt_);
   }
  */
}

void ggNtuplizer::fillMET(const edm::Event& e, const edm::EventSetup& es) {

  metFilters_ = 0;
  if (addFilterInfoAOD_ && e.isRealData()){
    
    edm::Handle<bool> hcalNoiseHandle;
    e.getByLabel("HBHENoiseFilterResultProducer", "HBHENoiseFilterResult", hcalNoiseHandle);
    bool HBHENoiseResult_ = *hcalNoiseHandle;

    edm::Handle<bool> hcalIsoNoiseHandle;
    e.getByLabel("HBHENoiseFilterResultProducer", "HBHEIsoNoiseFilterResult", hcalIsoNoiseHandle);
    bool HBHEIsoNoiseResult_ = *hcalIsoNoiseHandle;
    
    edm::Handle<bool> cSCHandle;
    e.getByLabel("CSCTightHaloFilter", "", cSCHandle);
    bool CSCHaloResult_ = *cSCHandle;
    
    edm::Handle<bool> eCALTPHandle;
    e.getByLabel("EcalDeadCellTriggerPrimitiveFilter", "", eCALTPHandle);
    bool EcalDeadCellTFResult_ = *eCALTPHandle;

    edm::Handle<bool> bADSCHandle;
    e.getByLabel("eeBadScFilter", "", bADSCHandle);
    bool EEBadSCResult_ = *bADSCHandle;

    //#edm::Handle<bool> gOODVertexHandle;
    //e.getByLabel("primaryVertexFilter", "GoodVertexFilter", gOODVertexHandle);
    //bool goodVertexResult_ = *gOODVertexHandle;
     
    if ( !HBHENoiseResult_      ) metFilters_ += 1; 
    if ( !HBHEIsoNoiseResult_   ) metFilters_ += 2; 
    if ( !CSCHaloResult_        ) metFilters_ += 4; 
    //if ( !goodVertexResult_     ) metFilters_ += 8; 
    if ( !EEBadSCResult_        ) metFilters_ += 16; 
    if ( !EcalDeadCellTFResult_ ) metFilters_ += 32; 
  }

  if (addFilterInfoMINIAOD_ && e.isRealData()){
    string filterNamesToCheck[7] = {
      "Flag_HBHENoiseFilter",
      "Flag_HBHENoiseIsoFilter", 
      "Flag_CSCTightHaloFilter",
      "Flag_goodVertices",
      "Flag_eeBadScFilter",
      "Flag_EcalDeadCellTriggerPrimitiveFilter",
      "Flag_EcalDeadCellBoundaryEnergyFilter"
    };

    edm::Handle<edm::TriggerResults> patFilterResultsHandle;
    e.getByToken(patTrgResultsLabel_, patFilterResultsHandle);
    edm::TriggerResults const& patFilterResults = *patFilterResultsHandle;
    
    auto&& filterNames = e.triggerNames(patFilterResults);
    for (unsigned iF = 0; iF != 7; ++iF) {
      unsigned index = filterNames.triggerIndex(filterNamesToCheck[iF]);
      if ( index == filterNames.size() ) 
	edm::LogError("Unknown MET filter label") 
	  << filterNamesToCheck[iF] << " is missing, exiting";
      else {
	if ( !patFilterResults.accept(index) ) {
	  metFilters_ += pow(2, iF);
	}
      }
    }
  }

  edm::Handle<edm::View<pat::MET> > pfMETHandle;
  e.getByToken(pfMETlabel_, pfMETHandle);

  genMET_      = -99;
  genMETPhi_   = -99;
  pfMET_       = -99;
  pfMETPhi_    = -99;
  pfMETsumEt_  = -99;
  pfMETmEtSig_ = -99;
  pfMETSig_    = -99;

  if (pfMETHandle.isValid()) {
    const pat::MET *pfMET = 0;
    pfMET = &(pfMETHandle->front());

    pfMET_       = pfMET->et();
    pfMETPhi_    = pfMET->phi();
    pfMETsumEt_  = pfMET->sumEt();
    pfMETmEtSig_ = (pfMET->mEtSig() < 1.e10) ? pfMET->mEtSig() : 0;
    pfMETSig_    = (pfMET->significance() < 1.e10) ? pfMET->significance() : 0;;

    if (!isAOD_) {
      // Type1MET uncertainties =======================================
      pfMET_T1JERUp_ = pfMET->shiftedPt(pat::MET::JetResUp);
      pfMET_T1JERDo_ = pfMET->shiftedPt(pat::MET::JetResDown);
      pfMET_T1JESUp_ = pfMET->shiftedPt(pat::MET::JetEnUp);
      pfMET_T1JESDo_ = pfMET->shiftedPt(pat::MET::JetEnDown);
      pfMET_T1MESUp_ = pfMET->shiftedPt(pat::MET::MuonEnUp);
      pfMET_T1MESDo_ = pfMET->shiftedPt(pat::MET::MuonEnDown);
      pfMET_T1EESUp_ = pfMET->shiftedPt(pat::MET::ElectronEnUp);
      pfMET_T1EESDo_ = pfMET->shiftedPt(pat::MET::ElectronEnDown);
      pfMET_T1PESUp_ = pfMET->shiftedPt(pat::MET::PhotonEnUp);
      pfMET_T1PESDo_ = pfMET->shiftedPt(pat::MET::PhotonEnDown);
      pfMET_T1TESUp_ = pfMET->shiftedPt(pat::MET::TauEnUp);
      pfMET_T1TESDo_ = pfMET->shiftedPt(pat::MET::TauEnDown);
      pfMET_T1UESUp_ = pfMET->shiftedPt(pat::MET::UnclusteredEnUp);
      pfMET_T1UESDo_ = pfMET->shiftedPt(pat::MET::UnclusteredEnDown);
    }

    if (!e.isRealData()) {
      genMET_    = pfMET->genMET()->et();
      genMETPhi_ = pfMET->genMET()->phi();
    }

  } 

  /*
    else
    edm::LogWarning("ggNtuplizer") << "MET info unavailable";

   if (doNoHFMET_){
    pat::METCollection newMet;
    pat::METCollection t1txyMet;

    edm::Handle < pat::METCollection > newMetHandle;
    edm::InputTag _newMetLabel("slimmedMETsNoHF","","ggKit");
    edm::InputTag _t1txyMetLabel("patPFMetT1TxyNoHF");
    e.getByLabel(_newMetLabel,newMetHandle);
    if (newMetHandle.isValid() ) newMet = *newMetHandle;

    edm::Handle < pat::METCollection > t1txyMetHandle;
    e.getByLabel(_t1txyMetLabel,t1txyMetHandle);
    if (t1txyMetHandle.isValid() ) t1txyMet = *t1txyMetHandle;
    if (t1txyMetHandle.isValid() && newMetHandle.isValid()){
      noHFMET_   = newMet[0].pt();
      noHFMETPhi_= newMet[0].phi();

      //alternate way to get the Type1 value
      noHFMETType1_  = newMet[0].shiftedPt(pat::MET::NoShift, pat::MET::Type1); //second argument is Type1 per default

      //Type1MET uncertainties =======================================
      noHF_T1JERUp_ = newMet[0].shiftedPt(pat::MET::JetResUp);
      noHF_T1JERDo_ = newMet[0].shiftedPt(pat::MET::JetResDown);
      noHF_T1JESUp_ = newMet[0].shiftedPt(pat::MET::JetEnUp);
      noHF_T1JESDo_ = newMet[0].shiftedPt(pat::MET::JetEnDown);
      noHF_T1MESUp_ = newMet[0].shiftedPt(pat::MET::MuonEnUp);
      noHF_T1MESDo_ = newMet[0].shiftedPt(pat::MET::MuonEnDown);
      noHF_T1EESUp_ = newMet[0].shiftedPt(pat::MET::ElectronEnUp);
      noHF_T1EESDo_ = newMet[0].shiftedPt(pat::MET::ElectronEnDown);
      noHF_T1TESUp_ = newMet[0].shiftedPt(pat::MET::TauEnUp);
      noHF_T1TESDo_ = newMet[0].shiftedPt(pat::MET::TauEnDown);
      noHF_T1UESUp_ = newMet[0].shiftedPt(pat::MET::UnclusteredEnUp);
      noHF_T1UESDo_ = newMet[0].shiftedPt(pat::MET::UnclusteredEnDown);

      //other functions to access the shifted MET variables =================
      noHF_T1Phi_ = newMet[0].shiftedPhi(pat::MET::NoShift);  //second argument is Type1 per default
      noHF_T1Px_  = newMet[0].shiftedPx(pat::MET::NoShift);  //second argument is Type1 per default
      noHF_T1Py_  = newMet[0].shiftedPy(pat::MET::NoShift);  //second argument is Type1 per default
      noHF_T1SumEt_ = newMet[0].shiftedSumEt(pat::MET::NoShift);  //second argument is Type1 per default


      //extra Txy stuff =======================================================
      noHF_T1TxyPhi_ = t1txyMet[0].phi();
      noHF_T1TxyPt_ = t1txyMet[0].pt();
    }else
    edm::LogWarning("ggNtuplizer") << "MET info unavailable";


  }
  */

}
