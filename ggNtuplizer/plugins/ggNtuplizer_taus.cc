// Updated by Abdollah Mohammadi (KSU)  [18 May 2015]
// abdollah.mohammadi@cern.ch

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "ggAnalysis/ggNtuplizer/interface/ggNtuplizer.h"

using namespace std;

// (local) variables associated with tree branches
Int_t          nTau_;

// decay mode discriminators

//Tau Id & Isolation
vector<bool>   pfTausDiscriminationByDecayModeFinding_;
vector<bool>   pfTausDiscriminationByDecayModeFindingNewDMs_;
//vector<bool>   tauByLooseElectronRejection_;
//vector<bool>   tauByMediumElectronRejection_;
//vector<bool>   tauByTightElectronRejection_;
vector<bool>   tauByMVA5LooseElectronRejection_;
vector<bool>   tauByMVA5MediumElectronRejection_;
vector<bool>   tauByMVA5TightElectronRejection_;
vector<bool>   tauByMVA5VTightElectronRejection_;
vector<bool>   tauByLooseMuonRejection3_;
vector<bool>   tauByTightMuonRejection3_;
//vector<bool>   tauByMVALooseMuonRejection_;
//vector<bool>   tauByMVAMediumMuonRejection_;
//vector<bool>   tauByMVATightMuonRejection_;
vector<bool>   tauByLooseCombinedIsolationDeltaBetaCorr3Hits_;
vector<bool>   tauByMediumCombinedIsolationDeltaBetaCorr3Hits_;
vector<bool>   tauByTightCombinedIsolationDeltaBetaCorr3Hits_;
vector<float>   tauCombinedIsolationDeltaBetaCorrRaw3Hits_;
//vector<bool>   tauByVLooseIsolationMVA3newDMwoLT_;
//vector<bool>   tauByLooseIsolationMVA3newDMwoLT_;
//vector<bool>   tauByMediumIsolationMVA3newDMwoLT_;
//vector<bool>   tauByTightIsolationMVA3newDMwoLT_;
//vector<bool>   tauByVTightIsolationMVA3newDMwoLT_;
//vector<bool>   tauByVVTightIsolationMVA3newDMwoLT_;
//vector<float>   tauByIsolationMVA3newDMwoLTraw_;
vector<bool>   tauByVLooseIsolationMVA3oldDMwLT_;
vector<bool>   tauByLooseIsolationMVA3oldDMwLT_;
vector<bool>   tauByMediumIsolationMVA3oldDMwLT_;
vector<bool>   tauByTightIsolationMVA3oldDMwLT_;
vector<bool>   tauByVTightIsolationMVA3oldDMwLT_;
vector<bool>   tauByVVTightIsolationMVA3oldDMwLT_;
vector<float>   tauByIsolationMVA3oldDMwLTraw_;
//vector<bool>   tauByVLooseIsolationMVA3oldDMwoLT_;
//vector<bool>   tauByLooseIsolationMVA3oldDMwoLT_;
//vector<bool>   tauByTightIsolationMVA3oldDMwoLT_;
//vector<bool>   tauByVTightIsolationMVA3oldDMwoLT_;
//vector<bool>   tauByVVTightIsolationMVA3oldDMwoLT_;
//vector<float>   tauByIsolationMVA3oldDMwoLTraw_;
vector<bool>   tauByVLooseIsolationMVA3newDMwLT_;
vector<bool>   tauByLooseIsolationMVA3newDMwLT_;
vector<bool>   tauByMediumIsolationMVA3newDMwLT_;
vector<bool>   tauByTightIsolationMVA3newDMwLT_;
vector<bool>   tauByVTightIsolationMVA3newDMwLT_;
vector<bool>   tauByVVTightIsolationMVA3newDMwLT_;
vector<float>   tauByIsolationMVA3newDMwLTraw_;

//Tau Kinematics
vector<float> tauEta_;
vector<float> tauPhi_;
vector<float> tauPt_;
vector<float> tauEt_;
vector<float> tauCharge_;
vector<int>   tauDecayMode_;
vector<float> tauP_;
vector<float> tauPx_;
vector<float> tauPy_;
vector<float> tauPz_;
vector<float> tauVz_;
vector<float> tauEnergy_;
vector<float> tauMass_;
vector<float> tauDxy_;
vector<float> tauZImpact_;

//Tau Ingredients
vector<float> tauChargedIsoPtSum_;
vector<float> tauNeutralIsoPtSum_;
vector<float> tauPuCorrPtSum_;
vector<float> tauNumSignalPFChargedHadrCands_;
vector<float> tauNumSignalPFNeutrHadrCands_;
vector<float> tauNumSignalPFGammaCands_;
vector<float> tauNumSignalPFCands_;
vector<float> tauNumIsolationPFChargedHadrCands_;
vector<float> tauNumIsolationPFNeutrHadrCands_;
vector<float> tauNumIsolationPFGammaCands_;
vector<float> tauNumIsolationPFCands_;
//vector<float> tauEMFraction_;
//vector<float> tauHCAL3x3OverPLead_;
//vector<float> tauHCALMaxOverPLead_;
//vector<float> tauHCALTotOverPLead_;
//vector<float> tauIsolationPFChargedHadrCandsPtSum_;
//vector<float> tauIsolationPFGammaCandsEtSum_;
//vector<float> tauLeadPFChargedHadrCandsignedSipt_;
vector<bool>  tauLeadChargedHadronExists_;
vector<float> tauLeadChargedHadronEta_;
vector<float> tauLeadChargedHadronPhi_;
vector<float> tauLeadChargedHadronPt_;

void ggNtuplizer::branchesTaus(TTree* tree)
{
    
    
    
    tree->Branch("nTau", &nTau_);
    
    //Tau Id & Isolation
    tree->Branch("pfTausDiscriminationByDecayModeFinding", &pfTausDiscriminationByDecayModeFinding_);
    tree->Branch("pfTausDiscriminationByDecayModeFindingNewDMs", &pfTausDiscriminationByDecayModeFindingNewDMs_);
//    tree->Branch("tauByLooseElectronRejection", &tauByLooseElectronRejection_);
//    tree->Branch("tauByMediumElectronRejection", &tauByMediumElectronRejection_);
//    tree->Branch("tauByTightElectronRejection", &tauByTightElectronRejection_);
    tree->Branch("tauByMVA5LooseElectronRejection", &tauByMVA5LooseElectronRejection_);
    tree->Branch("tauByMVA5MediumElectronRejection", &tauByMVA5MediumElectronRejection_);
    tree->Branch("tauByMVA5TightElectronRejection", &tauByMVA5TightElectronRejection_);
    tree->Branch("tauByMVA5VTightElectronRejection", &tauByMVA5VTightElectronRejection_);
    tree->Branch("tauByLooseMuonRejection3", &tauByLooseMuonRejection3_);
    tree->Branch("tauByTightMuonRejection3", &tauByTightMuonRejection3_);
//    tree->Branch("tauByMVALooseMuonRejection", &tauByMVALooseMuonRejection_);
//    tree->Branch("tauByMVAMediumMuonRejection", &tauByMVAMediumMuonRejection_);
//    tree->Branch("tauByMVATightMuonRejection", &tauByMVATightMuonRejection_);
    tree->Branch("tauByLooseCombinedIsolationDeltaBetaCorr3Hits", &tauByLooseCombinedIsolationDeltaBetaCorr3Hits_);
    tree->Branch("tauByMediumCombinedIsolationDeltaBetaCorr3Hits", &tauByMediumCombinedIsolationDeltaBetaCorr3Hits_);
    tree->Branch("tauByTightCombinedIsolationDeltaBetaCorr3Hits", &tauByTightCombinedIsolationDeltaBetaCorr3Hits_);
    tree->Branch("tauCombinedIsolationDeltaBetaCorrRaw3Hits", &tauCombinedIsolationDeltaBetaCorrRaw3Hits_);
//    tree->Branch("tauByVLooseIsolationMVA3newDMwoLT", &tauByVLooseIsolationMVA3newDMwoLT_);
//    tree->Branch("tauByLooseIsolationMVA3newDMwoLT", &tauByLooseIsolationMVA3newDMwoLT_);
//    tree->Branch("tauByMediumIsolationMVA3newDMwoLT", &tauByMediumIsolationMVA3newDMwoLT_);
//    tree->Branch("tauByTightIsolationMVA3newDMwoLT", &tauByTightIsolationMVA3newDMwoLT_);
//    tree->Branch("tauByVTightIsolationMVA3newDMwoLT", &tauByVTightIsolationMVA3newDMwoLT_);
//    tree->Branch("tauByVVTightIsolationMVA3newDMwoLT", &tauByVVTightIsolationMVA3newDMwoLT_);
//    tree->Branch("tauByIsolationMVA3newDMwoLTraw", &tauByIsolationMVA3newDMwoLTraw_);
    tree->Branch("tauByVLooseIsolationMVA3oldDMwLT", &tauByVLooseIsolationMVA3oldDMwLT_);
    tree->Branch("tauByLooseIsolationMVA3oldDMwLT", &tauByLooseIsolationMVA3oldDMwLT_);
    tree->Branch("tauByMediumIsolationMVA3oldDMwLT", &tauByMediumIsolationMVA3oldDMwLT_);
    tree->Branch("tauByTightIsolationMVA3oldDMwLT", &tauByTightIsolationMVA3oldDMwLT_);
    tree->Branch("tauByVTightIsolationMVA3oldDMwLT", &tauByVTightIsolationMVA3oldDMwLT_);
    tree->Branch("tauByVVTightIsolationMVA3oldDMwLT", &tauByVVTightIsolationMVA3oldDMwLT_);
    tree->Branch("tauByIsolationMVA3oldDMwLTraw", &tauByIsolationMVA3oldDMwLTraw_);
//    tree->Branch("tauByVLooseIsolationMVA3oldDMwoLT", &tauByVLooseIsolationMVA3oldDMwoLT_);
//    tree->Branch("tauByLooseIsolationMVA3oldDMwoLT", &tauByLooseIsolationMVA3oldDMwoLT_);
//    tree->Branch("tauByTightIsolationMVA3oldDMwoLT", &tauByTightIsolationMVA3oldDMwoLT_);
//    tree->Branch("tauByVTightIsolationMVA3oldDMwoLT", &tauByVTightIsolationMVA3oldDMwoLT_);
//    tree->Branch("tauByVVTightIsolationMVA3oldDMwoLT", &tauByVVTightIsolationMVA3oldDMwoLT_);
//    tree->Branch("tauByIsolationMVA3oldDMwoLTraw", &tauByIsolationMVA3oldDMwoLTraw_);
    tree->Branch("tauByLooseIsolationMVA3newDMwLT", &tauByLooseIsolationMVA3newDMwLT_);
    tree->Branch("tauByVLooseIsolationMVA3newDMwLT", &tauByVLooseIsolationMVA3newDMwLT_);
    tree->Branch("tauByMediumIsolationMVA3newDMwLT", &tauByMediumIsolationMVA3newDMwLT_);
    tree->Branch("tauByTightIsolationMVA3newDMwLT", &tauByTightIsolationMVA3newDMwLT_);
    tree->Branch("tauByVTightIsolationMVA3newDMwLT", &tauByVTightIsolationMVA3newDMwLT_);
    tree->Branch("tauByVVTightIsolationMVA3newDMwLT", &tauByVVTightIsolationMVA3newDMwLT_);
    tree->Branch("tauByIsolationMVA3newDMwLTraw", &tauByIsolationMVA3newDMwLTraw_);
    
    //Tau Kinematics
    tree->Branch("tauEta"  ,&tauEta_);
    tree->Branch("tauPhi"  ,&tauPhi_);
    tree->Branch("tauPt"  ,&tauPt_);
    tree->Branch("tauEt"  ,&tauEt_);
    tree->Branch("tauCharge"  ,&tauCharge_);
    tree->Branch("tauP"  ,&tauP_);
    tree->Branch("tauPx"  ,&tauPx_);
    tree->Branch("tauPy"  ,&tauPy_);
    tree->Branch("tauPz"  ,&tauPz_);
    tree->Branch("tauVz"  ,&tauVz_);
    tree->Branch("tauEnergy"  ,&tauEnergy_);
    tree->Branch("tauMass"  ,&tauMass_);
    tree->Branch("tauDxy"  ,&tauDxy_);
    tree->Branch("tauZImpact"  ,&tauZImpact_);
    
    // Tau Ingredients
    tree->Branch("tauDecayMode"  ,&tauDecayMode_);
    //    tree->Branch("tauEMFraction"  ,&tauEMFraction_);
    //    tree->Branch("tauHCAL3x3OverPLead"  ,&tauHCAL3x3OverPLead_);
    //    tree->Branch("tauHCALMaxOverPLead"  ,&tauHCALMaxOverPLead_);
    //    tree->Branch("tauHCALTotOverPLead"  ,&tauHCALTotOverPLead_);
    //    tree->Branch("tauIsolationPFChargedHadrCandsPtSum"  ,&tauIsolationPFChargedHadrCandsPtSum_);
    //    tree->Branch("tauIsolationPFGammaCandsEtSum"  ,&tauIsolationPFGammaCandsEtSum_);
    //    tree->Branch("tauLeadPFChargedHadrCandsignedSipt"  ,&tauLeadPFChargedHadrCandsignedSipt_);
    tree->Branch("tauLeadChargedHadronExists"  ,&tauLeadChargedHadronExists_);
    tree->Branch("tauLeadChargedHadronEta"  ,&tauLeadChargedHadronEta_);
    tree->Branch("tauLeadChargedHadronPhi"  ,&tauLeadChargedHadronPhi_);
    tree->Branch("tauLeadChargedHadronPt"  ,&tauLeadChargedHadronPt_);
    tree->Branch("tauChargedIsoPtSum"  ,&tauChargedIsoPtSum_);
    tree->Branch("tauNeutralIsoPtSum"  ,&tauNeutralIsoPtSum_);
    tree->Branch("tauPuCorrPtSum"  ,&tauPuCorrPtSum_);
    tree->Branch("tauNumSignalPFChargedHadrCands"  ,&tauNumSignalPFChargedHadrCands_);
    tree->Branch("tauNumSignalPFNeutrHadrCands"  ,&tauNumSignalPFNeutrHadrCands_);
    tree->Branch("tauNumSignalPFGammaCands"  ,&tauNumSignalPFGammaCands_);
    tree->Branch("tauNumSignalPFCands"  ,&tauNumSignalPFCands_);
    tree->Branch("tauNumIsolationPFChargedHadrCands"  ,&tauNumIsolationPFChargedHadrCands_);
    tree->Branch("tauNumIsolationPFNeutrHadrCands"  ,&tauNumIsolationPFNeutrHadrCands_);
    tree->Branch("tauNumIsolationPFGammaCands"  ,&tauNumIsolationPFGammaCands_);
    tree->Branch("tauNumIsolationPFCands"  ,&tauNumIsolationPFCands_);
    
    
    
    
    
    
    
    
    
}

void ggNtuplizer::fillTaus(const edm::Event& e)
{
    
    // Tau Id & Isolation
//    tauByLooseElectronRejection_.clear();
//    tauByMediumElectronRejection_.clear();
//    tauByTightElectronRejection_.clear();
    tauByMVA5LooseElectronRejection_.clear();
    tauByMVA5MediumElectronRejection_.clear();
    tauByMVA5TightElectronRejection_.clear();
    tauByMVA5VTightElectronRejection_.clear();
    tauByLooseMuonRejection3_.clear();
    tauByTightMuonRejection3_.clear();
//    tauByMVALooseMuonRejection_.clear();
//    tauByMVAMediumMuonRejection_.clear();
//    tauByMVATightMuonRejection_.clear();
    pfTausDiscriminationByDecayModeFinding_.clear();
    pfTausDiscriminationByDecayModeFindingNewDMs_.clear();
    tauByLooseCombinedIsolationDeltaBetaCorr3Hits_.clear();
    tauByMediumCombinedIsolationDeltaBetaCorr3Hits_.clear();
    tauByTightCombinedIsolationDeltaBetaCorr3Hits_.clear();
    tauCombinedIsolationDeltaBetaCorrRaw3Hits_.clear();
//    tauByVLooseIsolationMVA3newDMwoLT_.clear();
//    tauByLooseIsolationMVA3newDMwoLT_.clear();
//    tauByMediumIsolationMVA3newDMwoLT_.clear();
//    tauByTightIsolationMVA3newDMwoLT_.clear();
//    tauByVTightIsolationMVA3newDMwoLT_.clear();
//    tauByVVTightIsolationMVA3newDMwoLT_.clear();
//    tauByIsolationMVA3newDMwoLTraw_.clear();
    tauByVLooseIsolationMVA3oldDMwLT_.clear();
    tauByLooseIsolationMVA3oldDMwLT_.clear();
    tauByMediumIsolationMVA3oldDMwLT_.clear();
    tauByTightIsolationMVA3oldDMwLT_.clear();
    tauByVTightIsolationMVA3oldDMwLT_.clear();
    tauByVVTightIsolationMVA3oldDMwLT_.clear();
    tauByIsolationMVA3oldDMwLTraw_.clear();
//    tauByVLooseIsolationMVA3oldDMwoLT_.clear();
//    tauByLooseIsolationMVA3oldDMwoLT_.clear();
//    tauByTightIsolationMVA3oldDMwoLT_.clear();
//    tauByVTightIsolationMVA3oldDMwoLT_.clear();
//    tauByVVTightIsolationMVA3oldDMwoLT_.clear();
//    tauByIsolationMVA3oldDMwoLTraw_.clear();
    tauByLooseIsolationMVA3newDMwLT_.clear();
    tauByVLooseIsolationMVA3newDMwLT_.clear();
    tauByMediumIsolationMVA3newDMwLT_.clear();
    tauByTightIsolationMVA3newDMwLT_.clear();
    tauByVTightIsolationMVA3newDMwLT_.clear();
    tauByVVTightIsolationMVA3newDMwLT_.clear();
    tauByIsolationMVA3newDMwLTraw_.clear();
    
    //Tau Kinematics
    tauEta_.clear();
    tauPhi_.clear();
    tauPt_.clear();
    tauEt_.clear();
    tauCharge_.clear();
    tauP_.clear();
    tauPx_.clear();
    tauPy_.clear();
    tauPz_.clear();
    tauVz_.clear();
    tauEnergy_.clear();
    tauMass_.clear();
    tauDxy_.clear();
    tauZImpact_.clear();
    
    // Tau Ingredients
    tauDecayMode_.clear();
    //    tauEMFraction_.clear();
    //    tauHCAL3x3OverPLead_.clear();
    //    tauHCALMaxOverPLead_.clear();
    //    tauHCALTotOverPLead_.clear();
    //    tauIsolationPFChargedHadrCandsPtSum_.clear();
    //    tauIsolationPFGammaCandsEtSum_.clear();
    //    tauLeadPFChargedHadrCandsignedSipt_.clear();
    tauLeadChargedHadronExists_.clear();
    tauLeadChargedHadronEta_.clear();
    tauLeadChargedHadronPhi_.clear();
    tauLeadChargedHadronPt_.clear();
    tauChargedIsoPtSum_.clear();
    tauNeutralIsoPtSum_.clear();
    tauPuCorrPtSum_.clear();
    tauNumSignalPFChargedHadrCands_.clear();
    tauNumSignalPFNeutrHadrCands_.clear();
    tauNumSignalPFGammaCands_.clear();
    tauNumSignalPFCands_.clear();
    tauNumIsolationPFChargedHadrCands_.clear();
    tauNumIsolationPFNeutrHadrCands_.clear();
    tauNumIsolationPFGammaCands_.clear();
    tauNumIsolationPFCands_.clear();
    
    
    
    nTau_ = 0;
    
    edm::Handle<vector<pat::Tau> > tauHandle;
    e.getByToken(tauCollection_, tauHandle);
    
    if (!tauHandle.isValid()) {
        edm::LogWarning("ggNtuplizer") << "no pat::Tau in event";
        return;
    }
    
    //startTaus Lvdp
    for(vector<pat::Tau>::const_iterator itau = tauHandle->begin(); itau != tauHandle->end(); ++itau) {
        
        
        // Tau Id & Isolation
//        tauByLooseElectronRejection_.push_back(itau->tauID("againstElectronLoose"));
//        tauByMediumElectronRejection_.push_back(itau->tauID("againstElectronMedium"));
//        tauByTightElectronRejection_.push_back(itau->tauID("againstElectronTight"));
        tauByMVA5LooseElectronRejection_.push_back(itau->tauID("againstElectronLooseMVA5"));
        tauByMVA5MediumElectronRejection_.push_back(itau->tauID("againstElectronMediumMVA5"));
        tauByMVA5TightElectronRejection_.push_back(itau->tauID("againstElectronTightMVA5"));
        tauByMVA5VTightElectronRejection_.push_back(itau->tauID("againstElectronVTightMVA5"));
        tauByLooseMuonRejection3_.push_back(itau->tauID("againstMuonLoose3"));
        tauByTightMuonRejection3_.push_back(itau->tauID("againstMuonTight3"));
//        tauByMVALooseMuonRejection_.push_back(itau->tauID("againstMuonLooseMVA"));
//        tauByMVAMediumMuonRejection_.push_back(itau->tauID("againstMuonMediumMVA"));
//        tauByMVATightMuonRejection_.push_back(itau->tauID("againstMuonTightMVA"));
        pfTausDiscriminationByDecayModeFinding_.push_back(itau->tauID("decayModeFinding"));
        pfTausDiscriminationByDecayModeFindingNewDMs_.push_back(itau->tauID("decayModeFindingNewDMs"));
        tauByLooseCombinedIsolationDeltaBetaCorr3Hits_.push_back(itau->tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits"));
        tauByMediumCombinedIsolationDeltaBetaCorr3Hits_.push_back(itau->tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits"));
        tauByTightCombinedIsolationDeltaBetaCorr3Hits_.push_back(itau->tauID("byTightCombinedIsolationDeltaBetaCorr3Hits"));
        tauCombinedIsolationDeltaBetaCorrRaw3Hits_.push_back(itau->tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits"));
//        tauByVLooseIsolationMVA3newDMwoLT_.push_back(itau->tauID("byVLooseIsolationMVA3newDMwoLT"));
//        tauByLooseIsolationMVA3newDMwoLT_.push_back(itau->tauID("byLooseIsolationMVA3newDMwoLT"));
//        tauByMediumIsolationMVA3newDMwoLT_.push_back(itau->tauID("byMediumIsolationMVA3newDMwoLT"));
//        tauByTightIsolationMVA3newDMwoLT_.push_back(itau->tauID("byTightIsolationMVA3newDMwoLT"));
//        tauByVTightIsolationMVA3newDMwoLT_.push_back(itau->tauID("byVTightIsolationMVA3newDMwoLT"));
//        tauByVVTightIsolationMVA3newDMwoLT_.push_back(itau->tauID("byVVTightIsolationMVA3newDMwoLT"));
//        tauByIsolationMVA3newDMwoLTraw_.push_back(itau->tauID("byIsolationMVA3newDMwoLTraw"));
        tauByVLooseIsolationMVA3oldDMwLT_.push_back(itau->tauID("byVLooseIsolationMVA3oldDMwLT"));
        tauByLooseIsolationMVA3oldDMwLT_.push_back(itau->tauID("byLooseIsolationMVA3oldDMwLT"));
        tauByMediumIsolationMVA3oldDMwLT_.push_back(itau->tauID("byMediumIsolationMVA3oldDMwLT"));
        tauByTightIsolationMVA3oldDMwLT_.push_back(itau->tauID("byTightIsolationMVA3oldDMwLT"));
        tauByVTightIsolationMVA3oldDMwLT_.push_back(itau->tauID("byVTightIsolationMVA3oldDMwLT"));
        tauByVVTightIsolationMVA3oldDMwLT_.push_back(itau->tauID("byVVTightIsolationMVA3oldDMwLT"));
        tauByIsolationMVA3oldDMwLTraw_.push_back(itau->tauID("byIsolationMVA3oldDMwLTraw"));
//        tauByVLooseIsolationMVA3oldDMwoLT_.push_back(itau->tauID("byVLooseIsolationMVA3oldDMwoLT"));
//        tauByLooseIsolationMVA3oldDMwoLT_.push_back(itau->tauID("byLooseIsolationMVA3oldDMwoLT"));
//        tauByTightIsolationMVA3oldDMwoLT_.push_back(itau->tauID("byTightIsolationMVA3oldDMwoLT"));
//        tauByVTightIsolationMVA3oldDMwoLT_.push_back(itau->tauID("byVTightIsolationMVA3oldDMwoLT"));
//        tauByVVTightIsolationMVA3oldDMwoLT_.push_back(itau->tauID("byVVTightIsolationMVA3oldDMwoLT"));
//        tauByIsolationMVA3oldDMwoLTraw_.push_back(itau->tauID("byIsolationMVA3oldDMwoLTraw"));
        tauByLooseIsolationMVA3newDMwLT_.push_back(itau->tauID("byLooseIsolationMVA3newDMwLT"));
        tauByVLooseIsolationMVA3newDMwLT_.push_back(itau->tauID("byVLooseIsolationMVA3newDMwLT"));
        tauByMediumIsolationMVA3newDMwLT_.push_back(itau->tauID("byMediumIsolationMVA3newDMwLT"));
        tauByTightIsolationMVA3newDMwLT_.push_back(itau->tauID("byTightIsolationMVA3newDMwLT"));
        tauByVTightIsolationMVA3newDMwLT_.push_back(itau->tauID("byVTightIsolationMVA3newDMwLT"));
        tauByVVTightIsolationMVA3newDMwLT_.push_back(itau->tauID("byVVTightIsolationMVA3newDMwLT"));
        tauByIsolationMVA3newDMwLTraw_.push_back(itau->tauID("byIsolationMVA3newDMwLTraw"));
        
        //Tau Kinematics
        tauEta_.push_back(itau->eta());
        tauPhi_.push_back(itau->phi());
        tauPt_.push_back(itau->pt());
        tauEt_.push_back(itau->et());
        tauCharge_.push_back(itau->charge());
        tauP_.push_back(itau->p() );
        tauPx_.push_back(itau->px() );
        tauPy_.push_back(itau->py() );
        tauPz_.push_back(itau->pz() );
        tauVz_.push_back(itau->vz() );
        tauEnergy_.push_back(itau->energy() );
        tauMass_.push_back(itau->mass());
        tauDxy_.push_back(itau->dxy() );
        tauZImpact_.push_back(itau->vertex().z() + 130./tan(itau->theta()));
        
        
        // Tau Ingredients
        tauDecayMode_.push_back(itau->decayMode());
        //        tauEMFraction_.push_back(itau->emFraction());
        //        tauHCAL3x3OverPLead_.push_back(itau->hcal3x3OverPLead());
        //        tauHCALMaxOverPLead_.push_back(itau->hcalMaxOverPLead());
        //        tauHCALTotOverPLead_.push_back(itau->hcalTotOverPLead());
        //        tauIsolationPFChargedHadrCandsPtSum_.push_back(itau->isolationPFChargedHadrCandsPtSum());
        //        tauIsolationPFGammaCandsEtSum_.push_back(itau->isolationPFGammaCandsEtSum());
        //        tauLeadPFChargedHadrCandsignedSipt_.push_back(itau->leadPFChargedHadrCandsignedSipt());
        tauChargedIsoPtSum_.push_back(itau->tauID("chargedIsoPtSum") );
        tauNeutralIsoPtSum_.push_back(itau->tauID("neutralIsoPtSum")  );
        tauPuCorrPtSum_.push_back(itau->tauID("puCorrPtSum")  );
        tauNumSignalPFChargedHadrCands_.push_back(itau->signalPFChargedHadrCands().size());
        tauNumSignalPFNeutrHadrCands_.push_back(itau->signalPFNeutrHadrCands().size());
        tauNumSignalPFGammaCands_.push_back(itau->signalPFGammaCands().size());
        tauNumSignalPFCands_.push_back(itau->signalPFCands().size());
        tauNumIsolationPFChargedHadrCands_.push_back(itau->isolationPFChargedHadrCands().size());
        tauNumIsolationPFNeutrHadrCands_.push_back(itau->isolationPFNeutrHadrCands().size());
        tauNumIsolationPFGammaCands_.push_back(itau->isolationPFGammaCands().size());
        tauNumIsolationPFCands_.push_back(itau->isolationPFCands().size());
        
        const reco::PFCandidatePtr& leadPFChargedHadrCand_Ref = itau->leadPFChargedHadrCand();
        if(leadPFChargedHadrCand_Ref.isNonnull()){// this check is needed in case hpsTau fails decayModeFinding
            tauLeadChargedHadronExists_.push_back(true);
            tauLeadChargedHadronEta_.push_back(leadPFChargedHadrCand_Ref->eta());
            tauLeadChargedHadronPhi_.push_back(leadPFChargedHadrCand_Ref->phi());
            tauLeadChargedHadronPt_.push_back(leadPFChargedHadrCand_Ref->pt());
        } else {
            tauLeadChargedHadronExists_.push_back(false);
            tauLeadChargedHadronEta_.push_back(0);
            tauLeadChargedHadronPhi_.push_back(0);
            tauLeadChargedHadronPt_.push_back(0);
        }
        
        ++nTau_;
        
    } // loop over tau candidates
    
}
