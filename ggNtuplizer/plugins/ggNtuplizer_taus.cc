#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "ggAnalysis/ggNtuplizer/interface/ggNtuplizer.h"

using namespace std;

// (local) variables associated with tree branches
Int_t          nTau_;

// decay mode discriminators
vector<bool>   tauByLooseElectronRejection_;
vector<bool>   tauByMediumElectronRejection_;
vector<bool>   tauByTightElectronRejection_;
vector<bool>   tauByMVA5LooseElectronRejection_;
vector<bool>   tauByMVA5MediumElectronRejection_;
vector<bool>   tauByMVA5TightElectronRejection_;
vector<bool>   tauByMVA5VTightElectronRejection_;
vector<bool>   tauByLooseMuonRejection_;
vector<bool>   tauByMediumMuonRejection_;
vector<bool>   tauByTightMuonRejection_;
vector<bool>   tauByLooseMuonRejection3_;
vector<bool>   tauByTightMuonRejection3_;
vector<bool>   tauByMVALooseMuonRejection_;
vector<bool>   tauByMVAMediumMuonRejection_;
vector<bool>   tauByMVATightMuonRejection_;
vector<bool>   tauByMVArawMuonRejection_;
vector<bool>   pfTausDiscriminationByDecayModeFinding_;
vector<bool>   tauByVLooseIsolation_;
vector<bool>   tauByVLooseCombinedIsolationDBSumPtCorr_;
vector<bool>   tauByLooseCombinedIsolationDBSumPtCorr_;
vector<bool>   tauByMediumCombinedIsolationDBSumPtCorr_;
vector<bool>   tauByTightCombinedIsolationDBSumPtCorr_;
vector<bool>   tauByLooseCombinedIsolationDBSumPtCorr3Hits_;
vector<bool>   tauByMediumCombinedIsolationDBSumPtCorr3Hits_;
vector<bool>   tauByTightCombinedIsolationDBSumPtCorr3Hits_;
vector<bool>   tauByVLooseIsolationMVA3newDMwoLT_;
vector<bool>   tauByLooseIsolationMVA3newDMwoLT_;
vector<bool>   tauByMediumIsolationMVA3newDMwoLT_;
vector<bool>   tauByTightIsolationMVA3newDMwoLT_;
vector<bool>   tauByVTightIsolationMVA3newDMwoLT_;
vector<bool>   tauByVVTightIsolationMVA3newDMwoLT_;
vector<bool>   tauByIsolationMVA3newDMwoLTraw_;
vector<bool>   tauByVLooseIsolationMVA3oldDMwLT_;
vector<bool>   tauByLooseIsolationMVA3oldDMwLT_;
vector<bool>   tauByMediumIsolationMVA3oldDMwLT_;
vector<bool>   tauByTightIsolationMVA3oldDMwLT_;
vector<bool>   tauByVTightIsolationMVA3oldDMwLT_;
vector<bool>   tauByVVTightIsolationMVA3oldDMwLT_;
vector<bool>   tauByIsolationMVA3oldDMwLTraw_;
vector<bool>   tauByVLooseIsolationMVA3oldDMwoLT_;
vector<bool>   tauByLooseIsolationMVA3oldDMwoLT_;
vector<bool>   tauByTightIsolationMVA3oldDMwoLT_;
vector<bool>   tauByVTightIsolationMVA3oldDMwoLT_;
vector<bool>   tauByVVTightIsolationMVA3oldDMwoLT_;
vector<bool>   tauByIsolationMVA3oldDMwoLTraw_;
vector<bool>   tauByLooseIsolationMVA3newDMwLT_;
vector<bool>   tauByVLooseIsolationMVA3newDMwLT_;
vector<bool>   tauByMediumIsolationMVA3newDMwLT_;
vector<bool>   tauByTightIsolationMVA3newDMwLT_;
vector<bool>   tauByVTightIsolationMVA3newDMwLT_;
vector<bool>   tauByVVTightIsolationMVA3newDMwLT_;
vector<bool>   tauByIsolationMVA3newDMwLTraw_;

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

void ggNtuplizer::branchesTaus(TTree* tree)
{
  tree->Branch("nTau", &nTau_);
  tree->Branch("tauByLooseElectronRejection", &tauByLooseElectronRejection_);
  tree->Branch("tauByMediumElectronRejection", &tauByMediumElectronRejection_);
  tree->Branch("tauByTightElectronRejection", &tauByTightElectronRejection_);
  tree->Branch("tauByMVA5LooseElectronRejection", &tauByMVA5LooseElectronRejection_);
  tree->Branch("tauByMVA5MediumElectronRejection", &tauByMVA5MediumElectronRejection_);
  tree->Branch("tauByMVA5TightElectronRejection", &tauByMVA5TightElectronRejection_);
  tree->Branch("tauByMVA5VTightElectronRejection", &tauByMVA5VTightElectronRejection_);
  tree->Branch("tauByLooseMuonRejection", &tauByLooseMuonRejection_);
  tree->Branch("tauByMediumMuonRejection", &tauByMediumMuonRejection_);
  tree->Branch("tauByTightMuonRejection", &tauByTightMuonRejection_);
  tree->Branch("tauByLooseMuonRejection3", &tauByLooseMuonRejection3_);
  tree->Branch("tauByTightMuonRejection3", &tauByTightMuonRejection3_);
  tree->Branch("tauByMVALooseMuonRejection", &tauByMVALooseMuonRejection_);
  tree->Branch("tauByMVAMediumMuonRejection", &tauByMVAMediumMuonRejection_);
  tree->Branch("tauByMVATightMuonRejection", &tauByMVATightMuonRejection_);
  tree->Branch("tauByMVArawMuonRejection", &tauByMVArawMuonRejection_);
  //    tree->Branch("pfTausDiscriminationByDecayModeFinding", &pfTausDiscriminationByDecayModeFinding_);
  tree->Branch("tauByVLooseIsolation", &tauByVLooseIsolation_);
  tree->Branch("tauByVLooseCombinedIsolationDBSumPtCorr", &tauByVLooseCombinedIsolationDBSumPtCorr_);
  tree->Branch("tauByLooseCombinedIsolationDBSumPtCorr", &tauByLooseCombinedIsolationDBSumPtCorr_);
  tree->Branch("tauByMediumCombinedIsolationDBSumPtCorr", &tauByMediumCombinedIsolationDBSumPtCorr_);
  tree->Branch("tauByTightCombinedIsolationDBSumPtCorr", &tauByTightCombinedIsolationDBSumPtCorr_);
  tree->Branch("tauByLooseCombinedIsolationDBSumPtCorr3Hits", &tauByLooseCombinedIsolationDBSumPtCorr3Hits_);
  tree->Branch("tauByMediumCombinedIsolationDBSumPtCorr3Hits", &tauByMediumCombinedIsolationDBSumPtCorr3Hits_);
  tree->Branch("tauByTightCombinedIsolationDBSumPtCorr3Hits", &tauByTightCombinedIsolationDBSumPtCorr3Hits_);
  tree->Branch("tauByVLooseIsolationMVA3newDMwoLT", &tauByVLooseIsolationMVA3newDMwoLT_);
  tree->Branch("tauByLooseIsolationMVA3newDMwoLT", &tauByLooseIsolationMVA3newDMwoLT_);
  tree->Branch("tauByMediumIsolationMVA3newDMwoLT", &tauByMediumIsolationMVA3newDMwoLT_);
  tree->Branch("tauByTightIsolationMVA3newDMwoLT", &tauByTightIsolationMVA3newDMwoLT_);
  tree->Branch("tauByVTightIsolationMVA3newDMwoLT", &tauByVTightIsolationMVA3newDMwoLT_);
  tree->Branch("tauByVVTightIsolationMVA3newDMwoLT", &tauByVVTightIsolationMVA3newDMwoLT_);
  tree->Branch("tauByIsolationMVA3newDMwoLTraw", &tauByIsolationMVA3newDMwoLTraw_);
  tree->Branch("tauByVLooseIsolationMVA3oldDMwLT", &tauByVLooseIsolationMVA3oldDMwLT_);
  tree->Branch("tauByLooseIsolationMVA3oldDMwLT", &tauByLooseIsolationMVA3oldDMwLT_);
  tree->Branch("tauByMediumIsolationMVA3oldDMwLT", &tauByMediumIsolationMVA3oldDMwLT_);
  tree->Branch("tauByTightIsolationMVA3oldDMwLT", &tauByTightIsolationMVA3oldDMwLT_);
  tree->Branch("tauByVTightIsolationMVA3oldDMwLT", &tauByVTightIsolationMVA3oldDMwLT_);
  tree->Branch("tauByVVTightIsolationMVA3oldDMwLT", &tauByVVTightIsolationMVA3oldDMwLT_);
  tree->Branch("tauByIsolationMVA3oldDMwLTraw", &tauByIsolationMVA3oldDMwLTraw_);
  tree->Branch("tauByVLooseIsolationMVA3oldDMwoLT", &tauByVLooseIsolationMVA3oldDMwoLT_);
  tree->Branch("tauByLooseIsolationMVA3oldDMwoLT", &tauByLooseIsolationMVA3oldDMwoLT_);
  tree->Branch("tauByTightIsolationMVA3oldDMwoLT", &tauByTightIsolationMVA3oldDMwoLT_);
  tree->Branch("tauByVTightIsolationMVA3oldDMwoLT", &tauByVTightIsolationMVA3oldDMwoLT_);
  tree->Branch("tauByVVTightIsolationMVA3oldDMwoLT", &tauByVVTightIsolationMVA3oldDMwoLT_);
  tree->Branch("tauByIsolationMVA3newDMwoLTraw", &tauByIsolationMVA3newDMwoLTraw_);
  tree->Branch("tauByLooseIsolationMVA3newDMwLT", &tauByLooseIsolationMVA3newDMwLT_);
  tree->Branch("tauByVLooseIsolationMVA3newDMwLT", &tauByVLooseIsolationMVA3newDMwLT_);
  tree->Branch("tauByMediumIsolationMVA3newDMwLT", &tauByMediumIsolationMVA3newDMwLT_);
  tree->Branch("tauByTightIsolationMVA3newDMwLT", &tauByTightIsolationMVA3newDMwLT_);
  tree->Branch("tauByVTightIsolationMVA3newDMwLT", &tauByVTightIsolationMVA3newDMwLT_);
  tree->Branch("tauByVVTightIsolationMVA3newDMwLT", &tauByVVTightIsolationMVA3newDMwLT_);
  tree->Branch("tauByIsolationMVA3newDMwLTraw", &tauByIsolationMVA3newDMwLTraw_);

  tree->Branch("tauEta"  ,&tauEta_);
  tree->Branch("tauPhi"  ,&tauPhi_);
  tree->Branch("tauPt"  ,&tauPt_);
  tree->Branch("tauEt"  ,&tauEt_);
  tree->Branch("tauCharge"  ,&tauCharge_);
  tree->Branch("tauDecayMode"  ,&tauDecayMode_);
  tree->Branch("tauEMFraction"  ,&tauEMFraction_);
  tree->Branch("tauHCAL3x3OverPLead"  ,&tauHCAL3x3OverPLead_);
  tree->Branch("tauHCALMaxOverPLead"  ,&tauHCALMaxOverPLead_);
  tree->Branch("tauHCALTotOverPLead"  ,&tauHCALTotOverPLead_);
  tree->Branch("tauIsolationPFChargedHadrCandsPtSum"  ,&tauIsolationPFChargedHadrCandsPtSum_);
  tree->Branch("tauIsolationPFGammaCandsEtSum"  ,&tauIsolationPFGammaCandsEtSum_);
  tree->Branch("tauLeadPFChargedHadrCandsignedSipt"  ,&tauLeadPFChargedHadrCandsignedSipt_);
  tree->Branch("tauLeadChargedHadronExists"  ,&tauLeadChargedHadronExists_);
  tree->Branch("tauLeadChargedHadronEta"  ,&tauLeadChargedHadronEta_);
  tree->Branch("tauLeadChargedHadronPhi"  ,&tauLeadChargedHadronPhi_);
  tree->Branch("tauLeadChargedHadronPt"  ,&tauLeadChargedHadronPt_);

}

void ggNtuplizer::fillTaus(const edm::Event& e)
{

  // cleanup from previous execution
  //@Lovedeep: 09/14 (Updated TauId Info)
  tauByLooseElectronRejection_.clear();
  tauByMediumElectronRejection_.clear();
  tauByTightElectronRejection_.clear();
  tauByMVA5LooseElectronRejection_.clear();
  tauByMVA5MediumElectronRejection_.clear();
  tauByMVA5TightElectronRejection_.clear();
  tauByMVA5VTightElectronRejection_.clear();
  tauByLooseMuonRejection_.clear();
  tauByMediumMuonRejection_.clear();
  tauByTightMuonRejection_.clear();
  tauByLooseMuonRejection3_.clear();
  tauByTightMuonRejection3_.clear();
  tauByMVALooseMuonRejection_.clear();
  tauByMVAMediumMuonRejection_.clear();
  tauByMVATightMuonRejection_.clear();
  tauByMVArawMuonRejection_.clear();
  pfTausDiscriminationByDecayModeFinding_.clear();
  tauByVLooseIsolation_.clear();
  tauByVLooseCombinedIsolationDBSumPtCorr_.clear();
  tauByLooseCombinedIsolationDBSumPtCorr_.clear();
  tauByMediumCombinedIsolationDBSumPtCorr_.clear();
  tauByTightCombinedIsolationDBSumPtCorr_.clear();
  tauByLooseCombinedIsolationDBSumPtCorr3Hits_.clear();
  tauByMediumCombinedIsolationDBSumPtCorr3Hits_.clear();
  tauByTightCombinedIsolationDBSumPtCorr3Hits_.clear();
  tauByVLooseIsolationMVA3newDMwoLT_.clear();
  tauByLooseIsolationMVA3newDMwoLT_.clear();
  tauByMediumIsolationMVA3newDMwoLT_.clear();
  tauByTightIsolationMVA3newDMwoLT_.clear();
  tauByVTightIsolationMVA3newDMwoLT_.clear();
  tauByVVTightIsolationMVA3newDMwoLT_.clear();
  tauByIsolationMVA3newDMwoLTraw_.clear();
  tauByVLooseIsolationMVA3oldDMwLT_.clear();
  tauByLooseIsolationMVA3oldDMwLT_.clear();
  tauByMediumIsolationMVA3oldDMwLT_.clear();
  tauByTightIsolationMVA3oldDMwLT_.clear();
  tauByVTightIsolationMVA3oldDMwLT_.clear();
  tauByVVTightIsolationMVA3oldDMwLT_.clear();
  tauByIsolationMVA3oldDMwLTraw_.clear();
  tauByVLooseIsolationMVA3oldDMwoLT_.clear();
  tauByLooseIsolationMVA3oldDMwoLT_.clear();
  tauByTightIsolationMVA3oldDMwoLT_.clear();
  tauByVTightIsolationMVA3oldDMwoLT_.clear();
  tauByVVTightIsolationMVA3oldDMwoLT_.clear();
  tauByIsolationMVA3newDMwoLTraw_.clear();
  tauByLooseIsolationMVA3newDMwLT_.clear();
  tauByVLooseIsolationMVA3newDMwLT_.clear();
  tauByMediumIsolationMVA3newDMwLT_.clear();
  tauByTightIsolationMVA3newDMwLT_.clear();
  tauByVTightIsolationMVA3newDMwLT_.clear();
  tauByVVTightIsolationMVA3newDMwLT_.clear();
  tauByIsolationMVA3newDMwLTraw_.clear();

  tauEta_.clear();
  tauPhi_.clear();
  tauPt_.clear();
  tauEt_.clear();
  tauCharge_.clear();
  tauDecayMode_.clear();
  tauEMFraction_.clear();
  tauHCAL3x3OverPLead_.clear();
  tauHCALMaxOverPLead_.clear();
  tauHCALTotOverPLead_.clear();
  tauIsolationPFChargedHadrCandsPtSum_.clear();
  tauIsolationPFGammaCandsEtSum_.clear();
  tauLeadPFChargedHadrCandsignedSipt_.clear();
  tauLeadChargedHadronExists_.clear();
  tauLeadChargedHadronEta_.clear();
  tauLeadChargedHadronPhi_.clear();
  tauLeadChargedHadronPt_.clear();

  nTau_ = 0;

  edm::Handle<vector<pat::Tau> > tauHandle;
  e.getByToken(tauCollection_, tauHandle);

  if (!tauHandle.isValid()) {
    edm::LogWarning("ggNtuplizer") << "no pat::Tau in event";
    return;
  }

  //startTaus Lvdp
  for(vector<pat::Tau>::const_iterator itau = tauHandle->begin(); itau != tauHandle->end(); ++itau) {

    tauByLooseElectronRejection_.push_back(itau->tauID("byLooseElectronRejection"));
    tauByMediumElectronRejection_.push_back(itau->tauID("byMediumElectronRejection"));
    tauByTightElectronRejection_.push_back(itau->tauID("byTightElectronRejection"));
    tauByMVA5LooseElectronRejection_.push_back(itau->tauID("byMVA5LooseElectronRejection"));
    tauByMVA5MediumElectronRejection_.push_back(itau->tauID("byMVA5MediumElectronRejection"));
    tauByMVA5TightElectronRejection_.push_back(itau->tauID("byMVA5TightElectronRejection"));
    tauByMVA5VTightElectronRejection_.push_back(itau->tauID("byMVA5VTightElectronRejection"));
    tauByLooseMuonRejection_.push_back(itau->tauID("byLooseMuonRejection"));
    tauByMediumMuonRejection_.push_back(itau->tauID("byMediumMuonRejection"));
    tauByTightMuonRejection_.push_back(itau->tauID("byTightMuonRejection"));
    tauByLooseMuonRejection3_.push_back(itau->tauID("byLooseMuonRejection3"));
    tauByTightMuonRejection3_.push_back(itau->tauID("byTightMuonRejection3"));
    tauByMVALooseMuonRejection_.push_back(itau->tauID("byMVALooseMuonRejection"));
    tauByMVAMediumMuonRejection_.push_back(itau->tauID("byMVAMediumMuonRejection"));
    tauByMVATightMuonRejection_.push_back(itau->tauID("byMVATightMuonRejection"));
    tauByMVArawMuonRejection_.push_back(itau->tauID("byMVArawMuonRejection"));
    //  pfTausDiscriminationByDecayModeFinding_.push_back(itau->tauID("pfTausDiscriminationByDecayModeFinding"));
    tauByVLooseIsolation_.push_back(itau->tauID("byVLooseIsolation"));
    tauByVLooseCombinedIsolationDBSumPtCorr_.push_back(itau->tauID("byVLooseCombinedIsolationDBSumPtCorr"));
    tauByLooseCombinedIsolationDBSumPtCorr_.push_back(itau->tauID("byLooseCombinedIsolationDBSumPtCorr"));
    tauByMediumCombinedIsolationDBSumPtCorr_.push_back(itau->tauID("byMediumCombinedIsolationDBSumPtCorr"));
    tauByTightCombinedIsolationDBSumPtCorr_.push_back(itau->tauID("byTightCombinedIsolationDBSumPtCorr"));
    tauByLooseCombinedIsolationDBSumPtCorr3Hits_.push_back(itau->tauID("byLooseCombinedIsolationDBSumPtCorr3Hits"));
    tauByMediumCombinedIsolationDBSumPtCorr3Hits_.push_back(itau->tauID("byMediumCombinedIsolationDBSumPtCorr3Hits"));
    tauByTightCombinedIsolationDBSumPtCorr3Hits_.push_back(itau->tauID("byTightCombinedIsolationDBSumPtCorr3Hits"));
    tauByVLooseIsolationMVA3newDMwoLT_.push_back(itau->tauID("byVLooseIsolationMVA3newDMwoLT"));
    tauByLooseIsolationMVA3newDMwoLT_.push_back(itau->tauID("byLooseIsolationMVA3newDMwoLT"));
    tauByMediumIsolationMVA3newDMwoLT_.push_back(itau->tauID("byMediumIsolationMVA3newDMwoLT"));
    tauByTightIsolationMVA3newDMwoLT_.push_back(itau->tauID("byTightIsolationMVA3newDMwoLT"));
    tauByVTightIsolationMVA3newDMwoLT_.push_back(itau->tauID("byVTightIsolationMVA3newDMwoLT"));
    tauByVVTightIsolationMVA3newDMwoLT_.push_back(itau->tauID("byVVTightIsolationMVA3newDMwoLT"));
    tauByIsolationMVA3newDMwoLTraw_.push_back(itau->tauID("byIsolationMVA3newDMwoLTraw"));
    tauByVLooseIsolationMVA3oldDMwLT_.push_back(itau->tauID("byVLooseIsolationMVA3oldDMwLT"));
    tauByLooseIsolationMVA3oldDMwLT_.push_back(itau->tauID("byLooseIsolationMVA3oldDMwLT"));
    tauByMediumIsolationMVA3oldDMwLT_.push_back(itau->tauID("byMediumIsolationMVA3oldDMwLT"));
    tauByTightIsolationMVA3oldDMwLT_.push_back(itau->tauID("byTightIsolationMVA3oldDMwLT"));
    tauByVTightIsolationMVA3oldDMwLT_.push_back(itau->tauID("byVTightIsolationMVA3oldDMwLT"));
    tauByVVTightIsolationMVA3oldDMwLT_.push_back(itau->tauID("byVVTightIsolationMVA3oldDMwLT"));
    tauByIsolationMVA3oldDMwLTraw_.push_back(itau->tauID("byIsolationMVA3oldDMwLTraw"));
    tauByVLooseIsolationMVA3oldDMwoLT_.push_back(itau->tauID("byVLooseIsolationMVA3oldDMwoLT"));
    tauByLooseIsolationMVA3oldDMwoLT_.push_back(itau->tauID("byLooseIsolationMVA3oldDMwoLT"));
    tauByTightIsolationMVA3oldDMwoLT_.push_back(itau->tauID("byTightIsolationMVA3oldDMwoLT"));
    tauByVTightIsolationMVA3oldDMwoLT_.push_back(itau->tauID("byVTightIsolationMVA3oldDMwoLT"));
    tauByVVTightIsolationMVA3oldDMwoLT_.push_back(itau->tauID("byVVTightIsolationMVA3oldDMwoLT"));
    tauByIsolationMVA3newDMwoLTraw_.push_back(itau->tauID("byIsolationMVA3newDMwoLTraw"));
    tauByLooseIsolationMVA3newDMwLT_.push_back(itau->tauID("byLooseIsolationMVA3newDMwLT"));
    tauByVLooseIsolationMVA3newDMwLT_.push_back(itau->tauID("byVLooseIsolationMVA3newDMwLT"));
    tauByMediumIsolationMVA3newDMwLT_.push_back(itau->tauID("byMediumIsolationMVA3newDMwLT"));
    tauByTightIsolationMVA3newDMwLT_.push_back(itau->tauID("byTightIsolationMVA3newDMwLT"));
    tauByVTightIsolationMVA3newDMwLT_.push_back(itau->tauID("byVTightIsolationMVA3newDMwLT"));
    tauByVVTightIsolationMVA3newDMwLT_.push_back(itau->tauID("byVVTightIsolationMVA3newDMwLT"));
    tauByIsolationMVA3newDMwLTraw_.push_back(itau->tauID("byIsolationMVA3newDMwLTraw"));

    tauEta_.push_back(itau->eta());
    tauPhi_.push_back(itau->phi());
    tauPt_ .push_back(itau->pt());
    tauEt_ .push_back(itau->et());
    tauCharge_.push_back(itau->charge());
    tauDecayMode_.push_back(itau->decayMode());
    tauEMFraction_.push_back(itau->emFraction());
    tauHCAL3x3OverPLead_.push_back(itau->hcal3x3OverPLead());
    tauHCALMaxOverPLead_.push_back(itau->hcalMaxOverPLead());
    tauHCALTotOverPLead_.push_back(itau->hcalTotOverPLead());
    tauIsolationPFChargedHadrCandsPtSum_.push_back(itau->isolationPFChargedHadrCandsPtSum());
    tauIsolationPFGammaCandsEtSum_.push_back(itau->isolationPFGammaCandsEtSum());
    tauLeadPFChargedHadrCandsignedSipt_.push_back(itau->leadPFChargedHadrCandsignedSipt());
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
