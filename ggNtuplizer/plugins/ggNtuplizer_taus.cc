// Updated by Abdollah Mohammadi (KSU)  [18 May 2015]
// abdollah.mohammadi@cern.ch

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "ggAnalysis/ggNtuplizer/interface/ggNtuplizer.h"

using namespace std;

// (local) variables associated with tree branches
Int_t          nTau_;

// decay mode discriminators

//Tau Id & Isolation
vector<bool>   taupfTausDiscriminationByDecayModeFinding_;
vector<bool>   taupfTausDiscriminationByDecayModeFindingNewDMs_;
vector<bool>   tauByMVA5VLooseElectronRejection_;
vector<bool>   tauByMVA5LooseElectronRejection_;
vector<bool>   tauByMVA5MediumElectronRejection_;
vector<bool>   tauByMVA5TightElectronRejection_;
vector<bool>   tauByMVA5VTightElectronRejection_;
vector<bool>   tauByLooseMuonRejection3_;
vector<bool>   tauByTightMuonRejection3_;

vector<bool>   tauByLooseCombinedIsolationDeltaBetaCorr3Hits_;
vector<bool>   tauByMediumCombinedIsolationDeltaBetaCorr3Hits_;
vector<bool>   tauByTightCombinedIsolationDeltaBetaCorr3Hits_;
vector<float>   tauCombinedIsolationDeltaBetaCorrRaw3Hits_;

vector<bool>   tauByVLooseIsolationMVA3oldDMwLT_;
vector<bool>   tauByLooseIsolationMVA3oldDMwLT_;
vector<bool>   tauByMediumIsolationMVA3oldDMwLT_;
vector<bool>   tauByTightIsolationMVA3oldDMwLT_;
vector<bool>   tauByVTightIsolationMVA3oldDMwLT_;
vector<bool>   tauByVVTightIsolationMVA3oldDMwLT_;
vector<float>   tauByIsolationMVA3oldDMwLTraw_;

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
vector<int> tauNumSignalPFChargedHadrCands_;
vector<int> tauNumSignalPFNeutrHadrCands_;
vector<int> tauNumSignalPFGammaCands_;
vector<int> tauNumSignalPFCands_;
vector<int> tauNumIsolationPFChargedHadrCands_;
vector<int> tauNumIsolationPFNeutrHadrCands_;
vector<int> tauNumIsolationPFGammaCands_;
vector<int> tauNumIsolationPFCands_;
vector<bool>  tauLeadChargedHadronExists_;
vector<float> tauLeadChargedHadronEta_;
vector<float> tauLeadChargedHadronPhi_;
vector<float> tauLeadChargedHadronPt_;


vector<bool> taubyLoosePileupWeightedIsolation3Hits_;
vector<bool> taubyMediumPileupWeightedIsolation3Hits_;
vector<bool> taubyTightPileupWeightedIsolation3Hits_;
vector<bool> taubyPhotonPtSumOutsideSignalCone_;
vector<bool> taubyPileupWeightedIsolationRaw3Hits_;
vector<float> tauneutralIsoPtSumWeight_;
vector<float> taufootprintCorrection_;
vector<float> tauphotonPtSumOutsideSignalCone_;
vector<float> taudz_;
vector<float> taudxy_;



void ggNtuplizer::branchesTaus(TTree* tree)
{
    
    
    
    tree->Branch("nTau", &nTau_);
    
    //Tau Id & Isolation
    tree->Branch("taupfTausDiscriminationByDecayModeFinding", &taupfTausDiscriminationByDecayModeFinding_);
    tree->Branch("taupfTausDiscriminationByDecayModeFindingNewDMs", &taupfTausDiscriminationByDecayModeFindingNewDMs_);
    
    tree->Branch("tauByMVA5VLooseElectronRejection", &tauByMVA5VLooseElectronRejection_);
    tree->Branch("tauByMVA5LooseElectronRejection", &tauByMVA5LooseElectronRejection_);
    tree->Branch("tauByMVA5MediumElectronRejection", &tauByMVA5MediumElectronRejection_);
    tree->Branch("tauByMVA5TightElectronRejection", &tauByMVA5TightElectronRejection_);
    tree->Branch("tauByMVA5VTightElectronRejection", &tauByMVA5VTightElectronRejection_);
    
    tree->Branch("tauByLooseMuonRejection3", &tauByLooseMuonRejection3_);
    tree->Branch("tauByTightMuonRejection3", &tauByTightMuonRejection3_);
    
    tree->Branch("tauByLooseCombinedIsolationDeltaBetaCorr3Hits", &tauByLooseCombinedIsolationDeltaBetaCorr3Hits_);
    tree->Branch("tauByMediumCombinedIsolationDeltaBetaCorr3Hits", &tauByMediumCombinedIsolationDeltaBetaCorr3Hits_);
    tree->Branch("tauByTightCombinedIsolationDeltaBetaCorr3Hits", &tauByTightCombinedIsolationDeltaBetaCorr3Hits_);
    tree->Branch("tauCombinedIsolationDeltaBetaCorrRaw3Hits", &tauCombinedIsolationDeltaBetaCorrRaw3Hits_);
    
    
    tree->Branch("tauByVLooseIsolationMVA3oldDMwLT", &tauByVLooseIsolationMVA3oldDMwLT_);
    tree->Branch("tauByLooseIsolationMVA3oldDMwLT", &tauByLooseIsolationMVA3oldDMwLT_);
    tree->Branch("tauByMediumIsolationMVA3oldDMwLT", &tauByMediumIsolationMVA3oldDMwLT_);
    tree->Branch("tauByTightIsolationMVA3oldDMwLT", &tauByTightIsolationMVA3oldDMwLT_);
    tree->Branch("tauByVTightIsolationMVA3oldDMwLT", &tauByVTightIsolationMVA3oldDMwLT_);
    tree->Branch("tauByVVTightIsolationMVA3oldDMwLT", &tauByVVTightIsolationMVA3oldDMwLT_);
    tree->Branch("tauByIsolationMVA3oldDMwLTraw", &tauByIsolationMVA3oldDMwLTraw_);
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
    
    tree->Branch("taubyLoosePileupWeightedIsolation3Hits"  ,&taubyLoosePileupWeightedIsolation3Hits_);
    tree->Branch("taubyMediumPileupWeightedIsolation3Hits"  ,&taubyMediumPileupWeightedIsolation3Hits_);
    tree->Branch("taubyTightPileupWeightedIsolation3Hits"  ,&taubyTightPileupWeightedIsolation3Hits_);
    tree->Branch("taubyPhotonPtSumOutsideSignalCone"  ,&taubyPhotonPtSumOutsideSignalCone_);
    tree->Branch("taubyPileupWeightedIsolationRaw3Hits"  ,&taubyPileupWeightedIsolationRaw3Hits_);
    tree->Branch("tauneutralIsoPtSumWeight"  ,&tauneutralIsoPtSumWeight_);
    
    tree->Branch("taufootprintCorrection"  ,&taufootprintCorrection_);
    tree->Branch("tauphotonPtSumOutsideSignalCone"  ,&tauphotonPtSumOutsideSignalCone_);
    tree->Branch("taudz"  ,&taudz_);
    tree->Branch("taudxy"  ,&taudxy_);
    
    
    
    
    
    
    
    
    
}

void ggNtuplizer::fillTaus(const edm::Event& e)
{
    
    // Tau Id & Isolation
    
    tauByMVA5VLooseElectronRejection_.clear();
    tauByMVA5LooseElectronRejection_.clear();
    tauByMVA5MediumElectronRejection_.clear();
    tauByMVA5TightElectronRejection_.clear();
    tauByMVA5VTightElectronRejection_.clear();
    tauByLooseMuonRejection3_.clear();
    tauByTightMuonRejection3_.clear();
    taupfTausDiscriminationByDecayModeFinding_.clear();
    taupfTausDiscriminationByDecayModeFindingNewDMs_.clear();
    tauByLooseCombinedIsolationDeltaBetaCorr3Hits_.clear();
    tauByMediumCombinedIsolationDeltaBetaCorr3Hits_.clear();
    tauByTightCombinedIsolationDeltaBetaCorr3Hits_.clear();
    tauCombinedIsolationDeltaBetaCorrRaw3Hits_.clear();
    tauByVLooseIsolationMVA3oldDMwLT_.clear();
    tauByLooseIsolationMVA3oldDMwLT_.clear();
    tauByMediumIsolationMVA3oldDMwLT_.clear();
    tauByTightIsolationMVA3oldDMwLT_.clear();
    tauByVTightIsolationMVA3oldDMwLT_.clear();
    tauByVVTightIsolationMVA3oldDMwLT_.clear();
    tauByIsolationMVA3oldDMwLTraw_.clear();
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
    
    taubyLoosePileupWeightedIsolation3Hits_.clear();
    taubyMediumPileupWeightedIsolation3Hits_.clear();
    taubyTightPileupWeightedIsolation3Hits_.clear();
    taubyPhotonPtSumOutsideSignalCone_.clear();
    taubyPileupWeightedIsolationRaw3Hits_.clear();
    tauneutralIsoPtSumWeight_.clear();
    taufootprintCorrection_.clear();
    tauphotonPtSumOutsideSignalCone_.clear();
    taudz_.clear();
    taudxy_.clear();
    
    
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
        tauByMVA5VLooseElectronRejection_.push_back(itau->tauID("againstElectronVLooseMVA5"));
        tauByMVA5LooseElectronRejection_.push_back(itau->tauID("againstElectronLooseMVA5"));
        tauByMVA5MediumElectronRejection_.push_back(itau->tauID("againstElectronMediumMVA5"));
        tauByMVA5TightElectronRejection_.push_back(itau->tauID("againstElectronTightMVA5"));
        tauByMVA5VTightElectronRejection_.push_back(itau->tauID("againstElectronVTightMVA5"));
        
        tauByLooseMuonRejection3_.push_back(itau->tauID("againstMuonLoose3"));
        tauByTightMuonRejection3_.push_back(itau->tauID("againstMuonTight3"));
        
        taupfTausDiscriminationByDecayModeFinding_.push_back(itau->tauID("decayModeFinding"));
        taupfTausDiscriminationByDecayModeFindingNewDMs_.push_back(itau->tauID("decayModeFindingNewDMs"));
        
        tauByLooseCombinedIsolationDeltaBetaCorr3Hits_.push_back(itau->tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits"));
        tauByMediumCombinedIsolationDeltaBetaCorr3Hits_.push_back(itau->tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits"));
        tauByTightCombinedIsolationDeltaBetaCorr3Hits_.push_back(itau->tauID("byTightCombinedIsolationDeltaBetaCorr3Hits"));
        tauCombinedIsolationDeltaBetaCorrRaw3Hits_.push_back(itau->tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits"));
        
        tauByVLooseIsolationMVA3oldDMwLT_.push_back(itau->tauID("byVLooseIsolationMVA3oldDMwLT"));
        tauByLooseIsolationMVA3oldDMwLT_.push_back(itau->tauID("byLooseIsolationMVA3oldDMwLT"));
        tauByMediumIsolationMVA3oldDMwLT_.push_back(itau->tauID("byMediumIsolationMVA3oldDMwLT"));
        tauByTightIsolationMVA3oldDMwLT_.push_back(itau->tauID("byTightIsolationMVA3oldDMwLT"));
        tauByVTightIsolationMVA3oldDMwLT_.push_back(itau->tauID("byVTightIsolationMVA3oldDMwLT"));
        tauByVVTightIsolationMVA3oldDMwLT_.push_back(itau->tauID("byVVTightIsolationMVA3oldDMwLT"));
        tauByIsolationMVA3oldDMwLTraw_.push_back(itau->tauID("byIsolationMVA3oldDMwLTraw"));
        
        tauByLooseIsolationMVA3newDMwLT_.push_back(itau->tauID("byLooseIsolationMVA3newDMwLT"));
        tauByVLooseIsolationMVA3newDMwLT_.push_back(itau->tauID("byVLooseIsolationMVA3newDMwLT"));
        tauByMediumIsolationMVA3newDMwLT_.push_back(itau->tauID("byMediumIsolationMVA3newDMwLT"));
        tauByTightIsolationMVA3newDMwLT_.push_back(itau->tauID("byTightIsolationMVA3newDMwLT"));
        tauByVTightIsolationMVA3newDMwLT_.push_back(itau->tauID("byVTightIsolationMVA3newDMwLT"));
        tauByVVTightIsolationMVA3newDMwLT_.push_back(itau->tauID("byVVTightIsolationMVA3newDMwLT"));
        tauByIsolationMVA3newDMwLTraw_.push_back(itau->tauID("byIsolationMVA3newDMwLTraw"));

        taubyLoosePileupWeightedIsolation3Hits_.push_back(itau->tauID("byLoosePileupWeightedIsolation3Hits"));
        taubyMediumPileupWeightedIsolation3Hits_.push_back(itau->tauID("byMediumPileupWeightedIsolation3Hits"));
        taubyTightPileupWeightedIsolation3Hits_.push_back(itau->tauID("byTightPileupWeightedIsolation3Hits"));
        taubyPhotonPtSumOutsideSignalCone_.push_back(itau->tauID("byPhotonPtSumOutsideSignalCone"));
        taubyPileupWeightedIsolationRaw3Hits_.push_back(itau->tauID("byPileupWeightedIsolationRaw3Hits"));

        
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
        tauChargedIsoPtSum_.push_back(itau->tauID("chargedIsoPtSum") );
        tauNeutralIsoPtSum_.push_back(itau->tauID("neutralIsoPtSum")  );
        tauPuCorrPtSum_.push_back(itau->tauID("puCorrPtSum")  );
        tauneutralIsoPtSumWeight_.push_back(itau->tauID("neutralIsoPtSumWeight"));
        taufootprintCorrection_.push_back(itau->tauID("footprintCorrection"));
        tauphotonPtSumOutsideSignalCone_.push_back(itau->tauID("photonPtSumOutsideSignalCone"));
        
        tauNumSignalPFChargedHadrCands_.push_back(itau->signalChargedHadrCands().size());
        tauNumSignalPFNeutrHadrCands_.push_back(itau->signalNeutrHadrCands().size());
        tauNumSignalPFGammaCands_.push_back(itau->signalGammaCands().size());
        tauNumSignalPFCands_.push_back(itau->signalCands().size());
        
        tauNumIsolationPFChargedHadrCands_.push_back(itau->isolationChargedHadrCands().size());
        tauNumIsolationPFNeutrHadrCands_.push_back(itau->isolationNeutrHadrCands().size());
        tauNumIsolationPFGammaCands_.push_back(itau->isolationGammaCands().size());
        tauNumIsolationPFCands_.push_back(itau->isolationCands().size());
        
         edm::Handle<reco::VertexCollection> vertexs;
         e.getByToken(vtxLabel_, vertexs);

        if (vertexs->size()>0) {
            pat::PackedCandidate const* packedLeadTauCand = dynamic_cast<pat::PackedCandidate const*>(itau->leadChargedHadrCand().get());
            taudz_.push_back(packedLeadTauCand->dz());
            taudxy_.push_back(packedLeadTauCand->dxy());
            tauLeadChargedHadronExists_.push_back(true);
            tauLeadChargedHadronEta_.push_back(packedLeadTauCand->eta());
            tauLeadChargedHadronPhi_.push_back(packedLeadTauCand->phi());
            tauLeadChargedHadronPt_.push_back(packedLeadTauCand->pt());
        }
        
        
        
        
        ++nTau_;
        
    } // loop over tau candidates
    
}
