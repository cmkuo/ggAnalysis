#include "../interface/ggNtuplizer.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/ESGetToken.h"
//#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"

using namespace std;
using namespace edm;

/////EB                                                                                                                                   
Int_t nebRechit_;
vector<Float_t> ebRechitE_;
vector<Int_t> ebRechitiEta_;
vector<Int_t> ebRechitiPhi_;
vector<Int_t> ebRechitZside_;
vector<Float_t> ebRechitEta_;
vector<Float_t> ebRechitPhi_;
vector<Float_t> ebRechitX_;
vector<Float_t> ebRechitY_;
vector<Float_t> ebRechitZ_;
vector<Int_t> ebRechitDetID_;

/////EE                                                                                                                                   
Int_t neeRechit_;
vector<Float_t> eeRechitE_;
vector<Int_t> eeRechitiEta_;
vector<Int_t> eeRechitiPhi_;
vector<Int_t> eeRechitZside_;
vector<Float_t> eeRechitEta_;
vector<Float_t> eeRechitPhi_;
vector<Float_t> eeRechitX_;
vector<Float_t> eeRechitY_;
vector<Float_t> eeRechitZ_;
vector<Int_t> eeRechitDetID_;

/////ES                                                                                                                                   
Int_t nesRechit_;
vector<Float_t> esRechitE_;
vector<Int_t> esRechitiEta_;
vector<Int_t> esRechitiPhi_;
vector<Int_t> esRechitZside_;
vector<Int_t> esRechitPlane_;
vector<Int_t> esRechitStrip_;
vector<Float_t> esRechitEta_;
vector<Float_t> esRechitPhi_;
vector<Float_t> esRechitX_;
vector<Float_t> esRechitY_;
vector<Float_t> esRechitZ_;


void ggNtuplizer::branchesRechits(TTree* tree) {

  //EB
  if(dumpCrystalinfo_){
  tree->Branch("nebRechit",         &nebRechit_);
  tree->Branch("ebRechitE",         &ebRechitE_);
  tree->Branch("ebRechitiEta",         &ebRechitiEta_);
  tree->Branch("ebRechitiPhi",         &ebRechitiPhi_);
  tree->Branch("ebRechitZside",         &ebRechitZside_);
  tree->Branch("ebRechitEta",         &ebRechitEta_);
  tree->Branch("ebRechitPhi",         &ebRechitPhi_);
  tree->Branch("ebRechitX",         &ebRechitX_);
  tree->Branch("ebRechitY",         &ebRechitY_);
  tree->Branch("ebRechitZ",         &ebRechitZ_);
  tree->Branch("ebRechitDetID",         &ebRechitDetID_);

  ///EE
  tree->Branch("neeRechit",         &neeRechit_);
  tree->Branch("eeRechitE",         &eeRechitE_);
  tree->Branch("eeRechitiEta",         &eeRechitiEta_);
  tree->Branch("eeRechitiPhi",         &eeRechitiPhi_);
  tree->Branch("eeRechitZside",         &eeRechitZside_);
  tree->Branch("eeRechitEta",         &eeRechitEta_);
  tree->Branch("eeRechitPhi",         &eeRechitPhi_);
  tree->Branch("eeRechitX",         &eeRechitX_);
  tree->Branch("eeRechitY",         &eeRechitY_);
  tree->Branch("eeRechitZ",         &eeRechitZ_);
  tree->Branch("eeRechitDetID",         &eeRechitDetID_);

  ///ES
  tree->Branch("nesRechit",         &nesRechit_);
  tree->Branch("esRechitE",         &esRechitE_);
  tree->Branch("esRechitiEta",         &esRechitiEta_);
  tree->Branch("esRechitiPhi",         &esRechitiPhi_);
  tree->Branch("esRechitZside",         &esRechitZside_);
  tree->Branch("esRechitPlane",         &esRechitPlane_);
  tree->Branch("esRechitStrip",         &esRechitStrip_);
  tree->Branch("esRechitEta",         &esRechitEta_);
  tree->Branch("esRechitPhi",         &esRechitPhi_);
  tree->Branch("esRechitX",         &esRechitX_);
  tree->Branch("esRechitY",         &esRechitY_);
  tree->Branch("esRechitZ",         &esRechitZ_);}
}

//
// member functions
//

// ------------ method called for each event  ------------
void ggNtuplizer::fillRechits(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  //bool debugRH = true;
  bool debugRH = false;

  //////////////////////////////ECAL rechits/////////////////////////////////////////
  //EB
 
  nebRechit_ = 0;  

  ebRechitE_.clear();                 
  ebRechitiEta_.clear();
  ebRechitiPhi_.clear();
  ebRechitZside_.clear();

  ebRechitEta_.clear();
  ebRechitPhi_.clear();
  ebRechitX_.clear();
  ebRechitY_.clear();
  ebRechitZ_.clear();
  ebRechitDetID_.clear();
  
  //EE
  neeRechit_ = 0;  

  eeRechitE_.clear();                 
  eeRechitiEta_.clear();
  eeRechitiPhi_.clear();
  eeRechitZside_.clear();

  eeRechitEta_.clear();
  eeRechitPhi_.clear();
  eeRechitX_.clear();
  eeRechitY_.clear();
  eeRechitZ_.clear();
  eeRechitDetID_.clear();

  //ES
  nesRechit_ = 0;  

  esRechitE_.clear();                 
  esRechitiEta_.clear();
  esRechitiPhi_.clear();
  esRechitZside_.clear();
  esRechitPlane_.clear();
  esRechitStrip_.clear();

  esRechitEta_.clear();
  esRechitPhi_.clear();
  esRechitX_.clear();
  esRechitY_.clear();
  esRechitZ_.clear();

  iSetup.get<CaloGeometryRecord>().get(pG_);
  const CaloGeometry* geo = pG_.product();

  edm::Handle<EcalRecHitCollection> barrelRecHitsHandle;
  edm::Handle<EcalRecHitCollection> endcapRecHitsHandle;
  edm::Handle<EcalRecHitCollection> esRecHitsHandle;
  iEvent.getByToken(ebReducedRecHitCollection_,barrelRecHitsHandle);
  iEvent.getByToken(eeReducedRecHitCollection_,endcapRecHitsHandle);
  iEvent.getByToken(esReducedRecHitCollection_,esRecHitsHandle);

  const EcalRecHitCollection* EBRecHits = nullptr;
  const EcalRecHitCollection* EERecHits = nullptr;
  const EcalRecHitCollection* ESRecHits = nullptr;
  
  if(debugRH){
    std::cout<<"Getting ECAL RECHIT"<<std::endl;
  }

  if ( !barrelRecHitsHandle.isValid() ){
    LogDebug("") << "ggNtuplizer_rechits: Error! EB rechits can't get product!" << std::endl;
  } else{
    EBRecHits = barrelRecHitsHandle.product();
  }

  if ( !endcapRecHitsHandle.isValid() ){
    LogDebug("") << "ggNtuplizer_rechits: Error! EE rechits can't get product!" << std::endl;
  } else{
    EERecHits = endcapRecHitsHandle.product();
  }
  
  if ( !esRecHitsHandle.isValid() ){
    LogDebug("") << "ggNtuplizer_rechits: Error! ES rechits can't get product!" << std::endl;
  } else{
    ESRecHits = esRecHitsHandle.product();
  }

  if(debugRH){
    std::cout<<"Looping EB"<<std::endl;
  }

  
  //EB
  EcalRecHitCollection::const_iterator ebrechit;

  if(debugRH){
    std::cout<<"Start looping EB now"<<std::endl;
  }

  if(debugRH) std::cout<<"EB EH size "<<EBRecHits->size()<<std::endl;

  for( ebrechit = EBRecHits->begin(); ebrechit != EBRecHits->end(); ebrechit++ ){

    if(debugRH){
      std::cout<<"1st EB rechit"<<std::endl;
    }

    double Energy = ebrechit->energy();
    if(debugRH){
      std::cout<<"E "<<Energy<<std::endl;
    }

    EBDetId det = ebrechit->id();
    if(debugRH){
      std::cout<<"detID "<<det<<std::endl;
    }

    int ieta = det.ieta();
    if(debugRH){
      std::cout<<"ieta "<<ieta<<std::endl;
    }


    int iphi = det.iphi();
    int zside = det.zside();

    uint32_t detid = det.rawId();
    if(dumpCrystalinfo_){
    ebRechitE_.push_back(Energy);
    ebRechitiEta_.push_back(ieta);
    ebRechitiPhi_.push_back(iphi);
    ebRechitZside_.push_back(zside);
    
    ebRechitDetID_.push_back(detid);}
    
    const GlobalPoint & rechitPoint = geo->getPosition(det);

    if(dumpCrystalinfo_){
    ebRechitEta_.push_back(rechitPoint.eta());}

    if(debugRH){
      std::cout<<"eta "<<rechitPoint.eta()<<std::endl;
    }

    if(dumpCrystalinfo_){
    ebRechitPhi_.push_back(rechitPoint.phi());
    ebRechitX_.push_back(rechitPoint.x());
    ebRechitY_.push_back(rechitPoint.y());
    ebRechitZ_.push_back(rechitPoint.z());}

    nebRechit_++;
  }


  if(debugRH){
    std::cout<<"Looping EE"<<std::endl;
  }
  ///EE
  EcalRecHitCollection::const_iterator eerechit;
  for( eerechit = EERecHits->begin(); eerechit != EERecHits->end(); eerechit++ ){

    double Energy = eerechit->energy();
    EEDetId det = eerechit->id();
    int ieta = det.ix();
    int iphi = det.iy();
    int zside = det.zside();

    uint32_t detid = det.rawId();

    if(dumpCrystalinfo_){
    eeRechitE_.push_back(Energy);
    eeRechitiEta_.push_back(ieta);
    eeRechitiPhi_.push_back(iphi);
    eeRechitZside_.push_back(zside);
    eeRechitDetID_.push_back(detid);}
    
    const GlobalPoint & rechitPoint = geo->getPosition(det);
    if(dumpCrystalinfo_){
    eeRechitEta_.push_back(rechitPoint.eta());
    eeRechitPhi_.push_back(rechitPoint.phi());
    eeRechitX_.push_back(rechitPoint.x());
    eeRechitY_.push_back(rechitPoint.y());
    eeRechitZ_.push_back(rechitPoint.z());}

    neeRechit_++;
  }


  ///ES
  EcalRecHitCollection::const_iterator esrechit;
  for( esrechit = ESRecHits->begin(); esrechit != ESRecHits->end(); esrechit++ ){

    double Energy = esrechit->energy();
    ESDetId det = esrechit->id();
    int ieta = det.six();
    int iphi = det.siy();
    int zside = det.zside();
    int strip = det.strip();

    int plane = det.plane();

    if(dumpCrystalinfo_){
    esRechitE_.push_back(Energy);
    esRechitiEta_.push_back(ieta);
    esRechitiPhi_.push_back(iphi);
    esRechitPlane_.push_back(plane);
    esRechitStrip_.push_back(strip);
    esRechitZside_.push_back(zside);}

    const GlobalPoint & rechitPoint = geo->getPosition(det);
    if(dumpCrystalinfo_){
    esRechitEta_.push_back(rechitPoint.eta());
    esRechitPhi_.push_back(rechitPoint.phi());
    esRechitX_.push_back(rechitPoint.x());
    esRechitY_.push_back(rechitPoint.y());
    esRechitZ_.push_back(rechitPoint.z());}

    nesRechit_++;
  }
}

