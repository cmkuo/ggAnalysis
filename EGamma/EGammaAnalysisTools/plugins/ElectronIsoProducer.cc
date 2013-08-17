// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "EGamma/EGammaAnalysisTools/src/PFIsolationEstimator.cc"
//
// class declaration
//

using namespace std;
using namespace reco;
class ElectronIsoProducer : public edm::EDFilter {
      public:
         explicit ElectronIsoProducer(const edm::ParameterSet&);
         ~ElectronIsoProducer();
      private:
        virtual bool filter(edm::Event&, const edm::EventSetup&);
  
// ----------member data ---------------------------
        bool verbose_;
        edm::InputTag vertexTag_;
        edm::InputTag electronTag_;
        edm::InputTag particleFlowTag_;
        std::string nameIsoCh_;
        std::string nameIsoPh_;
        std::string nameIsoNh_;

        PFIsolationEstimator isolator;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
ElectronIsoProducer::ElectronIsoProducer(const edm::ParameterSet& iConfig) {
        verbose_ = iConfig.getUntrackedParameter<bool>("verbose", false);
        vertexTag_ = iConfig.getParameter<edm::InputTag>("vertexTag");
        electronTag_ = iConfig.getParameter<edm::InputTag>("electronTag");
        particleFlowTag_ = iConfig.getParameter<edm::InputTag>("particleFlowTag");
  
        nameIsoCh_ = iConfig.getParameter<std::string>("nameValueMapIsoCh");
        nameIsoPh_ = iConfig.getParameter<std::string>("nameValueMapIsoPh");
        nameIsoNh_ = iConfig.getParameter<std::string>("nameValueMapIsoNh");
  
  
        produces<edm::ValueMap<double> >(nameIsoCh_);
        produces<edm::ValueMap<double> >(nameIsoPh_);
        produces<edm::ValueMap<double> >(nameIsoNh_);
  
        isolator.initializeElectronIsolation(kTRUE); //NOTE: this automatically set all the correct defaul veto values 
	//IMPORTANT!!!!: 0.3 is the default value for the cut-based-id. But it is analysis dependent, eg HZZ is 0.4
        isolator.setConeSize(0.3); 
	// isolator.setApplyMissHitPhVeto(kTRUE); // NOTE: for the moment set to FALSE just to be in synch wit the isoDep: 26May become true by default in the PFIsolationEstimator.cc

	// below just examples on how to change parameters but don't uncomment unless you have a good motivation to change
	// the default values aumatically set by initializeElectronIsolation(kTRUE)
	//   isolator.setParticleType(-1);  //-1 electrons and 1 photons
	//   isolator.setApplyDzDxyVeto(kFALSE);
	//   isolator.setApplyPFPUVeto(kTRUE);
	//   isolator.setDeltaRVetoEndcapCharged(0.015);
	//   isolator.setDeltaRVetoEndcapPhotons(0.08);
  
}


ElectronIsoProducer::~ElectronIsoProducer()
{
  
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool ElectronIsoProducer::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {
	using namespace edm;

        std::auto_ptr<edm::ValueMap<double> > chIsoMap(new edm::ValueMap<double>() );
	edm::ValueMap<double>::Filler chFiller(*chIsoMap);

        std::auto_ptr<edm::ValueMap<double> > phIsoMap(new edm::ValueMap<double>() );
	edm::ValueMap<double>::Filler phFiller(*phIsoMap);

        std::auto_ptr<edm::ValueMap<double> > nhIsoMap(new edm::ValueMap<double>() );
	edm::ValueMap<double>::Filler nhFiller(*nhIsoMap);

	Handle<reco::VertexCollection>  vertexCollection;
	iEvent.getByLabel(vertexTag_, vertexCollection);

   
	Handle<reco::GsfElectronCollection> egCollection;
	iEvent.getByLabel(electronTag_,egCollection);
        const reco::GsfElectronCollection egCandidates = (*egCollection.product());


	// All PF Candidate for alternate isolation
	Handle<reco::PFCandidateCollection> pfCandidatesH;
	iEvent.getByLabel(particleFlowTag_, pfCandidatesH);
	const  PFCandidateCollection thePfColl = *(pfCandidatesH.product());


        std::vector<double> chIsoValues;
	std::vector<double> phIsoValues;
	std::vector<double> nhIsoValues;
        chIsoValues.reserve(egCollection->size());
        phIsoValues.reserve(egCollection->size());
        nhIsoValues.reserve(egCollection->size());
   
	unsigned int ivtx = 0;
	VertexRef myVtxRef(vertexCollection, ivtx);

        for ( reco::GsfElectronCollection::const_iterator egIter = egCandidates.begin(); egIter != egCandidates.end(); ++egIter) {
	  isolator.fGetIsolation(&*egIter,
				 &thePfColl,
				 myVtxRef,
				 vertexCollection);
	  
	  if(verbose_) {
	    std::cout << " run " << iEvent.id().run() << " lumi " << iEvent.id().luminosityBlock() << " event " << iEvent.id().event();
	    std::cout << " pt " <<  egIter->pt() << " eta " << egIter->eta() << " phi " << egIter->phi() 
		      << " charge " << egIter->charge()<< " : " << std::endl;;
	    
	    std::cout << " ChargedIso " << isolator.getIsolationCharged() << std::endl;
	    std::cout << " PhotonIso " << isolator.getIsolationPhoton() << std::endl;
	    std::cout << " NeutralHadron Iso " << isolator.getIsolationNeutral()  << std::endl;
	  }

	  chIsoValues.push_back(isolator.getIsolationCharged());
	  phIsoValues.push_back(isolator.getIsolationPhoton());
	  nhIsoValues.push_back(isolator.getIsolationNeutral());

	}

	  
	chFiller.insert(egCollection, chIsoValues.begin(), chIsoValues.end() );
	chFiller.fill();
	
	phFiller.insert(egCollection, phIsoValues.begin(), phIsoValues.end() );
	phFiller.fill();
	
	nhFiller.insert(egCollection, nhIsoValues.begin(), nhIsoValues.end() );
	nhFiller.fill();
	
	
	iEvent.put(chIsoMap,nameIsoCh_);
	iEvent.put(phIsoMap,nameIsoPh_);
	iEvent.put(nhIsoMap,nameIsoNh_);
	
	  
	return true;
}

//define this as a plug-in
DEFINE_FWK_MODULE(ElectronIsoProducer);
