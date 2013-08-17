// -*- C++ -*-
//
// Package:    Filter
// Class:      ElectronHLTMatching
// 
/**\class Filter ElectronHLTMatching.cc HLTMatching/Filter/src/ElectronHLTMatching.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Matteo Sani,40 3-A02,+41227671577,
//         Created:  Wed Mar 28 18:07:47 CEST 2012
// $Id: ElectronHLTMatching.cc,v 1.1 2013/05/25 13:00:08 cmkuo Exp $
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "CommonTools/Utils/interface/StringCutObjectSelector.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "TPRegexp.h"

class ElectronHLTMatching : public edm::EDProducer {
public:
  explicit ElectronHLTMatching(const edm::ParameterSet&);
  ~ElectronHLTMatching();
  
private:

  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void produce(edm::Event&, const edm::EventSetup&);

  template <class T1, class T2>
  std::vector<int> matchByDeltaR(const std::vector<T1> &, 
				 const std::vector<T2> &,
				 const double);

  trigger::TriggerObjectCollection selectTriggerObjects(const trigger::TriggerObjectCollection &,
							const trigger::TriggerEvent &,
							const std::string);

  reco::GsfElectronCollection selectElectrons(const reco::GsfElectronCollection&, 
					      const StringCutObjectSelector<reco::GsfElectron> &);
  edm::InputTag inputCollection;
  edm::InputTag triggerResultsLabel;
  edm::InputTag triggerSummaryLabel;
  std::vector<std::string> triggerPaths;
  std::vector<std::string> modules;
  std::string recoCuts;
  std::string hltCuts;
  bool doMatching;
  float dR;

  std::vector<edm::InputTag> moduleLabels;
  HLTConfigProvider hltConfig;
  std::vector<std::string> realHltPaths;
};      

ElectronHLTMatching::ElectronHLTMatching(const edm::ParameterSet& iConfig) {
  
  inputCollection     = iConfig.getParameter<edm::InputTag>("InputCollection");
  triggerResultsLabel = iConfig.getParameter<edm::InputTag>("TriggerResults");
  triggerSummaryLabel = iConfig.getParameter<edm::InputTag>("HLTTriggerSummaryAOD");
  modules             = iConfig.getParameter<std::vector<std::string> >("ModuleLabels");
  triggerPaths        = iConfig.getParameter<std::vector<std::string> >("TriggerPaths");

  recoCuts = iConfig.getParameter<std::string>("RecoCuts");
  hltCuts  = iConfig.getParameter<std::string>("HLTCuts");

  dR = iConfig.getParameter<double>("DeltaR");

  doMatching = iConfig.getParameter<bool>("DoMatching");

  if (modules.size() != triggerPaths.size()) {
    edm::LogWarning("HLTMatchingFilter") << "You are specifying a different number of modules and trigger paths !!"; 
  }

  produces<reco::GsfElectronCollection>();
}

ElectronHLTMatching::~ElectronHLTMatching() 
{}

void ElectronHLTMatching::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup) 
{}


void ElectronHLTMatching::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  // Update when necessary the trigger table
  bool changedConfig = false;
  if (!hltConfig.init(iEvent.getRun(), iSetup, triggerResultsLabel.process(), changedConfig)) {
    edm::LogError("HLTMatchingFilter") << "Initialization of HLTConfigProvider failed!!"; 
    return;
  }

  if (changedConfig or realHltPaths.size() == 0) {
    realHltPaths.clear();
    moduleLabels.clear();

    for (size_t i = 0; i < triggerPaths.size(); i++) {
      TPRegexp pattern(triggerPaths[i]);
      for (size_t j = 0; j < hltConfig.triggerNames().size(); j++) {
	if (TString(hltConfig.triggerNames()[j]).Contains(pattern))
	  realHltPaths.push_back(hltConfig.triggerNames()[j]);
      }
    }

    for (unsigned int j=0; j<modules.size(); j++)
      moduleLabels.push_back(edm::InputTag(modules[j], "", triggerResultsLabel.process()));
    /*
    for (size_t j=0; j<realHltPaths.size(); j++) {
      std::vector<std::string> temp = hltConfig.saveTagsModules(realHltPaths[j]);
      moduleLabels.push_back(edm::InputTag(temp[temp.size()-1], "", triggerResultsLabel.process()));
    }
    */
  }
  std::auto_ptr<reco::GsfElectronCollection> filteredElectrons(new reco::GsfElectronCollection() );

  // Get the collections
  edm::Handle<reco::GsfElectronCollection> elH;
  iEvent.getByLabel(inputCollection, elH);
  const reco::GsfElectronCollection electrons = (*elH.product());

  edm::Handle<trigger::TriggerEvent> triggerSummary;
  iEvent.getByLabel(triggerSummaryLabel, triggerSummary);

  edm::Handle<edm::TriggerResults> triggerResults;
  
  if(!triggerSummary.isValid()) {
    edm::LogError("HLTMatching/ElectronMatching") << "Missing triggerSummary with label " << triggerSummaryLabel <<std::endl;
    return;
  }

  iEvent.getByLabel(triggerResultsLabel, triggerResults);
  if(!triggerResults.isValid()) {
    edm::LogError("HLTMatching/ElectronMatching") << "Missing triggerResults with label " << triggerResultsLabel <<std::endl;
    return;
  }

  // Check if the event passes the specified trigger
  bool passTrigger = false;
  for (size_t i = 0; i < realHltPaths.size(); i++) {
    for (size_t j = 0; j < hltConfig.size(); j++) {
      if (hltConfig.triggerName(j) == realHltPaths[i]) {
	if (triggerResults->accept(j)) {
	  passTrigger = true;
	  break;
	}
      }
    }
  }

  if (!passTrigger) {
    iEvent.put(filteredElectrons);
    return;
  }

  // Select reco and HLT objects and do the matching
  trigger::TriggerObjectCollection allTriggerObjects = triggerSummary->getObjects();
  trigger::TriggerObjectCollection hltElectrons = selectTriggerObjects(allTriggerObjects, *triggerSummary, hltCuts);

  reco::GsfElectronCollection selectedElectrons = selectElectrons(electrons, recoCuts);

  std::vector<int> matches = matchByDeltaR(selectedElectrons, hltElectrons, dR);
  
  for (size_t i=0; i<selectedElectrons.size(); i++) {
    if ((doMatching && matches[i] != -1) ||
	(!doMatching && matches[i] == -1)) {
      filteredElectrons->push_back(selectedElectrons[i]);
    }
  }

  iEvent.put(filteredElectrons);
}

template <class T1, class T2> 
std::vector<int> ElectronHLTMatching::matchByDeltaR(const std::vector<T1> & collection1, 
						    const std::vector<T2> & collection2,
						    const double maxDeltaR) {

  const size_t n1 = collection1.size();
  const size_t n2 = collection2.size();

  std::vector<int> result(n1, -1);
  std::vector<std::vector<double> > deltaRMatrix(n1, std::vector<double>(n2, 999.));

  for (size_t i = 0; i < n1; i++)
    for (size_t j = 0; j < n2; j++) 
      deltaRMatrix[i][j] = deltaR(collection1[i], collection2[j]); 

  while (1) {
    double minDeltaR = maxDeltaR;
    int min_i = -1;
    int min_j = -1;
    for (size_t i=0; i<deltaRMatrix.size(); i++) {
      for (size_t j=0; j<deltaRMatrix[i].size(); j++) {
	if (deltaRMatrix[i][j] < minDeltaR) {
	  minDeltaR = deltaRMatrix[i][j];
	  min_i = i;
	  min_j = j;
	}
      }
    }

    if (min_i == -1)
      break;

    result[min_i] = min_j;
    deltaRMatrix[min_i] = std::vector<double>(n2, 999.);
    for (size_t i=0; i<n1; i++)
      deltaRMatrix[i][min_j] = 999.;
  }

  return result;
}


trigger::TriggerObjectCollection ElectronHLTMatching::selectTriggerObjects(const trigger::TriggerObjectCollection & triggerObjects,
									   const trigger::TriggerEvent & triggerSummary,
									   const std::string hltcuts) {
  
  trigger::TriggerObjectCollection selectedObjects;
  StringCutObjectSelector<trigger::TriggerObject> selector(hltcuts);

  for (size_t t=0; t<moduleLabels.size(); t++) {

    size_t filterIndex = triggerSummary.filterIndex(moduleLabels[t]);
   
    if (filterIndex < triggerSummary.sizeFilters()) {
      const trigger::Keys &keys = triggerSummary.filterKeys(filterIndex);
   
      for (size_t j = 0; j < keys.size(); j++) {
	trigger::TriggerObject foundObject = triggerObjects[keys[j]];
	if (selector(foundObject)) {
	  selectedObjects.push_back(foundObject);
	}
      }
    }
  }

  return selectedObjects;
}

reco::GsfElectronCollection ElectronHLTMatching::selectElectrons(const reco::GsfElectronCollection& allElectrons, 
								 const StringCutObjectSelector<reco::GsfElectron> &selector) {

  reco::GsfElectronCollection selected(allElectrons);
  reco::GsfElectronCollection::iterator iter = selected.begin();
  while (iter != selected.end()) {
    if (selector(*iter))
      ++iter;
    else 
      selected.erase(iter);
  }

  return selected;
}

DEFINE_FWK_MODULE(ElectronHLTMatching);
