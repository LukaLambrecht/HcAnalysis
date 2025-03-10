// -*- C++ -*-
//
// Package:    HcAnalysis/HcAnalysis
// Class:      HcAnalysis
//
/**\class HcAnalysis HcAnalysis.cc HcAnalysis/HcAnalysis/plugins/HcAnalysis.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Luka Lambrecht
//         Created:  Mon, 10 Mar 2025 13:20:43 GMT
//
//


// include header
#include "HcAnalysis/HcAnalysis/interface/HcAnalysis.h"


//
// constructors and destructor
//
HcAnalysis::HcAnalysis(const edm::ParameterSet& iConfig):
  // read tokens
  // note: the tokens must be initialized in the header file,
  //       and their value must be set in the config file.
  packedPFCandidatesToken(consumes<pat::PackedCandidateCollection>(
    iConfig.getUntrackedParameter<edm::InputTag>("packedPFCandidates"))),
  lostTracksToken(consumes<pat::PackedCandidateCollection>(
    iConfig.getUntrackedParameter<edm::InputTag>("lostTracks")))
{
  // initialize specific analyzers
  dsMesonAnalyzer = new DsMesonAnalyzer(iConfig, this);
  higgsAnalyzer = new HiggsAnalyzer(iConfig, this);
}

HcAnalysis::~HcAnalysis() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

  // delete specific analyzers
  delete dsMesonAnalyzer;
  delete higgsAnalyzer;
}

//
// member functions
//

// ------------ method called for each event  ------------
void HcAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  // set event number
  _event = (unsigned long) iEvent.id().event();

  // run specific analyzers
  dsMesonAnalyzer->analyze(iEvent);
  higgsAnalyzer->analyze(iEvent);

  // get objects from tokens
  //edm::Handle<pat::PackedCandidateCollection> packedPFCandidates;
  //iEvent.getByToken(packedPFCandidatesToken, packedPFCandidates);
  //edm::Handle<pat::PackedCandidateCollection> lostTracks;
  //iEvent.getByToken(lostTracksToken, lostTracks);

  // fill output tree
  outputTree->Fill();
}

// ------------ method called once each job just before starting event loop  ------------
void HcAnalysis::beginJob() {
  
  // initialize output tree
  outputTree = outputFileService->make<TTree>("Events", "Events");
  outputTree->Branch("run", &_run, "run/l");
  outputTree->Branch("luminosityBlock", &_luminosityBlock, "luminosityBlock/l");
  outputTree->Branch("event", &_event, "_event/l");

  // do begin job for specific analyzers
  dsMesonAnalyzer->beginJob(outputTree);
  higgsAnalyzer->beginJob(outputTree);

  // initialize run, lumiblock and event number
  _run = 0;
  _luminosityBlock = 0;
  _event = 0;
}

// ------------ method called for each lumi block ---------
void HcAnalysis::beginLuminosityBlock(
  const edm::LuminosityBlock& iLumi, const edm::EventSetup& iSetup){
  _luminosityBlock = (unsigned long) iLumi.id().luminosityBlock();
}

//------------- method called for each run -------------
void HcAnalysis::beginRun(
  const edm::Run& iRun, edm::EventSetup const& iSetup){
  _run = (unsigned long) iRun.id().run();
}

// define this as a plug-in
DEFINE_FWK_MODULE(HcAnalysis);
