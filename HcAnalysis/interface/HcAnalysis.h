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

#ifndef HcAnalysis_h
#define HcAnalysis_h

// system include files
#include <memory>

// framework include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

// specific data format include files
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/PATObject.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

// TFile service
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// ROOT
#include "TTree.h"
#include "math.h"

// specificd analyzers
#include "HcAnalysis/HcAnalysis/interface/DsMesonAnalyzer.h"
#include "HcAnalysis/HcAnalysis/interface/HiggsAnalyzer.h"

//
// class declaration
//

// declare specific analyzers
class DsMesonAnalyzer;
class HiggsAnalyzer;

// define total analyzer
class HcAnalysis : public edm::one::EDAnalyzer<edm::one::WatchLuminosityBlocks, edm::one::WatchRuns, edm::one::SharedResources> {

  // define specific analyzers as friends
  friend DsMesonAnalyzer;
  friend HiggsAnalyzer;

  public:

    // constructor and destructor
    explicit HcAnalysis(const edm::ParameterSet&);
    ~HcAnalysis();
    
    // public variables and member functions
    // ...

  private:
    
    // template functions
    virtual void beginJob() override;
    virtual void beginRun(const edm::Run&, edm::EventSetup const&) override;
    virtual void beginLuminosityBlock(const edm::LuminosityBlock&, const edm::EventSetup&) override;
    virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
    virtual void endLuminosityBlock(const edm::LuminosityBlock&, const edm::EventSetup&) override {};
    virtual void endRun(const edm::Run&, edm::EventSetup const&) override {};
    virtual void endJob() override {};

    // tokens
    edm::EDGetTokenT<pat::PackedCandidateCollection> packedPFCandidatesToken;
    edm::EDGetTokenT<pat::PackedCandidateCollection> lostTracksToken;

    // specific analyzers
    DsMesonAnalyzer* dsMesonAnalyzer;
    HiggsAnalyzer* higgsAnalyzer;

    // output file and tree
    edm::Service<TFileService> outputFileService;
    TTree* outputTree;

    // private variables and member functions
    unsigned long _run;
    unsigned long _luminosityBlock;
    unsigned long _event;

};

#endif
