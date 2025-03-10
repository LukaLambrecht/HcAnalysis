/*
Custom analyzer class for finding H bosons.
*/

#ifndef HIGGS_ANALYZER_H
#define HIGGS_ANALYZER_H

// include other parts of the framework
#include "HcAnalysis/HcAnalysis/interface/HcAnalysis.h"

// system include files
#include <memory>

// main include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

// general include files
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/LuminosityBlock.h"
#include "TTree.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Math/interface/deltaR.h"

class HcAnalysis;

class HiggsAnalyzer {
  friend class HcAnalysis;
  private:

    HcAnalysis* hcAnalyzer;

    // constants
    static constexpr double hmass = 125.;

    // ROOT tree variable declarations
    static const unsigned nHiggs_max = 5;
    unsigned _nHiggs = 0;
    double _HiggsInvMass[nHiggs_max];

  public:
    HiggsAnalyzer(const edm::ParameterSet& iConfig, HcAnalysis* vars);
    ~HiggsAnalyzer();
    
    // template member functions
    void beginJob(TTree*);
    void analyze(const edm::Event&);
    
};

#endif
