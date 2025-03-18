/*
Custom analyzer class for investigating c-quark fragmentation.
*/

#ifndef CFRAGMENTATION_ANALYZER_H
#define CFRAGMENTATION_ANALYZER_H

// include other parts of the framework
#include "HcAnalysis/HcAnalysis/interface/HcAnalysis.h"

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


class HcAnalysis;

class cFragmentationAnalyzer {
  friend class HcAnalysis;
  private:

    HcAnalysis* hcAnalyzer;

    int _cFragmentationPdgId = 0;
    int _cBarFragmentationPdgId = 0;

  public:
    cFragmentationAnalyzer(const edm::ParameterSet& iConfig, HcAnalysis* vars);
    ~cFragmentationAnalyzer();
    // template member functions
    void beginJob(TTree*);
    void analyze(const edm::Event&);
};

#endif
