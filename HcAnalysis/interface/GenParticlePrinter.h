/*
Simple analyzer to print some gen particle info.
For testing and debugging, not meant to be used in production.
*/

#ifndef GEN_PARTICLE_PRINTER_H
#define GEN_PARTICLE_PRINTER_H

// include other parts of the framework
#include "HcAnalysis/HcAnalysis/interface/HcAnalysis.h"

// system include files
#include <memory>
#include <unordered_map>

// root classes
#include <Math/Vector4D.h>

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

// dataformats include files
#include "DataFormats/Math/interface/LorentzVector.h"

// include other analyzers and tools
#include "HcAnalysis/HcAnalysis/interface/GenTools.h"


class HcAnalysis;

class GenParticlePrinter {
  friend class HcAnalysis;
  private:

    HcAnalysis* hcAnalyzer;

  public:
    GenParticlePrinter(const edm::ParameterSet& iConfig, HcAnalysis* vars);
    ~GenParticlePrinter();
    // template member functions
    void beginJob(TTree*);
    void analyze(const edm::Event&);
};

#endif
