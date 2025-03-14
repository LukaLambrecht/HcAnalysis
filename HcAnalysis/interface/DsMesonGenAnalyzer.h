/*
Custom analyzer class for investigating gen-level Ds meson decays.
*/

#ifndef DSGEN_ANALYZER_H
#define DSGEN_ANALYZER_H

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


class HcAnalysis;

class DsMesonGenAnalyzer {
  friend class HcAnalysis;
  private:

    HcAnalysis* hcAnalyzer;

    int _nGenDsMeson = 0;
    int _nGenDsMesonToKKPi = 0;
    int _genDsMesonDecayType = 0;

  public:
    DsMesonGenAnalyzer(const edm::ParameterSet& iConfig, HcAnalysis* vars);
    ~DsMesonGenAnalyzer();
    // template member functions
    void beginJob(TTree*);
    void analyze(const edm::Event&);

    static int find_Ds_decay_type(const std::vector<reco::GenParticle>&);
    static std::map< std::string, const reco::GenParticle* > find_Ds_to_PhiPi_to_KKPi(
      const std::vector<reco::GenParticle>&);
};

#endif
