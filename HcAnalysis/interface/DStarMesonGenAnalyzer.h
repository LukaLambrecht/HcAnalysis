/*
Custom analyzer class for investigating gen-level D* meson decays.
*/

#ifndef DSTARGEN_ANALYZER_H
#define DSTARGEN_ANALYZER_H

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

class DStarMesonGenAnalyzer {
  friend class HcAnalysis;
  private:

    HcAnalysis* hcAnalyzer;

    int _nGenDStarMeson = 0;
    int _nGenDStarMesonToKPiPi = 0;
    int _genDStarMesonDecayType = 0;

  public:
    DStarMesonGenAnalyzer(const edm::ParameterSet& iConfig, HcAnalysis* vars);
    ~DStarMesonGenAnalyzer();
    // template member functions
    void beginJob(TTree*);
    void analyze(const edm::Event&);

    static int find_DStar_decay_type(const std::vector<reco::GenParticle>&);
    static std::map< std::string, const reco::GenParticle* > find_DStar_to_DZeroPi_to_KPiPi(
      const std::vector<reco::GenParticle>&);
};

#endif
