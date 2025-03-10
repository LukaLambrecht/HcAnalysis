/*
Custom analyzer class for finding Ds mesons from triplets of tracks.

The targeted decay chain is the following:
Ds -> phi pi+ -> K+ K- pi+
*/

#ifndef DSMESON_ANALYZER_H
#define DSMESON_ANALYZER_H

// include other parts of the framework
#include "HcAnalysis/HcAnalysis/interface/HcAnalysis.h"

// system include files
#include <memory>
#include <unordered_map>
#include <Math/Vector4D.h>
#include <Math/SVector.h> // root high-performance vector class
#include <Math/SMatrix.h> // root high-performance matrix class
#include "TMatrixDSym.h" // for fixTrackCovariance
#include "TVectorD.h" // for fixTrackCovariance

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

// vertex fitter include files
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/ParametrizedEngine/src/OAEParametrizedMagneticField.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"

// data format include files
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Candidate/interface/VertexCompositePtrCandidate.h"

class HcAnalysis;

class DsMesonAnalyzer {
  friend class HcAnalysis;
  private:

    HcAnalysis* hcAnalyzer;

    // constants
    static constexpr double pimass = 0.13957;
    static constexpr double pimass2 = pimass*pimass;
    static constexpr double kmass = 0.493677;
    static constexpr double kmass2 = kmass*kmass;
    static constexpr double phimass = 1.019461;
    static constexpr double dsmass = 1.96847;

    // ROOT tree variable declarations
    static const unsigned nDsMeson_max = 50;
    unsigned _nDsMeson = 0;
    double _DsMeson_mass[nDsMeson_max];
    double _DsMeson_pt[nDsMeson_max];
    double _DsMeson_eta[nDsMeson_max];
    double _DsMeson_phi[nDsMeson_max];
    double _DsMeson_PhiMeson_mass[nDsMeson_max];
    double _DsMeson_PhiMeson_pt[nDsMeson_max];
    double _DsMeson_PhiMeson_eta[nDsMeson_max];
    double _DsMeson_PhiMeson_phi[nDsMeson_max];

  public:
    DsMesonAnalyzer(const edm::ParameterSet& iConfig, HcAnalysis* vars);
    ~DsMesonAnalyzer();
    
    // template member functions
    void beginJob(TTree*);
    void analyze(const edm::Event&);
    
    // helper functions
};

#endif
