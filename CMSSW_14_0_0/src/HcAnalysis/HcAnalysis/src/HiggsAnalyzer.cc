/*
Custom analyzer class for finding H bosons.
*/

// include header
// note: all other includes are defined in the header
#include "HcAnalysis/HcAnalysis/interface/HiggsAnalyzer.h"

// constructor //
HiggsAnalyzer::HiggsAnalyzer(
    const edm::ParameterSet& iConfig, HcAnalysis* hcAnalyzer):
    hcAnalyzer(hcAnalyzer){
};

// destructor //
HiggsAnalyzer::~HiggsAnalyzer(){
}

// beginJob //
void HiggsAnalyzer::beginJob(TTree* outputTree){
    // initialize branches in the output tree

    outputTree->Branch("_nHiggs", &_nHiggs, "_nHiggs/i");
    outputTree->Branch("_HiggsInvMass", &_HiggsInvMass, "_HiggsInvMass[_nHiggs]/D");
}

// analyze (main method) //
void HiggsAnalyzer::analyze(const edm::Event& iEvent){
    // to do
}
