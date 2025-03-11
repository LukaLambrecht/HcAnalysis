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

    outputTree->Branch("nHiggs", &_nHiggs, "nHiggs/i");
    outputTree->Branch("HiggsInvMass", &_HiggsInvMass, "HiggsInvMass[nHiggs]/D");
}

// analyze (main method) //
void HiggsAnalyzer::analyze(const edm::Event& iEvent){
    // to do
}
