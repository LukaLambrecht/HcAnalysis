/*
Simple analyzer to print some gen particle info.
For testing and debugging, not meant to be used in production.
*/

#include "HcAnalysis/HcAnalysis/interface/GenParticlePrinter.h"

// constructor //
GenParticlePrinter::GenParticlePrinter(const edm::ParameterSet& iConfig, HcAnalysis* hcAnalyzer):
    hcAnalyzer(hcAnalyzer){
};

// destructor //
GenParticlePrinter::~GenParticlePrinter(){
}

// beginJob //
void GenParticlePrinter::beginJob(TTree* outputTree){
    // nothing to be done here so far
}

// analyze (main method) //
void GenParticlePrinter::analyze(const edm::Event& iEvent){

    // get gen particles
    edm::Handle<std::vector<reco::GenParticle>> genParticles;
    iEvent.getByToken(hcAnalyzer->prunedGenParticlesToken, genParticles);
    if(!genParticles.isValid()){
        std::cout << "WARNING: genParticle collection not valid" << std::endl;
        return;
    }

    // print a few relevant gen particles in the event
    for( const reco::GenParticle& p : *genParticles ){
        if(!p.isLastCopy()) continue;
        int pdgid = p.pdgId();
        if(std::abs(pdgid) < 400 || std::abs(pdgid) > 500) continue;
        int mompdgid = GenTools::getMotherPdgId(p, *genParticles);
        std::cout << "Particle " << pdgid << std::endl;
        std::cout << "  kinematics: " << p.pt() << " " << p.eta() << " " << p.phi() << std::endl;
        std::cout << "  mass: " << p.mass() << std::endl;
        std::cout << "  mom pdg id: " << mompdgid << std::endl;
    }
}
