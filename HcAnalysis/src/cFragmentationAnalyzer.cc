/*
Custom analyzer class for investigating c-quark fragmentation.
*/

#include "HcAnalysis/HcAnalysis/interface/cFragmentationAnalyzer.h"

// constructor //
cFragmentationAnalyzer::cFragmentationAnalyzer(const edm::ParameterSet& iConfig, HcAnalysis* hcAnalyzer):
    hcAnalyzer(hcAnalyzer){
};

// destructor //
cFragmentationAnalyzer::~cFragmentationAnalyzer(){
}

// beginJob //
void cFragmentationAnalyzer::beginJob(TTree* outputTree){
    outputTree->Branch("cFragmentationPdgId", &_cFragmentationPdgId, "cFragmentationPdgId/I");
    outputTree->Branch("cBarFragmentationPdgId", &_cBarFragmentationPdgId, "cFragmentationPdgId/I");
}

// analyze (main method) //
void cFragmentationAnalyzer::analyze(const edm::Event& iEvent){

    // initialization
    _cFragmentationPdgId = 0;
    _cBarFragmentationPdgId = 0;

    // get gen particles
    edm::Handle<std::vector<reco::GenParticle>> genParticles;
    iEvent.getByToken(hcAnalyzer->prunedGenParticlesToken, genParticles);
    if(!genParticles.isValid()){
        std::cout << "WARNING: genParticle collection not valid" << std::endl;
        return;
    }

    // find all gen particles from the hard scattering
    // (implemented here as having a proton as their mother)
    std::vector<const reco::GenParticle*> hardScatterParticles;
    for( const reco::GenParticle& p : *genParticles ){
        if( !p.isLastCopy() ) continue;
        int mompdgid = GenTools::getMotherPdgId(p, *genParticles);
        if( std::abs(mompdgid)!=2212 ) continue;
        hardScatterParticles.push_back(&p);
    }
    if( hardScatterParticles.size() < 1 ) return;

    // find charmed hadrons
    for( const reco::GenParticle* p : hardScatterParticles ){
        int pdgid = p->pdgId();
        if( (std::abs(pdgid) > 400 && std::abs(pdgid) < 500)
            || (std::abs(pdgid) > 4000 && std::abs(pdgid) < 5000) ){
            if(pdgid > 0) _cFragmentationPdgId = pdgid;
            else _cBarFragmentationPdgId = pdgid;
        }
    }

    // printouts
    if( _cFragmentationPdgId==0 || _cBarFragmentationPdgId==0 ){
        std::cout << "WARNING in cFragmentationAnalyzer:";
        std::cout << " no c-meson and/or cbar-meson found." << std::endl;
        std::cout << "Pdgids of hard scattering particles are:" << std::endl;
        for( const reco::GenParticle* p : hardScatterParticles ){
            int pdgid = p->pdgId();
            std::cout << pdgid << " ";
        }
        std::cout << std::endl;
    }
}
