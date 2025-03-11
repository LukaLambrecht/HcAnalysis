/*
Custom analyzer class for investigating gen-level Ds meson decays.
*/

#include "HcAnalysis/HcAnalysis/interface/DsMesonGenAnalyzer.h"

// constructor //
DsMesonGenAnalyzer::DsMesonGenAnalyzer(const edm::ParameterSet& iConfig, HcAnalysis* hcAnalyzer):
    hcAnalyzer(hcAnalyzer){
};

// destructor //
DsMesonGenAnalyzer::~DsMesonGenAnalyzer(){
}

// beginJob //
void DsMesonGenAnalyzer::beginJob(TTree* outputTree){
    outputTree->Branch("nGenDsMeson", &_nGenDsMeson, "nGenDsMeson/i");
    outputTree->Branch("nGenDsMesonToKKPi", &_nGenDsMesonToKKPi, "nGenDsMesonToKKPi/i");
}

// analyze (main method) //
void DsMesonGenAnalyzer::analyze(const edm::Event& iEvent){

    // get gen particles
    edm::Handle<std::vector<reco::GenParticle>> genParticles;
    iEvent.getByToken(hcAnalyzer->prunedGenParticlesToken, genParticles);
    if(!genParticles.isValid()){
        std::cout << "WARNING: genParticle collection not valid" << std::endl;
        return;
    }

    // find Ds -> K K pi
    std::map< std::string, const reco::GenParticle* > DsGenParticles;
    DsGenParticles = find_Ds_To_PhiPi_To_KKPi( *genParticles );
    
    // set variables to write to the tree
    if( DsGenParticles["Ds"]==nullptr ) _nGenDsMeson = 0;
    else _nGenDsMeson = 1;
    if( DsGenParticles["KPlus"]==nullptr ) _nGenDsMesonToKKPi = 0;
    else _nGenDsMesonToKKPi = 1;
}

std::map< std::string, const reco::GenParticle* > DsMesonGenAnalyzer::find_Ds_To_PhiPi_To_KKPi(
        const std::vector<reco::GenParticle>& genParticles){
    // find Ds -> phi pi -> K K pi at GEN level
    // note: this function assumes maximum one such decay per event,
    //       to be extended to include multiple.

    // initialize output
    std::map< std::string, const reco::GenParticle* > res = {
        {"Ds", nullptr},
        {"Phi", nullptr},
        {"Pi", nullptr},
        {"KPlus", nullptr},
        {"KMinus", nullptr},
    };

    // find the Ds
    int dsindex = -1;
    for(unsigned i=0; i<genParticles.size(); ++i){
        int absId = std::abs( genParticles.at(i).pdgId() );
        if(absId == 431) dsindex = i;
    }
    if( dsindex<0 ){
        std::cout << "WARNING: no Ds found" << std::endl;
        return res;
    }
    const reco::GenParticle ds = genParticles.at(dsindex);
    res["Ds"] = &ds;

    // find its daughters
    std::vector<const reco::GenParticle*> dsdaughters;
    for(unsigned int i=0; i<ds.numberOfDaughters(); ++i){
        dsdaughters.push_back( &genParticles[ds.daughterRef(i).key()] );
    }

    // printouts
    std::cout << "ds daughters:" << std::endl;
    for( const reco::GenParticle* p: dsdaughters ){ std::cout << p->pdgId() << " "; }
    std::cout << std::endl;

    // find the pion and phi
    if( dsdaughters.size()!=2 ) return res;
    const reco::GenParticle* pi;
    const reco::GenParticle* phi;
    if( std::abs(dsdaughters.at(0)->pdgId())==333
        && std::abs(dsdaughters.at(1)->pdgId())==211 ){
        phi = dsdaughters.at(0);
        pi = dsdaughters.at(1);
    } else if( std::abs(dsdaughters.at(0)->pdgId())==211
        && std::abs(dsdaughters.at(1)->pdgId())==333 ){
        phi = dsdaughters.at(1);
        pi = dsdaughters.at(0);
    } else return res;
    res["Phi"] = phi;
    res["Pi"] = pi;

    // printouts
    std::cout << "  -> found Ds -> phi + pi" << std::endl;

    // find the daughters of the phi
    std::vector<const reco::GenParticle*> phidaughters;
    for(unsigned int i=0; i<phi->numberOfDaughters(); ++i){
        phidaughters.push_back( &genParticles[phi->daughterRef(i).key()] );
    }

    // printouts
    std::cout << "phi daughters" << std::endl;
    for( const reco::GenParticle* p: phidaughters ){ std::cout << p->pdgId() << " "; } 
    std::cout << std::endl;

    // find the kaons
    const reco::GenParticle* K1;
    const reco::GenParticle* K2;
    const reco::GenParticle* KPlus;
    const reco::GenParticle* KMinus;
    if( phidaughters.size()!=2 ) return res;
    if( std::abs(phidaughters.at(0)->pdgId())==321
        && std::abs(phidaughters.at(1)->pdgId())==321 ){
        K1 = phidaughters.at(0);
        K2 = phidaughters.at(1);
    } else return res;

    // find which one is the positive and which one the negative
    if( K1->charge() > 0 && K2->charge() < 0 ){
        KPlus = K1; 
        KMinus = K2;
    }
    else{
        KPlus = K2;
        KMinus = K1;
    }

    // print kinematics
    std::cout << "Ds kinematics:" << std::endl;
    std::cout << ds.pt() << " " << ds.eta() << " " << ds.phi() << std::endl;
    std::cout << "pion kinematics:" << std::endl;
    std::cout << pi->pt() << " " << pi->eta() << " " << pi->phi() << std::endl;
    std::cout << "phi kinematics:" << std::endl;
    std::cout << phi->pt() << " " << phi->eta() << " " << phi->phi() << std::endl;
    std::cout << "kaon1 kinematics:" << std::endl;
    std::cout << K1->pt() << " " << K1->eta() << " " << K1->phi() << std::endl;
    std::cout << "kaon2 kinematics:" << std::endl;
    std::cout << K2->pt() << " " << K2->eta() << " " << K2->phi() << std::endl;

    // set the particles
    res["Ds"] = &ds;
    res["Phi"] = phi;
    res["Pi"] = pi;
    res["KPlus"] = KPlus;
    res["KMinus"] = KMinus;
    return res;
}
