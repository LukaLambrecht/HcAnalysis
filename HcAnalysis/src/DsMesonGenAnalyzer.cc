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
    outputTree->Branch("genDsMesonDecayType", &_genDsMesonDecayType, "genDsMesonDecayType/i");
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

    // find decay type
    _genDsMesonDecayType = find_Ds_decay_type( *genParticles );

    // find Ds -> K K pi
    std::map< std::string, const reco::GenParticle* > DsGenParticles;
    DsGenParticles = find_Ds_to_PhiPi_to_KKPi( *genParticles );
    
    // set variables to write to the tree
    if( DsGenParticles["Ds"]==nullptr ) _nGenDsMeson = 0;
    else _nGenDsMeson = 1;
    if( DsGenParticles["KPlus"]==nullptr ) _nGenDsMesonToKKPi = 0;
    else _nGenDsMesonToKKPi = 1;
}

int DsMesonGenAnalyzer::find_Ds_decay_type(
        const std::vector<reco::GenParticle>& genParticles){
    // find what type of event this is concerning the production and decay of Ds mesons.
    // the numbering convention is as follows:
    // 0: undefined, none of the below
    // 1: Ds -> phi pi -> K K pi (i.e. the decay of interest)
    // 2: Ds -> phi pi -> other pi
    // 3: Ds -> other
    // 4: multiple Ds mesons, to be decided how to handle if significant
    // 10: no Ds meson, but another charmed meson
    
    // find all gen particles from the hard scattering
    // (implemented here as having a proton as their mother)
    std::vector<const reco::GenParticle*> hardScatterParticles;
    for( const reco::GenParticle& p : genParticles ){
        if( !p.isLastCopy() ) continue;
        int mompdgid = GenTools::getMotherPdgId(p, genParticles);
        if( std::abs(mompdgid)!=2212 ) continue;
        hardScatterParticles.push_back(&p);
    }
    if( hardScatterParticles.size() < 1 ) return 0;

    // find a Ds meson or other charmed mesons
    int nDsMesons = 0;
    reco::GenParticle ds;
    bool hasOtherCharmedMeson = false;
    for( const reco::GenParticle* p : hardScatterParticles ){
        int pdgid = p->pdgId();
        if(std::abs(pdgid) == 431){
            nDsMesons ++;
            ds = *p;
        }
        else if(std::abs(pdgid) > 400 || std::abs(pdgid) < 500){
            hasOtherCharmedMeson = true;
        }
    }

    // handle case of multiple Ds mesons
    if( nDsMesons > 1 ) return 4;
    
    // handle case of no Ds mesons
    if( nDsMesons == 0 ){
        if( hasOtherCharmedMeson ) return 10;
        else return 0;
    }

    // find the decay products of the Ds meson
    std::vector<const reco::GenParticle*> dsDaughters;
    for(unsigned int i=0; i < ds.numberOfDaughters(); ++i){
        dsDaughters.push_back( &genParticles[ds.daughterRef(i).key()] );
    }

    // find if they are a pion and a phi meson
    if( dsDaughters.size()!=2 ) return 3;
    const reco::GenParticle* phi;
    if( std::abs(dsDaughters.at(0)->pdgId())==333
        && std::abs(dsDaughters.at(1)->pdgId())==211 ){
        phi = dsDaughters.at(0);
    } else if( std::abs(dsDaughters.at(0)->pdgId())==211
        && std::abs(dsDaughters.at(1)->pdgId())==333 ){
        phi = dsDaughters.at(1);
    } else return 3;

    // find the daughters of the phi
    std::vector<const reco::GenParticle*> phiDaughters;
    for(unsigned int i=0; i < phi->numberOfDaughters(); ++i){
        phiDaughters.push_back( &genParticles[phi->daughterRef(i).key()] );
    }

    // find if they are kaons
    if( phiDaughters.size()!=2 ) return 2;
    if( std::abs(phiDaughters.at(0)->pdgId())==321
        && std::abs(phiDaughters.at(1)->pdgId())==321 ){
        // pass
    } else return 2;

    // if all checks above succeeded,
    // we have a genuine Ds -> phi pi -> K K pi event
    return 1;
}


std::map< std::string, const reco::GenParticle* > DsMesonGenAnalyzer::find_Ds_to_PhiPi_to_KKPi(
        const std::vector<reco::GenParticle>& genParticles){
    // find Ds -> phi pi -> K K pi at GEN level

    // initialize output
    std::map< std::string, const reco::GenParticle* > res = {
        {"Ds", nullptr},
        {"Phi", nullptr},
        {"Pi", nullptr},
        {"KPlus", nullptr},
        {"KMinus", nullptr},
    };

    // find all gen particles from the hard scattering
    // (implemented here as having a proton as their mother)
    std::vector<const reco::GenParticle*> hardScatterParticles;
    for( const reco::GenParticle& p : genParticles ){
        if( !p.isLastCopy() ) continue;
        int mompdgid = GenTools::getMotherPdgId(p, genParticles);
        if( std::abs(mompdgid)!=2212 ) continue;
        hardScatterParticles.push_back(&p);
    }
    if( hardScatterParticles.size() < 1 ) return res;

    // find the Ds meson
    int nDsMesons = 0;
    const reco::GenParticle* ds;
    for( const reco::GenParticle* p : hardScatterParticles ){
        int pdgid = p->pdgId();
        if(std::abs(pdgid) == 431){
            nDsMesons ++;
            ds = p;
        }
    }

    // handle the case of no or multiple Ds mesons found
    if( nDsMesons==0 ){
        //std::cout << "WARNING: no Ds found" << std::endl;
        return res;
    }
    else if( nDsMesons > 1 ){
        //std::cout << "WARNING: multiple Ds found" << std::endl;
        return res;
    }

    // store the Ds meson in the output map
    res["Ds"] = ds;

    // find its daughters
    std::vector<const reco::GenParticle*> dsDaughters;
    for(unsigned int i=0; i<ds->numberOfDaughters(); ++i){
        dsDaughters.push_back( &genParticles[ds->daughterRef(i).key()] );
    }

    // printouts
    //std::cout << "ds daughters:" << std::endl;
    //for( const reco::GenParticle* p: dsDaughters ){ std::cout << p->pdgId() << " "; }
    //std::cout << std::endl;

    // find the pion and phi
    if( dsDaughters.size()!=2 ) return res;
    const reco::GenParticle* pi;
    const reco::GenParticle* phi;
    if( std::abs(dsDaughters.at(0)->pdgId())==333
        && std::abs(dsDaughters.at(1)->pdgId())==211 ){
        phi = dsDaughters.at(0);
        pi = dsDaughters.at(1);
    } else if( std::abs(dsDaughters.at(0)->pdgId())==211
        && std::abs(dsDaughters.at(1)->pdgId())==333 ){
        phi = dsDaughters.at(1);
        pi = dsDaughters.at(0);
    } else return res;

    // store the particles in the output map
    res["Phi"] = phi;
    res["Pi"] = pi;

    // printouts
    //std::cout << "  -> found Ds -> phi + pi" << std::endl;

    // find the daughters of the phi
    std::vector<const reco::GenParticle*> phiDaughters;
    for(unsigned int i=0; i<phi->numberOfDaughters(); ++i){
        phiDaughters.push_back( &genParticles[phi->daughterRef(i).key()] );
    }

    // printouts
    //std::cout << "phi daughters" << std::endl;
    //for( const reco::GenParticle* p: phiDaughters ){ std::cout << p->pdgId() << " "; } 
    //std::cout << std::endl;

    // find the kaons
    const reco::GenParticle* K1;
    const reco::GenParticle* K2;
    const reco::GenParticle* KPlus;
    const reco::GenParticle* KMinus;
    if( phiDaughters.size()!=2 ) return res;
    if( std::abs(phiDaughters.at(0)->pdgId())==321
        && std::abs(phiDaughters.at(1)->pdgId())==321 ){
        K1 = phiDaughters.at(0);
        K2 = phiDaughters.at(1);
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
    /*std::cout << "Ds kinematics:" << std::endl;
    std::cout << ds.pt() << " " << ds.eta() << " " << ds.phi() << std::endl;
    std::cout << "pion kinematics:" << std::endl;
    std::cout << pi->pt() << " " << pi->eta() << " " << pi->phi() << std::endl;
    std::cout << "phi kinematics:" << std::endl;
    std::cout << phi->pt() << " " << phi->eta() << " " << phi->phi() << std::endl;
    std::cout << "kaon1 kinematics:" << std::endl;
    std::cout << K1->pt() << " " << K1->eta() << " " << K1->phi() << std::endl;
    std::cout << "kaon2 kinematics:" << std::endl;
    std::cout << K2->pt() << " " << K2->eta() << " " << K2->phi() << std::endl;*/

    // set the particles in the output map
    res["KPlus"] = KPlus;
    res["KMinus"] = KMinus;
    return res;
}
