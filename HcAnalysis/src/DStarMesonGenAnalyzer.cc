/*
Custom analyzer class for investigating gen-level D* meson decays.
*/

#include "HcAnalysis/HcAnalysis/interface/DStarMesonGenAnalyzer.h"

// constructor //
DStarMesonGenAnalyzer::DStarMesonGenAnalyzer(const edm::ParameterSet& iConfig, HcAnalysis* hcAnalyzer):
    hcAnalyzer(hcAnalyzer){
};

// destructor //
DStarMesonGenAnalyzer::~DStarMesonGenAnalyzer(){
}

// beginJob //
void DStarMesonGenAnalyzer::beginJob(TTree* outputTree){
    outputTree->Branch("nGenDStarMeson", &_nGenDStarMeson, "nGenDStarMeson/i");
    outputTree->Branch("nGenDStarMesonToKPiPi", &_nGenDStarMesonToKPiPi, "nGenDStarMesonToKPiPi/i");
    outputTree->Branch("genDStarMesonDecayType", &_genDStarMesonDecayType, "genDStarMesonDecayType/i");
}

// analyze (main method) //
void DStarMesonGenAnalyzer::analyze(const edm::Event& iEvent){

    // get gen particles
    edm::Handle<std::vector<reco::GenParticle>> genParticles;
    iEvent.getByToken(hcAnalyzer->prunedGenParticlesToken, genParticles);
    if(!genParticles.isValid()){
        std::cout << "WARNING: genParticle collection not valid" << std::endl;
        return;
    }

    // find decay type
    _genDStarMesonDecayType = find_DStar_decay_type( *genParticles );

    // find D* -> pi D0 -> pi K pi
    std::map< std::string, const reco::GenParticle* > DStarGenParticles;
    DStarGenParticles = find_DStar_to_DZeroPi_to_KPiPi( *genParticles );
    
    // set variables to write to the tree
    if( DStarGenParticles["DStar"]==nullptr ) _nGenDStarMeson = 0;
    else _nGenDStarMeson = 1;
    if( DStarGenParticles["K"]==nullptr ) _nGenDStarMesonToKPiPi = 0;
    else _nGenDStarMesonToKPiPi = 1;
}

int DStarMesonGenAnalyzer::find_DStar_decay_type(
        const std::vector<reco::GenParticle>& genParticles){
    // find what type of event this is concerning the production and decay of D* mesons.
    // the numbering convention is as follows:
    // 0: undefined, none of the below
    // 1: D* -> D0 pi -> K pi pi (i.e. the decay of interest)
    // 2: D* -> D0 pi -> other pi
    // 3: D* -> other
    // 4: multiple D* mesons, to be decided how to handle if significant
    // 10: no D* meson, but another charmed meson
    
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

    // find a D* meson or other charmed mesons
    int nDStarMesons = 0;
    reco::GenParticle dstar;
    bool hasOtherCharmedMeson = false;
    for( const reco::GenParticle* p : hardScatterParticles ){
        int pdgid = p->pdgId();
        if(std::abs(pdgid) == 413){
            nDStarMesons ++;
            dstar = *p;
        }
        else if(std::abs(pdgid) > 400 || std::abs(pdgid) < 500){
            hasOtherCharmedMeson = true;
        }
    }

    // handle case of multiple Ds mesons
    if( nDStarMesons > 1 ) return 4;
    
    // handle case of no Ds mesons
    if( nDStarMesons == 0 ){
        if( hasOtherCharmedMeson ) return 10;
        else return 0;
    }

    // find the decay products of the DStar meson
    std::vector<const reco::GenParticle*> dstarDaughters;
    for(unsigned int i=0; i < dstar.numberOfDaughters(); ++i){
        dstarDaughters.push_back( &genParticles[dstar.daughterRef(i).key()] );
    }

    // find if they are a D0 meson and a pion
    if( dstarDaughters.size()!=2 ) return 3;
    const reco::GenParticle* dzero;
    if( std::abs(dstarDaughters.at(0)->pdgId())==421
        && std::abs(dstarDaughters.at(1)->pdgId())==211 ){
        dzero = dstarDaughters.at(0);
    } else if( std::abs(dstarDaughters.at(0)->pdgId())==211
        && std::abs(dstarDaughters.at(1)->pdgId())==421 ){
        dzero = dstarDaughters.at(1);
    } else return 3;

    // find the daughters of the D0
    std::vector<const reco::GenParticle*> dzeroDaughters;
    for(unsigned int i=0; i < dzero->numberOfDaughters(); ++i){
        dzeroDaughters.push_back( &genParticles[dzero->daughterRef(i).key()] );
    }

    // find if they are a kaon and a pion
    if( dzeroDaughters.size()!=2 ) return 2;
    if( std::abs(dzeroDaughters.at(0)->pdgId())==321
        && std::abs(dzeroDaughters.at(1)->pdgId())==211 ){
        // pass
    } else if( std::abs(dzeroDaughters.at(0)->pdgId())==211
        && std::abs(dzeroDaughters.at(1)->pdgId())==321 ){
        // pass
    } else return 2;

    // if all checks above succeeded,
    // we have a genuine D* -> D0 pi -> K pi pi event
    return 1;
}


std::map< std::string, const reco::GenParticle* > DStarMesonGenAnalyzer::find_DStar_to_DZeroPi_to_KPiPi(
        const std::vector<reco::GenParticle>& genParticles){
    // find D* -> D0 pi -> K pi pi at GEN level

    // initialize output
    std::map< std::string, const reco::GenParticle* > res = {
        {"DStar", nullptr},
        {"DZero", nullptr},
        {"Pi1", nullptr},
        {"K", nullptr},
        {"Pi2", nullptr},
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
    int nDStarMesons = 0;
    const reco::GenParticle* dstar;
    for( const reco::GenParticle* p : hardScatterParticles ){
        int pdgid = p->pdgId();
        if(std::abs(pdgid) == 413){
            nDStarMesons ++;
            dstar = p;
        }
    }

    // handle the case of no or multiple D* mesons found
    if( nDStarMesons==0 ){
        //std::cout << "WARNING: no D* found" << std::endl;
        return res;
    }
    else if( nDStarMesons > 1 ){
        //std::cout << "WARNING: multiple D* found" << std::endl;
        return res;
    }

    // store the D* meson in the output map
    res["DStar"] = dstar;

    // find its daughters
    std::vector<const reco::GenParticle*> dstarDaughters;
    for(unsigned int i=0; i<dstar->numberOfDaughters(); ++i){
        dstarDaughters.push_back( &genParticles[dstar->daughterRef(i).key()] );
    }

    // printouts
    //std::cout << "D* daughters:" << std::endl;
    //for( const reco::GenParticle* p: dstarDaughters ){ std::cout << p->pdgId() << " "; }
    //std::cout << std::endl;

    // find the D0 meson and the pion
    if( dstarDaughters.size()!=2 ) return res;
    const reco::GenParticle* dzero;
    const reco::GenParticle* pi;
    if( std::abs(dstarDaughters.at(0)->pdgId())==421
        && std::abs(dstarDaughters.at(1)->pdgId())==211 ){
        dzero = dstarDaughters.at(0);
        pi = dstarDaughters.at(1);
    } else if( std::abs(dstarDaughters.at(0)->pdgId())==211
        && std::abs(dstarDaughters.at(1)->pdgId())==421 ){
        dzero = dstarDaughters.at(1);
        pi = dstarDaughters.at(0);
    } else return res;

    // store the particles in the output map
    res["DZero"] = dzero;
    res["Pi1"] = pi;

    // printouts
    //std::cout << "  -> found D* -> D0 + pi" << std::endl;

    // find the daughters of the D0
    std::vector<const reco::GenParticle*> dzeroDaughters;
    for(unsigned int i=0; i<dzero->numberOfDaughters(); ++i){
        dzeroDaughters.push_back( &genParticles[dzero->daughterRef(i).key()] );
    }

    // printouts
    //std::cout << "D0 daughters" << std::endl;
    //for( const reco::GenParticle* p: dzeroDaughters ){ std::cout << p->pdgId() << " "; } 
    //std::cout << std::endl;

    // find the kaon and the pion
    const reco::GenParticle* K;
    const reco::GenParticle* pi2;
    if( dzeroDaughters.size()!=2 ) return res;
    if( std::abs(dzeroDaughters.at(0)->pdgId())==321
        && std::abs(dzeroDaughters.at(1)->pdgId())==211 ){
        K = dzeroDaughters.at(0);
        pi2 = dzeroDaughters.at(1);
    } else if( std::abs(dzeroDaughters.at(0)->pdgId())==211
        && std::abs(dzeroDaughters.at(1)->pdgId())==321 ){
        K = dzeroDaughters.at(1);
        pi2 = dzeroDaughters.at(0);
    } else return res;

    // print kinematics
    /*std::cout << "D* kinematics:" << std::endl;
    std::cout << dstar.pt() << " " << dstar.eta() << " " << dstar.phi() << std::endl;
    std::cout << "pion kinematics:" << std::endl;
    std::cout << pi->pt() << " " << pi->eta() << " " << pi->phi() << std::endl;
    std::cout << "D0 kinematics:" << std::endl;
    std::cout << dzero->pt() << " " << dzero->eta() << " " << dzero->phi() << std::endl;
    std::cout << "kaon kinematics:" << std::endl;
    std::cout << K->pt() << " " << K->eta() << " " << K->phi() << std::endl;
    std::cout << "pion2 kinematics:" << std::endl;
    std::cout << pi2->pt() << " " << pi2->eta() << " " << pi2->phi() << std::endl;*/

    // set the particles in the output map
    res["K"] = K;
    res["Pi2"] = pi2;
    return res;
}
