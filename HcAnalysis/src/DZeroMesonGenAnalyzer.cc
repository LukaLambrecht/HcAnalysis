/*
Custom analyzer class for investigating gen-level D0 meson decays.
*/

#include "HcAnalysis/HcAnalysis/interface/DZeroMesonGenAnalyzer.h"

// constructor //
DZeroMesonGenAnalyzer::DZeroMesonGenAnalyzer(const edm::ParameterSet& iConfig, HcAnalysis* hcAnalyzer):
    hcAnalyzer(hcAnalyzer){
};

// destructor //
DZeroMesonGenAnalyzer::~DZeroMesonGenAnalyzer(){
}

// beginJob //
void DZeroMesonGenAnalyzer::beginJob(TTree* outputTree){
    outputTree->Branch("nGenDZeroMeson", &_nGenDZeroMeson, "nGenDZeroMeson/i");
    outputTree->Branch("nGenDZeroMesonToKPi", &_nGenDZeroMesonToKPi, "nGenDZeroMesonToKPi/i");
    outputTree->Branch("genDZeroMesonDecayType", &_genDZeroMesonDecayType, "genDZeroMesonDecayType/i");
}

// analyze (main method) //
void DZeroMesonGenAnalyzer::analyze(const edm::Event& iEvent){

    // get gen particles
    edm::Handle<std::vector<reco::GenParticle>> genParticles;
    iEvent.getByToken(hcAnalyzer->prunedGenParticlesToken, genParticles);
    if(!genParticles.isValid()){
        std::cout << "WARNING: genParticle collection not valid" << std::endl;
        return;
    }

    // find decay type
    _genDZeroMesonDecayType = find_DZero_decay_type( *genParticles );

    // find Ds -> K K pi
    std::map< std::string, const reco::GenParticle* > DZeroGenParticles;
    DZeroGenParticles = find_DZero_to_KPi( *genParticles );
    
    // set variables to write to the tree
    if( DZeroGenParticles["DZero"]==nullptr ) _nGenDZeroMeson = 0;
    else _nGenDZeroMeson = 1;
    if( DZeroGenParticles["K"]==nullptr ) _nGenDZeroMesonToKPi = 0;
    else _nGenDZeroMesonToKPi = 1;
}

int DZeroMesonGenAnalyzer::find_DZero_decay_type(
        const std::vector<reco::GenParticle>& genParticles){
    // find what type of event this is concerning the production and decay of D0 mesons.
    // the numbering convention is as follows:
    // 0: undefined, none of the below
    // 1: D0 -> K pi (i.e. the decay of interest)
    // 2: D0 -> other
    // 3: multiple D0 mesons, to be decided how to handle if significant
    // 10: no D0 meson, but another charmed meson
    
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

    // find a D0 meson or other charmed mesons
    int nDZeroMesons = 0;
    reco::GenParticle dzero;
    bool hasOtherCharmedMeson = false;
    for( const reco::GenParticle* p : hardScatterParticles ){
        int pdgid = p->pdgId();
        if(std::abs(pdgid) == 421){
            nDZeroMesons ++;
            dzero = *p;
        }
        else if(std::abs(pdgid) > 400 || std::abs(pdgid) < 500){
            hasOtherCharmedMeson = true;
        }
    }

    // handle case of multiple D0 mesons
    if( nDZeroMesons > 1 ) return 3;
    
    // handle case of no Ds mesons
    if( nDZeroMesons == 0 ){
        if( hasOtherCharmedMeson ) return 10;
        else return 0;
    }

    // find the decay products of the D0 meson
    std::vector<const reco::GenParticle*> dzeroDaughters;
    for(unsigned int i=0; i < dzero.numberOfDaughters(); ++i){
        dzeroDaughters.push_back( &genParticles[dzero.daughterRef(i).key()] );
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
    // we have a genuine D0 -> K pi event
    return 1;
}


std::map< std::string, const reco::GenParticle* > DZeroMesonGenAnalyzer::find_DZero_to_KPi(
        const std::vector<reco::GenParticle>& genParticles){
    // find D0 -> K pi at GEN level

    // initialize output
    std::map< std::string, const reco::GenParticle* > res = {
        {"DZero", nullptr},
        {"K", nullptr},
        {"Pi", nullptr}
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

    // find the D0 meson
    int nDZeroMesons = 0;
    const reco::GenParticle* dzero;
    for( const reco::GenParticle* p : hardScatterParticles ){
        int pdgid = p->pdgId();
        if(std::abs(pdgid) == 421){
            nDZeroMesons ++;
            dzero = p;
        }
    }

    // handle the case of no or multiple D0 mesons found
    if( nDZeroMesons==0 ){
        //std::cout << "WARNING: no D0 found" << std::endl;
        return res;
    }
    else if( nDZeroMesons > 1 ){
        //std::cout << "WARNING: multiple D0 found" << std::endl;
        return res;
    }

    // store the Ds meson in the output map
    res["DZero"] = dzero;

    // find its daughters
    std::vector<const reco::GenParticle*> dzeroDaughters;
    for(unsigned int i=0; i<dzero->numberOfDaughters(); ++i){
        dzeroDaughters.push_back( &genParticles[dzero->daughterRef(i).key()] );
    }

    //std::cout << "dzero daughters:" << std::endl;
    //for( const reco::GenParticle* p: dzeroDaughters ){ std::cout << p->pdgId() << " "; }
    //std::cout << std::endl;

    // find the kaon and the pion
    if( dzeroDaughters.size()!=2 ) return res;
    const reco::GenParticle* K;
    const reco::GenParticle* pi;
    if( std::abs(dzeroDaughters.at(0)->pdgId())==321
        && std::abs(dzeroDaughters.at(1)->pdgId())==211 ){
        K = dzeroDaughters.at(0);
        pi = dzeroDaughters.at(1);
    } else if( std::abs(dzeroDaughters.at(0)->pdgId())==211
        && std::abs(dzeroDaughters.at(1)->pdgId())==321 ){
        K = dzeroDaughters.at(1);
        pi = dzeroDaughters.at(0);
    } else return res;

    // store the particles in the output map
    res["K"] = K;
    res["Pi"] = pi;

    // printouts
    //std::cout << "  -> found D0 -> K + pi" << std::endl;

    /*std::cout << "D0 kinematics:" << std::endl;
    std::cout << dzero.pt() << " " << dzero.eta() << " " << dzero.phi() << std::endl;
    std::cout << "kaon kinematics:" << std::endl;
    std::cout << K->pt() << " " << K->eta() << " " << K->phi() << std::endl;
    std::cout << "pion kinematics:" << std::endl;
    std::cout << pi->pt() << " " << pi->eta() << " " << pi->phi() << std::endl;*/

    return res;
}
