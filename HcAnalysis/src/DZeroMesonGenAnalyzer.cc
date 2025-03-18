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
    std::vector< std::map< std::string, const reco::GenParticle* > > DZeroGenParticles;
    DZeroGenParticles = find_DZero_to_KPi( *genParticles );
    
    // set variables to write to the tree
    _nGenDZeroMesonToKPi = DZeroGenParticles.size();
}

int DZeroMesonGenAnalyzer::find_DZero_decay_type(
        const std::vector<reco::GenParticle>& genParticles){
    // find what type of event this is concerning the production and decay of D0 mesons.
    // the numbering convention is as follows:
    // 0: undefined, none of the below.
    // 1: at least one D0 -> K pi (i.e. the decay of interest).
    // 2: at least one D0, but excluding the above.
    // 3: at least one charmed hadron, but excluding the above.
    
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

    // initialize result
    int res = 99;

    // loop over hard scattering particles
    for( const reco::GenParticle* p : hardScatterParticles ){

        // check if it is a charmed hadron
        int pdgid = p->pdgId();
        bool cMeson = (std::abs(pdgid) > 400 && std::abs(pdgid) < 500);
        bool cBaryon = (std::abs(pdgid) > 4000 && std::abs(pdgid) < 5000);
        if( !(cMeson || cBaryon) ) continue;
        if(res > 3) res = 3;

        // check if it is a D0 meson
        if(std::abs(pdgid) != 421) continue;
        const reco::GenParticle* dzero = p;
        if(res > 2) res = 2;

        // find the decay products of the D0 meson
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
        // we have a genuine D0 -> K pi event
        if(res > 1) res = 1;
        break;
    }
    if(res > 3) res = 0;
    return res;
}


std::vector< std::map< std::string, const reco::GenParticle* > > DZeroMesonGenAnalyzer::find_DZero_to_KPi(
        const std::vector<reco::GenParticle>& genParticles){
    // find D0 -> K pi at GEN level

    // initialize output
    std::vector< std::map< std::string, const reco::GenParticle* > > res;

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

    // loop over hard scattering particles
    for( const reco::GenParticle* p : hardScatterParticles ){

        // find if it is a D0 meson
        int pdgid = p->pdgId();
        if(std::abs(pdgid) != 421) continue;
        const reco::GenParticle* dzero = p;

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
        std::map< std::string, const reco::GenParticle* > thisres = {
            {"DZero", dzero},
            {"K", K},
            {"Pi", pi}
        };
        res.push_back(thisres);
    
        // printouts
        //std::cout << "  -> found D0 -> K + pi" << std::endl;

        /*std::cout << "D0 kinematics:" << std::endl;
        std::cout << dzero.pt() << " " << dzero.eta() << " " << dzero.phi() << std::endl;
        std::cout << "kaon kinematics:" << std::endl;
        std::cout << K->pt() << " " << K->eta() << " " << K->phi() << std::endl;
        std::cout << "pion kinematics:" << std::endl;
        std::cout << pi->pt() << " " << pi->eta() << " " << pi->phi() << std::endl;*/
    }
    return res;
}
