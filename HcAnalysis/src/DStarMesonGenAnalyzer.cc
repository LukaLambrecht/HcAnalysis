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
    outputTree->Branch("genDStarMesonDecayType", &_genDStarMesonDecayType, "genDStarMesonDecayType/i");
    outputTree->Branch("nGenDStarMesonToKPiPi", &_nGenDStarMesonToKPiPi, "nGenDStarMesonToKPiPi/i");

    outputTree->Branch("genDStarMeson_DStar_pt", &_genDStarMeson_DStar_pt, "genDStarMeson_DStar_pt[nGenDStarMesonToKPiPi]/D");
    outputTree->Branch("genDStarMeson_DStar_eta", &_genDStarMeson_DStar_eta, "genDStarMeson_DStar_eta[nGenDStarMesonToKPiPi]/D");
    outputTree->Branch("genDStarMeson_DStar_phi", &_genDStarMeson_DStar_phi, "genDStarMeson_DStar_phi[nGenDStarMesonToKPiPi]/D");
    outputTree->Branch("genDStarMeson_DZero_pt", &_genDStarMeson_DZero_pt, "genDStarMeson_DZero_pt[nGenDStarMesonToKPiPi]/D");
    outputTree->Branch("genDStarMeson_DZero_eta", &_genDStarMeson_DZero_eta, "genDStarMeson_DZero_eta[nGenDStarMesonToKPiPi]/D");
    outputTree->Branch("genDStarMeson_DZero_phi", &_genDStarMeson_DZero_phi, "genDStarMeson_DZero_phi[nGenDStarMesonToKPiPi]/D");
    outputTree->Branch("genDStarMeson_Pi1_pt", &_genDStarMeson_Pi1_pt, "genDStarMeson_Pi1_pt[nGenDStarMesonToKPiPi]/D");
    outputTree->Branch("genDStarMeson_Pi1_eta", &_genDStarMeson_Pi1_eta, "genDStarMeson_Pi1_eta[nGenDStarMesonToKPiPi]/D");
    outputTree->Branch("genDStarMeson_Pi1_phi", &_genDStarMeson_Pi1_phi, "genDStarMeson_Pi1_phi[nGenDStarMesonToKPiPi]/D");
    outputTree->Branch("genDStarMeson_K_pt", &_genDStarMeson_K_pt, "genDStarMeson_K_pt[nGenDStarMesonToKPiPi]/D");
    outputTree->Branch("genDStarMeson_K_eta", &_genDStarMeson_K_eta, "genDStarMeson_K_eta[nGenDStarMesonToKPiPi]/D");
    outputTree->Branch("genDStarMeson_K_phi", &_genDStarMeson_K_phi, "genDStarMeson_K_phi[nGenDStarMesonToKPiPi]/D");
    outputTree->Branch("genDStarMeson_Pi2_pt", &_genDStarMeson_Pi2_pt, "genDStarMeson_Pi2_pt[nGenDStarMesonToKPiPi]/D");
    outputTree->Branch("genDStarMeson_Pi2_eta", &_genDStarMeson_Pi2_eta, "genDStarMeson_Pi2_eta[nGenDStarMesonToKPiPi]/D");
    outputTree->Branch("genDStarMeson_Pi2_phi", &_genDStarMeson_Pi2_phi, "genDStarMeson_Pi2_phi[nGenDStarMesonToKPiPi]/D");
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
    std::vector< std::map< std::string, const reco::GenParticle* > > DStarGenParticles;
    DStarGenParticles = find_DStar_to_DZeroPi_to_KPiPi( *genParticles );
    
    // set variables to write to the tree
    _nGenDStarMesonToKPiPi = DStarGenParticles.size();
    for(unsigned int idx=0; idx < DStarGenParticles.size(); idx++){
        _genDStarMeson_DStar_pt[idx] = DStarGenParticles[idx].at("DStar")->pt();
        _genDStarMeson_DStar_eta[idx] = DStarGenParticles[idx].at("DStar")->eta();
        _genDStarMeson_DStar_phi[idx] = DStarGenParticles[idx].at("DStar")->phi();
        _genDStarMeson_DZero_pt[idx] = DStarGenParticles[idx].at("DZero")->pt();
        _genDStarMeson_DZero_eta[idx] = DStarGenParticles[idx].at("DZero")->eta();
        _genDStarMeson_DZero_phi[idx] = DStarGenParticles[idx].at("DZero")->phi();
        _genDStarMeson_Pi1_pt[idx] = DStarGenParticles[idx].at("Pi1")->pt();
        _genDStarMeson_Pi1_eta[idx] = DStarGenParticles[idx].at("Pi1")->eta();
        _genDStarMeson_Pi1_phi[idx] = DStarGenParticles[idx].at("Pi1")->phi();
        _genDStarMeson_K_pt[idx] = DStarGenParticles[idx].at("K")->pt();
        _genDStarMeson_K_eta[idx] = DStarGenParticles[idx].at("K")->eta();
        _genDStarMeson_K_phi[idx] = DStarGenParticles[idx].at("K")->phi();
        _genDStarMeson_Pi2_pt[idx] = DStarGenParticles[idx].at("Pi2")->pt();
        _genDStarMeson_Pi2_eta[idx] = DStarGenParticles[idx].at("Pi2")->eta();
        _genDStarMeson_Pi2_phi[idx] = DStarGenParticles[idx].at("Pi2")->phi();
    }
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


std::vector< std::map< std::string, const reco::GenParticle* > > DStarMesonGenAnalyzer::find_DStar_to_DZeroPi_to_KPiPi(
        const std::vector<reco::GenParticle>& genParticles){
    // find D* -> D0 pi -> K pi pi at GEN level

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

    // loop over all hard scattering particles
    for( const reco::GenParticle* p : hardScatterParticles ){
    
        // check if it is a D* meson
        int pdgid = p->pdgId();
        if(std::abs(pdgid) != 413) continue;
        const reco::GenParticle* dstar = p;

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
        if( dstarDaughters.size()!=2 ) continue;
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
        } else continue;

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
        if( dzeroDaughters.size()!=2 ) continue;
        if( std::abs(dzeroDaughters.at(0)->pdgId())==321
            && std::abs(dzeroDaughters.at(1)->pdgId())==211 ){
            K = dzeroDaughters.at(0);
            pi2 = dzeroDaughters.at(1);
        } else if( std::abs(dzeroDaughters.at(0)->pdgId())==211
            && std::abs(dzeroDaughters.at(1)->pdgId())==321 ){
            K = dzeroDaughters.at(1);
            pi2 = dzeroDaughters.at(0);
        } else continue;

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
        std::map< std::string, const reco::GenParticle* > thisres = {
            {"DStar", dstar},
            {"DZero", dzero},
            {"Pi1", pi},
            {"K", K},
            {"Pi2", pi2}
        };
        res.push_back(thisres);
    }
    return res;
}
