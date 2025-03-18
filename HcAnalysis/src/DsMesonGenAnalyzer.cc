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
    outputTree->Branch("genDsMesonDecayType", &_genDsMesonDecayType, "genDsMesonDecayType/i");
    outputTree->Branch("nGenDsMesonToKKPi", &_nGenDsMesonToKKPi, "nGenDsMesonToKKPi/i");

    outputTree->Branch("genDsMeson_Ds_pt", &_genDsMeson_Ds_pt, "genDsMeson_Ds_pt[nGenDsMesonToKKPi]/D");
    outputTree->Branch("genDsMeson_Ds_eta", &_genDsMeson_Ds_eta, "genDsMeson_Ds_eta[nGenDsMesonToKKPi]/D");
    outputTree->Branch("genDsMeson_Ds_phi", &_genDsMeson_Ds_phi, "genDsMeson_Ds_phi[nGenDsMesonToKKPi]/D");
    outputTree->Branch("genDsMeson_Phi_pt", &_genDsMeson_Phi_pt, "genDsMeson_Phi_pt[nGenDsMesonToKKPi]/D");
    outputTree->Branch("genDsMeson_Phi_eta", &_genDsMeson_Phi_eta, "genDsMeson_Phi_eta[nGenDsMesonToKKPi]/D");
    outputTree->Branch("genDsMeson_Phi_phi", &_genDsMeson_Phi_phi, "genDsMeson_Phi_phi[nGenDsMesonToKKPi]/D");
    outputTree->Branch("genDsMeson_Pi_pt", &_genDsMeson_Pi_pt, "genDsMeson_Pi_pt[nGenDsMesonToKKPi]/D");
    outputTree->Branch("genDsMeson_Pi_eta", &_genDsMeson_Pi_eta, "genDsMeson_Pi_eta[nGenDsMesonToKKPi]/D");
    outputTree->Branch("genDsMeson_Pi_phi", &_genDsMeson_Pi_phi, "genDsMeson_Pi_phi[nGenDsMesonToKKPi]/D");
    outputTree->Branch("genDsMeson_KPlus_pt", &_genDsMeson_KPlus_pt, "genDsMeson_KPlus_pt[nGenDsMesonToKKPi]/D");
    outputTree->Branch("genDsMeson_KPlus_eta", &_genDsMeson_KPlus_eta, "genDsMeson_KPlus_eta[nGenDsMesonToKKPi]/D");
    outputTree->Branch("genDsMeson_KPlus_phi", &_genDsMeson_KPlus_phi, "genDsMeson_KPlus_phi[nGenDsMesonToKKPi]/D");
    outputTree->Branch("genDsMeson_KMinus_pt", &_genDsMeson_KMinus_pt, "genDsMeson_KMinus_pt[nGenDsMesonToKKPi]/D");
    outputTree->Branch("genDsMeson_KMinus_eta", &_genDsMeson_KMinus_eta, "genDsMeson_KMinus_eta[nGenDsMesonToKKPi]/D");
    outputTree->Branch("genDsMeson_KMinus_phi", &_genDsMeson_KMinus_phi, "genDsMeson_KMinus_phi[nGenDsMesonToKKPi]/D");
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

    // find Ds -> phi pi -> K K pi
    std::vector< std::map< std::string, const reco::GenParticle* > > DsGenParticles;
    DsGenParticles = find_Ds_to_PhiPi_to_KKPi( *genParticles );

    // set variables to write to the tree
    _nGenDsMesonToKKPi = DsGenParticles.size();
    for(unsigned int idx=0; idx < DsGenParticles.size(); idx++){
        _genDsMeson_Ds_pt[idx] = DsGenParticles[idx].at("Ds")->pt();
        _genDsMeson_Ds_eta[idx] = DsGenParticles[idx].at("Ds")->eta();
        _genDsMeson_Ds_phi[idx] = DsGenParticles[idx].at("Ds")->phi();
        _genDsMeson_Phi_pt[idx] = DsGenParticles[idx].at("Phi")->pt();
        _genDsMeson_Phi_eta[idx] = DsGenParticles[idx].at("Phi")->eta();
        _genDsMeson_Phi_phi[idx] = DsGenParticles[idx].at("Phi")->phi();
        _genDsMeson_Pi_pt[idx] = DsGenParticles[idx].at("Pi")->pt();
        _genDsMeson_Pi_eta[idx] = DsGenParticles[idx].at("Pi")->eta();
        _genDsMeson_Pi_phi[idx] = DsGenParticles[idx].at("Pi")->phi();
        _genDsMeson_KPlus_pt[idx] = DsGenParticles[idx].at("KPlus")->pt();
        _genDsMeson_KPlus_eta[idx] = DsGenParticles[idx].at("KPlus")->eta();
        _genDsMeson_KPlus_phi[idx] = DsGenParticles[idx].at("KPlus")->phi();
        _genDsMeson_KMinus_pt[idx] = DsGenParticles[idx].at("KMinus")->pt();
        _genDsMeson_KMinus_eta[idx] = DsGenParticles[idx].at("KMinus")->eta();
        _genDsMeson_KMinus_phi[idx] = DsGenParticles[idx].at("KMinus")->phi();
    }
}

int DsMesonGenAnalyzer::find_Ds_decay_type(
        const std::vector<reco::GenParticle>& genParticles){
    // find what type of event this is concerning the production and decay of Ds mesons.
    // the numbering convention is as follows:
    // 0: undefined, none of the below.
    // 1: at least one Ds -> phi pi -> K K pi (i.e. the decay of interest).
    // 2: at least one Ds -> ph pi, but excluding the above.
    // 3: at least one Ds, but excluding the above.
    // 4: at least one charmed hadron, but excluding the above.
    
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
        if(res > 4) res = 4;

        // check if it is a D* meson
        if(std::abs(pdgid) != 431) continue;
        const reco::GenParticle* ds = p;
        if(res > 3) res = 3;

        // find the decay products of the Ds meson
        std::vector<const reco::GenParticle*> dsDaughters;
        for(unsigned int i=0; i < ds->numberOfDaughters(); ++i){
            dsDaughters.push_back( &genParticles[ds->daughterRef(i).key()] );
        }

        // find if they are a pion and a phi meson
        if( dsDaughters.size()!=2 ) continue;
        const reco::GenParticle* phi;
        if( std::abs(dsDaughters.at(0)->pdgId())==333
            && std::abs(dsDaughters.at(1)->pdgId())==211 ){
            phi = dsDaughters.at(0);
        } else if( std::abs(dsDaughters.at(0)->pdgId())==211
            && std::abs(dsDaughters.at(1)->pdgId())==333 ){
            phi = dsDaughters.at(1);
        } else continue;
        if(res > 2) res = 2;

        // find the daughters of the phi
        std::vector<const reco::GenParticle*> phiDaughters;
        for(unsigned int i=0; i < phi->numberOfDaughters(); ++i){
            phiDaughters.push_back( &genParticles[phi->daughterRef(i).key()] );
        }

        // find if they are kaons
        if( phiDaughters.size()!=2 ) continue;
        if( std::abs(phiDaughters.at(0)->pdgId())==321
            && std::abs(phiDaughters.at(1)->pdgId())==321 ){
            // pass
        } else continue;

        // if all checks above succeeded,
        // we have a genuine Ds -> phi pi -> K K pi event
        if(res > 1) res = 1;
        break;
    }
    if(res > 4) res = 0;
    return res;
}


std::vector< std::map< std::string, const reco::GenParticle* > > DsMesonGenAnalyzer::find_Ds_to_PhiPi_to_KKPi(
        const std::vector<reco::GenParticle>& genParticles){
    // find Ds -> phi pi -> K K pi at GEN level

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

        // check if it is a Ds meson
        int pdgid = p->pdgId();
        if(std::abs(pdgid) != 431) continue;
        const reco::GenParticle* ds = p;

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
        if( dsDaughters.size()!=2 ) continue;
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
        } else continue;

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
        if( phiDaughters.size()!=2 ) continue;
        if( std::abs(phiDaughters.at(0)->pdgId())==321
            && std::abs(phiDaughters.at(1)->pdgId())==321 ){
            K1 = phiDaughters.at(0);
            K2 = phiDaughters.at(1);
        } else continue;

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
        std::map< std::string, const reco::GenParticle* > thisres = {
            {"Ds", ds},
            {"Phi", phi},
            {"Pi", pi},
            {"KPlus", KPlus},
            {"KMinus", KMinus}
        };
        res.push_back(thisres);
    }
    return res;
}
