/*
Tools for working with gen-level particles and decay chains
*/

#include "HcAnalysis/HcAnalysis/interface/GenTools.h"

const reco::GenParticle* GenTools::getFirstMother(
        const reco::GenParticle& gen,
        const std::vector<reco::GenParticle>& genParticles){
    if(gen.numberOfMothers() == 0) return nullptr;
    return &genParticles[gen.motherRef(0).key()];
}

const int GenTools::getFirstMotherIndex(
        const reco::GenParticle& gen,
        const std::vector<reco::GenParticle>& genParticles){
    if(gen.numberOfMothers() == 0) return -1;
    return gen.motherRef(0).key();
}

const reco::GenParticle* GenTools::getMother(
        const reco::GenParticle& gen,
        const std::vector<reco::GenParticle>& genParticles){
    const reco::GenParticle* mom = getFirstMother(gen, genParticles);
    if(!mom) return nullptr;
    else if(mom->pdgId() == gen.pdgId()) return getMother(*mom, genParticles);
    else return mom;
}

int GenTools::getMotherPdgId(
        const reco::GenParticle& gen,
        const std::vector<reco::GenParticle>& genParticles){
    const reco::GenParticle* mom = getMother(gen, genParticles);
    if(!mom) return 0;
    return mom->pdgId();
}
