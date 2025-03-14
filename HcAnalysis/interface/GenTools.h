/*
Tools for working with gen-level particles and decay chains
*/

#ifndef GenTools_H
#define GenTools_H

#include <set>

#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "TLorentzVector.h"

namespace GenTools{
    const reco::GenParticle* getFirstMother(const reco::GenParticle&, const std::vector<reco::GenParticle>&);
    const int getFirstMotherIndex(const reco::GenParticle&, const std::vector<reco::GenParticle>&);
    const reco::GenParticle* getMother(const reco::GenParticle&, const std::vector<reco::GenParticle>&);
    int getMotherPdgId(const reco::GenParticle&, const std::vector<reco::GenParticle>&);
}
#endif
