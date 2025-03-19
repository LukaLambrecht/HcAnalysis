/*
Custom analyzer class for finding D* mesons triplets of tracks.
*/

// include header
// note: all other includes are defined in the header
#include "HcAnalysis/HcAnalysis/interface/DStarMesonAnalyzer.h"

// constructor //
DStarMesonAnalyzer::DStarMesonAnalyzer(
    const edm::ParameterSet& iConfig, HcAnalysis* hcAnalyzer):
    hcAnalyzer(hcAnalyzer){
};

// destructor //
DStarMesonAnalyzer::~DStarMesonAnalyzer(){
}

// beginJob //
void DStarMesonAnalyzer::beginJob(TTree* outputTree){
    // initialize branches in the output tree
    outputTree->Branch("nDStarMeson", &_nDStarMeson, "nDStarMeson/i");
    outputTree->Branch("DStarMeson_mass", &_DStarMeson_mass, "DStarMeson_mass[nDStarMeson]/D");
    outputTree->Branch("DStarMeson_pt", &_DStarMeson_pt, "DStarMeson_pt[nDStarMeson]/D");
    outputTree->Branch("DStarMeson_eta", &_DStarMeson_eta, "DStarMeson_eta[nDStarMeson]/D");
    outputTree->Branch("DStarMeson_phi", &_DStarMeson_phi, "DStarMeson_phi[nDStarMeson]/D");
    outputTree->Branch("DStarMeson_hasFastGenMatch", &_DStarMeson_hasFastGenMatch, "DStarMeson_hasFastGenMatch[nDStarMeson]/O");
    outputTree->Branch("DStarMeson_hasFastPartialGenMatch", &_DStarMeson_hasFastPartialGenMatch, "DStarMeson_hasFastPartialGenMatch[nDStarMeson]/O");
    outputTree->Branch("DStarMeson_DZeroMeson_mass", &_DStarMeson_DZeroMeson_mass, "DStarMeson_DZeroMeson_mass[nDStarMeson]/D");
    outputTree->Branch("DStarMeson_DZeroMeson_pt", &_DStarMeson_DZeroMeson_pt, "DStarMeson_DZeroMeson_pt[nDStarMeson]/D");
    outputTree->Branch("DStarMeson_DZeroMeson_eta", &_DStarMeson_DZeroMeson_eta, "DStarMeson_DZeroMeson_eta[nDStarMeson]/D");
    outputTree->Branch("DStarMeson_DZeroMeson_phi", &_DStarMeson_DZeroMeson_phi, "DStarMeson_DZeroMeson_phi[nDStarMeson]/D");
    outputTree->Branch("DStarMeson_DZeroMeson_massDiff", &_DStarMeson_DZeroMeson_massDiff, "DStarMeson_DZeroMeson_massDiff[nDStarMeson]/D");
    outputTree->Branch("DStarMeson_Pi1_pt", &_DStarMeson_Pi1_pt, "DStarMeson_Pi1_pt[nDStarMeson]/D");
    outputTree->Branch("DStarMeson_Pi1_eta", &_DStarMeson_Pi1_eta, "DStarMeson_Pi1_eta[nDStarMeson]/D");
    outputTree->Branch("DStarMeson_Pi1_phi", &_DStarMeson_Pi1_phi, "DStarMeson_Pi1_phi[nDStarMeson]/D");
    outputTree->Branch("DStarMeson_K_pt", &_DStarMeson_K_pt, "DStarMeson_K_pt[nDStarMeson]/D");
    outputTree->Branch("DStarMeson_K_eta", &_DStarMeson_K_eta, "DStarMeson_K_eta[nDStarMeson]/D");
    outputTree->Branch("DStarMeson_K_phi", &_DStarMeson_K_phi, "DStarMeson_K_phi[nDStarMeson]/D");
    outputTree->Branch("DStarMeson_Pi2_pt", &_DStarMeson_Pi2_pt, "DStarMeson_Pi2_pt[nDStarMeson]/D");
    outputTree->Branch("DStarMeson_Pi2_eta", &_DStarMeson_Pi2_eta, "DStarMeson_Pi2_eta[nDStarMeson]/D");
    outputTree->Branch("DStarMeson_Pi2_phi", &_DStarMeson_Pi2_phi, "DStarMeson_Pi2_phi[nDStarMeson]/D");
    outputTree->Branch("DStarMeson_tr1tr2_deltaR", &_DStarMeson_tr1tr2_deltaR, "DStarMeson_tr1tr2_deltaR[nDStarMeson]/D");
    outputTree->Branch("DStarMeson_tr3d0_deltaR", &_DStarMeson_tr3d0_deltaR, "DStarMeson_tr3d0_deltaR[nDStarMeson]/D");
    outputTree->Branch("DStarMeson_d0vtx_normchi2", &_DStarMeson_d0vtx_normchi2, "DStarMeson_d0vtx_normchi2[nDStarMeson]/D");
    outputTree->Branch("DStarMeson_dstarvtx_normchi2", &_DStarMeson_dstarvtx_normchi2, "DStarMeson_dstarvtx_normchi2[nDStarMeson]/D");
    outputTree->Branch("DStarMeson_tr1tr2_sepx", &_DStarMeson_tr1tr2_sepx, "DStarMeson_tr1tr2_sepx[nDStarMeson]/D");
    outputTree->Branch("DStarMeson_tr1tr2_sepy", &_DStarMeson_tr1tr2_sepy, "DStarMeson_tr1tr2_sepy[nDStarMeson]/D");
    outputTree->Branch("DStarMeson_tr1tr2_sepz", &_DStarMeson_tr1tr2_sepz, "DStarMeson_tr1tr2_sepz[nDStarMeson]/D");
    outputTree->Branch("DStarMeson_tr3d0_sepx", &_DStarMeson_tr3d0_sepx, "DStarMeson_tr3d0_sepx[nDStarMeson]/D");
    outputTree->Branch("DStarMeson_tr3d0_sepy", &_DStarMeson_tr3d0_sepy, "DStarMeson_tr3d0_sepy[nDStarMeson]/D");
    outputTree->Branch("DStarMeson_tr3d0_sepz", &_DStarMeson_tr3d0_sepz, "DStarMeson_tr3d0_sepz[nDStarMeson]/D");
}

// analyze (main method) //
void DStarMesonAnalyzer::analyze(const edm::Event& iEvent){

    // get all required objects from tokens
    edm::Handle<std::vector<pat::PackedCandidate>> packedPFCandidates;
    iEvent.getByToken(hcAnalyzer->packedPFCandidatesToken, packedPFCandidates);
    edm::Handle<std::vector<pat::PackedCandidate>> lostTracks;
    iEvent.getByToken(hcAnalyzer->lostTracksToken, lostTracks);
    MagneticField* bfield = new OAEParametrizedMagneticField("3_8T");

    // settings for gen-matching
    edm::Handle<std::vector<reco::GenParticle>> genParticles;
    iEvent.getByToken(hcAnalyzer->prunedGenParticlesToken, genParticles);
    std::vector< std::map< std::string, const reco::GenParticle* > > DStarGenParticles;
    bool doMatching = true; // to do: explicitly disable for data
    if( !genParticles.isValid() ) doMatching = false;
    if( doMatching ){
        DStarGenParticles = DStarMesonGenAnalyzer::find_DStar_to_DZeroPi_to_KPiPi( *genParticles );
        if( DStarGenParticles.size()==0 ) doMatching = false;
    }

    // merge packed candidate tracks and lost tracks
    std::vector<reco::Track> allTracks;
    for(const pat::PackedCandidate& pc: *packedPFCandidates){
        if(pc.hasTrackDetails()){
            reco::Track track = *pc.bestTrack();
            allTracks.push_back(track);
        }
    }

    for(const pat::PackedCandidate& pc: *lostTracks){
        if(pc.hasTrackDetails()){
            reco::Track track = *pc.bestTrack();
            allTracks.push_back(track);
        }
    }

    // preselect tracks
    std::vector<reco::Track> selectedTracks;
    for(const reco::Track& track: allTracks){
        if(!track.quality(reco::TrackBase::qualityByName("highPurity"))) continue;
        if(track.pt() < 0.3) continue;
	    selectedTracks.push_back(track);
    }

    // initialize number of DStarMesons
    _nDStarMeson = 0;

    // loop over pairs of tracks
    for(unsigned i=0; i<selectedTracks.size(); i++){
      for(unsigned j=i+1; j<selectedTracks.size(); j++){
        const reco::Track tr1 = selectedTracks.at(i);
        const reco::Track tr2 = selectedTracks.at(j);

        // candidates must have opposite charge
        if(tr1.charge() * tr2.charge() > 0) continue;

        // candidates must point approximately in the same direction
        if( reco::deltaR(tr1, tr2) > 0.4 ) continue;

        // reference points of both tracks must be close together
        const math::XYZPoint tr1refpoint = tr1.referencePoint();
        const math::XYZPoint tr2refpoint = tr2.referencePoint();
        double twotracksepx = std::abs(tr1refpoint.x()-tr2refpoint.x());
        double twotracksepy = std::abs(tr1refpoint.y()-tr2refpoint.y());
        double twotracksepz = std::abs(tr1refpoint.z()-tr2refpoint.z());
        if( twotracksepx>0.1 || twotracksepy>0.1 || twotracksepz>0.1 ) continue;
       
        // find which track is positive and which is negative
        reco::Track postrack;
        reco::Track negtrack;
        if(tr1.charge()>0. and tr2.charge()<0){
            postrack = tr1;
            negtrack = tr2;
        } else if(tr1.charge()<0. and tr2.charge()>0){
	        postrack = tr2;
	        negtrack = tr1;
        } else continue; // should not normally happen but just for safety

        // make invariant mass
        // (note: we can safely assume the negative track is a kaon, and the positive one a pion,
        //  as the opposite case has a negligible branching fraction)
        ROOT::Math::PtEtaPhiMVector pi2P4(postrack.pt(), postrack.eta(), postrack.phi(), pimass);
        ROOT::Math::PtEtaPhiMVector KP4(negtrack.pt(), negtrack.eta(), negtrack.phi(), kmass);
        ROOT::Math::PtEtaPhiMVector dzeroP4 = pi2P4 + KP4;
        double dzeroInvMass = dzeroP4.M();
        
        // invariant mass must be close to resonance mass
        if( std::abs(dzeroInvMass - dzeromass) > 0.035 ) continue;

        // fit a vertex
        std::vector<reco::TransientTrack> transpair;
        transpair.push_back(reco::TransientTrack(tr1, bfield));
        transpair.push_back(reco::TransientTrack(tr2, bfield));
        KalmanVertexFitter vtxFitter(false);
        TransientVertex dzerovtx = vtxFitter.vertex(transpair);
        // vertex must be valid
        if(!dzerovtx.isValid()) continue;
        // chi squared of fit must be small
        if(dzerovtx.normalisedChiSquared()>5.) continue;
        if(dzerovtx.normalisedChiSquared()<0.) continue;
        
        // loop over third track
	    for(unsigned k=0; k<selectedTracks.size(); k++){
            if(k==i or k==j) continue;
            const reco::Track tr3 = selectedTracks.at(k);

            // candidates must point approximately in the same direction
            if( reco::deltaR(tr3, dzeroP4) > 0.4 ) continue;

            // reference point of third track must be close to phi vertex
            const math::XYZPoint tr3refpoint = tr3.referencePoint();
            double trackvtxsepx = std::abs(tr3refpoint.x()-dzerovtx.position().x());
            double trackvtxsepy = std::abs(tr3refpoint.y()-dzerovtx.position().y());
            double trackvtxsepz = std::abs(tr3refpoint.z()-dzerovtx.position().z());
            if( trackvtxsepx>0.1 || trackvtxsepy>0.1 || trackvtxsepz>0.1 ) continue;

            // make invariant mass (under the assumption of pi mass for the third track)
            ROOT::Math::PtEtaPhiMVector pi1P4(tr3.pt(), tr3.eta(), tr3.phi(), pimass);
            ROOT::Math::PtEtaPhiMVector dstarP4 = dzeroP4 + pi1P4;
            double dstarInvMass = dstarP4.M();

            // check if mass is close enough to D* mass
            if(std::abs(dstarInvMass - dstarmass) > 0.1) continue;

            // do a vertex fit
            std::vector<reco::TransientTrack> transtriplet;
            transtriplet.push_back(reco::TransientTrack(tr1, bfield));
            transtriplet.push_back(reco::TransientTrack(tr2, bfield));
            transtriplet.push_back(reco::TransientTrack(tr3, bfield));
            TransientVertex dstarvtx = vtxFitter.vertex(transtriplet);
            if(!dstarvtx.isValid()) continue;
            if(dstarvtx.normalisedChiSquared()>5.) continue;
            if(dstarvtx.normalisedChiSquared()<0.) continue;

            // set properties of the D* candidate
            _DStarMeson_mass[_nDStarMeson] = dstarP4.M();
            _DStarMeson_pt[_nDStarMeson] = dstarP4.pt();
            _DStarMeson_eta[_nDStarMeson] = dstarP4.eta();
            _DStarMeson_phi[_nDStarMeson] = dstarP4.phi();
            _DStarMeson_DZeroMeson_mass[_nDStarMeson] = dzeroP4.M();
            _DStarMeson_DZeroMeson_pt[_nDStarMeson] = dzeroP4.pt();
            _DStarMeson_DZeroMeson_eta[_nDStarMeson] = dzeroP4.eta();
            _DStarMeson_DZeroMeson_phi[_nDStarMeson] = dzeroP4.phi();
            _DStarMeson_DZeroMeson_massDiff[_nDStarMeson] = dstarP4.M() - dzeroP4.M();
            _DStarMeson_Pi1_pt[_nDStarMeson] = pi1P4.pt();
            _DStarMeson_Pi1_eta[_nDStarMeson] = pi1P4.eta();
            _DStarMeson_Pi1_phi[_nDStarMeson] = pi1P4.phi();
            _DStarMeson_K_pt[_nDStarMeson] = KP4.pt();
            _DStarMeson_K_eta[_nDStarMeson] = KP4.eta();
            _DStarMeson_K_phi[_nDStarMeson] = KP4.phi();
            _DStarMeson_Pi2_pt[_nDStarMeson] = pi2P4.pt();
            _DStarMeson_Pi2_eta[_nDStarMeson] = pi2P4.eta();
            _DStarMeson_Pi2_phi[_nDStarMeson] = pi2P4.phi();
            _DStarMeson_tr1tr2_deltaR[_nDStarMeson] = reco::deltaR(tr1, tr2);
            _DStarMeson_tr3d0_deltaR[_nDStarMeson] = reco::deltaR(tr3, dzeroP4);
            _DStarMeson_d0vtx_normchi2[_nDStarMeson] = dzerovtx.normalisedChiSquared();
            _DStarMeson_dstarvtx_normchi2[_nDStarMeson] = dstarvtx.normalisedChiSquared();
            _DStarMeson_tr1tr2_sepx[_nDStarMeson] = twotracksepx;
            _DStarMeson_tr1tr2_sepy[_nDStarMeson] = twotracksepy;
            _DStarMeson_tr1tr2_sepz[_nDStarMeson] = twotracksepz;
            _DStarMeson_tr3d0_sepx[_nDStarMeson] = trackvtxsepx;
            _DStarMeson_tr3d0_sepy[_nDStarMeson] = trackvtxsepy;
            _DStarMeson_tr3d0_sepz[_nDStarMeson] = trackvtxsepz;            

            // check if this candidate can be matched to gen-level
            _DStarMeson_hasFastGenMatch[_nDStarMeson] = false;
            _DStarMeson_hasFastGenMatch[_nDStarMeson] = false;
            if( doMatching ){
                for( const auto& pmap : DStarGenParticles){
                    double dRThreshold = 0.05;
                    if( GenTools::isGeometricTrackMatch( tr3, *pmap.at("Pi1"), dRThreshold )
                        && ( (GenTools::isGeometricTrackMatch( postrack, *pmap.at("K"), dRThreshold )
                              && GenTools::isGeometricTrackMatch( negtrack, *pmap.at("Pi2"), dRThreshold ) )
                             || (GenTools::isGeometricTrackMatch( postrack, *pmap.at("Pi2"), dRThreshold )
                              && GenTools::isGeometricTrackMatch( negtrack, *pmap.at("K"), dRThreshold ) ) ) ){
                        _DStarMeson_hasFastGenMatch[_nDStarMeson] = true;
                    }
                    if( GenTools::isGeometricTrackMatch( tr3, *pmap.at("Pi1"), dRThreshold )
                         || GenTools::isGeometricTrackMatch( postrack, *pmap.at("K"), dRThreshold )
                         || GenTools::isGeometricTrackMatch( negtrack, *pmap.at("Pi2"), dRThreshold )
                         || GenTools::isGeometricTrackMatch( postrack, *pmap.at("Pi2"), dRThreshold )
                         || GenTools::isGeometricTrackMatch( negtrack, *pmap.at("K"), dRThreshold ) ){
                        _DStarMeson_hasFastPartialGenMatch[_nDStarMeson] = true;
                    }
                }
            }

            // update counter
            _nDStarMeson++;
 
           if( _nDStarMeson == nDStarMeson_max ) break;

        } // end loop over third track
        if( _nDStarMeson == nDStarMeson_max ) break;
      }
      if( _nDStarMeson == nDStarMeson_max) break;
    } // end loop over first and second track
    delete bfield;
}
