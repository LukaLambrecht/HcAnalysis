/*
Custom analyzer class for finding Ds mesons triplets of tracks.
*/

// include header
// note: all other includes are defined in the header
#include "HcAnalysis/HcAnalysis/interface/DsMesonAnalyzer.h"

// constructor //
DsMesonAnalyzer::DsMesonAnalyzer(
    const edm::ParameterSet& iConfig, HcAnalysis* hcAnalyzer):
    hcAnalyzer(hcAnalyzer){
};

// destructor //
DsMesonAnalyzer::~DsMesonAnalyzer(){
}

// beginJob //
void DsMesonAnalyzer::beginJob(TTree* outputTree){
    // initialize branches in the output tree
    outputTree->Branch("nDsMeson", &_nDsMeson, "nDsMeson/i");
    outputTree->Branch("DsMeson_mass", &_DsMeson_mass, "DsMeson_mass[nDsMeson]/D");
    outputTree->Branch("DsMeson_pt", &_DsMeson_pt, "DsMeson_pt[nDsMeson]/D");
    outputTree->Branch("DsMeson_eta", &_DsMeson_eta, "DsMeson_eta[nDsMeson]/D");
    outputTree->Branch("DsMeson_phi", &_DsMeson_phi, "DsMeson_phi[nDsMeson]/D");
    outputTree->Branch("DsMeson_hasFastGenMatch", &_DsMeson_hasFastGenMatch, "DsMeson_hasFastGenMatch[nDsMeson]/O");
    outputTree->Branch("DsMeson_hasFastPartialGenMatch", &_DsMeson_hasFastPartialGenMatch, "DsMeson_hasFastPartialGenMatch[nDsMeson]/O");
    outputTree->Branch("DsMeson_PhiMeson_mass", &_DsMeson_PhiMeson_mass, "DsMeson_PhiMeson_mass[nDsMeson]/D");
    outputTree->Branch("DsMeson_PhiMeson_pt", &_DsMeson_PhiMeson_pt, "DsMeson_PhiMeson_pt[nDsMeson]/D");
    outputTree->Branch("DsMeson_PhiMeson_eta", &_DsMeson_PhiMeson_eta, "DsMeson_PhiMeson_eta[nDsMeson]/D");
    outputTree->Branch("DsMeson_PhiMeson_phi", &_DsMeson_PhiMeson_phi, "DsMeson_PhiMeson_phi[nDsMeson]/D");
    outputTree->Branch("DsMeson_PhiMeson_massDiff", &_DsMeson_PhiMeson_massDiff, "DsMeson_PhiMeson_massDiff[nDsMeson]/D");
    outputTree->Branch("DsMeson_Pi_pt", &_DsMeson_Pi_pt, "DsMeson_Pi_pt[nDsMeson]/D");
    outputTree->Branch("DsMeson_Pi_eta", &_DsMeson_Pi_eta, "DsMeson_Pi_eta[nDsMeson]/D");
    outputTree->Branch("DsMeson_Pi_phi", &_DsMeson_Pi_phi, "DsMeson_Pi_phi[nDsMeson]/D");
    outputTree->Branch("DsMeson_KPlus_pt", &_DsMeson_KPlus_pt, "DsMeson_KPlus_pt[nDsMeson]/D");
    outputTree->Branch("DsMeson_KPlus_eta", &_DsMeson_KPlus_eta, "DsMeson_KPlus_eta[nDsMeson]/D");
    outputTree->Branch("DsMeson_KPlus_phi", &_DsMeson_KPlus_phi, "DsMeson_KPlus_phi[nDsMeson]/D");
    outputTree->Branch("DsMeson_KMinus_pt", &_DsMeson_KMinus_pt, "DsMeson_KMinus_pt[nDsMeson]/D");
    outputTree->Branch("DsMeson_KMinus_eta", &_DsMeson_KMinus_eta, "DsMeson_KMinus_eta[nDsMeson]/D");
    outputTree->Branch("DsMeson_KMinus_phi", &_DsMeson_KMinus_phi, "DsMeson_KMinus_phi[nDsMeson]/D");
    outputTree->Branch("DsMeson_tr1tr2_deltaR", &_DsMeson_tr1tr2_deltaR, "DsMeson_tr1tr2_deltaR[nDsMeson]/D");
    outputTree->Branch("DsMeson_tr3phi_deltaR", &_DsMeson_tr3phi_deltaR, "DsMeson_tr3phi_deltaR[nDsMeson]/D");
    outputTree->Branch("DsMeson_phivtx_normchi2", &_DsMeson_phivtx_normchi2, "DsMeson_phivtx_normchi2[nDsMeson]/D");
    outputTree->Branch("DsMeson_dsvtx_normchi2", &_DsMeson_dsvtx_normchi2, "DsMeson_dsvtx_normchi2[nDsMeson]/D");
    outputTree->Branch("DsMeson_tr1tr2_sepx", &_DsMeson_tr1tr2_sepx, "DsMeson_tr1tr2_sepx[nDsMeson]/D");
    outputTree->Branch("DsMeson_tr1tr2_sepy", &_DsMeson_tr1tr2_sepy, "DsMeson_tr1tr2_sepy[nDsMeson]/D");
    outputTree->Branch("DsMeson_tr1tr2_sepz", &_DsMeson_tr1tr2_sepz, "DsMeson_tr1tr2_sepz[nDsMeson]/D");
    outputTree->Branch("DsMeson_tr3phi_sepx", &_DsMeson_tr3phi_sepx, "DsMeson_tr3phi_sepx[nDsMeson]/D");
    outputTree->Branch("DsMeson_tr3phi_sepy", &_DsMeson_tr3phi_sepy, "DsMeson_tr3phi_sepy[nDsMeson]/D");
    outputTree->Branch("DsMeson_tr3phi_sepz", &_DsMeson_tr3phi_sepz, "DsMeson_tr3phi_sepz[nDsMeson]/D");
}

// analyze (main method) //
void DsMesonAnalyzer::analyze(const edm::Event& iEvent){

    // get all required objects from tokens
    edm::Handle<std::vector<pat::PackedCandidate>> packedPFCandidates;
    iEvent.getByToken(hcAnalyzer->packedPFCandidatesToken, packedPFCandidates);
    edm::Handle<std::vector<pat::PackedCandidate>> lostTracks;
    iEvent.getByToken(hcAnalyzer->lostTracksToken, lostTracks);
    MagneticField* bfield = new OAEParametrizedMagneticField("3_8T");

    // settings for gen-matching
    edm::Handle<std::vector<reco::GenParticle>> genParticles;
    iEvent.getByToken(hcAnalyzer->prunedGenParticlesToken, genParticles);
    std::vector< std::map< std::string, const reco::GenParticle* > > DsGenParticles;
    bool doMatching = true; // to do: explicitly disable for data
    if( !genParticles.isValid() ) doMatching = false;
    if( doMatching ){
        DsGenParticles = DsMesonGenAnalyzer::find_Ds_to_PhiPi_to_KKPi( *genParticles );
        if( DsGenParticles.size()==0 ) doMatching = false;
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

    // initialize number of DsMesons
    _nDsMeson = 0;

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

        // make invariant mass (under the assumption of K mass for both tracks)
        ROOT::Math::PtEtaPhiMVector KPlusP4(postrack.pt(), postrack.eta(), postrack.phi(), kmass);
        ROOT::Math::PtEtaPhiMVector KMinusP4(negtrack.pt(), negtrack.eta(), negtrack.phi(), kmass);
        ROOT::Math::PtEtaPhiMVector phiP4 = KPlusP4 + KMinusP4;
        double phiInvMass = phiP4.M();

        // check if mass is close enough to phi mass
        if(std::abs(phiInvMass - phimass) > 0.07) continue;

        // fit a vertex
        std::vector<reco::TransientTrack> transpair;
        transpair.push_back(reco::TransientTrack(tr1, bfield));
        transpair.push_back(reco::TransientTrack(tr2, bfield));
        KalmanVertexFitter vtxFitter(false);
        TransientVertex phivtx = vtxFitter.vertex(transpair);
        // vertex must be valid
        if(!phivtx.isValid()) continue;
        // chi squared of fit must be small
        if(phivtx.normalisedChiSquared()>5.) continue;
        if(phivtx.normalisedChiSquared()<0.) continue;
        
        // loop over third track
	    for(unsigned k=0; k<selectedTracks.size(); k++){
            if(k==i or k==j) continue;
            const reco::Track tr3 = selectedTracks.at(k);

            // candidates must point approximately in the same direction
            if( reco::deltaR(tr3, phiP4) > 0.4 ) continue;

            // reference point of third track must be close to phi vertex
            const math::XYZPoint tr3refpoint = tr3.referencePoint();
            double trackvtxsepx = std::abs(tr3refpoint.x()-phivtx.position().x());
            double trackvtxsepy = std::abs(tr3refpoint.y()-phivtx.position().y());
            double trackvtxsepz = std::abs(tr3refpoint.z()-phivtx.position().z());
            if( trackvtxsepx>0.1 || trackvtxsepy>0.1 || trackvtxsepz>0.1 ) continue;

            // make invariant mass (under the assumption of pi mass for the third track)
            ROOT::Math::PtEtaPhiMVector piP4(tr3.pt(), tr3.eta(), tr3.phi(), pimass);
            ROOT::Math::PtEtaPhiMVector dsP4 = phiP4 + piP4;
            double dsInvMass = dsP4.M();

            // check if mass is close enough to Ds mass
            if(std::abs(dsInvMass - dsmass) > 0.1) continue;

            // do a vertex fit
            std::vector<reco::TransientTrack> transtriplet;
            transtriplet.push_back(reco::TransientTrack(tr1, bfield));
            transtriplet.push_back(reco::TransientTrack(tr2, bfield));
            transtriplet.push_back(reco::TransientTrack(tr3, bfield));
            TransientVertex dsvtx = vtxFitter.vertex(transtriplet);
            if(!dsvtx.isValid()) continue;
            if(dsvtx.normalisedChiSquared()>5.) continue;
            if(dsvtx.normalisedChiSquared()<0.) continue;

            // set properties of the Ds candidate
            _DsMeson_mass[_nDsMeson] = dsP4.M();
            _DsMeson_pt[_nDsMeson] = dsP4.pt();
            _DsMeson_eta[_nDsMeson] = dsP4.eta();
            _DsMeson_phi[_nDsMeson] = dsP4.phi();
            _DsMeson_PhiMeson_mass[_nDsMeson] = phiP4.M();
            _DsMeson_PhiMeson_pt[_nDsMeson] = phiP4.pt();
            _DsMeson_PhiMeson_eta[_nDsMeson] = phiP4.eta();
            _DsMeson_PhiMeson_phi[_nDsMeson] = phiP4.phi();
            _DsMeson_PhiMeson_massDiff[_nDsMeson] = dsP4.M() - phiP4.M();
            _DsMeson_Pi_pt[_nDsMeson] = piP4.pt();
            _DsMeson_Pi_eta[_nDsMeson] = piP4.eta();
            _DsMeson_Pi_phi[_nDsMeson] = piP4.phi();
            _DsMeson_KPlus_pt[_nDsMeson] = KPlusP4.pt();
            _DsMeson_KPlus_eta[_nDsMeson] = KPlusP4.eta();
            _DsMeson_KPlus_phi[_nDsMeson] = KPlusP4.phi();
            _DsMeson_KMinus_pt[_nDsMeson] = KMinusP4.pt();
            _DsMeson_KMinus_eta[_nDsMeson] = KMinusP4.eta();
            _DsMeson_KMinus_phi[_nDsMeson] = KMinusP4.phi();
            _DsMeson_tr1tr2_deltaR[_nDsMeson] = reco::deltaR(tr1, tr2);
            _DsMeson_tr3phi_deltaR[_nDsMeson] = reco::deltaR(tr3, phiP4);
            _DsMeson_phivtx_normchi2[_nDsMeson] = phivtx.normalisedChiSquared();
            _DsMeson_dsvtx_normchi2[_nDsMeson] = dsvtx.normalisedChiSquared();
            _DsMeson_tr1tr2_sepx[_nDsMeson] = twotracksepx;
            _DsMeson_tr1tr2_sepy[_nDsMeson] = twotracksepy;
            _DsMeson_tr1tr2_sepz[_nDsMeson] = twotracksepz;
            _DsMeson_tr3phi_sepx[_nDsMeson] = trackvtxsepx;
            _DsMeson_tr3phi_sepy[_nDsMeson] = trackvtxsepy;
            _DsMeson_tr3phi_sepz[_nDsMeson] = trackvtxsepz;

            // check if this candidate can be matched to gen-level
            _DsMeson_hasFastGenMatch[_nDsMeson] = false;
            _DsMeson_hasFastGenMatch[_nDsMeson] = false;
            if( doMatching ){
                for( const auto& pmap : DsGenParticles){
                    double dRThreshold = 0.05;
                    if( GenTools::isGeometricTrackMatch( tr3, *pmap.at("Pi"), dRThreshold )
                        && GenTools::isGeometricTrackMatch( postrack, *pmap.at("KPlus"), dRThreshold )
                        && GenTools::isGeometricTrackMatch( negtrack, *pmap.at("KMinus"), dRThreshold ) ){
                        _DsMeson_hasFastGenMatch[_nDsMeson] = true;
                    }
                    if( GenTools::isGeometricTrackMatch( tr3, *pmap.at("Pi"), dRThreshold )
                        || GenTools::isGeometricTrackMatch( postrack, *pmap.at("KPlus"), dRThreshold )
                        || GenTools::isGeometricTrackMatch( negtrack, *pmap.at("KMinus"), dRThreshold ) ){
                        _DsMeson_hasFastPartialGenMatch[_nDsMeson] = true;
                    }
                }
            }

            // update counter
            _nDsMeson++;

            // break loop over third track in case maximum number was reached
            if( _nDsMeson == nDsMeson_max ) break;

        } // end loop over third track
        if( _nDsMeson == nDsMeson_max ) break;
      }
      if( _nDsMeson == nDsMeson_max) break;
    } // end loop over first and second track
    delete bfield;
}
