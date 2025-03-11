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
    outputTree->Branch("DsMeson_PhiMeson_mass", &_DsMeson_PhiMeson_mass, "DsMeson_PhiMeson_mass[nDsMeson]/D");
    outputTree->Branch("DsMeson_PhiMeson_pt", &_DsMeson_PhiMeson_pt, "DsMeson_PhiMeson_pt[nDsMeson]/D");
    outputTree->Branch("DsMeson_PhiMeson_eta", &_DsMeson_PhiMeson_eta, "DsMeson_PhiMeson_eta[nDsMeson]/D");
    outputTree->Branch("DsMeson_PhiMeson_phi", &_DsMeson_PhiMeson_phi, "DsMeson_PhiMeson_phi[nDsMeson]/D");
}

// analyze (main method) //
void DsMesonAnalyzer::analyze(const edm::Event& iEvent){

    // get all required objects from tokens
    edm::Handle<std::vector<pat::PackedCandidate>> packedPFCandidates;
    iEvent.getByToken(hcAnalyzer->packedPFCandidatesToken, packedPFCandidates);
    edm::Handle<std::vector<pat::PackedCandidate>> lostTracks;
    iEvent.getByToken(hcAnalyzer->lostTracksToken, lostTracks);
    MagneticField* bfield = new OAEParametrizedMagneticField("3_8T");

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
        if(track.pt() < 0.5) continue;
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

        // the following fragment applies some conditions on the point of closest approach,
        // namely that the distance of closest approach must be reasonably small
        // and that the point of closest approach must be situated not too far from the center of CMS.
        // update: now replaced by a method using the track reference point below.
        /*reco::TransientTrack trtr1(tr1, bfield);
        reco::TransientTrack trtr2(tr2, bfield);
        if(!trtr1.impactPointTSCP().isValid() or !trtr2.impactPointTSCP().isValid()) continue;
        FreeTrajectoryState state1 = trtr1.impactPointTSCP().theState();
        FreeTrajectoryState state2 = trtr2.impactPointTSCP().theState();
        ClosestApproachInRPhi capp; capp.calculate(state1,state2);
        if(!capp.status()) continue;
        double dca = fabs(capp.distance());
        if(dca < 0.) continue;
        if(dca > 1.) continue;
        GlobalPoint cxpt = capp.crossingPoint();
        if(std::sqrt(cxpt.x()*cxpt.x() + cxpt.y()*cxpt.y())>10.
            or std::abs(cxpt.z())>10.) continue;*/

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
        ROOT::Math::PtEtaPhiMVector posP4(postrack.pt(), postrack.eta(), postrack.phi(), kmass);
        ROOT::Math::PtEtaPhiMVector negP4(negtrack.pt(), negtrack.eta(), negtrack.phi(), kmass);
        ROOT::Math::PtEtaPhiMVector phiP4 = posP4 + negP4;
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
        if(phivtx.normalisedChiSquared()>15.) continue;
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
            if( trackvtxsepx>0.5 || trackvtxsepy>0.5 || trackvtxsepz>0.5 ) continue;

            // make invariant mass (under the assumption of pi mass for the third track)
            ROOT::Math::PtEtaPhiMVector thirdP4(tr3.pt(), tr3.eta(), tr3.phi(), pimass);
            ROOT::Math::PtEtaPhiMVector dsP4 = phiP4 + thirdP4;
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
            _nDsMeson++;
            if( _nDsMeson == nDsMeson_max ) break;

        } // end loop over third track
        if( _nDsMeson == nDsMeson_max ) break;
      }
      if( _nDsMeson == nDsMeson_max) break;
    } // end loop over first and second track
    delete bfield;
}
