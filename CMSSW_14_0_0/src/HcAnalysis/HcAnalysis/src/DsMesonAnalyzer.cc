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

    outputTree->Branch("_nDsMesons", &_nDsMesons, "_nDsMesons/i");
    outputTree->Branch("_DsMesonInvMass", &_DsMesonInvMass, "_DsMesonInvMass[_nDsMesons]/D");
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
        if(!track.quality(reco::TrackBase::qualityByName("loose"))) continue;
	    selectedTracks.push_back(track);
    }

    // initialize number of DsMesons
    _nDsMesons = 0;

    // loop over pairs of tracks
    for(unsigned i=0; i<selectedTracks.size(); i++){
      for(unsigned j=i+1; j<selectedTracks.size(); j++){
        const reco::Track tr1 = selectedTracks.at(i);
        const reco::Track tr2 = selectedTracks.at(j);

        // candidates must have opposite charge
        if(tr1.charge()*tr2.charge()>0) continue;

        // the following fragment applies some conditions on the point of closest approach,
        // namely that the distance of closest approach must be reasonably small
        // and that the point of closest approach must be situated not too far from the center of CMS.
        reco::TransientTrack trtr1(tr1, bfield);
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
            or std::abs(cxpt.z())>10.) continue;
       
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
        // todo

        // check if mass is close enough to phi mass
        // todo
        
        // loop over third track
	    for(unsigned k=0; k<selectedTracks.size(); k++){
            if(k==i or k==j) continue;
            const reco::Track tr3 = selectedTracks.at(k);

            // make invariant mass (under the assumption of pi mass for the third track)
            // todo

            // check if mass is close enough to Ds mass
            // todo

        } // end loop over third track
    } } // end loop over first and second track
    delete bfield;
}
