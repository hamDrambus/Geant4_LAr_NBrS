#include <geant4/SteppingAction.hh>

SteppingAction::SteppingAction() : N_iterations(3), epsilon(1e-10), stuck_N_iterations(0)
{}

void SteppingAction::UserSteppingAction(const G4Step* theStep)
{
	G4Track* theTrack = theStep->GetTrack();
	G4int trackID = theTrack->GetTrackID();
	G4ParticleDefinition* particleType = theTrack->GetDefinition();

	G4StepPoint* thePrePoint = theStep->GetPreStepPoint();
	G4VPhysicalVolume* thePrePV = thePrePoint->GetPhysicalVolume();

	G4StepPoint* thePostPoint = theStep->GetPostStepPoint();
	G4VPhysicalVolume* thePostPV = thePostPoint->GetPhysicalVolume();

	G4ThreeVector pos;
	std::string postVolName;
	std::string prevVolName;
	G4ThreeVector MomentumDirection;
	if (thePostPV != NULL) {
		pos = theStep->GetPostStepPoint()->GetPosition();
		postVolName = theStep->GetPostStepPoint()->GetPhysicalVolume()->GetName();
		MomentumDirection = theStep->GetPostStepPoint()->GetMomentumDirection();
	}
	if (thePrePV != nullptr)
	  prevVolName = thePrePV->GetName();
	if (CheckStuckParticle(theStep)) {
	  if (1 == stuck_N_iterations && gPars::general.teleportation_verbosity > 0) {
	    std::cerr<<"SteppingAction::UserSteppingAction:Warning! Particle stuck for "<<N_iterations<<" iterations!"<<std::endl;
	  }
	}

//	std::cout<<"************"<<std::endl;
//	std::cout<<"PrevVolume: \""<<prevVolName<<"\""<<std::endl;
//	std::cout<<"PostVolume: \""<<postVolName<<"\""<<std::endl;
	return G4UserSteppingAction::UserSteppingAction(theStep);
}

bool SteppingAction::CheckStuckParticle(const G4Step* step)
{
  //G4StepPoint* thePostPoint = step->GetPostStepPoint();
  bool is_stuck = false;
  G4VPhysicalVolume* thePostPV = step->GetPostStepPoint()->GetPhysicalVolume();
  if (thePostPV != NULL) {
    G4ThreeVector pos = step->GetPostStepPoint()->GetPosition();
    G4ThreeVector mom = step->GetPostStepPoint()->GetMomentumDirection();
    std::string vol_name = step->GetPostStepPoint()->GetPhysicalVolume()->GetName();
    double pos_mag = pos.mag();
    double mom_mag = mom.mag();
    is_stuck = (positions.size() >= N_iterations);
    for (std::size_t i = 0, i_end_ = positions.size(); i!=i_end_ && is_stuck; ++i) {
      if (pos_mag < epsilon*epsilon) {
        if ((pos - positions[i]).mag() > epsilon) {
          is_stuck = false;
          break;
        }
      } else {
        if ((pos - positions[i]).mag()/pos_mag > epsilon) {
          is_stuck = false;
          break;
        }
      }
      if (mom_mag < epsilon*epsilon) {
        if ((mom - momentums[i]).mag() > epsilon) {
          is_stuck = false;
          break;
        }
      } else {
        if ((mom - momentums[i]).mag()/mom_mag > epsilon) {
          is_stuck = false;
          break;
        }
      }
    }
    positions.push_back(pos);
    momentums.push_back(mom);
    if (positions.size() > N_iterations) {
      positions.pop_front();
      momentums.pop_front();
    }
  }
  if (!is_stuck)
    stuck_N_iterations = 0;
  else
    ++stuck_N_iterations;
  return is_stuck;
}
