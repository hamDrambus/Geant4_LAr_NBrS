#include "SteppingAction.hh"

SteppingAction::SteppingAction()
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
		MomentumDirection = theStep->GetTrack()->GetMomentumDirection();
	}
	if (thePrePV != nullptr)
	  prevVolName = thePrePV->GetName();

	//std::cout<<"************"<<std::endl;
	//std::cout<<"PrevVolume: \""<<prevVolName<<"\""<<std::endl;
	//std::cout<<"PostVolume: \""<<postVolName<<"\""<<std::endl;
	return G4UserSteppingAction::UserSteppingAction(theStep);
}
