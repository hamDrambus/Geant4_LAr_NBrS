#include <time.h>
#include <iostream>
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#include "GlobalParameters.hh"
#include "Randomize.hh"
#include "PhysicsList.hh"
#include "DetectorConstruction.hh"
#include "RunAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "SteppingAction.hh"
#include "TrackingAction.hh"

int main(int argc, char** argv)
{
	// Initialize detector and run parameters
	gPars::InitGlobals();
	// Choose the Random engine
	CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
	// Seed the random number generator manually
	G4long myseed = 43;
	CLHEP::HepRandom::setTheSeed(myseed);

	// Construct the default run manager
	G4RunManager * runManager = new G4RunManager;
	// Set mandatory initialization classes
	DetectorConstruction *detector = new DetectorConstruction();
	runManager->SetUserInitialization(detector);
	runManager->SetUserInitialization(new PhysicsList);
	// Set user action classes
	runManager->SetUserAction(new RunAction);
	runManager->SetUserAction(new PrimaryGeneratorAction(2.0)); //gPars::source.energy_spectrum
	runManager->SetUserAction(new SteppingAction(detector));
	runManager->SetUserAction(new TrackingAction);

	// Initialize G4 kernel
	runManager->Initialize();
	// Visualization manager
#ifdef G4VIS_USE
	G4VisManager* visManager = new G4VisExecutive;
	visManager->SetVerboseLevel(0);
	visManager->Initialize();
#endif

	// Get the pointer to the User Interface manager
	G4UImanager* UI = G4UImanager::GetUIpointer();
	UI->SetVerboseLevel(0);
	UI->ApplyCommand("/run/verbose 0");
	UI->ApplyCommand("/event/verbose 0");
	UI->ApplyCommand("/tracking/verbose 0");
	if (gPars::doView)
		UI->ApplyCommand("/control/execute vis.mac");

	runManager->BeamOn(gPars::source.N_events);

	if (gPars::doView)
		UI->ApplyCommand("vis/viewer/update");

	// Termination
#ifdef G4VIS_USE
	delete visManager;
#endif
	delete runManager;
	return 0;
}
