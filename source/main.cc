#include <time.h>
#include <iostream>
#include <G4RunManager.hh>
#include <G4MTRunManager.hh>
#include <G4UImanager.hh>
#include <G4UIterminal.hh>

#ifdef G4VIS_USE
#include <G4VisExecutive.hh>
#endif

#include "GlobalParameters.hh"
#include "GlobalData.hh"
#include "Randomize.hh"
#include "PhysicsList.hh"
#include "DetectorConstruction.hh"
#include "DetectorConstructionTHGEM1.hh"
#include "DetectorConstructionTHGEM1Shading.hh"
#include "UserInitialization.hh"

int main(int argc, char** argv)
{
	// Initialize detector and run parameters
	gPars::InitGlobals();
	// Choose the Random engine
	G4Random::setTheEngine(new CLHEP::RanecuEngine);
	// Seed the random number generator manually
	//G4long myseed = 43;
	G4long myseed = time(NULL);
	G4Random::setTheSeed(myseed);

	// Construct the default run manager
	G4MTRunManager *runManager = new G4MTRunManager;
	runManager->SetNumberOfEventsToBeStored(0);
	runManager->SetNumberOfThreads(std::max(gPars::general.thread_number, 1));
	// Set mandatory initialization classes
	runManager->SetUserInitialization(new DetectorConstruction);
	//runManager->SetUserInitialization(new DetectorConstructionTHGEM1);
	//runManager->SetUserInitialization(new DetectorConstructionTHGEM1Shading);
	runManager->SetUserInitialization(new PhysicsList);
	runManager->SetUserInitialization(new UserWorkerThread); // Only responsible for clearing merged runs.
	runManager->SetUserInitialization(new UserInitialization); // Primary generator, stepping, tracking, event and run actions

	// Initialize G4 kernel
	runManager->Initialize();
	gData.Initialize();
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
	UI->ApplyCommand("/tracking/storeTrajectory 1");
	if (gPars::general.doView)
		UI->ApplyCommand("/control/execute vis.mac");

	double THGEM1_z = gPars::det_dims.z_bottom_THGEM1 + gPars::det_dims.THGEM1_width_total / 2.0;
  //gData.PlotField("v00.01.01/center_axis_field.txt", G4ThreeVector(0, 0, THGEM1_z * mm - 3 * mm), G4ThreeVector(0, 0, THGEM1_z * mm + 3 * mm), 3000);
  //gData.PlotField("v00.01.01/x_p0.05_axis_field.txt", G4ThreeVector(0.05 * mm, 0, THGEM1_z * mm - 3 * mm), G4ThreeVector(0.05 * mm, 0, THGEM1_z * mm + 3 * mm), 3000);
  //gData.PlotField("v00.01.01/y_p0.05_axis_field.txt", G4ThreeVector(0, 0.05 * mm, THGEM1_z * mm - 3 * mm), G4ThreeVector(0, 0.05 * mm, THGEM1_z * mm + 3 * mm), 3000);
	runManager->BeamOn(gPars::source.N_events);

	if (gPars::general.doView)
		UI->ApplyCommand("vis/viewer/update");

	// Termination
#ifdef G4VIS_USE
	delete visManager;
#endif
	delete runManager;
	return 0;
}
