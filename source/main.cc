#include <time.h>
#include <iostream>
#include <G4RunManager.hh>
#include <G4MTRunManager.hh>
#include <G4UImanager.hh>
#include <G4UIterminal.hh>

#ifdef G4VIS_USE
#include <G4VisExecutive.hh>
#endif

#include <GlobalParameters.hh>
#include <GlobalData.hh>
#include <Randomize.hh>
#include <geant4/UserInitialization.hh>

int main(int argc, char** argv)
{
  std::string settings_fname;
  if (argc != 1) {
    if (argc>2)
      std::cout << "Warning! Only single parameter (settings filename) is used." << std::endl;
    settings_fname = argv[1];
  } else {
    std::cout << "Using default settings file \"settings.xml\"" << std::endl;
    settings_fname = "settings.xml";
  }
	// Initialize detector and run parameters
	if (!gPars::InitGlobals(settings_fname)) {
	  std::cerr<<"Failed to initialize globals."<<std::endl;
	  return -1;
	}
	// Choose the Random engine
	G4Random::setTheEngine(new CLHEP::RanecuEngine);
	// Seed the random number generator manually
	G4Random::setTheSeed(gPars::general.initial_seed);

	// Construct the default run manager
	G4MTRunManager *runManager = new G4MTRunManager;
	// Changed standard approach a little so that everything for
	// run manager is set up either in UserInitialization constructor or its overloaded methods.
	UserInitialization* userInitialization = new UserInitialization(runManager);

	// Initialize G4 kernel
	runManager->Initialize();
	gData.Initialize();
	// Visualization manager

  G4VisManager* visManager = nullptr;
	if (gPars::general.doView) {
    visManager = new G4VisExecutive;
    visManager->SetVerboseLevel(0);
    visManager->Initialize();
	}

	// Get the pointer to the User Interface manager
	G4UImanager* UI = G4UImanager::GetUIpointer();
	UI->SetVerboseLevel(0);
	UI->ApplyCommand("/run/verbose 0");
	UI->ApplyCommand("/event/verbose 0");
	UI->ApplyCommand("/tracking/verbose 0");
	UI->ApplyCommand("/tracking/storeTrajectory 2");
	if (gPars::general.doView) {
		UI->ApplyCommand("/control/execute vis.mac");
	}

	double THGEM1_z = gPars::det_dims->THGEM1_center.z() - gPars::det_dims->THGEM1_width_total / 2.0;
	//gData.PlotField("v00.01.01/center_axis_field.txt", G4ThreeVector(0, 0, THGEM1_z * mm - 3 * mm), G4ThreeVector(0, 0, THGEM1_z * mm + 3 * mm), 3000);
	//gData.PlotField("v00.01.01/x_p0.05_axis_field.txt", G4ThreeVector(0.05 * mm, 0, THGEM1_z * mm - 3 * mm), G4ThreeVector(0.05 * mm, 0, THGEM1_z * mm + 3 * mm), 3000);
	//gData.PlotField("v00.01.01/y_p0.05_axis_field.txt", G4ThreeVector(0, 0.05 * mm, THGEM1_z * mm - 3 * mm), G4ThreeVector(0, 0.05 * mm, THGEM1_z * mm + 3 * mm), 3000);
	runManager->BeamOn(gPars::source->N_events);

	if (gPars::general.doView) {
		UI->ApplyCommand("vis/viewer/update");
		rename_file("g4_00.wrl", gPars::general.output_folder + "g4_00.wrl");
	}

	if (visManager)
	  delete visManager;
	delete runManager;
	return 0;
}
