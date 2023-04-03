#include <geant4/UserInitialization.hh>

UserInitialization::UserInitialization(G4RunManager *run_manager) :
  runManager(run_manager)
{
  if (nullptr != runManager) {
    runManager->SetNumberOfEventsToBeStored(0);
    runManager->SetNumberOfThreads(std::max(gPars::general.thread_number, 1));
    // Create detector depending on settings given
    VDetectorConstruction* detector = nullptr;
    switch(gPars::det_dims->detector_type) {
    case (VDetectorDimensions::Full_detector_y2022): {
				detector = new Detector_full_y2022;
				break;
			}
    case (VDetectorDimensions::Full_detector): {
      detector = new Detector_full;
      break;
    }
    case (VDetectorDimensions::THGEM1_detailed): {
      detector = new Detector_THGEM1_detailed;
      break;
    }
    case (VDetectorDimensions::THGEM1_SiPM_shading): {
      detector = new Detector_THGEM1_SiPM_shading;
      break;
    }
    default: {
      G4Exception("UserInitialization::UserInitialization: ",
          "InvalidSetup", FatalException, "Unimplemented detector type is used.");
      return;
    }
    }
    // Set mandatory initialization classes
    runManager->SetUserInitialization(detector);
    runManager->SetUserInitialization(new PhysicsList);
    runManager->SetUserInitialization(new UserWorkerThread); // Only responsible for clearing merged runs.
    runManager->SetUserInitialization(this);
  }
}

UserInitialization::~UserInitialization()
{}

void UserInitialization::BuildForMaster() const
{
  SetUserAction(new RunAction);
}

void UserInitialization::Build() const
{
  // Create generator depending on settings given
  VGeneratePrimaries* generator = nullptr;
  switch(gPars::source->generator_type) {
  case (VSourceSettings::NBrS): {
    SettingsNBrSGenerator *settings = static_cast<SettingsNBrSGenerator*>(gPars::source);
    GenNBrS_InTHGEM *gen = new GenNBrS_InTHGEM;
    generator = gen;
    gen->SetNBrSYieldFactor(settings->NBrS_yield_factor);
    gen->SetSourceXYprofile(settings->xy_radius, settings->xy_radius_smearing);
    break;
  }
  case (VSourceSettings::PhotonsDirectly): {
    SettingsDirectPhotons *settings = static_cast<SettingsDirectPhotons*>(gPars::source);
    GenPhotonsDirectly *gen;
    if (settings->energy_spectrum.isValid())
      gen = new GenPhotonsDirectly(settings->energy_spectrum, settings->pattern, settings->angle);
    else
      gen = new GenPhotonsDirectly(settings->energy, settings->pattern, settings->angle);
    generator = gen;
    break;
  }
  case (VSourceSettings::ElectronPatterns): {
    SettingsElectronPattern *settings = static_cast<SettingsElectronPattern*>(gPars::source);
    GenElectronsPatterns *gen = new GenElectronsPatterns(settings->pattern);
    generator = gen;
    break;
  }
  case (VSourceSettings::ElectronPatternsNBrS): {
		SettingsElectronPatternNBrS *settings = static_cast<SettingsElectronPatternNBrS*>(gPars::source);
		GenElectronsPatternsNBrS *gen = new GenElectronsPatternsNBrS();
		gen->SetNBrSYieldFactor(settings->NBrS_yield_factor);
		gen->SetPattern(settings->pattern);
		generator = gen;
		break;
	}
  default: {
    G4Exception("UserInitialization::UserInitialization: ",
        "InvalidSetup", FatalException, "Unimplemented generator type is used.");
    return;
  }
  }
  SetUserAction(new RunAction);
  SetUserAction(generator);
  SetUserAction(new SteppingAction);
  SetUserAction(new TrackingAction);
}

UserWorkerThread::UserWorkerThread()
{}

G4WorkerRunManager* UserWorkerThread::CreateWorkerRunManager() const
{
  return new CleanerWorkerRunManager;
}
