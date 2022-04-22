#include <geant4/UserInitialization.hh>

UserInitialization::UserInitialization(G4RunManager *run_manager) :
  runManager(run_manager)
{
  if (nullptr != runManager) {
    runManager->SetNumberOfEventsToBeStored(0);
    runManager->SetNumberOfThreads(std::max(gPars::general.thread_number, 1));
    // Set mandatory initialization classes
    runManager->SetUserInitialization(new Detector_full);
    //runManager->SetUserInitialization(new Detector_THGEM1_detailed);
    //runManager->SetUserInitialization(new Detector_THGEM1_SiPM_shading);
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
  SetUserAction(new RunAction);
  //SetUserAction(new GenElectronsPatterns(GenElectronsPatterns::PatternElectron::UniformLineX));
  //SetUserAction(new GenPhotonsDirectly(3.1, GenPhotonsDirectly::PatternPhoton::THGEM1_hole_center));
  //SetUserAction(new GenPhotonsDirectly(3.1, GenPhotonsDirectly::PatternPhoton::Cathode_14mm_coll));
  //SetUserAction(new GenPhotonsDirectly(3.1, GenPhotonsDirectly::PatternPhoton::SiPM_shading, 20.0));
  SetUserAction(new GenNBrS_InTHGEM);
  SetUserAction(new SteppingAction);
  SetUserAction(new TrackingAction);
}

UserWorkerThread::UserWorkerThread()
{}

G4WorkerRunManager* UserWorkerThread::CreateWorkerRunManager() const
{
  return new CleanerWorkerRunManager;
}
