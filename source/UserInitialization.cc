#include "UserInitialization.hh"

UserInitialization::UserInitialization()
{}

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
  SetUserAction(new GenPhotonsDirectly(3.1, GenPhotonsDirectly::PatternPhoton::Cathode_14mm_coll));
  //SetUserAction(new PrimaryGeneratorAction(2.0));
  SetUserAction(new SteppingAction);
  SetUserAction(new TrackingAction);
}

UserWorkerThread::UserWorkerThread()
{}

G4WorkerRunManager* UserWorkerThread::CreateWorkerRunManager() const
{
  return new CleanerWorkerRunManager;
}
