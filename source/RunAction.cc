#include "RunAction.hh"

RunAction::RunAction()
{}

RunAction::~RunAction()
{}

G4Run* RunAction::GenerateRun()
{
  return new Run;
}

void RunAction::BeginOfRunAction(const G4Run* aRun)
{
  const Run* run = static_cast<const Run*>(aRun);
  G4cout << "### Run " << run->GetRunID() << " started." << G4endl;
  // Inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);
}

void RunAction::EndOfRunAction(const G4Run* aRun)
{
  const Run* run = static_cast<const Run*>(aRun);
  if (IsMaster()) {
    gData.progress_bar.progress_bar.mark_as_completed();
    run->results.Print(G4cout);
    run->AddToFile(run->results.generated_photons, gPars::results.generated_filename);
    run->AddToFile(run->results.recorded_photons, gPars::results.recorded_filename);
    G4cout << G4endl
    << "--------------------End of Global Run-----------------------"
    << "  #" << run->GetRunID() << " Event#: "<< run->GetNumberOfEvent() << G4endl;
  }
}
