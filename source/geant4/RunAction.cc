#include <geant4/RunAction.hh>

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
    if (!gData.progress_bar.progress_bar.is_completed())
      gData.progress_bar.progress_bar.mark_as_completed();
    run->results.Print(G4cout);
    for (std::size_t e = 0, e_end_ = run->results.generated_photons.size(); e!=e_end_; ++e)
      run->results.generated_photons[e].electron.track.Draw();
    G4cout << G4endl
        << "-----------------Writing results to files--------------------" << G4endl;
    run->AddToFile(run->results.generated_photons, gPars::general.output_folder + gPars::results.generated_filename);
    run->AddToFile(run->results.recorded_photons, gPars::general.output_folder + gPars::results.recorded_filename);
    G4cout << G4endl
    << "--------------------End of Global Run-----------------------"
    << "  #" << run->GetRunID() << " Event#: "<< run->GetNumberOfEvent() << G4endl;
  }
}
