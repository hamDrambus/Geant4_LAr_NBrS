#include <geant4/CleanerWorkerRunManager.hh>

CleanerWorkerRunManager::CleanerWorkerRunManager()
{}
CleanerWorkerRunManager::~CleanerWorkerRunManager()
{}

void CleanerWorkerRunManager::MergePartialResults()
{
  G4WorkerRunManager::MergePartialResults();
  if (currentRun) {
    Run *run = dynamic_cast<Run*>(currentRun);
    if (run) {
      run->Merged();
    }
  }
}
