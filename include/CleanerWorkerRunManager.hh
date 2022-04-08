#ifndef CleanerWorkerRunManager_h
#define CleanerWorkerRunManager_h

#include <G4WorkerRunManager.hh>
#include <G4MTRunManager.hh>
#include <G4ScoringManager.hh>

#include "Run.hh"

class CleanerWorkerRunManager : public G4WorkerRunManager
{
public:
  CleanerWorkerRunManager();
 ~CleanerWorkerRunManager();

protected:
  virtual void MergePartialResults() override; // Calls Run->Merged() function to notify it that the local data may be freed
};

#endif  // CleanerWorkerRunManager_h
