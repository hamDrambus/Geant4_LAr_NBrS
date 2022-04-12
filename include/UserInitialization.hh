#ifndef UserInitialization_h
#define UserInitialization_h

#include <G4VUserActionInitialization.hh>
#include <G4UserWorkerInitialization.hh>
#include <G4UserWorkerThreadInitialization.hh>
#include <Randomize.hh>
#include <G4MTRunManager.hh>

#include "GlobalParameters.hh"
#include "GlobalData.hh"
#include "GlobalParameters.hh"
#include "RunAction.hh"
#include "GenNBrS_InTHGEM.hh"
#include "GenElectronsPatterns.hh"
#include "GenPhotonsDirectly.hh"
#include "SteppingAction.hh"
#include "TrackingAction.hh"
#include "CleanerWorkerRunManager.hh"
#include "Detector_full.hh"
#include "Detector_THGEM1_detailed.hh"
#include "Detector_THGEM1_SiPM_shading.hh"
#include "PhysicsList.hh"

class UserInitialization : public G4VUserActionInitialization
{
public:
  UserInitialization(G4RunManager *run_manager);
  ~UserInitialization();

  void BuildForMaster() const override;
  void Build() const override;
protected:
  G4RunManager *runManager;
};

class UserWorkerThread : public G4UserWorkerThreadInitialization
{
public:
  UserWorkerThread();
  virtual G4WorkerRunManager* CreateWorkerRunManager() const override;
};

#endif // UserInitialization_h
