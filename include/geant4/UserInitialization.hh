#ifndef UserInitialization_h
#define UserInitialization_h

#include <G4VUserActionInitialization.hh>
#include <G4UserWorkerInitialization.hh>
#include <G4UserWorkerThreadInitialization.hh>
#include <Randomize.hh>
#include <G4MTRunManager.hh>

#include <GlobalParameters.hh>
#include <GlobalData.hh>
#include "RunAction.hh"
#include "SteppingAction.hh"
#include "TrackingAction.hh"
#include "CleanerWorkerRunManager.hh"
#include "PhysicsList.hh"
#include "detector/Detector_full_y2022.hh"
#include "detector/Detector_full_y2024_NIR.hh"
#include "detector/Detector_full.hh"
#include "detector/Detector_THGEM1_detailed.hh"
#include "detector/Detector_THGEM1_SiPM_shading.hh"
#include "generator/GenNBrS_InTHGEM.hh"
#include "generator/GenElectronsPatterns.hh"
#include "generator/GenElectronsPatternsNBrS.hh"
#include "generator/GenPhotonsDirectly.hh"

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
