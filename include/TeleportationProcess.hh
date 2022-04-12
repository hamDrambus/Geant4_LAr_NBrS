#ifndef TELEPORTATION_PROCESS_H_
#define TELEPORTATION_PROCESS_H_

#include "G4VDiscreteProcess.hh"
#include "G4ParallelWorldProcess.hh"
#include "G4TransportationManager.hh"

#include "HexagonalMapping.hh"
#include "GlobalParameters.hh"
#include "TrackUserInfo.hh"
#include "GlobalData.hh"

enum TeleportProcessSubType
{
  fTeleportTHGEM = 100
};

// This process teleports particle according to HexagonalMapping.
// Simply changing particle position (e.g. in SteppingAction) confuses geant4 and G4Navigator
// so instead particle is killed and secondary one is created where necessary.
// Code snippets are taken from G4OpBoundaryProcess (framework) and G4OpWLS (secondary particle creation)
// TODO*: implement this process in general terms, i.e. via list of abstract mapping classes instead
// of concrete HexagonalMapping with fixed hard-coded trigger conditions and used volumes.
class TeleportationProcess : public G4VDiscreteProcess
{
 public:
  explicit TeleportationProcess(const G4String& processName = "TeleportBoundary",
                               G4ProcessType type = fUserDefined);
  virtual ~TeleportationProcess() override;

  virtual G4bool IsApplicable(const G4ParticleDefinition& aParticleType) override;

  virtual G4double GetMeanFreePath(const G4Track&, G4double, G4ForceCondition* condition) override;
  // Returns infinity; i. e. the process does not limit the step, but sets the
  // 'Forced' condition for the DoIt to be invoked at every step. However, only
  // at a boundary will any action be taken.

  G4VParticleChange* PostStepDoIt(const G4Track& aTrack, const G4Step& aStep) override;
  // This is the method implementing boundary processes.

  virtual void Initialise();

  void SetVerboseLevel(G4int verbosity);

 private:
  TeleportationProcess(const TeleportationProcess& right) = delete;
  TeleportationProcess& operator=(const TeleportationProcess& right) = delete;

  void TeleportationProcessVerbose(void) const;

  HexagonalMappingData fOldState;
  HexagonalMappingData fNewState;
  double fTolerance;
};

#endif // TELEPORTATION_PROCESS_H_
