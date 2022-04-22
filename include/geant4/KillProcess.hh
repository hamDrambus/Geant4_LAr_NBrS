#ifndef KILL_PROCESS_H_
#define KILL_PROCESS_H_

#include <G4VDiscreteProcess.hh>
#include <G4ParallelWorldProcess.hh>
#include <G4TransportationManager.hh>
#include <G4Track.hh>
#include <G4OpticalPhoton.hh>

#include <GlobalParameters.hh>

//  This process is required only to kill stuck photons.
class KillProcess : public G4VDiscreteProcess
{
 public:
  explicit KillProcess(const G4String& processName = "KillProcess",
                               G4ProcessType type = fUserDefined);
  virtual ~KillProcess();

  virtual G4bool IsApplicable(const G4ParticleDefinition& aParticleType) override;

  virtual G4double GetMeanFreePath(const G4Track&, G4double, G4ForceCondition* condition) override;
  // Returns infinity; i. e. the process does not limit the step, but sets the
  // 'Forced' condition for the DoIt to be invoked at every step. However, only
  // at a boundary will any action be taken.

  G4VParticleChange* PostStepDoIt(const G4Track& aTrack, const G4Step& aStep) override;
  // This is the method implementing boundary processes.
private:
  KillProcess(const KillProcess& right) = delete;
  KillProcess& operator=(const KillProcess& right) = delete;
};

inline G4bool KillProcess::IsApplicable(const G4ParticleDefinition& aParticleType)
{
  return true;
}

inline G4double KillProcess::GetMeanFreePath(const G4Track&, G4double, G4ForceCondition* condition)
{
  *condition = Forced;
  return DBL_MAX;
}

#endif //KILL_PROCESS_H_
