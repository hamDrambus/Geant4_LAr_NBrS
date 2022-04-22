#include <geant4/KillProcess.hh>

KillProcess::KillProcess(const G4String& processName, G4ProcessType type)
  : G4VDiscreteProcess(processName, type)
{}

KillProcess::~KillProcess() {}

G4VParticleChange* KillProcess::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep)
{
  aParticleChange.Initialize(aTrack);
  if(aTrack.GetDynamicParticle()->GetDefinition() == G4OpticalPhoton::OpticalPhoton()) {
    if(aTrack.GetLocalTime() > gPars::general.photon_max_time) {
      aParticleChange.ProposeTrackStatus(fStopAndKill);
      aParticleChange.ProposeLocalEnergyDeposit(0.0);
      G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
    }
  }
  return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
}

