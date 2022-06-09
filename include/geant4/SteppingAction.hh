#ifndef SteppingAction_h
#define SteppingAction_h


#include <G4Material.hh>
#include <G4LogicalVolume.hh>
#include <G4ThreeVector.hh>
#include <G4SteppingManager.hh>
#include <G4Step.hh>
#include <G4StepPoint.hh>
#include <G4ProcessManager.hh>
#include <G4SDManager.hh>
#include <G4ParticleDefinition.hh>
#include <G4ParticleTypes.hh>
#include <G4OpBoundaryProcess.hh>
#include <G4SystemOfUnits.hh>
#include <G4UserSteppingAction.hh>
#include <globals.hh>
#include <G4ThreeVector.hh>
#include <G4TransportationManager.hh>

#include <GlobalParameters.hh>
#include "PhotonHit.hh"
#include "DetectorSensor.hh"

class SteppingAction : public G4UserSteppingAction
{
public:
  SteppingAction();
  virtual void UserSteppingAction(const G4Step*);
  bool CheckStuckParticle(const G4Step*);
protected:
  const int N_iterations;
  const double epsilon;
  std::size_t stuck_N_iterations;
  std::deque<G4ThreeVector> positions;
  std::deque<G4ThreeVector> momentums;
};

#endif
