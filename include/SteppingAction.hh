#ifndef SteppingAction_h
#define SteppingAction_h

#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4SteppingManager.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4ProcessManager.hh"
#include "G4SDManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4SystemOfUnits.hh"
#include "G4UserSteppingAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4TransportationManager.hh"

#include "GlobalParameters.hh"
#include "PhotonHit.hh"
#include "DetectorSensor.hh"
#include "DetectorConstruction.hh"

class SteppingAction : public G4UserSteppingAction
{
public:
  SteppingAction(DetectorConstruction* myDC);
  virtual void UserSteppingAction(const G4Step*);

private:
  DetectorConstruction* myDetector;
  void PassThroughGEM(const G4Step*, G4double z_pos, G4double z_size);
  void ReturnFromSingleTHGEMHole(const G4Step* theStep);
  double delta_z;
  bool is_photon_inside_single_THGEM_hole;
  G4String Volume_init_THGEM;
  G4ThreeVector MomentumDirection_init_THGEM;
  G4ThreeVector Position_init_THGEM;
  G4ThreeVector Position_init_SingleTHGEMHole;
  G4ThreeVector Position_final_SingleTHGEMHole;
};

#endif
