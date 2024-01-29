/*  Construct only part of THGEM1 active area with mapping and large sensor (SiPM) above
 *  in order to test that THGEM1 shading at different angles is consistent with theory.
 *  All surfaces are absorbers.
 *
 *  If this class used is G4RunManager some of global parameters are overridden so
 *  that they are valid for this non-standard geometry
 */

#ifndef Detector_THGEM1_SiPM_shading_h
#define Detector_THGEM1_SiPM_shading_h

#include <string>
#include <sstream>
#include <limits>
#include <globals.hh>

#include <G4SystemOfUnits.hh>
#include <G4ThreeVector.hh>

#include <G4Material.hh>
#include <G4Box.hh>
#include <G4Trap.hh>
#include <G4Tubs.hh>
#include <G4Cons.hh>
#include <G4Ellipsoid.hh>
#include <G4SubtractionSolid.hh>
#include <G4UnionSolid.hh>
#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4PVParameterised.hh>
#include <G4SDManager.hh>
#include <G4GeometryTolerance.hh>
#include <G4GeometryManager.hh>
#include <G4NistManager.hh>
#include <G4RunManager.hh>
#include <G4PhysicalVolumeStore.hh>
#include <G4LogicalVolumeStore.hh>
#include <G4SolidStore.hh>
#include <G4LogicalSkinSurface.hh>
#include <G4LogicalBorderSurface.hh>
#include <G4SurfaceProperty.hh>
#include <G4UImanager.hh>
#include <G4OpticalSurface.hh>
#include <G4UserLimits.hh>
#include <G4VisAttributes.hh>
#include <G4Colour.hh>
#include <G4ios.hh>
#include <G4GeneralParticleSource.hh>
#include <G4RegionStore.hh>
#include <G4Trd.hh>
#include <G4PhysicalVolumesSearchScene.hh>

#include <geant4/DetectorSensor.hh>
#include "DetectorParameterisation.hh"
#include "VDetectorConstruction.hh"

class Detector_THGEM1_SiPM_shading : public VDetectorConstruction
{
public:
  Detector_THGEM1_SiPM_shading();
	~Detector_THGEM1_SiPM_shading() override;

	virtual G4VPhysicalVolume* Construct() override;
	virtual void ConstructSDandField() override;

protected:
	virtual void SetSizeAndPosition() override;

	G4LogicalVolume*   logic_sensor;

  //Sensor (SiPM)
  double x_size_sensor;
  double y_size_sensor;
  double z_size_sensor;
  G4ThreeVector position_sensor;
};

#endif // Detector_THGEM1_SiPM_shading_h
