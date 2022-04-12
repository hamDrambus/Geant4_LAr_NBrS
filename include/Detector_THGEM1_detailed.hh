/*  Construct only part of THGEM1 active area with both cell parameterisation
 *  and mapping in order to test and view electron drift. If this class used is G4RunManager
 *  some of global parameters are overridden so that they are valid for this non-standard geometry
 */

#ifndef Detector_THGEM1_detailed_h
#define Detector_THGEM1_detailed_h

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
#include <G4AssemblyVolume.hh>

#include "DetectorParameterisation.hh"
#include "GlobalParameters.hh"
#include "HexagonalMapping.hh"
#include "VDetectorConstruction.hh"

class Detector_THGEM1_detailed : public VDetectorConstruction
{
public:
  Detector_THGEM1_detailed();
	~Detector_THGEM1_detailed();
public:
	virtual G4VPhysicalVolume* Construct() override;

private:
	virtual void SetSizeAndPosition() override;
};

#endif // Detector_THGEM1_detailed_h
