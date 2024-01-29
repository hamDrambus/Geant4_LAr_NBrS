#ifndef Detector_full_y2024_NIR_hh
#define Detector_full_y2024_NIR_hh

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
#include <G4OpticalSurface.hh>
#include <G4UserLimits.hh>
#include <G4VisAttributes.hh>
#include <G4Colour.hh>
#include <G4ios.hh>
#include <G4GeneralParticleSource.hh>
#include <G4RegionStore.hh>
#include <G4Trd.hh>
#include <G4PhysicalVolumesSearchScene.hh>

#include "GlobalParameters.hh"
#include <geant4/DetectorSensor.hh>
#include "DetectorParameterisation.hh"
#include "Detector_full_y2022.hh"

// Model of the detector constructed in Jan 2024 for measurements with 715 nm NIR filter before SiPMs
// 75% Electroconnect's THGEM were used for cathode and THGEM0,
// Polish standard thin GEM (PL 2022 #1) was used as GEM1.
// No WLS.
class Detector_full_y2024_NIR : public Detector_full_y2022
{
public:
	Detector_full_y2024_NIR();
	~Detector_full_y2024_NIR() override;

public:
	virtual G4VPhysicalVolume* Construct() override;
	// virtual void ConstructSDandField() override;

protected:
	virtual void SetSizeAndPosition() override;
	// virtual void CreateTHGEM0Cell(G4Material* medium);
	// Finds THGEM1 and THGEM0 volumes and initializes mapping class gPars::THGEM1_mapping
	// virtual void SetupTHGEMsMapping() override;

	// These parameters of parent class are ignored in this one
	// PMMA plate, before SiPMs
	// double PMMA_plate_size_xy;
	// double PMMA_plate_size_z;
	// G4ThreeVector position_PMMA_plate;

	// FSQ RG715 filter plate (frame, support)
	double filter_frame_thickness;
	double filter_frame_size_xy; // Full size. Has hole of filter_size_xy x filter_size_xy size
	double filter_frame_z_bottom;
	// FSQ RG715 filter itself
	double filter_size_xy;
	double filter_size_z;
	G4ThreeVector position_filter;
	G4ThreeVector position_filter_frame;
	// If PMMA thickness differs from thickness of filter support frame, then everything above
	// (i.e. SiPMs) must be appropriately shifted.
	double SiPMs_z_offset;
};

#endif // Detector_full_y2024_NIR_hh
