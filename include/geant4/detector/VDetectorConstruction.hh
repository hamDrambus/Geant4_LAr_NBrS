/*  Virtual DetectorConstruction class with common functions and definitions
 *  for concrete geometries.
 */

#ifndef VDetectorConstruction_h
#define VDetectorConstruction_h

#include <string>
#include <sstream>
#include <limits>
#include <globals.hh>

#include <G4SystemOfUnits.hh>
#include <G4ThreeVector.hh>
#include <G4VUserDetectorConstruction.hh>
#include <G4Material.hh>
#include <G4LogicalVolume.hh>
#include <G4OpticalSurface.hh>
#include <G4VisAttributes.hh>
#include <G4Colour.hh>
#include <G4Polycone.hh>
#include <G4GenericPolycone.hh>
#include <G4RegionStore.hh>
#include <G4PVPlacement.hh>
#include <G4SubtractionSolid.hh>
#include <G4NistManager.hh>

#include <utilities/GlobalUtilities.hh>
#include "GlobalParameters.hh"
#include "GlobalData.hh"

class VDetectorConstruction : public G4VUserDetectorConstruction
{
public:

  VDetectorConstruction();
	~VDetectorConstruction() override;

public:
	virtual G4VPhysicalVolume* Construct() override = 0;

protected:
	DataVector convert_lambda_to_E(DataVector &in) const;

	virtual void defineMaterials();
	virtual void defineSurfaces();
	virtual void SetSizeAndPosition() = 0;
	virtual void CreateTHGEM1Cell();
	virtual void SetupTHGEMsMapping(); //Finds THGEM1 volumes and initializes mapping class gPars::THGEM1_mapping

	std::vector<G4PhysicalVolumesSearchScene::Findings> LocatePV(G4VPhysicalVolume* volume); //must be called after physiWorld is fully constructed

  G4RotationMatrix *rotX_90;
  G4RotationMatrix *rotY_90;
  G4RotationMatrix *rotZ_90;
  G4RotationMatrix *rotZ_180;
  G4RotationMatrix *rotZ_270;

  G4Material* matAl;
  G4Material* matFe;
	G4Material* matFR4;
  G4Material* matLAr;
  G4Material* matGas;
  G4Material* matPMMA;
  G4Material* matPMMA_UV;
  G4Material* matTPB;

  G4Box*             solidWorld;
  G4LogicalVolume*   logicWorld;
  G4VPhysicalVolume* physiWorld;
  G4LogicalVolume*   logic_THGEM1_cell;
  G4LogicalVolume*   logic_THGEM1_cell_LAr;
  G4LogicalVolume*   logic_THGEM1_cell_FR4;
  G4LogicalVolume*   logic_THGEM1_cell_copper;
  G4LogicalVolume*   logic_THGEM1_container;

  G4VPhysicalVolume* phys_THGEM1_cell_LAr;
	G4VPhysicalVolume* phys_THGEM1_container;

  // surfaces
  G4OpticalSurface *AbsorberMaterial;
  G4OpticalSurface *FR4_unified;
  G4OpticalSurface *Cu_THGEM;
  G4OpticalSurface *Cu_cathode; // Separate because it accounts for non 100% Cu coverage
  G4OpticalSurface *LAr_OpticalSurface;
  G4OpticalSurface *Anode_wire_unified;
  G4OpticalSurface *stainlessSteel;
  G4OpticalSurface *PMT_OpticalSurface;
  G4OpticalSurface *SiPM_OpticalSurface;

  G4bool  fCheckOverlaps; // option to activate checking of volumes overlaps

  G4double HalfWorldLength;
  //THGEM1
  double x_size_THGEM1_container; // area with electric field, Cu and holes
  double y_size_THGEM1_container;
  double z_size_THGEM1_container;

  G4ThreeVector position_SingleTHGEM1Cell; // in world volume
  G4ThreeVector position_THGEM1_container; // in LAr or world volume
};

#endif // VDetectorConstruction_h
