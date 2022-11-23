#ifndef Detector_full_y2022_h
#define Detector_full_y2022_h

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
#include "VDetectorConstruction.hh"

// Model of the detector during fall of 2022 measurements (in preparation as of 2022.10.12).
// Separate from the previous model (Detector_full), only materials and surfaces from VDetectorConstruction are shared.
class Detector_full_y2022 : public VDetectorConstruction
{
public:

	Detector_full_y2022();
	~Detector_full_y2022() override;

public:
	virtual G4VPhysicalVolume* Construct() override;
	virtual void ConstructSDandField() override;

private:
	virtual void SetSizeAndPosition() override;
	virtual void CreateTHGEM0Cell();
	virtual void SetupTHGEMsMapping() override; //Finds THGEM1 and THGEM0 volumes and initializes mapping class gPars::THGEM1_mapping

	G4LogicalVolume*   logic_SiPM;
	G4LogicalVolume*   logic_PMT;

	G4LogicalVolume*   logic_THGEM0_cell;
	G4LogicalVolume*   logic_THGEM0_cell_LAr;
	G4LogicalVolume*   logic_THGEM0_cell_FR4;
	G4LogicalVolume*   logic_THGEM0_cell_copper;
	G4LogicalVolume*   logic_THGEM0_container;

	G4VPhysicalVolume* phys_THGEM0_cell_LAr;
	G4VPhysicalVolume* phys_THGEM0_container;

	// SiPMs
	int Nx_SiPMs;
	int Ny_SiPMs;
	double SiPM_thickness;
	double SiPM_size;
	double SiPM_pitch;
	// SiPM container
	double SiPM_cont_size_xy;
	double SiPM_cont_size_z;
	double SiPM_FR4_size_z;
	G4ThreeVector position_SiPM_container;
	G4ThreeVector position_SiPMFR4;
	G4ThreeVector offset_SiPM_in_container; // So that top of SiPMs coincides with FR4 substrate's bottom

	// PMMA plate, before SiPMs
	double PMMA_plate_size_xy;
	double PMMA_plate_size_z;
	G4ThreeVector position_PMMA_plate;

	// Anode wires (ground before SiPMs)
	double anode_wire_radius;
	double anode_wire_length;
	double anode_wire_step;
	int anode_wire_N;
	// Anode grid
	double anode_grid_thickness;
	double anode_grid_size_xy; // Full size. Has hole of length_wire x length_wire size
	double anode_grid_z_bottom;
	// Anode grid's container
	double anode_grid_cont_size_xy;
	double anode_grid_cont_size_z;
	G4ThreeVector position_anode_grid_container;
	G4ThreeVector position_anode_grid_frame;

	// THGEM1 (above EL gap and below Anode wires)
	double THGEM1_size_xy; // Full real size, including dielectric, z size is set in DetectorDimsFullY2022
	double THGEM1_active_size_xy;
	G4ThreeVector position_THGEM1_frame;
	G4ThreeVector position_THGEM1_copper_plate;  // dummy inside container in case mapping is avoided

	// LAr inner (inside the insulator box (acrylic), stretches from EL gap to the cathode)
	double LAr_inner_size_xy;
	double LAr_inner_size_z;
	G4ThreeVector position_LAr_inner;

	// LAr outer (outside the insulator box (acrylic), stretches from EL gap to the cathode)
	double LAr_outer_in_size_xy; // Limited by the insulator box
	double LAr_outer_out_size_xy; // Limited by the steel box
	double LAr_outer_size_z;
	G4ThreeVector position_LAr_outer;

	// Interface THGEM (THGEM0)
	double THGEM0_size_xy; // Full real size, including dielectric, z size is set in DetectorDimsFullY2022
	double THGEM0_active_size_xy;
	G4ThreeVector position_THGEM0_frame;
	G4ThreeVector position_THGEM0_copper_plate;  // dummy inside container in case mapping is avoided
	G4ThreeVector position_THGEM0_container;
	G4ThreeVector position_SingleTHGEM0Cell;

	// Field-forming wires
	double FieldWire_radius;
	double FieldWire_length;
	double FieldWire_pos_x;
	double FieldWire_top_pos_z;
	double FieldWire_bot_pos_z;
	G4ThreeVector position_FieldWire_bottom1;
	G4ThreeVector position_FieldWire_bottom2;
	G4ThreeVector position_FieldWire_bottom3;
	G4ThreeVector position_FieldWire_bottom4;
	G4ThreeVector position_FieldWire_top1;
	G4ThreeVector position_FieldWire_top2;
	G4ThreeVector position_FieldWire_top3;
	G4ThreeVector position_FieldWire_top4;

	// Cathode
	double Cathode_size_xy;
	double Cathode_size_z;
	G4ThreeVector position_Cathode;

	// PMTs
	double PMT_radius;
	double PMT_size_z;
	G4ThreeVector position_PMT_0;
	G4ThreeVector position_PMT_1;
	G4ThreeVector position_PMT_2;
	G4ThreeVector position_PMT_3;

	// Steel box (made of 4 plates with hole for PMT)
	double steel_box_size_x;
	double steel_box_size_y;
	double steel_box_size_z;
	G4ThreeVector position_steel_box_0;
	G4ThreeVector position_steel_box_1;
	G4ThreeVector position_steel_box_2;
	G4ThreeVector position_steel_box_3;

	// PMT anode grid wire
	double PMT_grid_wire_radius;
	double PMT_grid_wire_pitch;
	// PMT anode grid's container
	double PMT_grid_cont_thickness;
	double PMT_grid_cont_gas_size_x;
	double PMT_grid_cont_gas_size_y;
	double PMT_grid_cont_liquid_size_x;
	double PMT_grid_cont_liquid_size_y;
	int PMT_grid_gas_Ncells;
	int PMT_grid_gas_inner_Ncells;
	int PMT_grid_liquid_Ncells;
	int PMT_grid_liquid_inner_Ncells;
	G4ThreeVector position_PMT_grid_cont_gas_1;
	G4ThreeVector position_PMT_grid_cont_gas_inner_1;
	G4ThreeVector position_PMT_grid_cont_liquid_1;
	G4ThreeVector position_PMT_grid_cont_liquid_inner_1;

	// Insulator box (acrylic)
	double insulator_box_inner_size_xy;
	double insulator_box_thickness;
	double insulator_box_outer_size_xy;
	double insulator_box_size_z;
	G4ThreeVector position_insulator_box;

	// TPB (WLS) on the inside of the insulator box.
	// In order to avoid splitting TPB cylinder into liquid and gas part, TPB is inserted into the insulator box
	double TPB_radius;
	double TPB_thickness;
	G4ThreeVector position_TPB_0;
	G4ThreeVector position_TPB_1;
	G4ThreeVector position_TPB_2;
	G4ThreeVector position_TPB_3;

};

#endif
