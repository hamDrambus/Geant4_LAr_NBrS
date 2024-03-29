#ifndef Detector_full_h
#define Detector_full_h

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

// Model of the detector during 2021.11.11 and 2022.02.03 measurements.
class Detector_full : public VDetectorConstruction
{
public:

  Detector_full();
	~Detector_full() override;

public:
	virtual G4VPhysicalVolume* Construct() override;
	virtual void ConstructSDandField() override;

protected:
	virtual void SetSizeAndPosition() override;

	G4LogicalVolume*    logic_SiPM;
	G4LogicalVolume*    logic_PMT;
	G4LogicalVolume*    logic_anode_grid;

	//anode wire
	double radius_wire;
	double length_wire;
	double step_wire;
	int N_wire;

	//Interface wire
	double radius_Interface_wire;
	double length_Interface_wire;
	double step_Interface_wire;
	int N_Interface_wire;

	//Anode_grid
	double thickness_anode_grid;
	double size_anode_grid;
	double size_anode_grid_hole;
	double z_anode_grid_bottom;

	//PMTGridWire
	double PMTGridWireRadius;
	double PMTGridWirePitch;

	//PMTAnodeGridTracker
	double PMTAnodeGridTrackerThickness;
	double PMTAnodeGridTrackerGasXSize;
	double PMTAnodeGridTrackerGasYSize;
	double PMTAnodeGridTrackerZbottom;
	int PMTAnodeGridNCellsGas;
	int PMTAnodeGridNCellsGasInner;
	int PMTAnodeGridNCellsLiquid;
	int PMTAnodeGridNCellsLiquidInner;
	double PMTAnodeGridTrackerLiquidXSize;
	double PMTAnodeGridTrackerLiquidYSize;

	//PMMA plate
	double x_size_PMMA_plate;
	double y_size_PMMA_plate;
	double z_size_PMMA_plate;

	//SiPMs
	int Nx_SiPMs;
	int Ny_SiPMs;
	double thickness_SiPM;
	double SiPM_size;
	double SiPM_pitch;
	double z_size_SiPMFR4;

	//tracker SiPM
	double x_size_tracker;
	double y_size_tracker;
	double z_size_tracker;


	//tracker Anode_grid
	double x_size_tracker_anode_grid;
	double y_size_tracker_anode_grid;
	double z_size_tracker_anode_grid;

	//tracker THGEM2 (active region with holes)
	double x_size_tracker_THGEM2;
	double y_size_tracker_THGEM2;
	double z_size_tracker_THGEM2;

	//Interface_grid
	double x_size_tracker_Interface_grid;
	double y_size_tracker_Interface_grid;
	//double z_size_tracker_Interface_grid;
	double x_size_Interface_grid_substrate;
	double y_size_Interface_grid_substrate;
	double z_size_Interface_grid_substrate;

	//THGEM1
	double x_size_THGEM1; //full real size, including dielectric
	double y_size_THGEM1;
	double z_size_THGEM1;

	//SteelBox
	double xSizeSteelBox;
	double ySizeSteelBox;
	double zSizeSteelBox;

	//THGEM_without_holes
	double x_size_THGEM_without_holes;
	double y_size_THGEM_without_holes;
	double z_size_THGEM_without_holes;

	//Insulator_box
	double x_size_Insulator_box_inner;
	double y_size_Insulator_box_inner;
	double thickness_Insulator_box;
	double x_size_Insulator_box_outer;
	double y_size_Insulator_box_outer;
	double z_size_Insulator_box;
	double z_Insulator_box_center;

	//PMTs
	double radius_PMT;
	double z_size_PMT;
	double x_pos_PMT;
	double y_pos_PMT;
	double z_pos_PMT;

	//LAr_inner
	double x_size_LAr_inner;
	double y_size_LAr_inner;
	double z_size_LAr_inner;

	//LArOuter
	double x_size_LAr_outer_in;
	double y_size_LAr_outer_in;
	double x_size_LAr_outer_out;
	double y_size_LAr_outer_out;
	double z_size_LAr_outer;

	//FieldTHGEM
	double x_size_FieldTHGEM;
	double y_size_FieldTHGEM;
	double z_size_FieldTHGEM;
	double z_center_FieldTHGEM_1;
	double z_center_FieldTHGEM_2;
	double hole_size_FieldTHGEM;

	//FieldWire
	double radius_FieldWire;
	double length_FieldWire;
	double x_pos_FieldWire;
	double z_pos_FieldWire_bottom;
	double z_pos_FieldWire_top;

	//Cathode
	double x_size_Cathode;
	double y_size_Cathode;
	double z_size_Cathode;

	//LArInactive
	double x_size_LArInactive;
	double y_size_LArInactive;
	double z_size_LArInactive;

	//TPB
	double radiusTPB;
	double z_size_TPB;

	//PMMA_bottom
	double x_size_PMMA_bottom;
	double y_size_PMMA_bottom;
	double z_size_PMMA_bottom;
	double PMMA_bottom_center;

	//Al_window
	double diameter_size_Al_window;
	double z_size_Al_window;
	double z_space_Al_window;
	double Al_window_top_center;
	double Al_window_bottom_center;

	//CryogenicChamberBottom
	double diameter_size_internal_CryogenicChamberBottom;
	double diameter_size_external_CryogenicChamberBottom;
	double z_size_CryogenicChamberBottom;
	double CryogenicChamberBottom_center;

	//ExternalColl
	double diameter_ExternalColl;
	double z_size_ExternalColl;
	double ExternalColl_center;
	double ExternalColl_bottom;

	//alpha
	double radiusAlphaFull;
	double z_size_Alpha;

	G4ThreeVector position_SiPM_container;
	G4ThreeVector position_SiPMFR4;
	G4ThreeVector position_PMMA_plate;
	G4ThreeVector position_interface_wire_container;
	G4ThreeVector position_interface_frame;
	G4ThreeVector position_Insulator_box;
	G4ThreeVector position_LAr_inner;
	G4ThreeVector position_LAr_outer;
	G4ThreeVector position_THGEM1_frame;
	G4ThreeVector position_THGEM1_copper_plate;  // inside container

	//
	G4ThreeVector positionSteelBox0;
	G4ThreeVector positionSteelBox1;
	G4ThreeVector positionSteelBox2;
	G4ThreeVector positionSteelBox3;

	//FieldWires
	G4ThreeVector position_FieldWire_bottom1;
	G4ThreeVector position_FieldWire_bottom2;
	G4ThreeVector position_FieldWire_bottom3;
	G4ThreeVector position_FieldWire_bottom4;
	G4ThreeVector position_FieldWire_top1;
	G4ThreeVector position_FieldWire_top2;
	G4ThreeVector position_FieldWire_top3;
	G4ThreeVector position_FieldWire_top4;
	G4ThreeVector position_Cathode;
	G4ThreeVector position_LArInactive;
	G4ThreeVector position_PMMA_bottom;
	G4ThreeVector position_Al_window_top;
	G4ThreeVector position_Al_window_bottom;
	G4ThreeVector position_CryogenicChamberBottom;
	G4ThreeVector position_ExternalColl;

	//PMT
	G4ThreeVector position_PMT_0;
	G4ThreeVector position_PMT_1;
	G4ThreeVector position_PMT_2;
	G4ThreeVector position_PMT_3;

	//anode_grid
	G4ThreeVector position_anode_grid;

	//PMTAnodeGridTracker
	G4ThreeVector position_PMTAnodeGridTrackerGas_1;
	G4ThreeVector position_PMTAnodeGridTrackerGasInner_1;
	G4ThreeVector position_PMTAnodeGridTrackerLiquid_1;
	G4ThreeVector position_PMTAnodeGridTrackerLiquidInner_1;
};

#endif
