#include "DetectorConstruction.hh"

DetectorConstruction::DetectorConstruction() : fCheckOverlaps(gPars::debugging.check_geometry_overlap)
{
	rotX_90 = new G4RotationMatrix();
	rotX_90->rotateX(90 * deg);
	rotY_90 = new G4RotationMatrix();
	rotY_90->rotateY(90 * deg);
	rotZ_90 = new G4RotationMatrix();
	rotZ_90->rotateZ(90 * deg);
	rotZ_180 = new G4RotationMatrix();
	rotZ_180->rotateZ(180 * deg);
	rotZ_270 = new G4RotationMatrix();
	rotZ_270->rotateZ(270 * deg);
}

DetectorConstruction::~DetectorConstruction()
{
	delete rotX_90;
	delete rotY_90;
	delete rotZ_90;
	delete rotZ_180;
	delete rotZ_270;
}

G4VPhysicalVolume * DetectorConstruction::Construct()
{
	defineMaterials();
	defineSurfaces();
	SetSizeAndPosition();

	const double diameter_inter_ExternalColl = gPars::det_dims.external_collimator_diameter * mm;
	//-------------------------------------------------------------------------------

	solidWorld = new G4Box("sworld", HalfWorldLength, HalfWorldLength, HalfWorldLength);
	logicWorld = new G4LogicalVolume(solidWorld, G4Material::GetMaterial("Air"), "lWorld", 0, 0, 0);
	//  Must place the World Physical volume unrotated at (0,0,0).
	physiWorld = new G4PVPlacement(0, G4ThreeVector(), logicWorld, "pWorld", 0, // its mother  volume
		false, 0);
	// Set user cuts to avoid deadlocks
	G4double maxStep = 10.0*m, maxLength = 10.0*m, maxTime = 100.0*ns, minEkin = 0.5*eV;
	logicWorld->SetUserLimits(new G4UserLimits(maxStep, maxLength, maxTime, minEkin));

	//-------------------------------------------------------------------------------
	//create LAr box contained inside the PMMA insulator
	G4Box* solid_LAr_inner = new G4Box("solid_LAr_inner", x_size_LAr_inner / 2.0, y_size_LAr_inner / 2.0, z_size_LAr_inner / 2.0);
	G4LogicalVolume* logic_LAr_inner = new G4LogicalVolume(solid_LAr_inner, G4Material::GetMaterial("LAr"), "logic_LAr_inner", 0, 0, 0);
	G4VPhysicalVolume* phys_LAr_inner = new G4PVPlacement(0, position_LAr_inner, logic_LAr_inner,  "phys_LAr_inner",
		logicWorld, false, 0, fCheckOverlaps);

	//--------------------------------------------------------------------------------
	//create LAr box contained outside the PMMA insulator
	G4Box* solid_LAr_outer_out = new G4Box("solid_LAr_outer_out", x_size_LAr_outer_out / 2.0, y_size_LAr_outer_out / 2.0, z_size_LAr_outer / 2.0);
	G4Box* solid_LAr_outer_in = new G4Box("solid_LAr_outer_in", x_size_LAr_outer_in / 2.0, y_size_LAr_outer_in / 2.0, z_size_LAr_outer / 2.0);
	G4SubtractionSolid* solid_LAr_outer = new G4SubtractionSolid("solid_LAr_outer", solid_LAr_outer_out, solid_LAr_outer_in);
	G4LogicalVolume* logic_LAr_outer = new G4LogicalVolume(solid_LAr_outer, G4Material::GetMaterial("LAr"), "logic_LAr_outer", 0, 0, 0);
	G4VPhysicalVolume* phys_LAr_outer = new G4PVPlacement(0, position_LAr_outer, logic_LAr_outer, "phys_LAr_outer",
		logicWorld, false, 0, fCheckOverlaps);

	//--------------------------------------------------------------------------------
	//create FieldWires
	G4Tubs* solid_FieldWire = new G4Tubs("solid_FieldWire", 0, radius_FieldWire, length_FieldWire / 2.0, 0.*deg, 360.*deg);
	G4LogicalVolume* logic_FieldWire = new G4LogicalVolume(solid_FieldWire, G4Material::GetMaterial("Al"), "logic_FieldWire", 0, 0, 0);

	G4VPhysicalVolume* phys_FieldWire_bottom1 = new G4PVPlacement(rotX_90,	position_FieldWire_bottom1, logic_FieldWire,
		"phys__FieldWire_bottom1", logic_LAr_inner, false, 0, fCheckOverlaps);
	G4VPhysicalVolume* phys_FieldWire_bottom2 = new G4PVPlacement(rotX_90, position_FieldWire_bottom2, logic_FieldWire,
		"phys__FieldWire_bottom2", logic_LAr_inner, false, 0, fCheckOverlaps);
	G4VPhysicalVolume* phys_FieldWire_bottom3 = new G4PVPlacement(rotY_90, position_FieldWire_bottom3, logic_FieldWire,
		"phys__FieldWire_bottom3", logic_LAr_inner, false, 0, fCheckOverlaps);
	G4VPhysicalVolume* phys_FieldWire_bottom4 = new G4PVPlacement(rotY_90, position_FieldWire_bottom4, logic_FieldWire,
		"phys__FieldWire_bottom4", logic_LAr_inner, false, 0, fCheckOverlaps);

	G4VPhysicalVolume* phys_FieldWire_top1 = new G4PVPlacement(rotX_90, position_FieldWire_top1, logic_FieldWire,
		"phys__FieldWire_top1", logic_LAr_inner, false, 0, fCheckOverlaps);
	G4VPhysicalVolume* phys_FieldWire_top2 = new G4PVPlacement(rotX_90, position_FieldWire_top2, logic_FieldWire,
		"phys__FieldWire_top2", logic_LAr_inner, false, 0, fCheckOverlaps);
	G4VPhysicalVolume* phys_FieldWire_top3 = new G4PVPlacement(rotY_90, position_FieldWire_top3, logic_FieldWire,
		"phys__FieldWire_top3", logic_LAr_inner, false, 0, fCheckOverlaps);
	G4VPhysicalVolume* phys_FieldWire_top4 = new G4PVPlacement(rotY_90, position_FieldWire_top4, logic_FieldWire,
		"phys__FieldWire_top4", logic_LAr_inner, false, 0, fCheckOverlaps);

	//-------------------------------------------------------------------------------
	CreateTHGEM1Cell();
	G4VPhysicalVolume* phys_THGEM1_cell = new G4PVPlacement(0, position_SingleTHGEMCell, logic_THGEM1_cell,
		"phys_THGEM1_cell_isolated", logicWorld, false, 0, fCheckOverlaps);

	//-------------------------------------------------------------------------------
	//create bCathode
	G4Box* solid_bCathode = new G4Box("solid_bCathode", x_size_Cathode / 2.0, y_size_Cathode / 2.0, z_size_Cathode / 2.0);
	G4LogicalVolume* logic_bCathode = new G4LogicalVolume(solid_bCathode, G4Material::GetMaterial("FR4"), "logic_bCathode", 0, 0, 0);
	G4VPhysicalVolume* phys_bCathode = new G4PVPlacement(0, position_Cathode, logic_bCathode, "phys_bCathode",
		logicWorld, false, 0, fCheckOverlaps);

	//-------------------------------------------------------------------------------
	//Inactive LAr layer below the cathode. It is considered that in reality there is gas there. No effect on light collection.
	G4Box* solid_LArInactive = new G4Box("solid_LArInactive", x_size_LArInactive / 2.0, y_size_LArInactive / 2.0, z_size_LArInactive / 2.0);
	G4LogicalVolume* logic_LArInactive = new G4LogicalVolume(solid_LArInactive, G4Material::GetMaterial(/*"LAr"*/"Air"), "logic_LArInactive", 0, 0, 0);
	G4VPhysicalVolume* phys_LArInactive = new G4PVPlacement(0, position_LArInactive, logic_LArInactive, "phys_LArInactive",
		logicWorld, false, 0, fCheckOverlaps);

	//-------------------------------------------------------------------------------
	// Create bottom PMMA plate below the cathode. No effect on light collection.
	G4Box* solid_PMMA_bottom = new G4Box("solid_PMMA_bottom", x_size_PMMA_bottom / 2.0, y_size_PMMA_bottom / 2.0, z_size_PMMA_bottom / 2.0);
	G4LogicalVolume* logic_PMMA_bottom = new G4LogicalVolume(solid_PMMA_bottom, G4Material::GetMaterial("PMMA"), "logic_PMMA_bottom", 0, 0, 0);
	G4VPhysicalVolume* phys_PMMA_bottom = new G4PVPlacement(0, position_PMMA_bottom, logic_PMMA_bottom, "phys_PMMA_bottom",
		logicWorld, false, 0, fCheckOverlaps);

	//-------------------------------------------------------------------------------
	// Create Al window. No effect on light collection.
	G4Tubs* solid_Al_window = new G4Tubs("solid_Al_window", 0, diameter_size_Al_window / 2.0, z_size_Al_window / 2.0, 0. *deg, 360.*deg);
	G4LogicalVolume* logic_Al_window = new G4LogicalVolume(solid_Al_window, G4Material::GetMaterial("Al"), "logic_Al_window", 0, 0, 0);
	G4VPhysicalVolume* phys_Al_window_top = new G4PVPlacement(0, position_Al_window_top, logic_Al_window, "phys_Al_window_top",
		logicWorld, false, 0, fCheckOverlaps);
	G4VPhysicalVolume* phys_Al_window_bottom = new G4PVPlacement(0, position_Al_window_bottom, logic_Al_window, "phys_Al_window_bottom",
		logicWorld, false, 0, fCheckOverlaps);

	//-------------------------------------------------------------------------------
	// Create chamber bottom. No effect on light collection.
	G4Tubs* solid_CryogenicChamberBottom = new G4Tubs("solid_CryogenicChamberBottom", diameter_size_internal_CryogenicChamberBottom / 2.0, diameter_size_external_CryogenicChamberBottom / 2.0, z_size_CryogenicChamberBottom / 2.0, 0. *deg, 360.*deg);
	G4LogicalVolume* logic_CryogenicChamberBottom = new G4LogicalVolume(solid_CryogenicChamberBottom, G4Material::GetMaterial("Fe"), "logic_CryogenicChamberBottom", 0, 0, 0);
	G4VPhysicalVolume* phys_CryogenicChamberBottom = new G4PVPlacement(0, position_CryogenicChamberBottom, logic_CryogenicChamberBottom, "phys_CryogenicChamberBottom",
		logicWorld, false, 0, fCheckOverlaps);

	//-------------------------------------------------------------------------------
	// Create external collimator. No effect on light collection.
	if (diameter_inter_ExternalColl < diameter_ExternalColl) {
		G4Tubs* solid_ExternalColl = new G4Tubs("solid_ExternalColl", diameter_inter_ExternalColl / 2.0,
			diameter_ExternalColl / 2.0, z_size_ExternalColl / 2.0, 0. *deg, 360.*deg);
		G4LogicalVolume* logic_ExternalColl	= new G4LogicalVolume(solid_ExternalColl, G4Material::GetMaterial("Fe"), "logic_ExternalColl", 0, 0, 0);
		G4VPhysicalVolume* phys_ExternalColl = new G4PVPlacement(0, // no rotation
			position_ExternalColl, logic_ExternalColl, "phys_ExternalColl", logicWorld, // its mother  volume
			false, 0, fCheckOverlaps);
	}

	//-------------------------------------------------------------------------------
	// Create interface grid. Parameterization is required.
	G4Box* solid_tracker_Interface_grid = new G4Box("solid_tracker_Interface_grid", x_size_tracker_Interface_grid / 2.0, y_size_tracker_Interface_grid / 2.0, radius_Interface_wire);
	G4LogicalVolume* logic_tracker_Interface_grid = new G4LogicalVolume(solid_tracker_Interface_grid, G4Material::GetMaterial("LAr"), "logic_tracker_Interface_grid", 0, 0, 0);
	G4VPhysicalVolume* phys_tracker_Interface_grid = new G4PVPlacement(0, position_interface_wire_container, logic_tracker_Interface_grid, "phys_tracker_Interface_grid",
		logic_LAr_inner, false, 0, fCheckOverlaps);
	// Interface wires
	G4Tubs* solid_Interface_wire = new G4Tubs("solid_Interface_wire", 0, radius_Interface_wire, length_Interface_wire / 2.0, 0.*deg, 360.*deg);
	G4LogicalVolume* logic_Interface_wire = new G4LogicalVolume(solid_Interface_wire, matAl, "logic_Interface_wire", 0, 0, 0);
	G4VPVParameterisation* param_Interface_wire = new DetectorParameterisation(N_Interface_wire, 1, 1, rotX_90, G4ThreeVector(0,0,0), G4ThreeVector(step_Interface_wire, 0, 0));
	G4VPhysicalVolume* phys_Interface_wire = new G4PVParameterised("phys_Interface_wire", logic_Interface_wire, logic_tracker_Interface_grid,
		kXAxis, N_Interface_wire, param_Interface_wire, fCheckOverlaps);

	//-------------------------------------------------------------------------------
	//FR4 substrate for interface grid
	G4Box* solid_interface_grid_substrate = new G4Box("solid_interface_grid_substrate", x_size_Interface_grid_substrate / 2.0, y_size_Interface_grid_substrate / 2.0, z_size_Interface_grid_substrate / 2.0);
	G4Box* solid_interface_grid_hole = new G4Box("solid_interface_grid_hole", x_size_tracker_Interface_grid / 2.0, y_size_tracker_Interface_grid / 2.0, z_size_Interface_grid_substrate / 2.0);
	G4SubtractionSolid* solid_interface_grid_subtraction = new G4SubtractionSolid("solid_interface_grid_subtraction", solid_interface_grid_substrate, solid_interface_grid_hole);
	G4LogicalVolume* logic_interface_grid = new G4LogicalVolume(solid_interface_grid_subtraction, matAl, "l_anode_grid", 0, 0, 0);
	G4VPhysicalVolume* phys_interface_grid = new G4PVPlacement(0, position_interface_frame, logic_interface_grid, "p_interface_grid",
		logic_LAr_inner, false, 0);

	// Create anode grid. Parameterisation is required.
	// Create tracker (container for parameterised volumes)
	G4Box* solid_tracker_anode_grid = new G4Box("solid_tracker_anode_grid", x_size_tracker_anode_grid / 2.0, y_size_tracker_anode_grid / 2.0, z_size_tracker_anode_grid / 2.0);
	G4LogicalVolume* logic_tracker_anode_grid = new G4LogicalVolume(solid_tracker_anode_grid, G4Material::GetMaterial("Air"), "logic_tracker_anode_grid", 0, 0, 0);
	G4VPhysicalVolume* phys_tracker_anode_grid = new G4PVPlacement(0, position_anode_grid, logic_tracker_anode_grid, "phys_tracker_anode_grid",
		logicWorld, false, 0, fCheckOverlaps);
	// Create anode grid frame
	G4Box* solid_anode_grid_substrate = new G4Box("anode_grid_substrate", size_anode_grid / 2.0, size_anode_grid / 2.0, thickness_anode_grid / 2.0);
	G4Box* solid_anode_grid_hole = new G4Box("anode_grid_hole", size_anode_grid_hole / 2.0, size_anode_grid_hole / 2.0, thickness_anode_grid / 2.0);
	G4SubtractionSolid* solid_anode_grid_frame = new G4SubtractionSolid("anode_grid__frame", solid_anode_grid_substrate, solid_anode_grid_hole);
	G4LogicalVolume* logic_anode_grid_frame = new G4LogicalVolume(solid_anode_grid_frame, matAl, "l_anode_grid", 0, 0, 0);
	G4VPhysicalVolume* phys_anode_grid_frame = new G4PVPlacement(0, position_anode_grid, logic_anode_grid_frame, "p_anode_grid",
		logicWorld, false, 0);
	// Create anode wire grid
	G4Tubs* solid_wire = new G4Tubs("solid_wire", 0, radius_wire, length_wire / 2.0, 0.*deg, 360.*deg);
	G4LogicalVolume* logic_wire = new G4LogicalVolume(solid_wire, matAl, "lwire", 0, 0, 0);
	G4VPVParameterisation* param_wire = new DetectorParameterisation(N_wire, 1, 1, rotX_90, G4ThreeVector(0,0,0), G4ThreeVector(step_wire, 0, 0));
	G4VPhysicalVolume* phys_wire = new G4PVParameterised("phys_wire", logic_wire,	logic_tracker_anode_grid,
		kXAxis, N_wire, param_wire, fCheckOverlaps);

	//--------------------------------------------------------------------------------
	// Create PMT's wire grid
	// Gas part
	G4Box* solid_PMTAnodeGridTrackerGas = new G4Box("solid_PMTAnodeGridTrackerGas", PMTAnodeGridTrackerGasXSize / 2.0, PMTAnodeGridTrackerGasYSize / 2.0, PMTAnodeGridTrackerThickness / 2.0);
	G4LogicalVolume* logic_PMTAnodeGridTrackerGas = new G4LogicalVolume(solid_PMTAnodeGridTrackerGas, G4Material::GetMaterial("Air"), "logic_PMTAnodeGridTrackerGas", 0, 0, 0);
	G4LogicalVolume* logic_PMTAnodeGridTrackerGasInner = new G4LogicalVolume(solid_PMTAnodeGridTrackerGas, G4Material::GetMaterial("Air"), "logic_PMTAnodeGridTrackerGasInner", 0, 0, 0);

	G4VPhysicalVolume* phys_PMTAnodeGridTrackerGas_1 = new G4PVPlacement(rotY_90, position_PMTAnodeGridTrackerGas_1,
		logic_PMTAnodeGridTrackerGas,	"phys_PMTAnodeGridTrackerGas_1", logicWorld, false, 0, fCheckOverlaps);

	G4VPhysicalVolume* phys_PMTAnodeGridTrackerGasInner_1 = new G4PVPlacement(rotY_90, position_PMTAnodeGridTrackerGasInner_1,
		logic_PMTAnodeGridTrackerGasInner, "phys_PMTAnodeGridTrackerGasInner_1", logicWorld, false, 0, fCheckOverlaps);

	// Liquid part
	G4Box* solid_PMTAnodeGridTrackerLiquid = new G4Box("solid_PMTAnodeGridTrackerLiquid", PMTAnodeGridTrackerLiquidXSize / 2.0, PMTAnodeGridTrackerLiquidYSize / 2.0, PMTAnodeGridTrackerThickness / 2.0);
	G4LogicalVolume* logic_PMTAnodeGridTrackerLiquid = new G4LogicalVolume(solid_PMTAnodeGridTrackerLiquid, G4Material::GetMaterial("LAr"), "logic_PMTAnodeGridTrackerLiquid", 0, 0, 0);
	G4LogicalVolume* logic_PMTAnodeGridTrackerLiquidInner = new G4LogicalVolume(solid_PMTAnodeGridTrackerLiquid, G4Material::GetMaterial("LAr"), "logic_PMTAnodeGridTrackerLiquidInner", 0, 0, 0);
	G4VPhysicalVolume* phys_PMTAnodeGridTrackerLiquid_1 = new G4PVPlacement(rotY_90, position_PMTAnodeGridTrackerLiquid_1,
		logic_PMTAnodeGridTrackerLiquid, "phys_PMTAnodeGridTrackerLiquid_1", logic_LAr_outer, false, 0, fCheckOverlaps);
	G4VPhysicalVolume* phys_PMTAnodeGridTrackerLiquidInner_1 = new G4PVPlacement(rotY_90, position_PMTAnodeGridTrackerLiquidInner_1,
		logic_PMTAnodeGridTrackerLiquidInner, "phys_PMTAnodeGridTrackerLiquidInner_1", logic_LAr_outer, false, 0, fCheckOverlaps);

	//--------------------------------------------------------------------------------
	//create PMTGridWire
	G4Tubs* solid_PMTGridWire = new G4Tubs("solid_PMTGridWire", 0, PMTGridWireRadius, PMTAnodeGridTrackerGasYSize / 2.0, 0.*deg, 360.*deg);
	G4LogicalVolume* logic_PMTGridWire = new G4LogicalVolume(solid_PMTGridWire, matAl, "lwire", 0, 0, 0);
	G4VPVParameterisation* param_PMTGridWireGas = new DetectorParameterisation(PMTAnodeGridNCellsGas, 1, 1, rotX_90, G4ThreeVector(0,0,0), G4ThreeVector(PMTGridWirePitch, 0, 0));
	G4VPhysicalVolume* phys_PMTGridWireGas_0 = new G4PVParameterised("phys_PMTGridWireGas_0", logic_PMTGridWire, logic_PMTAnodeGridTrackerGas, kXAxis,
		PMTAnodeGridNCellsGas, param_PMTGridWireGas, fCheckOverlaps);

	G4Tubs* solid_PMTGridWireGasInner = new G4Tubs("solid_PMTGridWireGasInner", 0, PMTGridWireRadius, PMTAnodeGridTrackerGasXSize / 2.0, 0.*deg, 360.*deg);
	G4LogicalVolume* logic_PMTGridWireGasInner = new G4LogicalVolume(solid_PMTGridWireGasInner, matAl, "lwire", 0, 0, 0);
	G4VPVParameterisation* param_PMTGridWireGasInner = new DetectorParameterisation(1, PMTAnodeGridNCellsGasInner, 1, rotY_90, G4ThreeVector(0,0,0), G4ThreeVector(0, PMTGridWirePitch, 0));
	G4VPhysicalVolume* phys_PMTGridWireGasInner_0 = new G4PVParameterised("phys_PMTGridWireGasInner_0", logic_PMTGridWireGasInner,
		logic_PMTAnodeGridTrackerGasInner, kXAxis, PMTAnodeGridNCellsGasInner, param_PMTGridWireGasInner, fCheckOverlaps);

	G4VPVParameterisation* param_PMTGridWireLiquid = new DetectorParameterisation(PMTAnodeGridNCellsLiquid, 1, 1, rotX_90, G4ThreeVector(0,0,0), G4ThreeVector(PMTGridWirePitch, 0, 0));
	G4VPhysicalVolume* phys_PMTGridWireLiquid_0 = new G4PVParameterised("phys_PMTGridWireLiquid_0", logic_PMTGridWire,
		logic_PMTAnodeGridTrackerLiquid, kXAxis, PMTAnodeGridNCellsLiquid, param_PMTGridWireLiquid, fCheckOverlaps);

	G4Tubs* solid_PMTGridWireLiquidInner = new G4Tubs("solid_PMTGridWireLiquidInner", 0, PMTGridWireRadius, PMTAnodeGridTrackerLiquidXSize / 2.0, 0.*deg, 360.*deg);
	G4LogicalVolume* logic_PMTGridWireLiquidInner = new G4LogicalVolume(solid_PMTGridWireLiquidInner, matAl, "lwire", 0, 0, 0);
	G4VPVParameterisation* param_PMTGridWireLiquidInner = new DetectorParameterisation(1, PMTAnodeGridNCellsLiquidInner, 1, rotY_90, G4ThreeVector(0,0,0), G4ThreeVector(0, PMTGridWirePitch, 0));
	G4VPhysicalVolume* phys_PMTGridWireLiquidInner_0 = new G4PVParameterised("phys_PMTGridWireLiquidInner_0", logic_PMTGridWireLiquidInner,
		logic_PMTAnodeGridTrackerLiquidInner, kXAxis, PMTAnodeGridNCellsLiquidInner, param_PMTGridWireLiquidInner, fCheckOverlaps);

	//--------------------------------------------------------------------------------
	// Create SiPM matrix
	// create tracker (this is required for SiPM_matrix parameterisation)
	G4Box* trackerS = new G4Box("tracker", x_size_tracker / 2.0, y_size_tracker / 2.0, z_size_tracker / 2.0);
	G4LogicalVolume* trackerLV = new G4LogicalVolume(trackerS, G4Material::GetMaterial("Air"), "Tracker", 0, 0, 0);
	G4VPhysicalVolume* phys_tracker = new G4PVPlacement(0, position_SiPM_container, trackerLV, "Tracker", logicWorld,
	 	false, 0, fCheckOverlaps);

	solid_SiPM = new G4Box("sscintillator", size_SiPM / 2.0, size_SiPM / 2.0, thickness_SiPM / 2.0);
	logic_SiPM = new G4LogicalVolume(solid_SiPM, matAl, "lSiPM", 0, 0, 0);
	G4VPVParameterisation* chamberParam = new DetectorParameterisation(Nx_SiPMs, Ny_SiPMs, 1, NULL, G4ThreeVector(0,0,0), G4ThreeVector(chamberSpacing, chamberSpacing, 0));
	G4VPhysicalVolume* phys_SiPM = new G4PVParameterised(gPars::det_dims.SiPM_device_name, logic_SiPM, trackerLV,
		kXAxis, Nx_SiPMs * Ny_SiPMs, chamberParam, fCheckOverlaps);

	G4Box* solid_SiPMFR4 = new G4Box("solid_SiPMFR4", x_size_tracker / 2.0, y_size_tracker / 2.0, z_size_SiPMFR4 / 2.0);
	G4LogicalVolume* logic_SiPMFR4 = new G4LogicalVolume(solid_SiPMFR4, G4Material::GetMaterial("Air"), "logic_SiPMFR4", 0, 0, 0);
	G4VPhysicalVolume* phys_SiPMFR4 = new G4PVPlacement(0, position_SiPMFR4, logic_SiPMFR4, "phys_SiPMFR4",
		logicWorld, false, 0, fCheckOverlaps);

	//--------------------------------------------------------------------------------
	// Create PMMA plate
	G4Box* solid_PMMA_plate = new G4Box("solid_tracker_anode_grid", x_size_PMMA_plate / 2.0, y_size_PMMA_plate / 2.0, z_size_PMMA_plate / 2.0);
	G4LogicalVolume* logic_PMMA_plate = new G4LogicalVolume(solid_PMMA_plate, G4Material::GetMaterial("PMMA"), "logic_PMMA_plate", 0, 0, 0);
	G4VPhysicalVolume* phys_PMMA_plate = new G4PVPlacement(0, position_PMMA_plate, logic_PMMA_plate, "phys_PMMA_plate",
		logicWorld, false, 0, fCheckOverlaps);

	//--------------------------------------------------------------------------------
	// Create THGEM1
	G4Box* solid_whole_THGEM1 = new G4Box("solid_whole_THGEM1", x_size_THGEM1 / 2.0, y_size_THGEM1 / 2.0, z_size_THGEM1 / 2.0);
	G4Box* solid_active_THGEM1 = new G4Box("solid_whole_THGEM1", x_size_THGEM1_container / 2.0, y_size_THGEM1_container / 2.0, z_size_THGEM1 / 2.0);
	G4SubtractionSolid* solid_THGEM1_frame = new G4SubtractionSolid("solid_THGEM1_frame", solid_whole_THGEM1, solid_active_THGEM1);
  G4LogicalVolume* logic_THGEM1_frame = new G4LogicalVolume(solid_THGEM1_frame, G4Material::GetMaterial("FR4"), "logic_THGEM1_frame", 0, 0, 0);
  G4VPhysicalVolume* phys_THGEM1_frame = new G4PVPlacement(0, position_THGEM1_frame, logic_THGEM1_frame, "phys_THGEM1_frame",
      logic_LAr_inner, false, 0, fCheckOverlaps);

  G4Box* solid_THGEM1_container = new G4Box("solid_THGEM1_container", x_size_THGEM1_container / 2.0, y_size_THGEM1_container / 2.0, z_size_THGEM1_container / 2.0);
  G4LogicalVolume* logic_THGEM1_container = new G4LogicalVolume(solid_THGEM1_container, G4Material::GetMaterial("LAr"), "logic_THGEM1_container", 0, 0, 0);
  G4VPhysicalVolume* phys_THGEM1_container = new G4PVPlacement(0, position_THGEM1_container, logic_THGEM1_container, gPars::det_dims.THGEM1_cell_container_name,
        logic_LAr_inner, false, 0, fCheckOverlaps);

  G4LogicalVolume* logic_THGEM1_copper = new G4LogicalVolume(solid_active_THGEM1, G4Material::GetMaterial("FR4"), "logic_THGEM1_copper", 0, 0, 0);
  G4VPhysicalVolume* phys_THGEM1_copper = new G4PVPlacement(0, position_THGEM1_copper_plate, logic_THGEM1_copper, "phys_THGEM1_copper",
      logic_THGEM1_container, false, 0, fCheckOverlaps);

	//--------------------------------------------------------------------------------
	// Create insulator box
	G4Box* solid_Insulator_box_inner = new G4Box("solid_Insulator_box_inner", x_size_Insulator_box_inner / 2.0, y_size_Insulator_box_inner / 2.0, z_size_Insulator_box / 2.0);
	G4Box* solid_Insulator_box_outer = new G4Box("solid_Insulator_box_outer", x_size_Insulator_box_outer / 2.0, y_size_Insulator_box_outer / 2.0, z_size_Insulator_box / 2.0);
	G4SubtractionSolid* solid_Insulator_box_subtraction = new G4SubtractionSolid("solid_Insulator_box_subtraction", solid_Insulator_box_outer, solid_Insulator_box_inner);
	G4LogicalVolume* logic_Insulator_box = new G4LogicalVolume(solid_Insulator_box_subtraction, G4Material::GetMaterial("PMMA_UV"), "logic_Insulator_box", 0, 0, 0);
	G4VPhysicalVolume* phys_Insulator_box = new G4PVPlacement(0, position_Insulator_box, logic_Insulator_box,
		"phys_Insulator_box", logicWorld, false, 0, fCheckOverlaps);

	//--------------------------------------------------------------------------------
	// Create PMTs
	G4Tubs* solid_PMT = new G4Tubs("solid_PMT", 0, radius_PMT, z_size_PMT / 2.0, 0.*deg, 360.*deg);
	logic_PMT = new G4LogicalVolume(solid_PMT, matAl, "logic_PMT", 0, 0, 0);

	G4VPhysicalVolume* phys_PMT0 = new G4PVPlacement(rotY_90, position_PMT_0, logic_PMT, gPars::det_dims.PMT_device_name,
		logicWorld, false, 0, fCheckOverlaps);
	G4VPhysicalVolume* phys_PMT1 = new G4PVPlacement(rotY_90, position_PMT_1, logic_PMT, gPars::det_dims.PMT_device_name,
		logicWorld, false, 1, fCheckOverlaps);
	G4VPhysicalVolume* phys_PMT2 = new G4PVPlacement(rotX_90, position_PMT_2, logic_PMT, gPars::det_dims.PMT_device_name,
		logicWorld, false, 2, fCheckOverlaps);
	G4VPhysicalVolume* phys_PMT3 = new G4PVPlacement(rotX_90, position_PMT_3, logic_PMT, gPars::det_dims.PMT_device_name,
		logicWorld, false, 3, fCheckOverlaps);

	//--------------------------------------------------------------------------------
	// Create PMTs' steel box
	G4Box* solidSteelBox = new G4Box("solidSteelBox", xSizeSteelBox / 2.0, ySizeSteelBox / 2.0, zSizeSteelBox / 2.0);
	G4VSolid* solidBoxSubtractPMT = new G4SubtractionSolid("solidBoxSubtractPMT", solidSteelBox, solid_PMT, rotY_90, G4ThreeVector((xSizeSteelBox - z_size_PMT)/2, 0., 0.));
	G4LogicalVolume* logicSteelBox = new G4LogicalVolume(solidBoxSubtractPMT, G4Material::GetMaterial("Air"), "logicSteelBox", 0, 0, 0);

	G4VPhysicalVolume* physSteelBox0 = new G4PVPlacement(0, positionSteelBox0, logicSteelBox, "physSteelBox0",
		logicWorld, false, 0, fCheckOverlaps);
	G4VPhysicalVolume* physSteelBox1 = new G4PVPlacement(rotZ_180, positionSteelBox1, logicSteelBox, "physSteelBox1",
		logicWorld, false, 0, fCheckOverlaps);
	G4VPhysicalVolume* physSteelBox2 = new G4PVPlacement(rotZ_270, positionSteelBox2, logicSteelBox, "physSteelBox2",
		logicWorld, false, 0, fCheckOverlaps);
	G4VPhysicalVolume* physSteelBox3 = new G4PVPlacement(rotZ_90, positionSteelBox3, logicSteelBox, "physSteelBox3",
		logicWorld, false, 0, fCheckOverlaps);

	//--------------------------------------------------------------------------------
	// Setting surfaces //adsf
	G4LogicalBorderSurface* LAr_inner2physiWorld = new G4LogicalBorderSurface("LAr_inner2physiWorld", phys_LAr_inner, physiWorld, LAr_OpticalSurface);
	G4LogicalBorderSurface* tracker2SiPM = new G4LogicalBorderSurface("tracker2SiPM", phys_tracker, phys_SiPM, SiPM_OpticalSurface);
	G4LogicalBorderSurface* trackerSiPM_SiPMFR4 = new G4LogicalBorderSurface("trackerSiPM_SiPMFR4", phys_tracker, phys_SiPMFR4, FR4_unified);

	// PMT photocathode has same surface for LAr and gas.
	G4LogicalSkinSurface* sur_PMT_cathode = new G4LogicalSkinSurface("PMT_cathode", logic_PMT, PMT_cathode);
	G4LogicalSkinSurface* sur_BottomCathode = new G4LogicalSkinSurface("BottomCathode", logic_bCathode, Cu_Cathode/*AbsorberMaterial*/);
	G4LogicalSkinSurface* CuReflector_THGEM0_surface = new G4LogicalSkinSurface("CuReflector_THGEM_surface", logic_THGEM1_copper, Cu_THGEM);
	G4LogicalSkinSurface* sur_SteelBox = new G4LogicalSkinSurface("SteelBox", logicSteelBox, stainlessSteel);

	// PMT grid
	G4LogicalBorderSurface* PMTAnodeGridTrackerLiquidGas = new G4LogicalBorderSurface("PMTAnodeGridTrackerLiquid-Gas", phys_PMTAnodeGridTrackerLiquid_1, phys_PMTAnodeGridTrackerGas_1, LAr_OpticalSurface);
	G4LogicalBorderSurface* PMTAnodeGridTrackerLiquidGasInner = new G4LogicalBorderSurface("PMTAnodeGridTrackerGas-Liquid", phys_PMTAnodeGridTrackerLiquidInner_1, phys_PMTAnodeGridTrackerGasInner_1, LAr_OpticalSurface);
	G4LogicalSkinSurface* PMTGridWire0 = new G4LogicalSkinSurface("PMTGridWire_surface0", logic_PMTGridWire, stainlessSteel);
	G4LogicalSkinSurface* PMTGridWire1 = new G4LogicalSkinSurface("PMTGridWire_surface1", logic_PMTGridWireGasInner, stainlessSteel);
	G4LogicalSkinSurface* PMTGridWire2 = new G4LogicalSkinSurface("PMTGridWire_surface2", logic_PMTGridWireLiquidInner, stainlessSteel);

	// SiPM (anode) grid
	G4LogicalBorderSurface* physiWorld2anode_grid = new G4LogicalBorderSurface("physiWorld2anode_grid", physiWorld, phys_anode_grid_frame, /*AbsorberMaterial*/ FR4_unified);
	G4LogicalBorderSurface* tracker_anode_grid2wire = new G4LogicalBorderSurface("tracker_anode_grid2wire", phys_tracker_anode_grid, phys_wire, /*AbsorberMaterial*/ Anode_wire_unified);
	G4LogicalBorderSurface* tracker_anode_grid2anode_grid = new G4LogicalBorderSurface("tracker_anode_grid2anode_grid", phys_tracker_anode_grid, phys_anode_grid_frame, /*AbsorberMaterial*/ FR4_unified);

	// Interface grid
	G4LogicalBorderSurface* phys_LAr_inner2Interface_grid = new G4LogicalBorderSurface("phys_LAr_inner2Interface_grid", phys_LAr_inner, phys_interface_grid, /*AbsorberMaterial*/ FR4_unified);
	G4LogicalBorderSurface* tracker_Interface_grid2wire = new G4LogicalBorderSurface("tracker_anode_grid2wire", phys_tracker_Interface_grid, phys_Interface_wire, /*AbsorberMaterial*/ Anode_wire_unified);
	G4LogicalBorderSurface* tracker_Interface_grid2Interface_grid = new G4LogicalBorderSurface("tracker_Interface_grid2Interface_grid", phys_tracker_Interface_grid, phys_interface_grid, /*AbsorberMaterial*/ FR4_unified);

	//G4LogicalBorderSurface* PMMA_plate2anode_grid = new G4LogicalBorderSurface("PMMA_plate2anode_grid", phys_PMMA_plate, phys_anode_grid_frame, FR4_unified);
	G4LogicalSkinSurface* AlWindow = new G4LogicalSkinSurface("AlWindow_surface", logic_Al_window, stainlessSteel);

	G4LogicalSkinSurface* THGEM1_cell_cu = new G4LogicalSkinSurface("THGEM1_cell_cu_surface", logic_THGEM1_cell_copper, Cu_THGEM);
  G4LogicalSkinSurface* THGEM1_cell_FR4 = new G4LogicalSkinSurface("THGEM1_cell_FR4_surface", logic_THGEM1_cell_FR4, FR4_unified);
  G4LogicalBorderSurface* THGEM1_cell_isolation = new G4LogicalBorderSurface("THGEM1_cell_isolation_surface", physiWorld, phys_THGEM1_cell, AbsorberMaterial);
  G4LogicalBorderSurface* THGEM1_cell_isolation1 = new G4LogicalBorderSurface("THGEM1_cell_isolation_surface1", phys_THGEM1_cell, physiWorld, AbsorberMaterial);

	G4LogicalBorderSurface* World2THGEM1_without_hole = new G4LogicalBorderSurface("World2THGEM1_without_hole", physiWorld, phys_THGEM1_frame, FR4_unified);
	G4LogicalBorderSurface* phys_LAr_inner_2_THGEM1_without_hole = new G4LogicalBorderSurface("phys_LAr_inner_2_THGEM1_without_hole", phys_LAr_inner, phys_THGEM1_frame, FR4_unified);

	//FieldWires
	G4LogicalSkinSurface* sur_FieldWire = new G4LogicalSkinSurface("sur_FieldWire", logic_FieldWire, stainlessSteel);

	//--------------------------------------------------------------------------------
	// Setting visualization
	G4VisAttributes* LAr_VisAtt = new G4VisAttributes(G4Colour(0.6, 0.6, 1.0, 0.0));
	G4VisAttributes* gas_VisAtt = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0, 0.0));
	G4VisAttributes* PMMA_VisAtt = new G4VisAttributes(G4Colour(0.2, 0.2, 0.8, 0.0));
	G4VisAttributes* Wires_VisAtt = new G4VisAttributes(G4Colour(1.0, 0.1, 0.1, 0.4));
	G4VisAttributes* FR4_VisAtt = new G4VisAttributes(G4Colour(1.0, 0.5, 0.2, 0.4));
	G4VisAttributes* Sensor_VisAtt = new G4VisAttributes(G4Colour(0.8, 0.7, 0.2, 0.8));
	G4VisAttributes* Cu_VisAtt = new G4VisAttributes(G4Colour(0.0, 1.0, 1.0, 0.4));
	G4VisAttributes* Steel_VisAtt = new G4VisAttributes(G4Colour(0.2, 0.2, 0.2, 0.1));

	// SiPMs
	logic_wire->SetVisAttributes(Wires_VisAtt);
	logic_anode_grid_frame->SetVisAttributes(FR4_VisAtt);
	logic_SiPM->SetVisAttributes(Sensor_VisAtt);
  trackerLV->SetVisAttributes(gas_VisAtt);
  logic_SiPMFR4->SetVisAttributes(FR4_VisAtt);

	// Interface grid
	logic_Interface_wire->SetVisAttributes(Wires_VisAtt);
	logic_interface_grid->SetVisAttributes(FR4_VisAtt);

	// PMTs
	logic_PMT->SetVisAttributes(Sensor_VisAtt);
	logic_PMTGridWire->SetVisAttributes(Wires_VisAtt);
	logic_PMTGridWireGasInner->SetVisAttributes(gas_VisAtt);
	logic_PMTGridWireLiquidInner->SetVisAttributes(LAr_VisAtt);
	logicSteelBox->SetVisAttributes(Steel_VisAtt);

	// Separate THGEM hole
	logic_THGEM1_cell_copper->SetVisAttributes(Cu_VisAtt);
	logic_THGEM1_cell_LAr->SetVisAttributes(LAr_VisAtt);
	logic_THGEM1_cell->SetVisAttributes(LAr_VisAtt);
	logic_THGEM1_cell_FR4->SetVisAttributes(FR4_VisAtt);

	// THGEM1
	logic_THGEM1_copper->SetVisAttributes(Cu_VisAtt);
	logic_THGEM1_container->SetVisAttributes(LAr_VisAtt);
	logic_THGEM1_frame->SetVisAttributes(FR4_VisAtt);

	// LAr
	logic_LAr_inner->SetVisAttributes(LAr_VisAtt);
	logic_PMTAnodeGridTrackerLiquid->SetVisAttributes(LAr_VisAtt);
	logic_PMTAnodeGridTrackerLiquidInner->SetVisAttributes(LAr_VisAtt);
	logic_LAr_outer->SetVisAttributes(LAr_VisAtt);
	logic_LArInactive->SetVisAttributes(LAr_VisAtt);

	logic_bCathode->SetVisAttributes(FR4_VisAtt);
	logic_Al_window->SetVisAttributes(new G4VisAttributes(G4Colour(1.0, 1.0, 0.0, 0.0)));
	logic_CryogenicChamberBottom->SetVisAttributes(Steel_VisAtt);
	logic_PMMA_bottom->SetVisAttributes(PMMA_VisAtt);
	logic_PMMA_plate->SetVisAttributes(PMMA_VisAtt);
	logic_Insulator_box->SetVisAttributes(PMMA_VisAtt);

	logic_FieldWire->SetVisAttributes(Wires_VisAtt);

	logicWorld->SetVisAttributes(G4Colour(1.0, 1.0, 1.0, 0.0));
	return physiWorld;
}

void DetectorConstruction::ConstructSDandField()
{
  thePhotoDetector = new DetectorSensor("/detector/sensitiveDetector");
  G4SDManager::GetSDMpointer()->AddNewDetector(thePhotoDetector);
  SetSensitiveDetector(logic_SiPM, thePhotoDetector);
  SetSensitiveDetector(logic_PMT, thePhotoDetector);
  //G4SDManager::GetSDMpointer()->Activate("/detector/sensitiveDetector", true);
}

void DetectorConstruction::CreateTHGEM1Cell() //Must be the same as in gmsh-elmer simulation
{
  double cell_size_x = gPars::det_dims.THGEM1_hole_pitch / 2.0 * mm;
  double cell_size_y = cell_size_x * std::sqrt(3.0);
  double cell_size_z = z_size_THGEM1_container;
  double diel_size_z = gPars::det_dims.THGEM1_dielectric_thickness * mm;
  double radius = gPars::det_dims.THGEM1_hole_radius * mm;
  double radius_cu = radius + gPars::det_dims.THGEM1_hole_rim * mm;
  double cu_size_z = gPars::det_dims.THGEM1_copper_thickness * mm;

  G4ThreeVector zero(0.0, 0.0, 0.0);
  G4ThreeVector hole_1_pos(cell_size_x / 2.0, cell_size_y / 2.0, 0.0);
  G4ThreeVector hole_2_pos(-cell_size_x / 2.0, -cell_size_y / 2.0, 0.0);
  G4ThreeVector cu_top_pos(0.0, 0.0, diel_size_z / 2.0 + cu_size_z / 2.0);
  G4ThreeVector cu_bot_pos(0.0, 0.0, -diel_size_z / 2.0 - cu_size_z / 2.0);

  G4Box* solid_THGEM1_cell_isolation = new G4Box("solid_THGEM1_cell_isolation", cell_size_x, cell_size_y, cell_size_z); //so that cell is in the same material
  logic_THGEM1_cell = new G4LogicalVolume(solid_THGEM1_cell_isolation, G4Material::GetMaterial("LAr"), "logic_THGEM1_cell_isolation", 0, 0, 0);

  G4Box* solid_THGEM1_cell_LAr = new G4Box("solid_THGEM1_cell_LAr", cell_size_x / 2.0, cell_size_y / 2.0, cell_size_z / 2.0);
  logic_THGEM1_cell_LAr = new G4LogicalVolume(solid_THGEM1_cell_LAr, G4Material::GetMaterial("LAr"), "logic_THGEM1_cell_LAr", 0, 0, 0);
  G4VPhysicalVolume* phys_THGEM1_cell_LAr = new G4PVPlacement(0, zero, logic_THGEM1_cell_LAr, gPars::det_dims.THGEM1_cell_name,
      logic_THGEM1_cell, false, 0, fCheckOverlaps);

  G4Box* solid_THGEM1_diel_box = new G4Box("solid_THGEM1_diel_box", cell_size_x / 2.0, cell_size_y / 2.0, diel_size_z / 2.0);
  G4Tubs* solid_THGEM1_diel_hole = new G4Tubs("solid_THGEM1_diel_hole", 0, radius, diel_size_z / 2.0, 0.*deg, 360.*deg);
  G4SubtractionSolid* solid_THGEM1_diel_tmp = new G4SubtractionSolid("solid_THGEM1_diel_tmp", solid_THGEM1_diel_box, solid_THGEM1_diel_hole, 0, hole_1_pos);
  G4SubtractionSolid* solid_THGEM1_diel = new G4SubtractionSolid("solid_THGEM1_diel", solid_THGEM1_diel_tmp, solid_THGEM1_diel_hole, 0, hole_2_pos);
  logic_THGEM1_cell_FR4 = new G4LogicalVolume(solid_THGEM1_diel, G4Material::GetMaterial("FR4"), "logic_THGEM1_cell_FR4", 0, 0, 0);
  G4VPhysicalVolume* phys_THGEM1_cell_FR4 = new G4PVPlacement(0, zero, logic_THGEM1_cell_FR4, "phys_THGEM1_cell_FR4",
      logic_THGEM1_cell_LAr, false, 0, fCheckOverlaps);

  G4Box* solid_THGEM1_cu_box = new G4Box("solid_THGEM1_cu_box", cell_size_x / 2.0, cell_size_y / 2.0, cu_size_z / 2.0);
  G4Tubs* solid_THGEM1_cu_hole = new G4Tubs("solid_THGEM1_cu_hole", 0, radius_cu, cu_size_z / 2.0, 0.*deg, 360.*deg);
  G4SubtractionSolid* solid_THGEM1_cu_tmp = new G4SubtractionSolid("solid_THGEM1_cu_tmp", solid_THGEM1_cu_box, solid_THGEM1_cu_hole, 0, hole_1_pos);
  G4SubtractionSolid* solid_THGEM1_cu = new G4SubtractionSolid("solid_THGEM1_cu", solid_THGEM1_cu_tmp, solid_THGEM1_cu_hole, 0, hole_2_pos);
  logic_THGEM1_cell_copper = new G4LogicalVolume(solid_THGEM1_cu, G4Material::GetMaterial("FR4"), "logic_THGEM1_cell_copper", 0, 0, 0);
  G4VPhysicalVolume* phys_THGEM1_cell_copper_top = new G4PVPlacement(0, cu_top_pos, logic_THGEM1_cell_copper, "phys_THGEM1_cell_copper_top",
      logic_THGEM1_cell_LAr, false, 0, fCheckOverlaps);
  G4VPhysicalVolume* phys_THGEM1_cell_copper_bot = new G4PVPlacement(0, cu_bot_pos, logic_THGEM1_cell_copper, "phys_THGEM1_cell_copper_bot",
      logic_THGEM1_cell_LAr, false, 0, fCheckOverlaps);
}

void DetectorConstruction::defineSurfaces()
{
	//*********************************************************************************
	//unified model
	G4double ener[2] = {.1*eV, 10.*eV};
	G4double zero[2] = {0.0, 0.0};
	//-------------------------------------------------------------------------------
	LAr_OpticalSurface = new G4OpticalSurface("LAr_OpticalSurface");
	LAr_OpticalSurface->SetType(dielectric_dielectric);
	LAr_OpticalSurface->SetModel(unified);
	LAr_OpticalSurface->SetFinish(polished);
	LAr_OpticalSurface->SetSigmaAlpha(0);//alpha in degrees, from 0 to 90.
	G4MaterialPropertiesTable *LAr_MaterialProperty = new G4MaterialPropertiesTable();
	G4double LAr_Materialeff[2] = { 0, 0 };
	LAr_MaterialProperty->AddProperty("EFFICIENCY", ener, LAr_Materialeff, 2);
	LAr_OpticalSurface->SetMaterialPropertiesTable(LAr_MaterialProperty);
	//-------------------------------------------------------------------------------
	//FR4_unified
	FR4_unified = new G4OpticalSurface("FR4_unified");
	FR4_unified->SetType(dielectric_metal);
	FR4_unified->SetModel(unified);
	FR4_unified->SetFinish(ground);
	FR4_unified->SetSigmaAlpha(gPars::det_opt.FR4_SigmaAlpha);//alpha in degrees, from 0 to 90.
	G4MaterialPropertiesTable *FR4_MaterialProperty = new G4MaterialPropertiesTable();
	G4double FR4_Materialrefl[2] = {0.2, 0.2};//model_25b
	//G4double FR4_Materialrefl[2] = { 0.05, 0.05 };//https://www.cetem.gov.br/images/congressos/2008/CAC00560008.pdf
	G4double FR4_Materialeff[2] = { 0, 0 };
	FR4_MaterialProperty->AddProperty("REFLECTIVITY", ener, gPars::debugging.no_reflections ? zero : FR4_Materialrefl, 2);
	FR4_MaterialProperty->AddProperty("EFFICIENCY", ener, FR4_Materialeff, 2);
	FR4_unified->SetMaterialPropertiesTable(FR4_MaterialProperty);
	//-------------------------------------------------------------------------------
	//Anode_Wire (Bronze-Berillium)
	Anode_wire_unified = new G4OpticalSurface("Anode_wire_unified");
	Anode_wire_unified->SetType(dielectric_metal);
	Anode_wire_unified->SetModel(unified);
	Anode_wire_unified->SetFinish(polished);
	Anode_wire_unified->SetSigmaAlpha(gPars::det_opt.Wire_SigmaAlpha * degree);//alpha in degrees, from 0 to 90.
	G4MaterialPropertiesTable *Anode_wire_MaterialProperty = new G4MaterialPropertiesTable();
	G4double Anode_wire_Materialrefl[2] = { 0.5, 0.5 };//approximately https://nvlpubs.nist.gov/nistpubs/bulletin/07/nbsbulletinv7n2p197_A2b.pdf
	G4double Anode_wire_Materialeff[2] = { 0, 0 };
	Anode_wire_MaterialProperty->AddProperty("REFLECTIVITY", ener, gPars::debugging.no_reflections ? zero : Anode_wire_Materialrefl, 2);
	Anode_wire_MaterialProperty->AddProperty("EFFICIENCY", ener, Anode_wire_Materialeff, 2);
	Anode_wire_unified->SetMaterialPropertiesTable(Anode_wire_MaterialProperty);
	//-------------------------------------------------------------------------------
	//Cu_THGEM
	Cu_THGEM = new G4OpticalSurface("Cu_THGEM");
	Cu_THGEM->SetType(dielectric_metal);
	Cu_THGEM->SetModel(unified);
	Cu_THGEM->SetFinish(ground);
	Cu_THGEM->SetSigmaAlpha(gPars::det_opt.Cu_SigmaAlpha * degree);//alpha in degrees, from 0 to 90.
	G4MaterialPropertiesTable *Cu_THGEM_MaterialProperty = new G4MaterialPropertiesTable();
	G4double Cu_THGEM_Materialrefl[2] = {0.36, 0.36}; // Bass M. Handbook of optics, Vol.4 Edition3
	G4double Cu_THGEM_Materialeff[2] = {0, 0};
	Cu_THGEM_MaterialProperty->AddProperty("REFLECTIVITY", ener, gPars::debugging.no_reflections ? zero : Cu_THGEM_Materialrefl, 2);
	Cu_THGEM_MaterialProperty->AddProperty("EFFICIENCY", ener, Cu_THGEM_Materialeff, 2);
	Cu_THGEM->SetMaterialPropertiesTable(Cu_THGEM_MaterialProperty);
	//-------------------------------------------------------------------------------
	//Cu_Cathode
	Cu_Cathode = new G4OpticalSurface("Cu_Cathode");
	Cu_Cathode->SetType(dielectric_metal);
	Cu_Cathode->SetModel(unified);
	Cu_Cathode->SetFinish(polished);
	Cu_Cathode->SetSigmaAlpha(gPars::det_opt.Cu_SigmaAlpha * degree);//alpha in degrees, from 0 to 90.
	G4MaterialPropertiesTable *Cu_Cathode_MaterialProperty = new G4MaterialPropertiesTable();
	G4double Cu_Cathode_Materialrefl[2] = {0.36, 0.36}; // Bass M. Handbook of optics, Vol.4 Edition3
	G4double Cu_Cathode_Materialeff[2] = {0, 0};
	Cu_Cathode_MaterialProperty->AddProperty("REFLECTIVITY", ener, gPars::debugging.no_reflections ? zero : Cu_Cathode_Materialrefl, 2);
	Cu_Cathode_MaterialProperty->AddProperty("EFFICIENCY", ener, Cu_Cathode_Materialeff, 2);
	Cu_Cathode->SetMaterialPropertiesTable(Cu_Cathode_MaterialProperty);
	//-------------------------------------------------------------------------------
	//stainlessSteel
	stainlessSteel = new G4OpticalSurface("stainlessSteel");
	stainlessSteel->SetType(dielectric_metal);
	stainlessSteel->SetModel(unified);
	stainlessSteel->SetFinish(polished);
	stainlessSteel->SetSigmaAlpha(gPars::det_opt.StainlessSteel_SigmaAlpha * degree);//alpha in degrees, from 0 to 90.
	G4MaterialPropertiesTable *stainlessSteelMaterialProperty = new G4MaterialPropertiesTable();
	G4double stainlessSteelMaterialrefl[2] = {0.5, 0.5}; // doi:10.1063/1.2202915
	G4double stainlessSteelMaterialeff[2] = {0, 0};
	stainlessSteelMaterialProperty->AddProperty("REFLECTIVITY", ener, gPars::debugging.no_reflections ? zero : stainlessSteelMaterialrefl, 2);
	stainlessSteelMaterialProperty->AddProperty("EFFICIENCY", ener, stainlessSteelMaterialeff, 2);
	stainlessSteel->SetMaterialPropertiesTable(stainlessSteelMaterialProperty);
	//-----------------------------------------------------------------------------
	// PMT cathode
	PMT_cathode = new G4OpticalSurface("PMT_cathode", unified);
	PMT_cathode->SetType(dielectric_metal);
	PMT_cathode->SetModel(unified);
	PMT_cathode->SetFinish(polished);
	PMT_cathode->SetSigmaAlpha(0.);
	G4MaterialPropertiesTable* PMT_cathodeMaterialProperty = new G4MaterialPropertiesTable();
	G4double cathoderefl[2] = {0, 0};
	G4double cathodeeff[2] = {1, 1};
	PMT_cathodeMaterialProperty->AddProperty("EFFICIENCY", ener, cathodeeff, 2);//dummy
	PMT_cathodeMaterialProperty->AddProperty("REFLECTIVITY", ener, cathoderefl, 2);
	PMT_cathode->SetMaterialPropertiesTable(PMT_cathodeMaterialProperty);
	//-----------------------------------------------------------------------------
	// SiPM_OpticalSurface
	SiPM_OpticalSurface = new G4OpticalSurface("SiPM_OpticalSurface", unified);
	SiPM_OpticalSurface->SetType(dielectric_metal);
	SiPM_OpticalSurface->SetModel(unified);
	SiPM_OpticalSurface->SetFinish(polished);
	SiPM_OpticalSurface->SetSigmaAlpha(0.);
	G4MaterialPropertiesTable* SiPM_MaterialProperty = new G4MaterialPropertiesTable();
	G4double SiPM_refl[2] = {0, 0};
	G4double SiPMeff[2] = {1, 1};
	SiPM_MaterialProperty->AddProperty("REFLECTIVITY", ener, SiPM_refl, 2);
	SiPM_MaterialProperty->AddProperty("EFFICIENCY", ener, SiPMeff, 2);
	SiPM_OpticalSurface->SetMaterialPropertiesTable(SiPM_MaterialProperty);
	//-----------------------------------------------------------------------------
	// Chamber walls
	AbsorberMaterial = new G4OpticalSurface("Absorber", unified);
	AbsorberMaterial->SetType(dielectric_metal);
	AbsorberMaterial->SetModel(unified);
	AbsorberMaterial->SetFinish(polished);
	AbsorberMaterial->SetSigmaAlpha(0.);
	G4MaterialPropertiesTable *AbsorberMaterialProperty = new G4MaterialPropertiesTable();
	G4double AbsorberMaterialrefl[2] = {0, 0};
	G4double AbsorberMaterialeff[2] = {0, 0};
	AbsorberMaterialProperty->AddProperty("REFLECTIVITY", ener, AbsorberMaterialrefl, 2);
	AbsorberMaterialProperty->AddProperty("EFFICIENCY", ener, AbsorberMaterialeff, 2);
	AbsorberMaterial->SetMaterialPropertiesTable(AbsorberMaterialProperty);
}

void DetectorConstruction::SetSizeAndPosition()
{
	HalfWorldLength = 17 * cm;

	//PMTs
	radius_PMT = 45 * mm / 2.0;
	z_size_PMT = 2 * mm; //1 * um;
	x_pos_PMT = 152 * mm / 2.0 + z_size_PMT / 2;
	y_pos_PMT = x_pos_PMT;
	z_pos_PMT = 27.2 * mm + 63 * mm / 2.0;

	//anode wire
	radius_wire = 100 / 2.0 * um;//you can understand this from photo
	length_wire = 60 * mm; // 108*mm future case
	step_wire = 1 * mm;
	N_wire = length_wire / step_wire - 1; // 107 in future case

	//interface wire
	radius_Interface_wire = 100 / 2.0 * um; //from Chegodaev
	length_Interface_wire = 100 * mm; // active region
	step_Interface_wire = 1 * mm; //from Chegodaev
	N_Interface_wire = length_Interface_wire / step_wire - 1;

	//Anode_grid
	thickness_anode_grid = 0.5 * mm;
	size_anode_grid = 127*mm ;//see Download:\DetectorPhotos\2021\THGEM_Electroconnect
	size_anode_grid_hole = length_wire;
	z_anode_grid_bottom = 78.2 * mm; //78.2 in case of one THGEM, 82.7*mm in case of two THGEMs;
	double z_anode_grid_center = z_anode_grid_bottom + thickness_anode_grid / 2.0;

	//PMMA plate
	x_size_PMMA_plate = size_anode_grid;
	y_size_PMMA_plate = size_anode_grid;
	z_size_PMMA_plate = 1.5 * mm;
	double z_PMMA_plate_center = z_anode_grid_center + thickness_anode_grid / 2.0 + z_size_PMMA_plate / 2.0;

	//SiPMs
	Nx_SiPMs = gPars::det_dims.n_SiPMs_rows;
	Ny_SiPMs = gPars::det_dims.n_SiPMs_rows;
	thickness_SiPM = 1 * nm;
	size_SiPM = 6.0 * mm;
	chamberSpacing = 10 * mm;
	double z_SiPM_bottom = z_anode_grid_bottom + thickness_anode_grid + z_size_PMMA_plate + (0.1*mm /*small gap between PMMA and SiPM*/);// 85.7*mm in case of two THGEMs
	double z_SiPM_center = z_SiPM_bottom + thickness_SiPM / 2.0;
	z_size_SiPMFR4 = 2*mm;

	//tracker SiPM
	x_size_tracker = Nx_SiPMs * chamberSpacing + size_SiPM / 2.0;
	y_size_tracker = Ny_SiPMs * chamberSpacing + size_SiPM / 2.0;
	z_size_tracker = 0.1 * mm;

	//tracker Anode_grid
	x_size_tracker_anode_grid = length_wire;
	y_size_tracker_anode_grid = x_size_tracker_anode_grid;
	z_size_tracker_anode_grid = thickness_anode_grid;

	//tracker THGEM2 (active region with holes)
	x_size_tracker_THGEM2 = 100 * mm;
	y_size_tracker_THGEM2 = 100 * mm;
	z_size_tracker_THGEM2 = gPars::det_dims.THGEM1_width_total * mm;

	// THGEM1
	x_size_THGEM1 = size_anode_grid;
  y_size_THGEM1 = size_anode_grid;
  z_size_THGEM1 = gPars::det_dims.THGEM1_width_total * mm;
  x_size_THGEM1_container = gPars::det_dims.THGEM1_active_area_size * mm;
  y_size_THGEM1_container = gPars::det_dims.THGEM1_active_area_size * mm;
  z_size_THGEM1_container = gPars::det_dims.THGEM1_container_width * mm;

	//Interface_grid
	x_size_tracker_Interface_grid = x_size_tracker_THGEM2;
	y_size_tracker_Interface_grid = y_size_tracker_THGEM2;
	x_size_Interface_grid_substrate = size_anode_grid;
	y_size_Interface_grid_substrate = size_anode_grid;
	z_size_Interface_grid_substrate = gPars::det_dims.width_interface_grid_frame * mm;

	//THGEM_without_holes
	x_size_THGEM_without_holes = size_anode_grid;
	y_size_THGEM_without_holes = size_anode_grid;
	z_size_THGEM_without_holes = z_size_tracker_THGEM2;

	//Insulator_box
	x_size_Insulator_box_inner = 143 * mm;
	y_size_Insulator_box_inner = x_size_Insulator_box_inner;
	thickness_Insulator_box = 4 * mm;
	x_size_Insulator_box_outer = x_size_Insulator_box_inner + thickness_Insulator_box * 2;
	y_size_Insulator_box_outer = x_size_Insulator_box_outer;
	z_size_Insulator_box = 150 * mm;
	z_Insulator_box_center = z_size_Insulator_box / 2.0;

	//WLS
	/*const double radius_WLS = 70 * mm / 2.0;
	const double z_size_WLS = 100 * um;
	const double x_pos_WLS = 152 * mm / 2.0 + z_size_PMT / 2;
	const double y_pos_WLS = x_pos_WLS;
	const double z_pos_WLS = z_pos_PMT;*/

	//LAr_inner
	x_size_LAr_inner = x_size_Insulator_box_inner;
	y_size_LAr_inner = y_size_Insulator_box_inner;
	z_size_LAr_inner = gPars::det_dims.z_top_interface_grid * mm + (22.0 - gPars::det_dims.EL_gap_thickness);

	//LArOuter
	x_size_LAr_outer_in = x_size_Insulator_box_outer;
	y_size_LAr_outer_in = y_size_Insulator_box_outer;
	x_size_LAr_outer_out = 152 * mm;
	y_size_LAr_outer_out = x_size_LAr_outer_out;
	z_size_LAr_outer = z_size_LAr_inner;

	//FieldTHGEM
	x_size_FieldTHGEM = size_anode_grid;
	y_size_FieldTHGEM = size_anode_grid;
	z_size_FieldTHGEM = z_size_tracker_THGEM2;
	z_center_FieldTHGEM_1 = 18.2*mm + z_size_FieldTHGEM / 2;
	z_center_FieldTHGEM_2 = 34.2*mm + z_size_FieldTHGEM / 2;
	hole_size_FieldTHGEM = 88 * mm;

	//FieldWire
	radius_FieldWire = 1.5 * mm / 2.0;
	length_FieldWire = 95 * mm;
	x_pos_FieldWire = x_size_tracker_THGEM2 / 2.0;
	z_pos_FieldWire_bottom = 20 /*18.2*mm - radius_FieldWire*/;//see 210415 1618834738620-1618834738657 photos
	z_pos_FieldWire_top = 33 /*34.2*mm + radius_FieldWire*/;//see 210415 1618834738620-1618834738657 photos

	//Cathode
	x_size_Cathode = x_size_LAr_outer_out;
	y_size_Cathode = y_size_LAr_outer_out;
	z_size_Cathode = 0.5 * mm;

	//LArInactive
	x_size_LArInactive = x_size_LAr_outer_out;
	y_size_LArInactive = y_size_LAr_outer_out;
	z_size_LArInactive = 2 * mm;

	//PMMA_bottom
	x_size_PMMA_bottom = x_size_LAr_outer_out;
	y_size_PMMA_bottom = y_size_LAr_outer_out;
	z_size_PMMA_bottom = 3 * mm;
	PMMA_bottom_center = -z_size_Cathode - z_size_LArInactive - z_size_PMMA_bottom / 2.0;

	//Al_window
	diameter_size_Al_window = 50 * mm;
	z_size_Al_window = 1.0 * mm;
	z_space_Al_window = 21 * mm;
	Al_window_top_center = PMMA_bottom_center - z_size_PMMA_bottom / 2.0 - z_size_Al_window / 2.0;
	Al_window_bottom_center = Al_window_top_center - z_size_Al_window / 2.0 - z_space_Al_window - z_size_Al_window / 2.0;

	//CryogenicChamberBottom
	diameter_size_internal_CryogenicChamberBottom = diameter_size_Al_window;
	diameter_size_external_CryogenicChamberBottom = x_size_LAr_inner * sqrt(2);
	z_size_CryogenicChamberBottom = 10 * mm;
	CryogenicChamberBottom_center = PMMA_bottom_center - z_size_PMMA_bottom / 2.0 - z_size_CryogenicChamberBottom / 2.0;

	//ExternalColl
	diameter_ExternalColl = diameter_size_Al_window;
	z_size_ExternalColl = 12 * mm;
	ExternalColl_center = Al_window_bottom_center - z_size_Al_window / 2.0 - z_size_ExternalColl / 2.0;
	ExternalColl_bottom = ExternalColl_center - z_size_ExternalColl / 2.0;

	if (thickness_anode_grid < 2 * radius_wire)
	{
		G4Exception("DetectorConstruction::SetSizeAndPosition: ",
			"InvalidSetup", FatalException, "thickness_anode_grid < 2*radius_wire");
		return;
	}

	//PMTGridWire
	PMTGridWireRadius = 150 / 2.0 * um;
	PMTGridWirePitch = 1.2 * mm;

	//PMTAnodeGridTracker
	PMTAnodeGridTrackerThickness = PMTGridWireRadius * 2;
	PMTAnodeGridTrackerGasYSize = 50 * mm;
	PMTAnodeGridTrackerGasXSize = radius_PMT - (z_size_LAr_inner - z_pos_PMT);
	PMTAnodeGridTrackerGasXSize = std::max(PMTAnodeGridTrackerGasXSize, 0.0);
	PMTAnodeGridTrackerZbottom = z_pos_PMT;
	PMTAnodeGridNCellsGas = PMTAnodeGridTrackerGasXSize / PMTGridWirePitch;
	PMTAnodeGridNCellsGasInner = PMTAnodeGridTrackerGasYSize / PMTGridWirePitch;

	PMTAnodeGridTrackerLiquidXSize = 2*radius_PMT - PMTAnodeGridTrackerGasXSize;
	PMTAnodeGridTrackerLiquidYSize = PMTAnodeGridTrackerGasYSize;
	PMTAnodeGridNCellsLiquid = PMTAnodeGridTrackerLiquidXSize / PMTGridWirePitch;
	PMTAnodeGridNCellsLiquidInner = PMTAnodeGridTrackerLiquidYSize / PMTGridWirePitch;

	//SteelBox
	xSizeSteelBox = 3;
	ySizeSteelBox = x_size_LAr_outer_out;
	zSizeSteelBox = 70;

	//Alpha
	radiusAlphaFull = 12;
	z_size_Alpha = 2;

	//TPB
	radiusTPB = 35;
	z_size_TPB = 0.2;

	position_SingleTHGEMCell = gPars::det_dims.THGEM1_single_cell_position;

	position_anode_grid = G4ThreeVector(0, 0, z_anode_grid_center);
	position_SiPM_container = G4ThreeVector(0, 0, z_SiPM_center);
	position_PMMA_plate = G4ThreeVector(0, 0, z_PMMA_plate_center);
	position_SiPMFR4 = G4ThreeVector(0, 0, z_SiPM_center + z_size_tracker /2.0 + z_size_SiPMFR4/2.0);

	position_THGEM1_frame = G4ThreeVector(0, 0, gPars::det_dims.z_bottom_THGEM1 * mm + z_size_THGEM1 / 2.0 - z_size_LAr_inner / 2.0);
  position_THGEM1_container = G4ThreeVector(0, 0, gPars::det_dims.z_bottom_THGEM1 * mm + z_size_THGEM1 / 2.0 - z_size_LAr_inner / 2.0);
  position_THGEM1_copper_plate = G4ThreeVector(0, 0, 0); //is inside THGEM1_container

  gPars::debugging.THGEM1_hole_center = //x!=0 because x=0 is just across anode wire before SiPM.
      G4ThreeVector(gPars::det_dims.THGEM1_hole_pitch * mm, 0, gPars::det_dims.z_bottom_THGEM1 * mm + gPars::det_dims.THGEM1_width_total / 2.0 * mm);
  gPars::debugging.EL_gap_center =
      G4ThreeVector(gPars::det_dims.THGEM1_hole_pitch * mm, 0, (gPars::det_dims.z_top_interface_grid + gPars::det_dims.z_bottom_THGEM1) / 2.0 * mm);

	position_interface_wire_container = G4ThreeVector(0, 0, gPars::det_dims.z_top_interface_grid * mm - radius_Interface_wire - z_size_LAr_inner / 2.0);
	position_interface_frame = G4ThreeVector(0, 0, gPars::det_dims.z_top_interface_grid * mm - z_size_Interface_grid_substrate / 2.0 - z_size_LAr_inner / 2.0);


	position_Insulator_box = G4ThreeVector(0, 0, z_Insulator_box_center);

	position_LAr_inner = G4ThreeVector(0, 0, z_size_LAr_inner / 2.0);
	position_LAr_outer = G4ThreeVector(0, 0, z_size_LAr_inner / 2.0);

	//FieldWires
	position_FieldWire_bottom1 = G4ThreeVector(x_pos_FieldWire, 0, z_pos_FieldWire_bottom - z_size_LAr_inner / 2.0);
	position_FieldWire_bottom2 = G4ThreeVector(-x_pos_FieldWire, 0, z_pos_FieldWire_bottom - z_size_LAr_inner / 2.0);
	position_FieldWire_bottom3 = G4ThreeVector(0, x_pos_FieldWire, z_pos_FieldWire_bottom - z_size_LAr_inner / 2.0);
	position_FieldWire_bottom4 = G4ThreeVector(0, -x_pos_FieldWire, z_pos_FieldWire_bottom - z_size_LAr_inner / 2.0);
	position_FieldWire_top1 = G4ThreeVector(x_pos_FieldWire, 0, z_pos_FieldWire_top - z_size_LAr_inner / 2.0);
	position_FieldWire_top2 = G4ThreeVector(-x_pos_FieldWire, 0, z_pos_FieldWire_top - z_size_LAr_inner / 2.0);
	position_FieldWire_top3 = G4ThreeVector(0, x_pos_FieldWire, z_pos_FieldWire_top - z_size_LAr_inner / 2.0);
	position_FieldWire_top4 = G4ThreeVector(0, -x_pos_FieldWire, z_pos_FieldWire_top - z_size_LAr_inner / 2.0);

	position_Cathode = G4ThreeVector(0, 0, -z_size_Cathode / 2.0);
	position_LArInactive = G4ThreeVector(0, 0, -z_size_Cathode - z_size_LArInactive / 2.0);
	position_PMMA_bottom = G4ThreeVector(0, 0, PMMA_bottom_center);
	position_Al_window_top = G4ThreeVector(0, 0, Al_window_top_center);
	position_Al_window_bottom = G4ThreeVector(0, 0, Al_window_bottom_center);
	position_CryogenicChamberBottom = G4ThreeVector(0, 0, CryogenicChamberBottom_center);
	position_ExternalColl = G4ThreeVector(0, 0, ExternalColl_center);

	//PMT
	position_PMT_0 = G4ThreeVector(-x_pos_PMT, 0, z_pos_PMT);
	position_PMT_1 = G4ThreeVector(x_pos_PMT, 0, z_pos_PMT);
	position_PMT_2 = G4ThreeVector(0, -y_pos_PMT, z_pos_PMT);
	position_PMT_3 = G4ThreeVector(0, y_pos_PMT, z_pos_PMT);

	//anode_grid
	position_anode_grid = G4ThreeVector(0, 0, z_anode_grid_center);

	//PMTAnodeGridTracker
	const double x_pos_PMTAnodeGridTracker = (x_pos_PMT - z_size_PMT / 2.0 - PMTAnodeGridTrackerThickness / 2.0);
	const double y_pos_PMTAnodeGridTracker = 0;
	const double z_pos_PMTAnodeGridTracker = PMTAnodeGridTrackerGasXSize /2.0 + z_size_LAr_inner;
	position_PMTAnodeGridTrackerGas_1 = G4ThreeVector(x_pos_PMTAnodeGridTracker, y_pos_PMTAnodeGridTracker, z_pos_PMTAnodeGridTracker);
	position_PMTAnodeGridTrackerGasInner_1 = G4ThreeVector(x_pos_PMTAnodeGridTracker - PMTAnodeGridTrackerThickness, y_pos_PMTAnodeGridTracker, z_pos_PMTAnodeGridTracker);

	const double x_pos_PMTAnodeGridTrackerLiquid = x_pos_PMTAnodeGridTracker;
	const double y_pos_PMTAnodeGridTrackerLiquid = 0;
	const double z_pos_PMTAnodeGridTrackerLiquid = z_size_LAr_inner - PMTAnodeGridTrackerLiquidXSize/2.0 - z_size_LAr_inner / 2.0;
	position_PMTAnodeGridTrackerLiquid_1 = G4ThreeVector(x_pos_PMTAnodeGridTrackerLiquid, y_pos_PMTAnodeGridTrackerLiquid, z_pos_PMTAnodeGridTrackerLiquid);
	position_PMTAnodeGridTrackerLiquidInner_1 = G4ThreeVector(x_pos_PMTAnodeGridTrackerLiquid - PMTAnodeGridTrackerThickness, y_pos_PMTAnodeGridTrackerLiquid, z_pos_PMTAnodeGridTrackerLiquid);

	//SteelBox0
	const double xPosSteelBox = x_pos_PMT - z_size_PMT / 2.0 + xSizeSteelBox/2.0;
	positionSteelBox0 = G4ThreeVector(-xPosSteelBox, 0, z_pos_PMT);
	positionSteelBox1 = G4ThreeVector(xPosSteelBox, 0, z_pos_PMT);
	positionSteelBox2 = G4ThreeVector(0, -xPosSteelBox, z_pos_PMT);
	positionSteelBox3 = G4ThreeVector(0, xPosSteelBox, z_pos_PMT);
}

void DetectorConstruction::defineMaterials()
{
	G4NistManager* man = G4NistManager::Instance();

	G4Element *C = man->FindOrBuildElement("C");
	G4Element *H = man->FindOrBuildElement("H");
	G4Element *Si = man->FindOrBuildElement("Si");
	G4Element *O = man->FindOrBuildElement("O");
	G4Element *Sb = man->FindOrBuildElement("Sb");
	G4Element *Rb = man->FindOrBuildElement("Rb");
	G4Element *Cs = man->FindOrBuildElement("Cs");
	G4Element *Lu = man->FindOrBuildElement("Lu");
	G4Element *Y = man->FindOrBuildElement("Y");
	G4Element *Ce = man->FindOrBuildElement("Ce");
	G4Element *La = man->FindOrBuildElement("La");
	G4Element *Br = man->FindOrBuildElement("Br");
	G4Element *Na = man->FindOrBuildElement("Na");
	G4Element *I = man->FindOrBuildElement("I");
	G4Element *Tl = man->FindOrBuildElement("Tl");
	G4Element *Gd = man->FindOrBuildElement("Gd");
	G4Element *Al = man->FindOrBuildElement("Al");
	G4Element *Pr = man->FindOrBuildElement("Pr");
	G4Element *Ar = man->FindOrBuildElement("Ar");
	G4Element *W = man->FindOrBuildElement("W");
	G4Element *Cu = man->FindOrBuildElement("Cu");
	G4Element *Zn = man->FindOrBuildElement("Zn");
	G4Element *Fe = man->FindOrBuildElement("Fe");

	G4double a;  // atomic mass
	G4double z;  // atomic number
	G4double density;
	matAl = new G4Material("Al", z = 13., a = 26.98*g / mole, density = 2.7*g / cm3);
	matFe = new G4Material("Fe", z = 26., a = 55.85*g / mole, density = 7.874 *g / cm3);

	// Air
	G4Material* Air = man->FindOrBuildMaterial("G4_AIR");
	Air->SetName("Air");
	const G4int numentries = 2;
	G4double energies[numentries] = { 0.1*eV, 10.0*eV };
	G4double vacrindices[numentries] = { 1., 1. };
	//G4double airabsorpti[numentries] = { 10*m, 10*m }; // avoid infinite light-paths
	G4double airabsorpti[numentries] = { 10000 * m, 10000 * m };
	G4MaterialPropertiesTable* airprop = new G4MaterialPropertiesTable();
	airprop->AddProperty("ABSLENGTH", energies, airabsorpti, numentries);
	airprop->AddProperty("RINDEX", energies, vacrindices, numentries);
	Air->SetMaterialPropertiesTable(airprop);
	//------------------------------
	// LAr
	G4Material* LAr = new G4Material("LAr", 1.400*g / cm3, 1);
	LAr->AddElement(Ar, 1);
	const G4int LAr_numentries = 12;
	G4double LAr_energies[LAr_numentries] = {1*eV, 2.95*eV, 3.2*eV, 3.4*eV, 3.8*eV, 4*eV, 5.63*eV, 6.89*eV, 7.75*eV, 8.86*eV, 9.69*eV, 10.33*eV};
	G4double LAr_rindices[LAr_numentries] = {1.23, 1.23, 1.23, 1.23, 1.23, 1.23, 1.26, 1.29, 1.31, 1.34, 1.36, 1.45}; //doi:10.1088/1748-0221/15/09/P09009
	//G4double LAr_absorpti[LAr_numentries] = { 2 * m, 2 * m }; // avoid infinite light-paths
	G4double LAr_absorpti[LAr_numentries] = { 20000 * m, 20000 * m };
	G4MaterialPropertiesTable* LAr_prop = new G4MaterialPropertiesTable();
	LAr_prop->AddProperty("ABSLENGTH", energies, LAr_absorpti, numentries);
	LAr_prop->AddProperty("RINDEX", LAr_energies, LAr_rindices, LAr_numentries);
	LAr->SetMaterialPropertiesTable(LAr_prop);
	//------------------------------
	//FR4
	//I don't know a real chemical composition. So it's a dummy
	G4Material* FR4 = new G4Material("FR4", 1.850*g / cm3 /* * (1 - 0.28)*/ /*THGEM transparence*/, 3);
	FR4->AddElement(C, 5); FR4->AddElement(O, 2); FR4->AddElement(H, 8);
	const G4int numentries_FR4 = 2;
	G4double energies_FR4[numentries_FR4] = { 0.1*eV, 10.0*eV };
	G4double rindices_FR4[numentries_FR4] = { 1.5, 1.5 };
	G4double absorpti_FR4[numentries_FR4] = { 1 * um, 1 * um };
	G4MaterialPropertiesTable* prop_FR4 = new G4MaterialPropertiesTable();
	prop_FR4->AddProperty("ABSLENGTH", energies_FR4, absorpti_FR4, numentries_FR4);
	prop_FR4->AddProperty("RINDEX", energies_FR4, rindices_FR4, numentries_FR4);
	FR4->SetMaterialPropertiesTable(prop_FR4);
	//------------------------------
	// PMMA
	G4Material* PMMA = new G4Material("PMMA", 1.18*g / cm3, 3);
	PMMA->AddElement(C, 5);	PMMA->AddElement(O, 2); PMMA->AddElement(H, 8);
	DataVector PMMA_ABSLENGTH(gPars::det_opt.pmma_absorption_length_filename);
	DataVector PMMA_RINDEX(gPars::det_opt.pmma_rindex_filename);
	G4MaterialPropertiesTable* prop_PMMA = new G4MaterialPropertiesTable();
	prop_PMMA->AddProperty("ABSLENGTH", PMMA_ABSLENGTH.get_Xs(eV), PMMA_ABSLENGTH.get_Ys(mm));
	prop_PMMA->AddProperty("RINDEX", PMMA_RINDEX.get_Xs(eV), PMMA_RINDEX.get_Ys());
	PMMA->SetMaterialPropertiesTable(prop_PMMA);
	//prop_PMMA->DumpTable();
	//------------------------------
	// TPB
	G4Material* materialTPB = new G4Material("TPB", 1.18*g / cm3, 3);
	materialTPB->AddElement(C, 5); materialTPB->AddElement(O, 2); materialTPB->AddElement(H, 8);
	const int numentries_TPB = 6;
	G4double photonEnergyTPB[] = { 0.1*eV, 1240.0 / 411.0 *eV, 1240.0 / 409.0 *eV, 1240.0 / 401.0 *eV, 1240.0 / 399.0 *eV, 10.0*eV };
	G4double RIndexTPB[] = { 1.23, 1.23, 1.23, 1.23, 1.23, 1.23 };
	G4double AbsTPB[] = { 10 * m, 10.0 * m, 10.0 * m, 0.001 * mm, 0.001 * mm, 10 * m };
	G4double EmissionTPB[] = { 0.0, 0.0, 1.0, 0.0, 0.0, 0.0 };
	G4MaterialPropertiesTable* tableTPB = new G4MaterialPropertiesTable();
	tableTPB->AddProperty("RINDEX", photonEnergyTPB, RIndexTPB, numentries_TPB);
	tableTPB->AddProperty("WLSABSLENGTH", photonEnergyTPB, AbsTPB, numentries_TPB);
	tableTPB->AddProperty("WLSCOMPONENT", photonEnergyTPB, EmissionTPB, numentries_TPB);
	tableTPB->AddConstProperty("WLSTIMECONSTANT", 0.5*ns);
	materialTPB->SetMaterialPropertiesTable(tableTPB);
	//------------------------------
	// PMMA_UV
	G4Material* PMMA_UV = new G4Material("PMMA_UV", 1.18*g / cm3, 3);
	PMMA_UV->AddElement(C, 5);
	PMMA_UV->AddElement(O, 2);
	PMMA_UV->AddElement(H, 8);
	DataVector PMMA_UV_ABSLENGTH(gPars::det_opt.pmma_uv_absorption_length_filename);
	const G4int numentries_PMMA_UV = 2;
	G4double energies_PMMA_UV[numentries_PMMA_UV] = { 0.1*eV, 10.0*eV };
	G4double rindices_PMMA_UV[numentries_PMMA_UV] = { 1.5, 1.5 };
	G4MaterialPropertiesTable* prop_PMMA_UV = new G4MaterialPropertiesTable();
	prop_PMMA_UV->AddProperty("ABSLENGTH", PMMA_UV_ABSLENGTH.get_Xs(eV), PMMA_UV_ABSLENGTH.get_Ys(mm));
	prop_PMMA_UV->AddProperty("RINDEX", energies_PMMA_UV, rindices_PMMA_UV, numentries_PMMA_UV);
	PMMA_UV->SetMaterialPropertiesTable(prop_PMMA_UV);
}
