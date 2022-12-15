#include <geant4/detector/Detector_full.hh>

Detector_full::Detector_full() :
  VDetectorConstruction()
{}

Detector_full::~Detector_full()
{}

G4VPhysicalVolume * Detector_full::Construct()
{
  SetSizeAndPosition();
  defineMaterials();
  defineSurfaces();
	//-------------------------------------------------------------------------------

	solidWorld = new G4Box("sworld", HalfWorldLength, HalfWorldLength, HalfWorldLength);
	logicWorld = new G4LogicalVolume(solidWorld, matGas, "lWorld", 0, 0, 0);
	//  Must place the World Physical volume unrotated at (0,0,0).
	physiWorld = new G4PVPlacement(0, G4ThreeVector(), logicWorld, "pWorld", 0, // its mother  volume
		false, 0);
	// Set user cuts to avoid deadlocks
	G4double maxStep = 10.0*m, maxLength = 10.0*m, maxTime = 100.0*ns, minEkin = 0.5*eV;
	logicWorld->SetUserLimits(new G4UserLimits(maxStep, maxLength, maxTime, minEkin));

	//-------------------------------------------------------------------------------
	//create LAr box contained inside the PMMA insulator
	G4Box* solid_LAr_inner = new G4Box("solid_LAr_inner", x_size_LAr_inner / 2.0, y_size_LAr_inner / 2.0, z_size_LAr_inner / 2.0);
	G4LogicalVolume* logic_LAr_inner = new G4LogicalVolume(solid_LAr_inner, matLAr, "logic_LAr_inner", 0, 0, 0);
	G4VPhysicalVolume* phys_LAr_inner = new G4PVPlacement(0, position_LAr_inner, logic_LAr_inner,  "phys_LAr_inner",
		logicWorld, false, 0, fCheckOverlaps);

	//--------------------------------------------------------------------------------
	//create LAr box contained outside the PMMA insulator
	G4Box* solid_LAr_outer_out = new G4Box("solid_LAr_outer_out", x_size_LAr_outer_out / 2.0, y_size_LAr_outer_out / 2.0, z_size_LAr_outer / 2.0);
	G4Box* solid_LAr_outer_in = new G4Box("solid_LAr_outer_in", x_size_LAr_outer_in / 2.0, y_size_LAr_outer_in / 2.0, z_size_LAr_outer / 2.0);
	G4SubtractionSolid* solid_LAr_outer = new G4SubtractionSolid("solid_LAr_outer", solid_LAr_outer_out, solid_LAr_outer_in);
	G4LogicalVolume* logic_LAr_outer = new G4LogicalVolume(solid_LAr_outer, matLAr, "logic_LAr_outer", 0, 0, 0);
	G4VPhysicalVolume* phys_LAr_outer = new G4PVPlacement(0, position_LAr_outer, logic_LAr_outer, "phys_LAr_outer",
		logicWorld, false, 0, fCheckOverlaps);

	//--------------------------------------------------------------------------------
	//create FieldWires
	G4Tubs* solid_FieldWire = new G4Tubs("solid_FieldWire", 0, radius_FieldWire, length_FieldWire / 2.0, 0.*deg, 360.*deg);
	G4LogicalVolume* logic_FieldWire = new G4LogicalVolume(solid_FieldWire, matAl, "logic_FieldWire", 0, 0, 0);

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
	G4VPhysicalVolume* phys_THGEM1_cell = new G4PVPlacement(0, position_SingleTHGEM1Cell, logic_THGEM1_cell,
		"phys_THGEM1_cell_isolated", logicWorld, false, 0, fCheckOverlaps);

	//-------------------------------------------------------------------------------
	//create bCathode
	G4Box* solid_bCathode = new G4Box("solid_bCathode", x_size_Cathode / 2.0, y_size_Cathode / 2.0, z_size_Cathode / 2.0);
	G4LogicalVolume* logic_bCathode = new G4LogicalVolume(solid_bCathode, matFR4, "logic_bCathode", 0, 0, 0);
	G4VPhysicalVolume* phys_bCathode = new G4PVPlacement(0, position_Cathode, logic_bCathode, "phys_bCathode",
		logicWorld, false, 0, fCheckOverlaps);

	//-------------------------------------------------------------------------------
	//Inactive LAr layer below the cathode. It is considered that in reality there is gas there. No effect on light collection.
	G4Box* solid_LArInactive = new G4Box("solid_LArInactive", x_size_LArInactive / 2.0, y_size_LArInactive / 2.0, z_size_LArInactive / 2.0);
	G4LogicalVolume* logic_LArInactive = new G4LogicalVolume(solid_LArInactive, matGas/*"matLAr"*/, "logic_LArInactive", 0, 0, 0);
	G4VPhysicalVolume* phys_LArInactive = new G4PVPlacement(0, position_LArInactive, logic_LArInactive, "phys_LArInactive",
		logicWorld, false, 0, fCheckOverlaps);

	//-------------------------------------------------------------------------------
	// Create bottom PMMA plate below the cathode. No effect on light collection.
	G4Box* solid_PMMA_bottom = new G4Box("solid_PMMA_bottom", x_size_PMMA_bottom / 2.0, y_size_PMMA_bottom / 2.0, z_size_PMMA_bottom / 2.0);
	G4LogicalVolume* logic_PMMA_bottom = new G4LogicalVolume(solid_PMMA_bottom, matPMMA, "logic_PMMA_bottom", 0, 0, 0);
	G4VPhysicalVolume* phys_PMMA_bottom = new G4PVPlacement(0, position_PMMA_bottom, logic_PMMA_bottom, "phys_PMMA_bottom",
		logicWorld, false, 0, fCheckOverlaps);

	//-------------------------------------------------------------------------------
	// Create Al window. No effect on light collection.
	G4Tubs* solid_Al_window = new G4Tubs("solid_Al_window", 0, diameter_size_Al_window / 2.0, z_size_Al_window / 2.0, 0. *deg, 360.*deg);
	G4LogicalVolume* logic_Al_window = new G4LogicalVolume(solid_Al_window, matAl, "logic_Al_window", 0, 0, 0);
	G4VPhysicalVolume* phys_Al_window_top = new G4PVPlacement(0, position_Al_window_top, logic_Al_window, "phys_Al_window_top",
		logicWorld, false, 0, fCheckOverlaps);
	G4VPhysicalVolume* phys_Al_window_bottom = new G4PVPlacement(0, position_Al_window_bottom, logic_Al_window, "phys_Al_window_bottom",
		logicWorld, false, 0, fCheckOverlaps);

	//-------------------------------------------------------------------------------
	// Create chamber bottom. No effect on light collection.
	G4Tubs* solid_CryogenicChamberBottom = new G4Tubs("solid_CryogenicChamberBottom", diameter_size_internal_CryogenicChamberBottom / 2.0, diameter_size_external_CryogenicChamberBottom / 2.0, z_size_CryogenicChamberBottom / 2.0, 0. *deg, 360.*deg);
	G4LogicalVolume* logic_CryogenicChamberBottom = new G4LogicalVolume(solid_CryogenicChamberBottom, matFe, "logic_CryogenicChamberBottom", 0, 0, 0);
	G4VPhysicalVolume* phys_CryogenicChamberBottom = new G4PVPlacement(0, position_CryogenicChamberBottom, logic_CryogenicChamberBottom, "phys_CryogenicChamberBottom",
		logicWorld, false, 0, fCheckOverlaps);

		//-------------------------------------------------------------------------------
	// Create interface grid. Parameterization is required.
	G4Box* solid_tracker_Interface_grid = new G4Box("solid_tracker_Interface_grid", x_size_tracker_Interface_grid / 2.0, y_size_tracker_Interface_grid / 2.0, radius_Interface_wire);
	G4LogicalVolume* logic_tracker_Interface_grid = new G4LogicalVolume(solid_tracker_Interface_grid, matLAr, "logic_tracker_Interface_grid", 0, 0, 0);
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
	G4LogicalVolume* logic_interface_grid = new G4LogicalVolume(solid_interface_grid_subtraction, matFR4, "l_anode_grid", 0, 0, 0);
	G4VPhysicalVolume* phys_interface_grid = new G4PVPlacement(0, position_interface_frame, logic_interface_grid, "p_interface_grid",
		logic_LAr_inner, false, 0);

	// Create anode grid. Parameterisation is required.
	// Create tracker (container for parameterised volumes)
	G4Box* solid_tracker_anode_grid = new G4Box("solid_tracker_anode_grid", x_size_tracker_anode_grid / 2.0, y_size_tracker_anode_grid / 2.0, z_size_tracker_anode_grid / 2.0);
	G4LogicalVolume* logic_tracker_anode_grid = new G4LogicalVolume(solid_tracker_anode_grid, matGas, "logic_tracker_anode_grid", 0, 0, 0);
	G4VPhysicalVolume* phys_tracker_anode_grid = new G4PVPlacement(0, position_anode_grid, logic_tracker_anode_grid, "phys_tracker_anode_grid",
		logicWorld, false, 0, fCheckOverlaps);
	// Create anode grid frame
	G4Box* solid_anode_grid_substrate = new G4Box("anode_grid_substrate", size_anode_grid / 2.0, size_anode_grid / 2.0, thickness_anode_grid / 2.0);
	G4Box* solid_anode_grid_hole = new G4Box("anode_grid_hole", size_anode_grid_hole / 2.0, size_anode_grid_hole / 2.0, thickness_anode_grid / 2.0);
	G4SubtractionSolid* solid_anode_grid_frame = new G4SubtractionSolid("anode_grid__frame", solid_anode_grid_substrate, solid_anode_grid_hole);
	G4LogicalVolume* logic_anode_grid_frame = new G4LogicalVolume(solid_anode_grid_frame, matFR4, "l_anode_grid", 0, 0, 0);
	G4VPhysicalVolume* phys_anode_grid_frame = new G4PVPlacement(0, position_anode_grid, logic_anode_grid_frame, "p_anode_grid",
		logicWorld, false, 0);
	// Create anode wire grid
	G4Tubs* solid_wire = new G4Tubs("solid_wire", 0, radius_wire, length_wire / 2.0, 0.*deg, 360.*deg);
	G4LogicalVolume* logic_anode_wire = new G4LogicalVolume(solid_wire, matAl, "lwire", 0, 0, 0);
	G4VPVParameterisation* param_wire = new DetectorParameterisation(N_wire, 1, 1, rotX_90, G4ThreeVector(0,0,0), G4ThreeVector(step_wire, 0, 0));
	G4VPhysicalVolume* phys_wire = new G4PVParameterised("phys_wire", logic_anode_wire,	logic_tracker_anode_grid,
		kXAxis, N_wire, param_wire, fCheckOverlaps);

	//--------------------------------------------------------------------------------
	// Create PMT's wire grid
	// Gas part
	G4Box* solid_PMTAnodeGridTrackerGas = new G4Box("solid_PMTAnodeGridTrackerGas", PMTAnodeGridTrackerGasXSize / 2.0, PMTAnodeGridTrackerGasYSize / 2.0, PMTAnodeGridTrackerThickness / 2.0);
	G4LogicalVolume* logic_PMTAnodeGridTrackerGas = new G4LogicalVolume(solid_PMTAnodeGridTrackerGas, matGas, "logic_PMTAnodeGridTrackerGas", 0, 0, 0);
	G4LogicalVolume* logic_PMTAnodeGridTrackerGasInner = new G4LogicalVolume(solid_PMTAnodeGridTrackerGas, matGas, "logic_PMTAnodeGridTrackerGasInner", 0, 0, 0);

	G4VPhysicalVolume* phys_PMTAnodeGridTrackerGas_1 = new G4PVPlacement(rotY_90, position_PMTAnodeGridTrackerGas_1,
		logic_PMTAnodeGridTrackerGas,	"phys_PMTAnodeGridTrackerGas_1", logicWorld, false, 0, fCheckOverlaps);

	G4VPhysicalVolume* phys_PMTAnodeGridTrackerGasInner_1 = new G4PVPlacement(rotY_90, position_PMTAnodeGridTrackerGasInner_1,
		logic_PMTAnodeGridTrackerGasInner, "phys_PMTAnodeGridTrackerGasInner_1", logicWorld, false, 0, fCheckOverlaps);

	// Liquid part
	G4Box* solid_PMTAnodeGridTrackerLiquid = new G4Box("solid_PMTAnodeGridTrackerLiquid", PMTAnodeGridTrackerLiquidXSize / 2.0, PMTAnodeGridTrackerLiquidYSize / 2.0, PMTAnodeGridTrackerThickness / 2.0);
	G4LogicalVolume* logic_PMTAnodeGridTrackerLiquid = new G4LogicalVolume(solid_PMTAnodeGridTrackerLiquid, matLAr, "logic_PMTAnodeGridTrackerLiquid", 0, 0, 0);
	G4LogicalVolume* logic_PMTAnodeGridTrackerLiquidInner = new G4LogicalVolume(solid_PMTAnodeGridTrackerLiquid, matLAr, "logic_PMTAnodeGridTrackerLiquidInner", 0, 0, 0);
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
	G4LogicalVolume* trackerLV = new G4LogicalVolume(trackerS, matGas, "Tracker", 0, 0, 0);
	G4VPhysicalVolume* phys_tracker = new G4PVPlacement(0, position_SiPM_container, trackerLV, "Tracker", logicWorld,
	 	false, 0, fCheckOverlaps);

	G4Box* solid_SiPM = new G4Box("sscintillator", SiPM_size / 2.0, SiPM_size / 2.0, thickness_SiPM / 2.0);
	logic_SiPM = new G4LogicalVolume(solid_SiPM, matAl, "lSiPM", 0, 0, 0);
	G4VPVParameterisation* chamberParam = new DetectorParameterisation(Nx_SiPMs, Ny_SiPMs, 1, NULL, G4ThreeVector(0,0,0), G4ThreeVector(SiPM_pitch, SiPM_pitch, 0));
	G4VPhysicalVolume* phys_SiPM = new G4PVParameterised(gPars::det_dims->SiPM_device_name, logic_SiPM, trackerLV,
		kXAxis, Nx_SiPMs * Ny_SiPMs, chamberParam, fCheckOverlaps);

	G4Box* solid_SiPMFR4 = new G4Box("solid_SiPMFR4", x_size_tracker / 2.0, y_size_tracker / 2.0, z_size_SiPMFR4 / 2.0);
	G4LogicalVolume* logic_SiPMFR4 = new G4LogicalVolume(solid_SiPMFR4, matGas, "logic_SiPMFR4", 0, 0, 0);
	G4VPhysicalVolume* phys_SiPMFR4 = new G4PVPlacement(0, position_SiPMFR4, logic_SiPMFR4, "phys_SiPMFR4",
		logicWorld, false, 0, fCheckOverlaps);

	//--------------------------------------------------------------------------------
	// Create PMMA plate
	G4Box* solid_PMMA_plate = new G4Box("solid_tracker_anode_grid", x_size_PMMA_plate / 2.0, y_size_PMMA_plate / 2.0, z_size_PMMA_plate / 2.0);
	G4LogicalVolume* logic_PMMA_plate = new G4LogicalVolume(solid_PMMA_plate, matPMMA, "logic_PMMA_plate", 0, 0, 0);
	G4VPhysicalVolume* phys_PMMA_plate = new G4PVPlacement(0, position_PMMA_plate, logic_PMMA_plate, "phys_PMMA_plate",
		logicWorld, false, 0, fCheckOverlaps);

	//--------------------------------------------------------------------------------
	// Create THGEM1
	G4Box* solid_whole_THGEM1 = new G4Box("solid_whole_THGEM1", x_size_THGEM1 / 2.0, y_size_THGEM1 / 2.0, z_size_THGEM1 / 2.0);
	G4Box* solid_active_THGEM1 = new G4Box("solid_whole_THGEM1", x_size_THGEM1_container / 2.0, y_size_THGEM1_container / 2.0, z_size_THGEM1 / 2.0);
	G4SubtractionSolid* solid_THGEM1_frame = new G4SubtractionSolid("solid_THGEM1_frame", solid_whole_THGEM1, solid_active_THGEM1);
  G4LogicalVolume* logic_THGEM1_frame = new G4LogicalVolume(solid_THGEM1_frame, matFR4, "logic_THGEM1_frame", 0, 0, 0);
  G4VPhysicalVolume* phys_THGEM1_frame = new G4PVPlacement(0, position_THGEM1_frame, logic_THGEM1_frame, "phys_THGEM1_frame",
      logic_LAr_inner, false, 0, fCheckOverlaps);

  G4Box* solid_THGEM1_container = new G4Box("solid_THGEM1_container", x_size_THGEM1_container / 2.0, y_size_THGEM1_container / 2.0, z_size_THGEM1_container / 2.0);
  G4LogicalVolume* logic_THGEM1_container = new G4LogicalVolume(solid_THGEM1_container, matLAr, "logic_THGEM1_container", 0, 0, 0);
  phys_THGEM1_container = new G4PVPlacement(0, position_THGEM1_container, logic_THGEM1_container, "phys_THGEM1_cell_container",
        logic_LAr_inner, false, 0, fCheckOverlaps);

  G4LogicalVolume* logic_THGEM1_copper = new G4LogicalVolume(solid_active_THGEM1, matFR4, "logic_THGEM1_copper", 0, 0, 0);
  G4VPhysicalVolume* phys_THGEM1_copper = new G4PVPlacement(0, position_THGEM1_copper_plate, logic_THGEM1_copper, "phys_THGEM1_copper",
      logic_THGEM1_container, false, 0, fCheckOverlaps);

	//--------------------------------------------------------------------------------
	// Create insulator box
	G4Box* solid_Insulator_box_inner = new G4Box("solid_Insulator_box_inner", x_size_Insulator_box_inner / 2.0, y_size_Insulator_box_inner / 2.0, z_size_Insulator_box / 2.0);
	G4Box* solid_Insulator_box_outer = new G4Box("solid_Insulator_box_outer", x_size_Insulator_box_outer / 2.0, y_size_Insulator_box_outer / 2.0, z_size_Insulator_box / 2.0);
	G4SubtractionSolid* solid_Insulator_box_subtraction = new G4SubtractionSolid("solid_Insulator_box_subtraction", solid_Insulator_box_outer, solid_Insulator_box_inner);
	G4LogicalVolume* logic_Insulator_box = new G4LogicalVolume(solid_Insulator_box_subtraction, matPMMA_UV, "logic_Insulator_box", 0, 0, 0);
	G4VPhysicalVolume* phys_Insulator_box = new G4PVPlacement(0, position_Insulator_box, logic_Insulator_box,
		"phys_Insulator_box", logicWorld, false, 0, fCheckOverlaps);

	//--------------------------------------------------------------------------------
	// Create PMTs
	G4Tubs* solid_PMT = new G4Tubs("solid_PMT", 0, radius_PMT, z_size_PMT / 2.0, 0.*deg, 360.*deg);
	logic_PMT = new G4LogicalVolume(solid_PMT, matAl, "logic_PMT", 0, 0, 0);

	G4VPhysicalVolume* phys_PMT0 = new G4PVPlacement(rotY_90, position_PMT_0, logic_PMT, gPars::det_dims->PMT_device_name,
		logicWorld, false, 0, fCheckOverlaps);
	G4VPhysicalVolume* phys_PMT1 = new G4PVPlacement(rotY_90, position_PMT_1, logic_PMT, gPars::det_dims->PMT_device_name,
		logicWorld, false, 1, fCheckOverlaps);
	G4VPhysicalVolume* phys_PMT2 = new G4PVPlacement(rotX_90, position_PMT_2, logic_PMT, gPars::det_dims->PMT_device_name,
		logicWorld, false, 2, fCheckOverlaps);
	G4VPhysicalVolume* phys_PMT3 = new G4PVPlacement(rotX_90, position_PMT_3, logic_PMT, gPars::det_dims->PMT_device_name,
		logicWorld, false, 3, fCheckOverlaps);

	//--------------------------------------------------------------------------------
	// Create PMTs' steel box
	G4Box* solidSteelBox = new G4Box("solidSteelBox", xSizeSteelBox / 2.0, ySizeSteelBox / 2.0, zSizeSteelBox / 2.0);
	G4VSolid* solidBoxSubtractPMT = new G4SubtractionSolid("solidBoxSubtractPMT", solidSteelBox, solid_PMT, rotY_90, G4ThreeVector((xSizeSteelBox - z_size_PMT)/2, 0., 0.));
	G4LogicalVolume* logicSteelBox = new G4LogicalVolume(solidBoxSubtractPMT, matGas, "logicSteelBox", 0, 0, 0);

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
	G4LogicalBorderSurface* LAr_physiWorld2inner = new G4LogicalBorderSurface("LAr_physiWorld2inner", physiWorld, phys_LAr_inner, LAr_OpticalSurface);
	G4LogicalSkinSurface* SiPM = new G4LogicalSkinSurface("SiPM_surf", logic_SiPM, SiPM_OpticalSurface);
	G4LogicalSkinSurface* SiPMFR4 = new G4LogicalSkinSurface("SiPMFR4", logic_SiPMFR4, FR4_unified);

	// PMT photocathode has same surface for LAr and gas.
	G4LogicalSkinSurface* sur_PMT_cathode = new G4LogicalSkinSurface("PMT_cathode", logic_PMT, PMT_OpticalSurface);
	G4LogicalSkinSurface* sur_BottomCathode = new G4LogicalSkinSurface("BottomCathode", logic_bCathode, Cu_cathode);
	G4LogicalSkinSurface* CuReflector_THGEM0_surface = new G4LogicalSkinSurface("CuReflector_THGEM_surface", logic_THGEM1_copper, Cu_THGEM);
	G4LogicalSkinSurface* sur_SteelBox = new G4LogicalSkinSurface("SteelBox", logicSteelBox, stainlessSteel);

	// PMT grid
	G4LogicalBorderSurface* PMTAnodeGridTrackerLiquidGas = new G4LogicalBorderSurface("PMTAnodeGridTrackerLiquid-Gas", phys_PMTAnodeGridTrackerLiquid_1, phys_PMTAnodeGridTrackerGas_1, LAr_OpticalSurface);
	G4LogicalBorderSurface* PMTAnodeGridTrackerLiquidGas1 = new G4LogicalBorderSurface("PMTAnodeGridTrackerLiquid-Gas1", phys_PMTAnodeGridTrackerGas_1, phys_PMTAnodeGridTrackerLiquid_1, LAr_OpticalSurface);
	G4LogicalBorderSurface* PMTAnodeGridTrackerLiquidGasInner = new G4LogicalBorderSurface("PMTAnodeGridTrackerGas-Liquid", phys_PMTAnodeGridTrackerLiquidInner_1, phys_PMTAnodeGridTrackerGasInner_1, LAr_OpticalSurface);
  G4LogicalBorderSurface* PMTAnodeGridTrackerLiquidGasInner1 = new G4LogicalBorderSurface("PMTAnodeGridTrackerGas-Liquid1", phys_PMTAnodeGridTrackerGasInner_1, phys_PMTAnodeGridTrackerLiquidInner_1, LAr_OpticalSurface);
	G4LogicalSkinSurface* PMTGridWire0 = new G4LogicalSkinSurface("PMTGridWire_surface0", logic_PMTGridWire, stainlessSteel);
	G4LogicalSkinSurface* PMTGridWire1 = new G4LogicalSkinSurface("PMTGridWire_surface1", logic_PMTGridWireGasInner, stainlessSteel);
	G4LogicalSkinSurface* PMTGridWire2 = new G4LogicalSkinSurface("PMTGridWire_surface2", logic_PMTGridWireLiquidInner, stainlessSteel);

	// SiPM (anode) grid
	G4LogicalSkinSurface* anode_grid_frame = new G4LogicalSkinSurface("anode_grid_frame_sur", logic_anode_grid_frame, FR4_unified);
	G4LogicalSkinSurface* anode_grid_wire = new G4LogicalSkinSurface("anode_grid_wire_sur", logic_anode_wire, Anode_wire_unified);

	G4LogicalSkinSurface* PMMA_plate = new G4LogicalSkinSurface("PMMA_plate", logic_PMMA_plate, PMMA_OpticalSurface);
	G4LogicalSkinSurface* PMMA_box = new G4LogicalSkinSurface("PMMA_box", logic_Insulator_box, PMMA_OpticalSurface);
	// Interface grid
	G4LogicalBorderSurface* phys_LAr_inner2Interface_grid = new G4LogicalBorderSurface("phys_LAr_inner2Interface_grid", phys_LAr_inner, phys_interface_grid, /*AbsorberMaterial*/ FR4_unified);
	G4LogicalBorderSurface* tracker_Interface_grid2wire = new G4LogicalBorderSurface("tracker_anode_grid2wire", phys_tracker_Interface_grid, phys_Interface_wire, /*AbsorberMaterial*/ Anode_wire_unified);
	G4LogicalBorderSurface* tracker_Interface_grid2Interface_grid = new G4LogicalBorderSurface("tracker_Interface_grid2Interface_grid", phys_tracker_Interface_grid, phys_interface_grid, /*AbsorberMaterial*/ FR4_unified);

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
	G4VisAttributes LAr_VisAtt(G4Colour(0.6, 0.6, 1.0, 0.0));
	G4VisAttributes gas_VisAtt(G4Colour(1.0, 1.0, 1.0, 0.0));
	gas_VisAtt.SetVisibility(false);
	G4VisAttributes PMMA_VisAtt(G4Colour(0.2, 0.2, 0.8, 0.1));
	G4VisAttributes Wires_VisAtt(G4Colour(1.0, 0.1, 0.1, 0.4));
	G4VisAttributes FR4_VisAtt(G4Colour(1.0, 0.5, 0.2, 0.4));
	G4VisAttributes Sensor_VisAtt(G4Colour(0.8, 0.7, 0.2, 0.8));
	G4VisAttributes Cu_VisAtt(G4Colour(0.0, 1.0, 1.0, 0.4));
	G4VisAttributes Steel_VisAtt(G4Colour(0.2, 0.2, 0.2, 0.1));
	G4VisAttributes Al_VisAtt(G4Colour(1.0, 1.0, 0.0, 0.0));
	Al_VisAtt.SetVisibility(false);

	// SiPMs
	logic_anode_wire->SetVisAttributes(Wires_VisAtt);
	logic_anode_grid_frame->SetVisAttributes(FR4_VisAtt);
	logic_SiPM->SetVisAttributes(Sensor_VisAtt);
  trackerLV->SetVisAttributes(gas_VisAtt);
  FR4_VisAtt.SetVisibility(false);
  logic_SiPMFR4->SetVisAttributes(FR4_VisAtt);
  FR4_VisAtt.SetVisibility(true);

	// Interface grid
	logic_Interface_wire->SetVisAttributes(Wires_VisAtt);
	logic_interface_grid->SetVisAttributes(FR4_VisAtt);

	// PMTs
	logic_PMT->SetVisAttributes(Sensor_VisAtt);
	logic_PMTGridWire->SetVisAttributes(Wires_VisAtt);
	logic_PMTGridWireGasInner->SetVisAttributes(gas_VisAtt);
	LAr_VisAtt.SetVisibility(false);
	logic_PMTGridWireLiquidInner->SetVisAttributes(LAr_VisAtt);
	LAr_VisAtt.SetVisibility(true);
	Steel_VisAtt.SetVisibility(false);
	logicSteelBox->SetVisAttributes(Steel_VisAtt);
	Steel_VisAtt.SetVisibility(true);

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
	logic_Al_window->SetVisAttributes(Al_VisAtt);
	Steel_VisAtt.SetVisibility(false);
	logic_CryogenicChamberBottom->SetVisAttributes(Steel_VisAtt);
	Steel_VisAtt.SetVisibility(true);
	PMMA_VisAtt.SetVisibility(false);
	logic_PMMA_bottom->SetVisAttributes(PMMA_VisAtt);
	PMMA_VisAtt.SetVisibility(true);
	logic_PMMA_plate->SetVisAttributes(PMMA_VisAtt);
	logic_Insulator_box->SetVisAttributes(PMMA_VisAtt);

	logic_FieldWire->SetVisAttributes(Wires_VisAtt);

	logicWorld->SetVisAttributes(gas_VisAtt);

	SetupTHGEMsMapping();
	return physiWorld;
}

void Detector_full::ConstructSDandField()
{
  DetectorSensor *thePhotoDetector = new DetectorSensor("/detector/sensitiveDetector");
  G4SDManager::GetSDMpointer()->AddNewDetector(thePhotoDetector);
  SetSensitiveDetector(logic_SiPM, thePhotoDetector); // takes ownership of DetectorSensor*
  SetSensitiveDetector(logic_PMT, thePhotoDetector);
}

void Detector_full::SetSizeAndPosition()
{
	HalfWorldLength = 17 * cm;

	double width_interface_grid_support = 1.4 * mm; // Interface grid frame width at support pillars
  double width_interface_grid_frame = 5.0 * mm; // Interface grid frame full width (there are holes drilled for pillars)
  double LAr_drift_width = 48.0 * mm; // Distance between cathode top and interface grid bottom
  double max_EL_gap_thickness = 22.0 * mm; // Distance between interface grid top and THGEM1 bottom
  double THGEM_cathode_width = 0.5 * mm;
  double THGEM1_active_area_size = 100 * mm;

  int n_SiPMs_rows = 5; // total number = n_SiPMs_rows^2
  SiPM_size = 6 * mm; // Full size of SiPM active area
  SiPM_pitch = 10 * mm; // Distance between SiPM centers
  double EL_gap_thickness = -4 * mm; // Must be negative (single phase). From THGEM1 real bottom
  double z_top_interface_grid = LAr_drift_width + width_interface_grid_support;
  double z_bottom_THGEM1 =z_top_interface_grid + max_EL_gap_thickness; //=71.4

	//PMTs
	radius_PMT = 45 * mm / 2.0;
	z_size_PMT = 2 * mm;
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
	size_anode_grid = 127 * mm ;//see Download:\DetectorPhotos\2021\THGEM_Electroconnect
	size_anode_grid_hole = length_wire;
	z_anode_grid_bottom = z_bottom_THGEM1 + gPars::det_dims->THGEM1_width_total + 5 * mm;
	double z_anode_grid_center = z_anode_grid_bottom + thickness_anode_grid / 2.0;

	//PMMA plate
	x_size_PMMA_plate = size_anode_grid;
	y_size_PMMA_plate = size_anode_grid;
	z_size_PMMA_plate = 1.5 * mm;
	double z_PMMA_plate_center = z_anode_grid_center + thickness_anode_grid / 2.0 + z_size_PMMA_plate / 2.0;

	//SiPMs
	Nx_SiPMs = n_SiPMs_rows;
	Ny_SiPMs = n_SiPMs_rows;
	thickness_SiPM = 1 * nm;
	double z_SiPM_bottom = z_anode_grid_bottom + thickness_anode_grid + z_size_PMMA_plate + (0.1*mm /*small gap between PMMA and SiPM*/);
	double z_SiPM_center = z_SiPM_bottom + thickness_SiPM / 2.0;
	z_size_SiPMFR4 = 2*mm;

	//tracker SiPM
	x_size_tracker = Nx_SiPMs * SiPM_pitch + SiPM_size / 2.0;
	y_size_tracker = Ny_SiPMs * SiPM_pitch + SiPM_size / 2.0;
	z_size_tracker = 0.1 * mm;

	//tracker Anode_grid
	x_size_tracker_anode_grid = length_wire;
	y_size_tracker_anode_grid = x_size_tracker_anode_grid;
	z_size_tracker_anode_grid = thickness_anode_grid;

	//tracker THGEM2 (active region with holes)
	x_size_tracker_THGEM2 = 100 * mm;
	y_size_tracker_THGEM2 = 100 * mm;
	z_size_tracker_THGEM2 = gPars::det_dims->THGEM1_width_total;

	// THGEM1
	x_size_THGEM1 = size_anode_grid;
  y_size_THGEM1 = size_anode_grid;
  z_size_THGEM1 = gPars::det_dims->THGEM1_width_total;
  x_size_THGEM1_container = THGEM1_active_area_size;
  y_size_THGEM1_container = THGEM1_active_area_size;
  z_size_THGEM1_container = gPars::det_dims->THGEM1_container_width;

	//Interface_grid
	x_size_tracker_Interface_grid = x_size_tracker_THGEM2;
	y_size_tracker_Interface_grid = y_size_tracker_THGEM2;
	x_size_Interface_grid_substrate = size_anode_grid;
	y_size_Interface_grid_substrate = size_anode_grid;
	z_size_Interface_grid_substrate = width_interface_grid_frame;

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
	z_size_LAr_inner = z_top_interface_grid + (max_EL_gap_thickness - EL_gap_thickness);

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
	z_size_Cathode = THGEM_cathode_width;

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
	z_size_TPB = 0.4 * mm; // See Borisova's dissertation

	position_SingleTHGEM1Cell = G4ThreeVector(150 * mm, 150 * mm, 150 * mm);

	position_anode_grid = G4ThreeVector(0, 0, z_anode_grid_center);
	position_SiPM_container = G4ThreeVector(0, 0, z_SiPM_center);
	position_PMMA_plate = G4ThreeVector(0, 0, z_PMMA_plate_center);
	position_SiPMFR4 = G4ThreeVector(0, 0, z_SiPM_center + z_size_tracker /2.0 + z_size_SiPMFR4/2.0);

	position_THGEM1_frame = G4ThreeVector(0, 0, z_bottom_THGEM1 + z_size_THGEM1 / 2.0 - z_size_LAr_inner / 2.0);
  position_THGEM1_container = G4ThreeVector(0, 0, z_bottom_THGEM1 + z_size_THGEM1 / 2.0 - z_size_LAr_inner / 2.0);
  position_THGEM1_copper_plate = G4ThreeVector(0, 0, 0); //is inside THGEM1_container

	position_interface_wire_container = G4ThreeVector(0, 0, z_top_interface_grid - radius_Interface_wire - z_size_LAr_inner / 2.0);
	position_interface_frame = G4ThreeVector(0, 0, z_top_interface_grid - z_size_Interface_grid_substrate / 2.0 - z_size_LAr_inner / 2.0);


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
	const double z_pos_PMTAnodeGridTracker = PMTAnodeGridTrackerGasXSize / 2.0 + z_size_LAr_inner;
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

	if ((-EL_gap_thickness < (gPars::det_dims->THGEM1_width_total + (gPars::det_dims->THGEM1_container_width - gPars::det_dims->THGEM1_width_total)*0.5)) &&
      (-EL_gap_thickness > (-(gPars::det_dims->THGEM1_container_width - gPars::det_dims->THGEM1_width_total)*0.5))) {
    G4Exception("gPars::InitGlobals(): ",
          "InvalidSetup", FatalException, "LAr level intersects with THGEM1 (or its container)! This case is not supported.");
    return;
  }

	gPars::det_dims->THGEM_hole_center = //x!=0 because x=0 is just across anode wire before SiPM.
      G4ThreeVector(gPars::det_dims->THGEM1_hole_pitch, 0, z_bottom_THGEM1 + gPars::det_dims->THGEM1_width_total / 2.0);
  gPars::det_dims->EL_gap_center =
      G4ThreeVector(gPars::det_dims->THGEM1_hole_pitch, 0, (z_top_interface_grid + z_bottom_THGEM1) / 2.0);
  gPars::det_dims->Cathode_top_center =
      G4ThreeVector(0, 0, 0);
  gPars::det_dims->THGEM1_center = G4ThreeVector(0, 0, z_bottom_THGEM1 + gPars::det_dims->THGEM1_width_total / 2.0);
  gPars::det_dims->THGEM1_single_cell_position = position_SingleTHGEM1Cell;
  gPars::det_dims->n_PMTs = 4;
  gPars::det_dims->n_SiPMs = Nx_SiPMs * Ny_SiPMs;
}

