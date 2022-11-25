#include <geant4/detector/Detector_full_y2022.hh>

Detector_full_y2022::Detector_full_y2022() :
  VDetectorConstruction()
{}

Detector_full_y2022::~Detector_full_y2022()
{}

void Detector_full_y2022::CreateTHGEM0Cell()
{
	DetectorDimsFullY2022 *dims = (DetectorDimsFullY2022*) gPars::det_dims;
	double cell_size_x = dims->THGEM0_hole_pitch / 2.0;
	double cell_size_y = cell_size_x * std::sqrt(3.0);
	double cell_size_z = dims->THGEM0_container_width;
	double diel_size_z = dims->THGEM0_dielectric_thickness;
	double radius = dims->THGEM0_hole_radius;
	double radius_center = gPars::det_dims->THGEM1_dielectric_radius;
	double radius_cu = radius + dims->THGEM0_hole_rim;
	double cu_size_z = dims->THGEM0_copper_thickness;

	G4ThreeVector zero(0.0, 0.0, 0.0);
	G4ThreeVector hole_1_pos(cell_size_x / 2.0, cell_size_y / 2.0, 0.0);
	G4ThreeVector hole_2_pos(-cell_size_x / 2.0, -cell_size_y / 2.0, 0.0);
	G4ThreeVector cu_top_pos(0.0, 0.0, diel_size_z / 2.0 + cu_size_z / 2.0);
	G4ThreeVector cu_bot_pos(0.0, 0.0, -diel_size_z / 2.0 - cu_size_z / 2.0);

	G4Box* solid_THGEM0_cell_isolation = new G4Box("solid_THGEM0_cell_isolation", cell_size_x, cell_size_y, cell_size_z); //so that cell is in the same material
	logic_THGEM0_cell = new G4LogicalVolume(solid_THGEM0_cell_isolation, matLAr, "logic_THGEM0_cell_isolation", 0, 0, 0);

	G4Box* solid_THGEM0_cell_LAr = new G4Box("solid_THGEM0_cell_LAr", cell_size_x / 2.0, cell_size_y / 2.0, cell_size_z / 2.0);
	logic_THGEM0_cell_LAr = new G4LogicalVolume(solid_THGEM0_cell_LAr, matLAr, "logic_THGEM0_cell_LAr", 0, 0, 0);
	phys_THGEM0_cell_LAr = new G4PVPlacement(0, zero, logic_THGEM0_cell_LAr, "phys_THGEM0_cell",
			logic_THGEM0_cell, false, 0, fCheckOverlaps);

	G4Box* solid_THGEM0_diel_box = new G4Box("solid_THGEM0_diel_box", cell_size_x / 2.0, cell_size_y / 2.0, diel_size_z / 2.0);
	double z_epsilon = 1.01; //solid which is subtracted must be larger than (not have coincidental faces with) parent solid.
	// Otherwise G4SubtractionSolid's visualization breaks. It may break with this approach anyway. It also depends on the absolute scale.
	// This is ONLY visualization issue. Error looks like:
	// ERROR: G4VSceneHandler::RequestPrimitives
	// Polyhedron not available for solid_THGEM1_diel
	//
	double r_epsilon = radius_center + (radius - radius_center) * z_epsilon;
	if (r_epsilon < 0) { // Although radius must always be > radius_center in GEMs simply due to their manufacture process.
		r_epsilon = 0;
		z_epsilon = radius_center / (radius_center - radius);
	}
	double Zs[] = {-diel_size_z / 2.0 * z_epsilon, 0, diel_size_z / 2.0 * z_epsilon};
	double Rs[] = {r_epsilon, radius_center, r_epsilon};
	G4GenericPolycone* solid_diel_hole1 = new G4GenericPolycone("solid_diel_hole1", 177.*deg, 273.*deg, 3, Rs, Zs);
	G4GenericPolycone* solid_diel_hole2 = new G4GenericPolycone("solid_diel_hole2", -3.*deg, 93.*deg, 3, Rs, Zs);
	G4SubtractionSolid* solid_THGEM0_diel_tmp = new G4SubtractionSolid("solid_THGEM0_diel_tmp", solid_THGEM0_diel_box, solid_diel_hole1, 0, hole_1_pos);
	G4SubtractionSolid* solid_THGEM0_diel = new G4SubtractionSolid("solid_THGEM0_diel", solid_THGEM0_diel_tmp, solid_diel_hole2, 0, hole_2_pos);
	logic_THGEM0_cell_FR4 = new G4LogicalVolume(solid_THGEM0_diel, matFR4, "logic_THGEM0_cell_FR4", 0, 0, 0);
	G4VPhysicalVolume* phys_THGEM0_cell_FR4 = new G4PVPlacement(0, zero, logic_THGEM0_cell_FR4, "phys_THGEM0_cell_FR4",
			logic_THGEM0_cell_LAr, false, 0, fCheckOverlaps);

	G4Box* solid_THGEM0_cu_box = new G4Box("solid_THGEM0_cu_box", cell_size_x / 2.0, cell_size_y / 2.0, cu_size_z / 2.0);
	G4Tubs* solid_THGEM0_cu_hole1 = new G4Tubs("solid_THGEM0_cu_hole", 0, radius_cu, cu_size_z / 1.9, 177.*deg, 273.*deg);
	G4Tubs* solid_THGEM0_cu_hole2 = new G4Tubs("solid_THGEM0_cu_hole", 0, radius_cu, cu_size_z / 1.9, -3.*deg, 93.*deg);
	G4SubtractionSolid* solid_THGEM0_cu_tmp = new G4SubtractionSolid("solid_THGEM0_cu_tmp", solid_THGEM0_cu_box, solid_THGEM0_cu_hole1, 0, hole_1_pos);
	G4SubtractionSolid* solid_THGEM0_cu = new G4SubtractionSolid("solid_THGEM0_cu", solid_THGEM0_cu_tmp, solid_THGEM0_cu_hole2, 0, hole_2_pos);
	logic_THGEM0_cell_copper = new G4LogicalVolume(solid_THGEM0_cu, matFR4, "logic_THGEM0_cell_copper", 0, 0, 0);
	G4VPhysicalVolume* phys_THGEM0_cell_copper_top = new G4PVPlacement(0, cu_top_pos, logic_THGEM0_cell_copper, "phys_THGEM0_cell_copper_top",
			logic_THGEM0_cell_LAr, false, 0, fCheckOverlaps);
	G4VPhysicalVolume* phys_THGEM0_cell_copper_bot = new G4PVPlacement(0, cu_bot_pos, logic_THGEM0_cell_copper, "phys_THGEM0_cell_copper_bot",
			logic_THGEM0_cell_LAr, false, 0, fCheckOverlaps);

	//--------------------------------------------------------------------------------
	// Setting visualization
	G4VisAttributes Invisible(G4Colour(1, 1, 1, 0.0));
	Invisible.SetVisibility(false);
	G4VisAttributes AlmostInvisible(G4Colour(0.6, 0.6, 1.0, 0.05));
	G4VisAttributes LAr_VisAtt(G4Colour(0.6, 0.6, 1.0, 0.0));
	G4VisAttributes FR4_VisAtt(G4Colour(0.8, 0.85, 0.11, 0.8));
	G4VisAttributes Cu_VisAtt(G4Colour(0.8, 0.45, 0.2, 0.9));

	Invisible.SetForceWireframe(true);
	Cu_VisAtt.SetForceWireframe(true);
	FR4_VisAtt.SetForceWireframe(true);
	//FR4_VisAtt.SetForceLineSegmentsPerCircle(200);
	AlmostInvisible.SetForceWireframe(true);
	// Separate THGEM hole
	logic_THGEM0_cell_copper->SetVisAttributes(Cu_VisAtt);
	logic_THGEM0_cell_LAr->SetVisAttributes(Invisible);
	logic_THGEM0_cell_FR4->SetVisAttributes(FR4_VisAtt);
	logic_THGEM0_cell->SetVisAttributes(AlmostInvisible);
}

G4VPhysicalVolume * Detector_full_y2022::Construct()
{
  SetSizeAndPosition();
  defineMaterials();
  defineSurfaces();
	//-------------------------------------------------------------------------------
  const DetectorDimsFullY2022 *dims = (DetectorDimsFullY2022*) gPars::det_dims;

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
	G4Box* solid_LAr_inner = new G4Box("solid_LAr_inner", LAr_inner_size_xy / 2.0, LAr_inner_size_xy / 2.0, LAr_inner_size_z / 2.0);
	G4LogicalVolume* logic_LAr_inner = new G4LogicalVolume(solid_LAr_inner, matLAr, "logic_LAr_inner", 0, 0, 0);
	G4VPhysicalVolume* phys_LAr_inner = new G4PVPlacement(0, position_LAr_inner, logic_LAr_inner,  "phys_LAr_inner",
		logicWorld, false, 0, fCheckOverlaps);

	//--------------------------------------------------------------------------------
	//create LAr box contained outside the PMMA insulator
	G4Box* solid_LAr_outer_out = new G4Box("solid_LAr_outer_out", LAr_outer_out_size_xy / 2.0, LAr_outer_out_size_xy / 2.0, LAr_outer_size_z / 2.0);
	G4Box* solid_LAr_outer_in = new G4Box("solid_LAr_outer_in", LAr_outer_in_size_xy / 2.0, LAr_outer_in_size_xy / 2.0, LAr_outer_size_z / 2.0);
	G4SubtractionSolid* solid_LAr_outer = new G4SubtractionSolid("solid_LAr_outer", solid_LAr_outer_out, solid_LAr_outer_in);
	G4LogicalVolume* logic_LAr_outer = new G4LogicalVolume(solid_LAr_outer, matLAr, "logic_LAr_outer", 0, 0, 0);
	G4VPhysicalVolume* phys_LAr_outer = new G4PVPlacement(0, position_LAr_outer, logic_LAr_outer, "phys_LAr_outer",
		logicWorld, false, 0, fCheckOverlaps);

	//--------------------------------------------------------------------------------
	//create FieldWires
	G4Tubs* solid_FieldWire = new G4Tubs("solid_FieldWire", 0, FieldWire_radius, FieldWire_length / 2.0, 0.*deg, 360.*deg);
	G4LogicalVolume* logic_FieldWire = new G4LogicalVolume(solid_FieldWire, matAl, "logic_FieldWire", 0, 0, 0);

	G4VPhysicalVolume* phys_FieldWire_bottom1 = new G4PVPlacement(rotX_90, position_FieldWire_bottom1, logic_FieldWire,
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
	// Create Cathode
	G4Box* solid_bCathode = new G4Box("solid_bCathode", Cathode_size_xy / 2.0, Cathode_size_xy / 2.0, Cathode_size_z / 2.0);
	G4LogicalVolume* logic_bCathode = new G4LogicalVolume(solid_bCathode, matFR4, "logic_bCathode", 0, 0, 0);
	G4VPhysicalVolume* phys_bCathode = new G4PVPlacement(0, position_Cathode, logic_bCathode, "phys_bCathode",
		logicWorld, false, 0, fCheckOverlaps);

	//--------------------------------------------------------------------------------
	// Create THGEM0
	CreateTHGEM0Cell();
	G4VPhysicalVolume* phys_THGEM0_cell = new G4PVPlacement(0, position_SingleTHGEM0Cell, logic_THGEM0_cell,
		"phys_THGEM0_cell_isolated", logicWorld, false, 0, fCheckOverlaps);
	G4Box* solid_whole_THGEM0 = new G4Box("solid_whole_THGEM0", THGEM0_size_xy / 2.0, THGEM0_size_xy / 2.0, dims->THGEM0_width_total / 2.0);
	G4Box* solid_active_THGEM0 = new G4Box("solid_whole_THGEM0", THGEM0_active_size_xy / 2.0, THGEM0_active_size_xy / 2.0, dims->THGEM0_width_total / 2.0);
	G4SubtractionSolid* solid_THGEM0_frame = new G4SubtractionSolid("solid_THGEM0_frame", solid_whole_THGEM0, solid_active_THGEM0);
	G4LogicalVolume* logic_THGEM0_frame = new G4LogicalVolume(solid_THGEM0_frame, matFR4, "logic_THGEM0_frame", 0, 0, 0);
	G4VPhysicalVolume* phys_THGEM0_frame = new G4PVPlacement(0, position_THGEM0_frame, logic_THGEM0_frame, "phys_THGEM0_frame",
			logic_LAr_inner, false, 0, fCheckOverlaps);

	G4Box* solid_THGEM0_container = new G4Box("solid_THGEM0_container", THGEM0_active_size_xy / 2.0, THGEM0_active_size_xy / 2.0, dims->THGEM0_container_width / 2.0);
	G4LogicalVolume* logic_THGEM0_container = new G4LogicalVolume(solid_THGEM0_container, matLAr, "logic_THGEM0_container", 0, 0, 0);
	phys_THGEM0_container = new G4PVPlacement(0, position_THGEM0_container, logic_THGEM0_container, "phys_THGEM0_cell_container",
				logic_LAr_inner, false, 0, fCheckOverlaps);

	G4LogicalVolume* logic_THGEM0_copper = new G4LogicalVolume(solid_active_THGEM0, matFR4, "logic_THGEM0_copper", 0, 0, 0);
	G4VPhysicalVolume* phys_THGEM0_copper = new G4PVPlacement(0, position_THGEM0_copper_plate, logic_THGEM0_copper, "phys_THGEM0_copper",
			logic_THGEM0_container, false, 0, fCheckOverlaps);

	//--------------------------------------------------------------------------------
	// Create THGEM1
	CreateTHGEM1Cell();
	G4VPhysicalVolume* phys_THGEM1_cell = new G4PVPlacement(0, position_SingleTHGEM1Cell, logic_THGEM1_cell,
		"phys_THGEM1_cell_isolated", logicWorld, false, 0, fCheckOverlaps);
	G4Box* solid_whole_THGEM1 = new G4Box("solid_whole_THGEM1", THGEM1_size_xy / 2.0, THGEM1_size_xy / 2.0, dims->THGEM1_width_total / 2.0);
	G4Box* solid_active_THGEM1 = new G4Box("solid_whole_THGEM1", THGEM1_active_size_xy / 2.0, THGEM1_active_size_xy / 2.0, dims->THGEM1_width_total / 2.0);
	G4SubtractionSolid* solid_THGEM1_frame = new G4SubtractionSolid("solid_THGEM1_frame", solid_whole_THGEM1, solid_active_THGEM1);
	G4LogicalVolume* logic_THGEM1_frame = new G4LogicalVolume(solid_THGEM1_frame, matFR4, "logic_THGEM1_frame", 0, 0, 0);
	G4VPhysicalVolume* phys_THGEM1_frame = new G4PVPlacement(0, position_THGEM1_frame, logic_THGEM0_frame, "phys_THGEM1_frame",
			logicWorld, false, 0, fCheckOverlaps);

	G4Box* solid_THGEM1_container = new G4Box("solid_THGEM1_container", THGEM1_active_size_xy / 2.0, THGEM1_active_size_xy / 2.0, dims->THGEM1_container_width / 2.0);
	G4LogicalVolume* logic_THGEM1_container = new G4LogicalVolume(solid_THGEM1_container, matLAr, "logic_THGEM1_container", 0, 0, 0);
	phys_THGEM1_container = new G4PVPlacement(0, position_THGEM1_container, logic_THGEM1_container, "phys_THGEM1_cell_container",
			logicWorld, false, 0, fCheckOverlaps);

	G4LogicalVolume* logic_THGEM1_copper = new G4LogicalVolume(solid_active_THGEM1, matFR4, "logic_THGEM1_copper", 0, 0, 0);
	G4VPhysicalVolume* phys_THGEM1_copper = new G4PVPlacement(0, position_THGEM1_copper_plate, logic_THGEM1_copper, "phys_THGEM1_copper",
			logic_THGEM1_container, false, 0, fCheckOverlaps);

	// Create anode grid. Parameterisation is required.
	// Create tracker (container for parameterised volumes)
	G4Box* solid_tracker_anode_grid = new G4Box("solid_tracker_anode_grid", anode_grid_cont_size_xy / 2.0, anode_grid_cont_size_xy / 2.0, anode_grid_cont_size_z / 2.0);
	G4LogicalVolume* logic_tracker_anode_grid = new G4LogicalVolume(solid_tracker_anode_grid, matGas, "logic_tracker_anode_grid", 0, 0, 0);
	G4VPhysicalVolume* phys_tracker_anode_grid = new G4PVPlacement(0, position_anode_grid_container, logic_tracker_anode_grid, "phys_tracker_anode_grid",
		logicWorld, false, 0, fCheckOverlaps);
	// Create anode grid frame
	G4Box* solid_anode_grid_substrate = new G4Box("anode_grid_substrate", anode_grid_size_xy / 2.0, anode_grid_size_xy / 2.0, anode_grid_thickness / 2.0);
	G4Box* solid_anode_grid_hole = new G4Box("anode_grid_hole", anode_grid_cont_size_xy / 2.0, anode_grid_cont_size_xy / 2.0, anode_grid_thickness / 2.0);
	G4SubtractionSolid* solid_anode_grid_frame = new G4SubtractionSolid("anode_grid_frame", solid_anode_grid_substrate, solid_anode_grid_hole);
	G4LogicalVolume* logic_anode_grid_frame = new G4LogicalVolume(solid_anode_grid_frame, matFR4, "l_anode_grid", 0, 0, 0);
	G4VPhysicalVolume* phys_anode_grid_frame = new G4PVPlacement(0, position_anode_grid_frame, logic_anode_grid_frame, "p_anode_grid",
		logicWorld, false, 0);
	// Create anode wire grid
	G4Tubs* solid_wire = new G4Tubs("solid_wire", 0, anode_wire_radius, anode_wire_length / 2.0, 0.*deg, 360.*deg);
	G4LogicalVolume* logic_anode_wire = new G4LogicalVolume(solid_wire, matAl, "lwire", 0, 0, 0);
	G4VPVParameterisation* param_wire = new DetectorParameterisation(anode_wire_N, 1, 1, rotX_90, G4ThreeVector(0,0,0), G4ThreeVector(anode_wire_step, 0, 0));
	G4VPhysicalVolume* phys_wire = new G4PVParameterised("phys_wire", logic_anode_wire,	logic_tracker_anode_grid,
		kXAxis, anode_wire_N, param_wire, fCheckOverlaps);

	//--------------------------------------------------------------------------------
	// Create PMMA plate
	G4Box* solid_PMMA_plate = new G4Box("solid_tracker_anode_grid", PMMA_plate_size_xy / 2.0, PMMA_plate_size_xy / 2.0, PMMA_plate_size_z / 2.0);
	G4LogicalVolume* logic_PMMA_plate = new G4LogicalVolume(solid_PMMA_plate, matPMMA, "logic_PMMA_plate", 0, 0, 0);
	G4VPhysicalVolume* phys_PMMA_plate = new G4PVPlacement(0, position_PMMA_plate, logic_PMMA_plate, "phys_PMMA_plate",
		logicWorld, false, 0, fCheckOverlaps);

	//--------------------------------------------------------------------------------
	// Create SiPM matrix
	// create container (this is required for SiPM_matrix parameterisation)
	G4Box* solid_SiPM_conntainer = new G4Box("solid_SiPM_conntainer", SiPM_cont_size_xy / 2.0, SiPM_cont_size_xy / 2.0, SiPM_cont_size_z / 2.0);
	G4LogicalVolume* logic_SiPM_container = new G4LogicalVolume(solid_SiPM_conntainer, matGas, "logic_SiPM_conntainer", 0, 0, 0);
	G4VPhysicalVolume* phys_SiPM_container = new G4PVPlacement(0, position_SiPM_container, logic_SiPM_container, "phys_SiPM_conntainer", logicWorld,
		false, 0, fCheckOverlaps);

	G4Box* solid_SiPM = new G4Box("solid_SiPM", SiPM_size / 2.0, SiPM_size / 2.0, SiPM_thickness / 2.0);
	logic_SiPM = new G4LogicalVolume(solid_SiPM, matAl, "logic_SiPM", 0, 0, 0);
	G4VPVParameterisation* chamberParam = new DetectorParameterisation(Nx_SiPMs, Ny_SiPMs, 1, NULL, offset_SiPM_in_container, G4ThreeVector(SiPM_pitch, SiPM_pitch, 0));
	G4VPhysicalVolume* phys_SiPM = new G4PVParameterised(gPars::det_dims->SiPM_device_name, logic_SiPM, logic_SiPM_container,
		kXAxis, Nx_SiPMs * Ny_SiPMs, chamberParam, fCheckOverlaps);

	G4Box* solid_SiPMFR4 = new G4Box("solid_SiPMFR4", SiPM_cont_size_xy / 2.0, SiPM_cont_size_xy / 2.0, SiPM_FR4_size_z / 2.0);
	G4LogicalVolume* logic_SiPMFR4 = new G4LogicalVolume(solid_SiPMFR4, matGas, "logic_SiPMFR4", 0, 0, 0);
	G4VPhysicalVolume* phys_SiPMFR4 = new G4PVPlacement(0, position_SiPMFR4, logic_SiPMFR4, "phys_SiPMFR4",
		logicWorld, false, 0, fCheckOverlaps);

	//--------------------------------------------------------------------------------
	// Create PMT's wire grid
	// Gas part
	G4Box* solid_PMTAnodeGridTrackerGas = new G4Box("solid_PMTAnodeGridTrackerGas", PMT_grid_cont_gas_size_x / 2.0, PMT_grid_cont_gas_size_y / 2.0, PMT_grid_cont_thickness / 2.0);
	G4LogicalVolume* logic_PMTAnodeGridTrackerGas = new G4LogicalVolume(solid_PMTAnodeGridTrackerGas, matGas, "logic_PMTAnodeGridTrackerGas", 0, 0, 0);
	G4LogicalVolume* logic_PMTAnodeGridTrackerGasInner = new G4LogicalVolume(solid_PMTAnodeGridTrackerGas, matGas, "logic_PMTAnodeGridTrackerGasInner", 0, 0, 0);

	G4VPhysicalVolume* phys_PMTAnodeGridTrackerGas_1 = new G4PVPlacement(rotY_90, position_PMT_grid_cont_gas_1,
		logic_PMTAnodeGridTrackerGas,	"phys_PMTAnodeGridTrackerGas_1", logicWorld, false, 0, fCheckOverlaps);

	G4VPhysicalVolume* phys_PMTAnodeGridTrackerGasInner_1 = new G4PVPlacement(rotY_90, position_PMT_grid_cont_gas_inner_1,
		logic_PMTAnodeGridTrackerGasInner, "phys_PMTAnodeGridTrackerGasInner_1", logicWorld, false, 0, fCheckOverlaps);

	// Liquid part
	G4Box* solid_PMTAnodeGridTrackerLiquid = new G4Box("solid_PMTAnodeGridTrackerLiquid", PMT_grid_cont_liquid_size_x / 2.0, PMT_grid_cont_liquid_size_y / 2.0, PMT_grid_cont_thickness / 2.0);
	G4LogicalVolume* logic_PMTAnodeGridTrackerLiquid = new G4LogicalVolume(solid_PMTAnodeGridTrackerLiquid, matLAr, "logic_PMTAnodeGridTrackerLiquid", 0, 0, 0);
	G4LogicalVolume* logic_PMTAnodeGridTrackerLiquidInner = new G4LogicalVolume(solid_PMTAnodeGridTrackerLiquid, matLAr, "logic_PMTAnodeGridTrackerLiquidInner", 0, 0, 0);
	G4VPhysicalVolume* phys_PMTAnodeGridTrackerLiquid_1 = new G4PVPlacement(rotY_90, position_PMT_grid_cont_liquid_1,
		logic_PMTAnodeGridTrackerLiquid, "phys_PMTAnodeGridTrackerLiquid_1", logic_LAr_outer, false, 0, fCheckOverlaps);
	G4VPhysicalVolume* phys_PMTAnodeGridTrackerLiquidInner_1 = new G4PVPlacement(rotY_90, position_PMT_grid_cont_liquid_inner_1,
		logic_PMTAnodeGridTrackerLiquidInner, "phys_PMTAnodeGridTrackerLiquidInner_1", logic_LAr_outer, false, 0, fCheckOverlaps);

	//--------------------------------------------------------------------------------
	//create PMTGridWire
	G4Tubs* solid_PMTGridWire = new G4Tubs("solid_PMTGridWire", 0, PMT_grid_wire_radius, PMT_grid_cont_gas_size_y / 2.0, 0.*deg, 360.*deg);
	G4LogicalVolume* logic_PMTGridWire = new G4LogicalVolume(solid_PMTGridWire, matAl, "lwire", 0, 0, 0);
	G4VPVParameterisation* param_PMTGridWireGas = new DetectorParameterisation(PMT_grid_gas_Ncells, 1, 1, rotX_90, G4ThreeVector(0,0,0), G4ThreeVector(PMT_grid_wire_pitch, 0, 0));
	G4VPhysicalVolume* phys_PMTGridWireGas_0 = new G4PVParameterised("phys_PMTGridWireGas_0", logic_PMTGridWire, logic_PMTAnodeGridTrackerGas, kXAxis,
			PMT_grid_gas_Ncells, param_PMTGridWireGas, fCheckOverlaps);

	G4Tubs* solid_PMTGridWireGasInner = new G4Tubs("solid_PMTGridWireGasInner", 0, PMT_grid_wire_radius, PMT_grid_cont_gas_size_x / 2.0, 0.*deg, 360.*deg);
	G4LogicalVolume* logic_PMTGridWireGasInner = new G4LogicalVolume(solid_PMTGridWireGasInner, matAl, "lwire", 0, 0, 0);
	G4VPVParameterisation* param_PMTGridWireGasInner = new DetectorParameterisation(1, PMT_grid_gas_inner_Ncells, 1, rotY_90, G4ThreeVector(0,0,0), G4ThreeVector(0, PMT_grid_wire_pitch, 0));
	G4VPhysicalVolume* phys_PMTGridWireGasInner_0 = new G4PVParameterised("phys_PMTGridWireGasInner_0", logic_PMTGridWireGasInner,
		logic_PMTAnodeGridTrackerGasInner, kXAxis, PMT_grid_gas_inner_Ncells, param_PMTGridWireGasInner, fCheckOverlaps);

	G4VPVParameterisation* param_PMTGridWireLiquid = new DetectorParameterisation(PMT_grid_liquid_Ncells, 1, 1, rotX_90, G4ThreeVector(0,0,0), G4ThreeVector(PMT_grid_wire_pitch, 0, 0));
	G4VPhysicalVolume* phys_PMTGridWireLiquid_0 = new G4PVParameterised("phys_PMTGridWireLiquid_0", logic_PMTGridWire,
		logic_PMTAnodeGridTrackerLiquid, kXAxis, PMT_grid_liquid_Ncells, param_PMTGridWireLiquid, fCheckOverlaps);

	G4Tubs* solid_PMTGridWireLiquidInner = new G4Tubs("solid_PMTGridWireLiquidInner", 0, PMT_grid_wire_radius, PMT_grid_cont_liquid_size_x / 2.0, 0.*deg, 360.*deg);
	G4LogicalVolume* logic_PMTGridWireLiquidInner = new G4LogicalVolume(solid_PMTGridWireLiquidInner, matAl, "lwire", 0, 0, 0);
	G4VPVParameterisation* param_PMTGridWireLiquidInner = new DetectorParameterisation(1, PMT_grid_liquid_inner_Ncells, 1, rotY_90, G4ThreeVector(0,0,0), G4ThreeVector(0, PMT_grid_wire_pitch, 0));
	G4VPhysicalVolume* phys_PMTGridWireLiquidInner_0 = new G4PVParameterised("phys_PMTGridWireLiquidInner_0", logic_PMTGridWireLiquidInner,
		logic_PMTAnodeGridTrackerLiquidInner, kXAxis, PMT_grid_liquid_inner_Ncells, param_PMTGridWireLiquidInner, fCheckOverlaps);


	//--------------------------------------------------------------------------------
	// Create insulator box
	G4Box* solid_Insulator_box_inner = new G4Box("solid_Insulator_box_inner", insulator_box_inner_size_xy / 2.0, insulator_box_inner_size_xy / 2.0, insulator_box_size_z / 2.0);
	G4Box* solid_Insulator_box_outer = new G4Box("solid_Insulator_box_outer", insulator_box_outer_size_xy / 2.0, insulator_box_outer_size_xy / 2.0, insulator_box_size_z / 2.0);
	G4SubtractionSolid* solid_Insulator_box_subtraction = new G4SubtractionSolid("solid_Insulator_box_subtraction", solid_Insulator_box_outer, solid_Insulator_box_inner);
	G4LogicalVolume* logic_Insulator_box = new G4LogicalVolume(solid_Insulator_box_subtraction, matPMMA_UV, "logic_Insulator_box", 0, 0, 0);
	G4VPhysicalVolume* phys_Insulator_box = new G4PVPlacement(0, position_insulator_box, logic_Insulator_box,
		"phys_Insulator_box", logicWorld, false, 0, fCheckOverlaps);

	// Create TPB inside the insulator box
	G4Tubs* solid_TPB = new G4Tubs("solid_TPB", 0, TPB_radius, TPB_thickness / 2.0, 0.*deg, 360.*deg);
	G4LogicalVolume* logic_TPB = new G4LogicalVolume(solid_TPB, matTPB, "logic_TPB", 0, 0, 0);
	if (dims->has_WLS) {
		G4VPhysicalVolume* phys_TPB0 = new G4PVPlacement(rotY_90, position_TPB_0, logic_TPB, "phys_TPB0",
			logic_Insulator_box, false, 0, fCheckOverlaps);
		G4VPhysicalVolume* phys_TPB1 = new G4PVPlacement(rotY_90, position_TPB_1, logic_TPB, "phys_TPB1",
			logic_Insulator_box, false, 1, fCheckOverlaps);
		G4VPhysicalVolume* phys_TPB2 = new G4PVPlacement(rotX_90, position_TPB_2, logic_TPB, "phys_TPB2",
			logic_Insulator_box, false, 2, fCheckOverlaps);
		G4VPhysicalVolume* phys_TPB3 = new G4PVPlacement(rotX_90, position_TPB_3, logic_TPB, "phys_TPB3",
			logic_Insulator_box, false, 3, fCheckOverlaps);
	}

	//--------------------------------------------------------------------------------
	// Create PMTs
	G4Tubs* solid_PMT = new G4Tubs("solid_PMT", 0, PMT_radius, PMT_size_z / 2.0, 0.*deg, 360.*deg);
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
	G4Box* solidSteelBox = new G4Box("solidSteelBox", steel_box_size_x / 2.0, steel_box_size_y / 2.0, steel_box_size_z / 2.0);
	G4VSolid* solidBoxSubtractPMT = new G4SubtractionSolid("solidBoxSubtractPMT", solidSteelBox, solid_PMT, rotY_90, G4ThreeVector((steel_box_size_x - PMT_size_z) / 2.0, 0., 0.));
	G4LogicalVolume* logicSteelBox = new G4LogicalVolume(solidBoxSubtractPMT, matGas, "logicSteelBox", 0, 0, 0);

	G4VPhysicalVolume* physSteelBox0 = new G4PVPlacement(0, position_steel_box_0, logicSteelBox, "physSteelBox0",
		logicWorld, false, 0, fCheckOverlaps);
	G4VPhysicalVolume* physSteelBox1 = new G4PVPlacement(rotZ_180, position_steel_box_1, logicSteelBox, "physSteelBox1",
		logicWorld, false, 0, fCheckOverlaps);
	G4VPhysicalVolume* physSteelBox2 = new G4PVPlacement(rotZ_270, position_steel_box_2, logicSteelBox, "physSteelBox2",
		logicWorld, false, 0, fCheckOverlaps);
	G4VPhysicalVolume* physSteelBox3 = new G4PVPlacement(rotZ_90, position_steel_box_3, logicSteelBox, "physSteelBox3",
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

	// THGEM1
	G4LogicalSkinSurface* THGEM1_cell_cu = new G4LogicalSkinSurface("THGEM1_cell_cu_surface", logic_THGEM1_cell_copper, Cu_THGEM);
  G4LogicalSkinSurface* THGEM1_cell_FR4 = new G4LogicalSkinSurface("THGEM1_cell_FR4_surface", logic_THGEM1_cell_FR4, FR4_unified);
  G4LogicalBorderSurface* THGEM1_cell_isolation = new G4LogicalBorderSurface("THGEM1_cell_isolation_surface", physiWorld, phys_THGEM1_cell, AbsorberMaterial);
  G4LogicalBorderSurface* THGEM1_cell_isolation1 = new G4LogicalBorderSurface("THGEM1_cell_isolation_surface1", phys_THGEM1_cell, physiWorld, AbsorberMaterial);
  G4LogicalSkinSurface* THGEM1_frame = new G4LogicalSkinSurface("THGEM1_frame", logic_THGEM1_frame, FR4_unified);
  G4LogicalSkinSurface* THGEM1_cu = new G4LogicalSkinSurface("THGEM1_frame", logic_THGEM1_copper, Cu_THGEM);

  // THGEM0
	G4LogicalSkinSurface* THGEM0_cell_cu = new G4LogicalSkinSurface("THGEM0_cell_cu_surface", logic_THGEM0_cell_copper, Cu_THGEM);
	G4LogicalSkinSurface* THGEM0_cell_FR4 = new G4LogicalSkinSurface("THGEM0_cell_FR4_surface", logic_THGEM0_cell_FR4, FR4_unified);
	G4LogicalBorderSurface* THGEM0_cell_isolation = new G4LogicalBorderSurface("THGEM0_cell_isolation_surface", physiWorld, phys_THGEM0_cell, AbsorberMaterial);
	G4LogicalBorderSurface* THGEM0_cell_isolation1 = new G4LogicalBorderSurface("THGEM0_cell_isolation_surface1", phys_THGEM0_cell, physiWorld, AbsorberMaterial);
	G4LogicalSkinSurface* THGEM0_frame = new G4LogicalSkinSurface("THGEM0_frame", logic_THGEM0_frame, FR4_unified);
	G4LogicalSkinSurface* THGEM0_cu = new G4LogicalSkinSurface("THGEM0_frame", logic_THGEM0_copper, Cu_THGEM);

	//FieldWires
	G4LogicalSkinSurface* sur_FieldWire = new G4LogicalSkinSurface("sur_FieldWire", logic_FieldWire, stainlessSteel);

	//--------------------------------------------------------------------------------
	// Setting visualization
	G4VisAttributes LAr_VisAtt(G4Colour(0.6, 0.6, 1.0, 0.0));
	G4VisAttributes gas_VisAtt(G4Colour(1.0, 1.0, 1.0, 0.0));
	gas_VisAtt.SetVisibility(false);
	G4VisAttributes PMMA_VisAtt(G4Colour(0.2, 0.2, 0.8, 0.1));
	G4VisAttributes TPB_VisAtt(G4Colour(0.4, 0.8, 0.4, 0.2));
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
  logic_SiPM_container->SetVisAttributes(gas_VisAtt);
  FR4_VisAtt.SetVisibility(false);
  logic_SiPMFR4->SetVisAttributes(FR4_VisAtt);
  FR4_VisAtt.SetVisibility(true);

	// Interface grid (THGEM0)
  logic_THGEM0_frame->SetVisAttributes(FR4_VisAtt);
  logic_THGEM0_copper->SetVisAttributes(Cu_VisAtt);
  logic_THGEM0_container->SetVisAttributes(LAr_VisAtt);

  // THGEM1
	logic_THGEM1_copper->SetVisAttributes(Cu_VisAtt);
	logic_THGEM1_container->SetVisAttributes(gas_VisAtt);
	logic_THGEM1_frame->SetVisAttributes(FR4_VisAtt);

	// PMTs
	logic_PMT->SetVisAttributes(Sensor_VisAtt);
	logic_PMTGridWire->SetVisAttributes(Wires_VisAtt);
	logic_PMTGridWireGasInner->SetVisAttributes(Wires_VisAtt);
	logic_PMTGridWireLiquidInner->SetVisAttributes(Wires_VisAtt);
	logic_PMTAnodeGridTrackerLiquidInner->SetVisAttributes(LAr_VisAtt);
	logic_PMTAnodeGridTrackerLiquid->SetVisAttributes(LAr_VisAtt);
	logic_PMTAnodeGridTrackerLiquidInner->SetVisAttributes(gas_VisAtt);
	logic_PMTAnodeGridTrackerLiquid->SetVisAttributes(gas_VisAtt);

	Steel_VisAtt.SetVisibility(false);
	logicSteelBox->SetVisAttributes(Steel_VisAtt);
	Steel_VisAtt.SetVisibility(true);

	// LAr
	logic_LAr_inner->SetVisAttributes(LAr_VisAtt);
	logic_PMTAnodeGridTrackerLiquid->SetVisAttributes(LAr_VisAtt);
	logic_PMTAnodeGridTrackerLiquidInner->SetVisAttributes(LAr_VisAtt);
	logic_LAr_outer->SetVisAttributes(LAr_VisAtt);

	logic_bCathode->SetVisAttributes(FR4_VisAtt);
	Steel_VisAtt.SetVisibility(true);
	PMMA_VisAtt.SetVisibility(true);
	logic_PMMA_plate->SetVisAttributes(PMMA_VisAtt);
	logic_Insulator_box->SetVisAttributes(PMMA_VisAtt);
	logic_TPB->SetVisAttributes(TPB_VisAtt);

	logic_FieldWire->SetVisAttributes(Wires_VisAtt);

	logicWorld->SetVisAttributes(gas_VisAtt);

	SetupTHGEMsMapping();
	return physiWorld;
}

void Detector_full_y2022::ConstructSDandField()
{
  DetectorSensor *thePhotoDetector = new DetectorSensor("/detector/sensitiveDetector");
  G4SDManager::GetSDMpointer()->AddNewDetector(thePhotoDetector);
  SetSensitiveDetector(logic_SiPM, thePhotoDetector); // takes ownership of DetectorSensor*
  SetSensitiveDetector(logic_PMT, thePhotoDetector);
}

void Detector_full_y2022::SetSizeAndPosition()
{
	// x = 0, y = 0 is on the detector axis.
	// z = 0 (global) is the cathode's top.

	DetectorDimsFullY2022 *dims = (DetectorDimsFullY2022*) gPars::det_dims;

	HalfWorldLength = 17 * cm;

	Cathode_size_z = 0.5 * mm;

  double LAr_drift_width = 48.0 * mm; // Distance between cathode top and interface grid bottom
  double max_EL_gap_thickness = 22.0 * mm; // Distance between interface grid top and THGEM1 bottom
  double THGEMs_size_xy = 127 * mm; //see Download:\DetectorPhotos\2021\THGEM_Electroconnect
  double THGEMs_active_size_xy = 100 * mm; //see Download:\DetectorPhotos\2021\THGEM_Electroconnect
  double EL_gap_thickness = 6 * mm; // Must be positive (double phase). Counted from THGEM1's real bottom
	double interface_grid_top_z = LAr_drift_width + dims->THGEM0_width_total;
	double THGEM1_bottom_z = interface_grid_top_z + max_EL_gap_thickness; //=71.4
	LAr_inner_size_z = interface_grid_top_z + (max_EL_gap_thickness - EL_gap_thickness) + Cathode_size_z;
	double LAr_inner_ref_z = -LAr_inner_size_z / 2.0 + Cathode_size_z; // global 0 in LAr_innner reference frame.
	double LAr_level = interface_grid_top_z + (max_EL_gap_thickness - EL_gap_thickness);

	// Cathode
	Cathode_size_xy = THGEMs_size_xy;
	position_Cathode = G4ThreeVector(0, 0, -Cathode_size_z / 2.0 + LAr_inner_ref_z); // LAr is parent

	// Insulator box (acrylic)
	insulator_box_inner_size_xy = 143 * mm;
	insulator_box_thickness = 2 * mm;
	insulator_box_outer_size_xy = insulator_box_inner_size_xy + 2 * insulator_box_thickness;
	insulator_box_size_z = 160 * mm;
	position_insulator_box = G4ThreeVector(0, 0, 0.5 * insulator_box_size_z - Cathode_size_z); // Insulator box bottom is at the same level as cathode's bottom

	// LAr inner (inside the insulator box (acrylic), stretches from EL gap to the cathode)
	LAr_inner_size_xy = insulator_box_inner_size_xy;
	position_LAr_inner = G4ThreeVector(0, 0, 0.5 * LAr_inner_size_z - Cathode_size_z);

	// LAr outer (outside the insulator box (acrylic), stretches from EL gap to the cathode)
	LAr_outer_in_size_xy = insulator_box_outer_size_xy; // Limited by the insulator box
	LAr_outer_out_size_xy = 152 * mm; // Limited by the steel box
	LAr_outer_size_z = LAr_inner_size_z;
	position_LAr_outer = G4ThreeVector(0, 0, 0.5 * LAr_outer_size_z - Cathode_size_z);

	// Field-forming wires
	FieldWire_radius = 1.5 * mm / 2.0;
	FieldWire_length = 98 * mm;
	FieldWire_pos_x = THGEMs_active_size_xy / 2.0;
	FieldWire_bot_pos_z = 20 * mm; //18.2*mm - radius_FieldWire, see 210415 1618834738620-1618834738657 photos
	FieldWire_top_pos_z = 33 * mm; //34.2*mm - radius_FieldWire, see 210415 1618834738620-1618834738657 photos;
	position_FieldWire_bottom1 = G4ThreeVector(FieldWire_pos_x, 0, FieldWire_bot_pos_z + LAr_inner_ref_z);
	position_FieldWire_bottom2 = G4ThreeVector(-FieldWire_pos_x, 0, FieldWire_bot_pos_z + LAr_inner_ref_z);
	position_FieldWire_bottom3 = G4ThreeVector(0, FieldWire_pos_x, FieldWire_bot_pos_z + LAr_inner_ref_z);
	position_FieldWire_bottom4 = G4ThreeVector(0, -FieldWire_pos_x, FieldWire_bot_pos_z + LAr_inner_ref_z);
	position_FieldWire_top1 = G4ThreeVector(FieldWire_pos_x, 0, FieldWire_top_pos_z + LAr_inner_ref_z);
	position_FieldWire_top2 = G4ThreeVector(-FieldWire_pos_x, 0, FieldWire_top_pos_z + LAr_inner_ref_z);
	position_FieldWire_top3 = G4ThreeVector(0, FieldWire_pos_x, FieldWire_top_pos_z + LAr_inner_ref_z);
	position_FieldWire_top4 = G4ThreeVector(0, -FieldWire_pos_x, FieldWire_top_pos_z + LAr_inner_ref_z);

	// Interface THGEM (THGEM0)
	THGEM0_size_xy = THGEMs_size_xy; // Full real size, including dielectric, z size is set in DetectorDimsFullY2022
	THGEM0_active_size_xy = THGEMs_active_size_xy;
	position_THGEM0_container = position_THGEM0_frame = G4ThreeVector(0, 0, interface_grid_top_z - dims->THGEM0_width_total / 2.0 + LAr_inner_ref_z);
	position_THGEM0_copper_plate = G4ThreeVector(0, 0, 0); // dummy inside container in case mapping is avoided

	// THGEM1 (above EL gap and below Anode wires)
	THGEM1_size_xy = THGEMs_size_xy; // Full real size, including dielectric, z size is set in DetectorDimsFullY2022
	THGEM1_active_size_xy = THGEMs_active_size_xy;
	position_THGEM1_container = position_THGEM1_frame = G4ThreeVector(0, 0, THGEM1_bottom_z + dims->THGEM1_width_total / 2.0);
	position_THGEM1_copper_plate = G4ThreeVector(0, 0, 0); // dummy inside container in case mapping is avoided

	// Anode wires (ground before SiPMs)
	anode_wire_radius = 100 / 2.0 * um; //you can understand this from photo;
	anode_wire_length = 60 * mm;
	anode_wire_step = 1 * mm;
	anode_wire_N = anode_wire_length / anode_wire_step - 1;
	// Anode grid
	anode_grid_thickness = 0.5 * mm;
	anode_grid_size_xy = THGEMs_size_xy; // Full size. Has hole of length_wire x length_wire size
	anode_grid_z_bottom = THGEM1_bottom_z + dims->THGEM1_width_total + 5 * mm;
	// Anode grid's container
	anode_grid_cont_size_xy = anode_wire_length;
	anode_grid_cont_size_z = std::max(anode_grid_thickness, anode_wire_radius * 2.1);
	position_anode_grid_container = G4ThreeVector(0, 0, anode_grid_z_bottom + 0.5 * anode_grid_cont_size_z);
	position_anode_grid_frame = G4ThreeVector(0, 0, anode_grid_z_bottom + 0.5 * anode_grid_thickness);

	// PMMA plate, before SiPMs
	PMMA_plate_size_xy = THGEMs_size_xy;
	PMMA_plate_size_z = 1.5 * mm;
	position_PMMA_plate = G4ThreeVector(0, 0, anode_grid_z_bottom + anode_grid_cont_size_z + PMMA_plate_size_z / 2.0);

	// SiPMs
	Nx_SiPMs = 5;
	Ny_SiPMs = 5;
	SiPM_thickness = 1 * nm;
	SiPM_size = 6 * mm;
	SiPM_pitch = 10 * mm;
	// SiPM container
	SiPM_cont_size_xy = std::max(Nx_SiPMs, Ny_SiPMs) * SiPM_pitch + SiPM_size / 2.0;
	SiPM_cont_size_z = 2 * nm;
	SiPM_FR4_size_z = 2 * mm;
	double z_SiPM_top = position_PMMA_plate.z() + PMMA_plate_size_z / 2.0 + SiPM_cont_size_z + (0.2*mm /*small gap between PMMA and SiPM*/);
	position_SiPM_container = G4ThreeVector(0, 0, z_SiPM_top - SiPM_cont_size_z / 2.0);
	position_SiPMFR4 = G4ThreeVector(0, 0, z_SiPM_top + SiPM_FR4_size_z / 2.0);
	offset_SiPM_in_container = G4ThreeVector(0, 0, SiPM_cont_size_z / 2.0 - SiPM_thickness / 2.0); // So that top of SiPMs coincides with FR4 substrate's bottom

	// PMTs
	PMT_radius = 45 * mm / 2.0;
	PMT_size_z = 2 * mm;
	double PMT_x_pos = LAr_outer_out_size_xy / 2.0 + PMT_size_z / 2;
	double PMT_z_pos = 27.2 * mm + 63 * mm / 2.0;
	position_PMT_0 = G4ThreeVector(-PMT_x_pos, 0, PMT_z_pos);
	position_PMT_1 = G4ThreeVector(PMT_x_pos, 0, PMT_z_pos);
	position_PMT_2 = G4ThreeVector(0, -PMT_x_pos, PMT_z_pos);
	position_PMT_3 = G4ThreeVector(0, PMT_x_pos, PMT_z_pos);

	// Steel box (made of 4 plates with hole for PMT)
	steel_box_size_x = 3 * mm;
	steel_box_size_y = LAr_outer_out_size_xy;
	steel_box_size_z = 70 * mm;
	const double steel_box_x_pos = PMT_x_pos - PMT_size_z / 2.0 + steel_box_size_x/2.0;
	position_steel_box_0 = G4ThreeVector(-steel_box_x_pos, 0, PMT_z_pos);
	position_steel_box_1 = G4ThreeVector(steel_box_x_pos, 0, PMT_z_pos);
	position_steel_box_2 = G4ThreeVector(0, -steel_box_x_pos, PMT_z_pos);
	position_steel_box_3 = G4ThreeVector(0, steel_box_x_pos, PMT_z_pos);

	// PMT anode grid wire
	PMT_grid_wire_radius = 150 / 2.0 * um;
	PMT_grid_wire_pitch = 1.2 * mm;
	// PMT anode grid's container
	PMT_grid_cont_thickness = PMT_grid_wire_radius * 2;
	PMT_grid_cont_gas_size_x =  std::max(PMT_radius - (LAr_level - PMT_z_pos), 0.0);
	PMT_grid_cont_gas_size_y = 50 * mm;
	PMT_grid_cont_liquid_size_x = std::max(2*PMT_radius - PMT_grid_cont_gas_size_x, 0.0);
	PMT_grid_cont_liquid_size_y = PMT_grid_cont_gas_size_y;
	PMT_grid_gas_Ncells = PMT_grid_cont_gas_size_x / PMT_grid_wire_pitch;
	PMT_grid_gas_inner_Ncells = PMT_grid_cont_gas_size_y / PMT_grid_wire_pitch;
	PMT_grid_liquid_Ncells = PMT_grid_cont_liquid_size_x / PMT_grid_wire_pitch;
	PMT_grid_liquid_inner_Ncells = PMT_grid_cont_liquid_size_y / PMT_grid_wire_pitch;
	position_PMT_grid_cont_gas_1 = G4ThreeVector(PMT_x_pos - PMT_size_z / 2.0 - PMT_grid_cont_thickness / 2.0,
			0, LAr_level + PMT_grid_cont_gas_size_x / 2.0);
	position_PMT_grid_cont_gas_inner_1 = position_PMT_grid_cont_gas_1 + G4ThreeVector(-PMT_grid_cont_thickness, 0 , 0);
	position_PMT_grid_cont_liquid_1 = G4ThreeVector(PMT_x_pos - PMT_size_z / 2.0 - PMT_grid_cont_thickness / 2.0,
			0, LAr_level - PMT_grid_cont_liquid_size_x / 2.0 + LAr_inner_ref_z); // + LAr_inner_ref_z because LAr_inner is parent volume.
	position_PMT_grid_cont_liquid_inner_1 = position_PMT_grid_cont_liquid_1 + G4ThreeVector(-PMT_grid_cont_thickness, 0, 0);

	// TPB (WLS) on the inside of the insulator box.
	// To avoid splitting TPB cylinder into liquid and gas part, TPB is inserted into the insulator box
	TPB_radius = 70 * mm / 2.0;
	TPB_thickness = 0.4 * mm; // See Borisova's dissertation
	double TPB_x_pos = insulator_box_inner_size_xy / 2.0 + TPB_thickness / 2.0;
	position_TPB_0 = G4ThreeVector(-TPB_x_pos, 0, PMT_z_pos - insulator_box_size_z / 2.0 + Cathode_size_z);
	position_TPB_1 = G4ThreeVector(TPB_x_pos, 0, PMT_z_pos - insulator_box_size_z / 2.0 + Cathode_size_z);
	position_TPB_2 = G4ThreeVector(0, -TPB_x_pos, PMT_z_pos - insulator_box_size_z / 2.0 + Cathode_size_z);
	position_TPB_3 = G4ThreeVector(0, TPB_x_pos, PMT_z_pos - insulator_box_size_z / 2.0 + Cathode_size_z);


	if (anode_grid_thickness < 2 * anode_wire_radius) {
		G4Exception("DetectorConstruction::SetSizeAndPosition: ",
			"InvalidSetup", JustWarning, "anode_grid_thickness < 2*anode_wire_radius");
	}
	if (EL_gap_thickness < (dims->THGEM1_container_width - dims->THGEM1_width_total) / 2.0) {
		G4Exception("DetectorConstruction::SetSizeAndPosition: ",
			"InvalidSetup", FatalException, "EL_gap_thickness is too small. Single phase mode is not supported.");
		return;
	}
	if (EL_gap_thickness > (max_EL_gap_thickness - (dims->THGEM0_container_width - dims->THGEM0_width_total) / 2.0)) {
		G4Exception("DetectorConstruction::SetSizeAndPosition: ",
			"InvalidSetup", FatalException, "EL_gap_thickness is too large.");
		return;
	}

	position_SingleTHGEM1Cell = G4ThreeVector(150 * mm, 150 * mm, 150 * mm);
	position_SingleTHGEM0Cell = G4ThreeVector(-150 * mm, 150 * mm, 150 * mm);


  dims->THGEM1_center = G4ThreeVector(0, 0, THGEM1_bottom_z + gPars::det_dims->THGEM1_width_total / 2.0);

  if (dims->is_NBrS_in_THGEM0) {
		dims->THGEM_hole_center = //x!=0 because x=0 is just across anode wire before SiPM.
				G4ThreeVector(dims->THGEM0_hole_pitch, 0, interface_grid_top_z - dims->THGEM0_width_total / 2.0);
  } else {
		dims->THGEM_hole_center = //x!=0 because x=0 is just across anode wire before SiPM.
				G4ThreeVector(dims->THGEM1_hole_pitch, 0, THGEM1_bottom_z + gPars::det_dims->THGEM1_width_total / 2.0);
  }
  dims->THGEM0_center = G4ThreeVector(0, 0, interface_grid_top_z - dims->THGEM0_width_total / 2.0);

	dims->EL_gap_center = G4ThreeVector(dims->THGEM1_hole_pitch, 0, (interface_grid_top_z + THGEM1_bottom_z) / 2.0);
  dims->Cathode_top_center = G4ThreeVector(0, 0, 0);
  dims->THGEM1_single_cell_position = position_SingleTHGEM1Cell;
  dims->THGEM0_single_cell_position = position_SingleTHGEM0Cell;
  dims->n_PMTs = 4;
  dims->n_SiPMs = Nx_SiPMs * Ny_SiPMs;
}

void Detector_full_y2022::SetupTHGEMsMapping()
{
	DetectorDimsFullY2022 *dims = (DetectorDimsFullY2022*) gPars::det_dims;
	std::vector<G4PhysicalVolumesSearchScene::Findings> cells = LocatePV(phys_THGEM1_cell_LAr);
	std::vector<G4PhysicalVolumesSearchScene::Findings> containers = LocatePV(phys_THGEM1_container);
	if (cells.size() == 0) {
		std::cout<<"********************************"<<std::endl;
		std::cout<<"Detector_full_y2022::SetupTHGEM1Mapping:Warning:\n"
				"\tTHGEM1 cell is not found. Geometry without THGEM1 mapping is used."<<std::endl;
		goto second;
	}
	if (cells.size() > 1) {
		std::cout<<"********************************"<<std::endl;
		std::cout<<"Detector_full_y2022::SetupTHGEM1Mapping:Warning:\n"
				"\tfound several cell physical volumes (name \""<<phys_THGEM1_cell_LAr->GetName()<<")\"!"<<std::endl;
		std::cout<<"\tUsing only the first one."<<std::endl;
	}
	if (containers.size() == 0) {
		std::cout<<"********************************"<<std::endl;
		std::cout<<"Detector_full_y2022::SetupTHGEM1Mapping:Warning:\n"
				"\tTHGEM1 container is not found. Geometry without THGEM1 mapping is used."<<std::endl;
		goto second;
	}
	if (containers.size() > 1) {
		std::cout<<"********************************"<<std::endl;
		std::cout<<"Detector_full_y2022::SetupTHGEM1Mapping:Warning:\n"
				"\tfound several THGEM1 containers' physical volumes (name \""<<phys_THGEM1_container->GetName()<<")\"!"<<std::endl;
		std::cout<<"\tUsing only the first one."<<std::endl;
	}

	{
		G4ThreeVector cell_pos = cells[0].fFoundObjectTransformation.getTranslation();
		G4ThreeVector thgem1_pos = containers[0].fFoundObjectTransformation.getTranslation();
		G4VSolid* cell_box = phys_THGEM1_cell_LAr->GetLogicalVolume()->GetSolid();
		G4VSolid* thgem1_box = phys_THGEM1_container->GetLogicalVolume()->GetSolid();
		G4ThreeVector bMin, bMax;
		cell_box->BoundingLimits(bMin, bMax);
		G4ThreeVector cell_sizes = bMax - bMin;
		thgem1_box->BoundingLimits(bMin, bMax);
		G4ThreeVector thgem1_sizes = bMax - bMin;
		HexagonalMapping thgem1_map("THGEM1", thgem1_pos, cell_pos, thgem1_sizes, cell_sizes, !dims->is_NBrS_in_THGEM0);
		thgem1_map.AddTrigger(MappingTrigger(phys_THGEM1_cell_LAr, false, "THGEM1_leaving_cell"));
		thgem1_map.AddTrigger(MappingTrigger(phys_THGEM1_container, true, "THGEM1_entering_container"));
		gData.mapping_manager.AddMapping(thgem1_map);
	}

second:
	cells = LocatePV(phys_THGEM0_cell_LAr);
	containers = LocatePV(phys_THGEM0_container);
	if (cells.size() == 0) {
		std::cout<<"********************************"<<std::endl;
		std::cout<<"Detector_full_y2022::SetupTHGEM1Mapping:Warning:\n"
				"\tTHGEM0 cell is not found. Geometry without THGEM0 mapping is used."<<std::endl;
		return;
	}
	if (cells.size() > 1) {
		std::cout<<"********************************"<<std::endl;
		std::cout<<"Detector_full_y2022::SetupTHGEM1Mapping:Warning:\n"
				"\tfound several cell physical volumes (name \""<<phys_THGEM0_cell_LAr->GetName()<<")\"!"<<std::endl;
		std::cout<<"\tUsing only the first one."<<std::endl;
	}
	if (containers.size() == 0) {
		std::cout<<"********************************"<<std::endl;
		std::cout<<"Detector_full_y2022::SetupTHGEM1Mapping:Warning:\n"
				"\tTHGEM0 container is not found. Geometry without THGEM0 mapping is used."<<std::endl;
		return;
	}
	if (containers.size() > 1) {
		std::cout<<"********************************"<<std::endl;
		std::cout<<"Detector_full_y2022::SetupTHGEM1Mapping:Warning:\n"
				"\tfound several THGEM0 containers' physical volumes (name \""<<phys_THGEM0_container->GetName()<<")\"!"<<std::endl;
		std::cout<<"\tUsing only the first one."<<std::endl;
	}
	{
		G4ThreeVector cell_pos = cells[0].fFoundObjectTransformation.getTranslation();
		G4ThreeVector thgem0_pos = containers[0].fFoundObjectTransformation.getTranslation();
		G4VSolid* cell_box = phys_THGEM0_cell_LAr->GetLogicalVolume()->GetSolid();
		G4VSolid* thgem0_box = phys_THGEM0_container->GetLogicalVolume()->GetSolid();
		G4ThreeVector bMin, bMax;
		cell_box->BoundingLimits(bMin, bMax);
		G4ThreeVector cell_sizes = bMax - bMin;
		thgem0_box->BoundingLimits(bMin, bMax);
		G4ThreeVector thgem0_sizes = bMax - bMin;
		HexagonalMapping thgem0_map("THGEM0", thgem0_pos, cell_pos, thgem0_sizes, cell_sizes, dims->is_NBrS_in_THGEM0);
		thgem0_map.AddTrigger(MappingTrigger(phys_THGEM0_cell_LAr, false, "THGEM0_leaving_cell"));
		thgem0_map.AddTrigger(MappingTrigger(phys_THGEM0_container, true, "THGEM0_entering_container"));
		gData.mapping_manager.AddMapping(thgem0_map);
	}
}

