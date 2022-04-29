#include <geant4/detector/Detector_THGEM1_SiPM_shading.hh>

Detector_THGEM1_SiPM_shading::Detector_THGEM1_SiPM_shading() :
  VDetectorConstruction()
{}

Detector_THGEM1_SiPM_shading::~Detector_THGEM1_SiPM_shading()
{}

G4VPhysicalVolume * Detector_THGEM1_SiPM_shading::Construct()
{
  SetSizeAndPosition();
  defineMaterials();
  defineSurfaces();

	solidWorld = new G4Box("sworld", HalfWorldLength, HalfWorldLength, HalfWorldLength);
	logicWorld = new G4LogicalVolume(solidWorld, matLAr, "lWorld", 0, 0, 0);
	//  Must place the World Physical volume unrotated at (0,0,0).
	physiWorld = new G4PVPlacement(0, G4ThreeVector(), logicWorld, "pWorld", 0, // its mother  volume
		false, 0);
	// Set user cuts to avoid deadlocks
	G4double maxStep = 1.0*m, maxLength = 1.0*m, maxTime = 10.0*ns, minEkin = 0.2*eV;
	logicWorld->SetUserLimits(new G4UserLimits(maxStep, maxLength, maxTime, minEkin));

	//-------------------------------------------------------------------------------
	CreateTHGEM1Cell();
	G4VPhysicalVolume* phys_THGEM1_cell = new G4PVPlacement(0, position_SingleTHGEMCell, logic_THGEM1_cell,
		"phys_THGEM1_cell_isolated", logicWorld, false, 0, fCheckOverlaps);

	//--------------------------------------------------------------------------------
	// Create THGEM1
  G4Box* solid_THGEM1_container = new G4Box("solid_THGEM1_container", x_size_THGEM1_container / 2.0, y_size_THGEM1_container / 2.0, z_size_THGEM1_container / 2.0);
  G4LogicalVolume* logic_THGEM1_container = new G4LogicalVolume(solid_THGEM1_container, matLAr, "logic_THGEM1_container", 0, 0, 0);
  G4VPhysicalVolume* phys_THGEM1_container = new G4PVPlacement(0, position_THGEM1_container, logic_THGEM1_container, gPars::det_dims->THGEM1_cell_container_name,
      logicWorld, false, 0, fCheckOverlaps);

  G4Box* solid_THGEM1_dummy = new G4Box("solid_THGEM1_dummy", x_size_THGEM1_container / 2.0, y_size_THGEM1_container / 2.0,
      gPars::det_dims->THGEM1_width_total / 2.0);
  G4LogicalVolume* logic_THGEM1_dummy = new G4LogicalVolume(solid_THGEM1_dummy, matFR4, "logic_THGEM1_dummy", 0, 0, 0);
  G4VPhysicalVolume* phys_THGEM1_dummy = new G4PVPlacement(0, G4ThreeVector(0,0,0), logic_THGEM1_dummy, "phys_THGEM1_dummy",
      logic_THGEM1_container, false, 0, fCheckOverlaps);

  //--------------------------------------------------------------------------------
  // Create sensor
  G4Box* solid_sensor = new G4Box("solid_sensor", x_size_sensor / 2.0, y_size_sensor / 2.0, z_size_sensor / 2.0);
  logic_sensor = new G4LogicalVolume(solid_sensor, matLAr, "logic_sensor", 0, 0, 0);
  G4VPhysicalVolume* phys_sensor = new G4PVPlacement(0, position_sensor, logic_sensor, gPars::det_dims->SiPM_device_name,
      logicWorld, false, 0, fCheckOverlaps);

  //--------------------------------------------------------------------------------
  // Setting surfaces
	G4LogicalBorderSurface* THGEM1_cell_isolation = new G4LogicalBorderSurface("THGEM1_cell_isolation_surface", physiWorld, phys_THGEM1_cell, AbsorberMaterial);
  G4LogicalBorderSurface* THGEM1_cell_isolation1 = new G4LogicalBorderSurface("THGEM1_cell_isolation_surface1", phys_THGEM1_cell, physiWorld, AbsorberMaterial);

  G4LogicalSkinSurface* sur_Cu = new G4LogicalSkinSurface("Copper_surface", logic_THGEM1_cell_copper, AbsorberMaterial);
  G4LogicalSkinSurface* sur_FR4 = new G4LogicalSkinSurface("FR4_surface", logic_THGEM1_cell_FR4, AbsorberMaterial);
  G4LogicalSkinSurface* sur_FR4_dummy = new G4LogicalSkinSurface("FR4_dummy_surface", logic_THGEM1_dummy, AbsorberMaterial);

  G4LogicalSkinSurface* sur_sensor = new G4LogicalSkinSurface("sensor_surface", logic_sensor, SiPM_OpticalSurface);

  //--------------------------------------------------------------------------------
  // Setting visualization
  G4VisAttributes Invisible(G4Colour(1, 1, 1, 0.0));
  Invisible.SetVisibility(false);
  G4VisAttributes AlmostInvisible(G4Colour(1, 1, 1, 0.05));
  G4VisAttributes LAr_VisAtt(G4Colour(0.6, 0.6, 1.0, 0.0));
  G4VisAttributes FR4_VisAtt(G4Colour(0.8, 0.85, 0.11, 0.8));
  G4VisAttributes Cu_VisAtt(G4Colour(0.8, 0.45, 0.2, 0.9));
  G4VisAttributes Sensor_VisAtt(G4Colour(0.4, 0.4, 0.4, 0.7));

  logic_THGEM1_dummy->SetVisAttributes(FR4_VisAtt);
  logic_sensor->SetVisAttributes(Sensor_VisAtt);

  // world
  logic_THGEM1_container->SetVisAttributes(Invisible);
  logicWorld->SetVisAttributes(Invisible);

	return physiWorld;
}

void Detector_THGEM1_SiPM_shading::ConstructSDandField()
{
  DetectorSensor *thePhotoDetector = new DetectorSensor("/detector/sensitiveDetector");
  G4SDManager::GetSDMpointer()->AddNewDetector(thePhotoDetector);
  SetSensitiveDetector(logic_sensor, thePhotoDetector); // takes ownership of DetectorSensor*
}

void Detector_THGEM1_SiPM_shading::SetSizeAndPosition()
{
	HalfWorldLength = 10 * cm / 2.0;

	x_size_THGEM1_container = 5 * cm;
  y_size_THGEM1_container = 5 * cm;
  z_size_THGEM1_container = gPars::det_dims->THGEM1_container_width;

  position_SingleTHGEMCell = G4ThreeVector(4 * cm, 4 * cm, -4 * cm);
  position_THGEM1_container = G4ThreeVector(0 * mm, 0 * mm, 0 * mm);

  x_size_sensor = 9 * cm;
  y_size_sensor = 9 * cm;
  z_size_sensor = 1 * mm;
  position_sensor = G4ThreeVector(0 * mm, 0 * mm, z_size_THGEM1_container / 2.0 + 1 * mm);

  double z_bottom_THGEM1 = position_THGEM1_container.x() - gPars::det_dims->THGEM1_width_total / 2.0;
  gPars::det_dims->THGEM1_single_cell_position = position_SingleTHGEMCell;
  gPars::det_dims->n_PMTs = 0;
  gPars::det_dims->n_SiPMs = 0;
  gPars::det_dims->THGEM1_hole_center = G4ThreeVector(0, 0, 0);
  gPars::det_dims->THGEM1_center = G4ThreeVector(0, 0, 0);
  gPars::det_dims->EL_gap_center = G4ThreeVector(DBL_MAX, DBL_MAX, DBL_MAX);
  gPars::det_dims->Cathode_top_center = G4ThreeVector(DBL_MAX, DBL_MAX, DBL_MAX);
}
