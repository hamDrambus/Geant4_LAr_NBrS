#include <geant4/detector/Detector_THGEM1_detailed.hh>

Detector_THGEM1_detailed::Detector_THGEM1_detailed() :
  VDetectorConstruction()
{}

Detector_THGEM1_detailed::~Detector_THGEM1_detailed()
{}

G4VPhysicalVolume * Detector_THGEM1_detailed::Construct()
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
  phys_THGEM1_container = new G4PVPlacement(0, position_THGEM1_container, logic_THGEM1_container, "phys_THGEM1_cell_container",
      logicWorld, false, 0, fCheckOverlaps);

  // Parameterisation does not work with arbitrary transformations.
  // Basically only reflections cannot be done which is exactly what is
  // required for THGEM parameterisation (can't set G4Transform3D for G4VPhysicalVolume
  // in G4VPVParameterisation::ComputeTransformation).
  // So G4AssemblyVolume is used instead.

  // Define one layer as one assembly volume
  G4AssemblyVolume* assemblyTHGEM = new G4AssemblyVolume();
  HexagonalMapping mapping(position_THGEM1_container, position_SingleTHGEMCell,
      G4ThreeVector(x_size_THGEM1_container, y_size_THGEM1_container, z_size_THGEM1_container),
      G4ThreeVector(gPars::det_dims->THGEM1_hole_pitch / 2.0, gPars::det_dims->THGEM1_hole_pitch * sqrt(3) / 2.0, z_size_THGEM1_container));

  for (int i = 0, i_end_ = mapping.GetNcells(); i!=i_end_; ++i) {
    std::pair<int, int> inds = mapping.GetIndices(i);
    G4Transform3D Tr = mapping.GetCellRelativePointTransform(inds.first, inds.second);
    assemblyTHGEM->AddPlacedVolume(logic_THGEM1_cell_LAr, Tr);
  }
  G4Transform3D Tr = G4Translate3D(G4ThreeVector(0, 0, 0)); //relative position
  assemblyTHGEM->MakeImprint(logic_THGEM1_container, Tr);

  //--------------------------------------------------------------------------------
  // Setting surfaces (not used)
	G4LogicalBorderSurface* THGEM1_cell_isolation = new G4LogicalBorderSurface("THGEM1_cell_isolation_surface", physiWorld, phys_THGEM1_cell, AbsorberMaterial);
  G4LogicalBorderSurface* THGEM1_cell_isolation1 = new G4LogicalBorderSurface("THGEM1_cell_isolation_surface1", phys_THGEM1_cell, physiWorld, AbsorberMaterial);

  G4LogicalSkinSurface* sur_Cu = new G4LogicalSkinSurface("Copper_surface", logic_THGEM1_cell_copper, Cu_THGEM);
  G4LogicalSkinSurface* sur_FR4 = new G4LogicalSkinSurface("FR4_surface", logic_THGEM1_cell_FR4, FR4_unified);

  //--------------------------------------------------------------------------------
  // Setting visualization
  G4VisAttributes Invisible(G4Colour(1, 1, 1, 0.0));
  Invisible.SetVisibility(false);
  G4VisAttributes AlmostInvisible(G4Colour(1, 1, 1, 0.05));
  G4VisAttributes LAr_VisAtt(G4Colour(0.6, 0.6, 1.0, 0.0));
  G4VisAttributes FR4_VisAtt(G4Colour(0.8, 0.85, 0.11, 0.8));
  G4VisAttributes Cu_VisAtt(G4Colour(0.8, 0.45, 0.2, 0.9));

  // world
  logic_THGEM1_container->SetVisAttributes(Invisible);
  logicWorld->SetVisAttributes(Invisible);

  SetupTHGEM1Mapping();
	return physiWorld;
}

void Detector_THGEM1_detailed::SetSizeAndPosition()
{
	HalfWorldLength = 2 * cm / 2.0;

	x_size_THGEM1_container = 3 * mm;
  y_size_THGEM1_container = 5.5 * mm;
  z_size_THGEM1_container = gPars::det_dims->THGEM1_container_width;

  position_SingleTHGEMCell = G4ThreeVector(8 * mm, 8 * mm, 8 * mm);
  position_THGEM1_container = G4ThreeVector(0 * mm, 0 * mm, 0 * mm);

  gPars::det_dims->THGEM1_single_cell_position = position_SingleTHGEMCell;
  gPars::det_dims->n_PMTs = 0;
  gPars::det_dims->n_SiPMs = 0;
  gPars::det_dims->THGEM1_hole_center = G4ThreeVector(0, 0, 0);
  gPars::det_dims->THGEM1_center = G4ThreeVector(0, 0, 0);
  gPars::det_dims->EL_gap_center = G4ThreeVector(DBL_MAX, DBL_MAX, DBL_MAX);
  gPars::det_dims->Cathode_top_center = G4ThreeVector(DBL_MAX, DBL_MAX, DBL_MAX);
}
