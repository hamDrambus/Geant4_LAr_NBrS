#include "DetectorConstructionTHGEM1.hh"

DetectorConstructionTHGEM1::DetectorConstructionTHGEM1() : fCheckOverlaps(gPars::general.check_geometry_overlap)
{
	SetSizeAndPosition();
}

DetectorConstructionTHGEM1::~DetectorConstructionTHGEM1()
{}

G4VPhysicalVolume * DetectorConstructionTHGEM1::Construct()
{
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
  G4VPhysicalVolume* phys_THGEM1_container = new G4PVPlacement(0, position_THGEM1_container, logic_THGEM1_container, gPars::det_dims.THGEM1_cell_container_name,
      logicWorld, false, 0, fCheckOverlaps);

  //--------------------------------------------------------------------------------
  // Setting visualization
  G4VisAttributes Invisible(G4Colour(1, 1, 1, 0.0));
  Invisible.SetVisibility(false);
  G4VisAttributes AlmostInvisible(G4Colour(1, 1, 1, 0.05));
  G4VisAttributes LAr_VisAtt(G4Colour(0.6, 0.6, 1.0, 0.0));
  G4VisAttributes FR4_VisAtt(G4Colour(0.8, 0.85, 0.11, 0.8));
  G4VisAttributes Cu_VisAtt(G4Colour(0.8, 0.45, 0.2, 0.9));

  // Separate THGEM hole
  logic_THGEM1_cell_copper->SetVisAttributes(Cu_VisAtt);
  logic_THGEM1_cell_LAr->SetVisAttributes(Invisible);
  logic_THGEM1_cell_FR4->SetVisAttributes(FR4_VisAtt);

  // world
  logic_THGEM1_cell->SetVisAttributes(LAr_VisAtt);
  logic_THGEM1_container->SetVisAttributes(Invisible);
  logicWorld->SetVisAttributes(Invisible);

  // Parameterisation does not work with arbitrary transformations.
  // Basically only reflections cannot be done which is exactly what is
  // required for THGEM parameterisation (can't set G4Transform3D for G4VPhysicalVolume
  // in G4VPVParameterisation::ComputeTransformation).
  // So G4AssemblyVolume is used instead.

  // Define one layer as one assembly volume
  G4AssemblyVolume* assemblyTHGEM = new G4AssemblyVolume();
  HexagonalMapping mapping(position_THGEM1_container, position_SingleTHGEMCell,
      G4ThreeVector(x_size_THGEM1_container, y_size_THGEM1_container, z_size_THGEM1_container),
      G4ThreeVector(gPars::det_dims.THGEM1_hole_pitch * mm / 2.0, gPars::det_dims.THGEM1_hole_pitch * sqrt(3) * mm / 2.0, z_size_THGEM1_container));

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

	return physiWorld;
}

void DetectorConstructionTHGEM1::CreateTHGEM1Cell() //Must be the same as in gmsh-elmer simulation
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
  G4Tubs* solid_THGEM1_diel_hole1 = new G4Tubs("solid_THGEM1_diel_hole", 0, radius, diel_size_z / 1.9, 179.*deg, 271.*deg);
  G4Tubs* solid_THGEM1_diel_hole2 = new G4Tubs("solid_THGEM1_diel_hole", 0, radius, diel_size_z / 1.9, -1.*deg, 91.*deg);
  G4SubtractionSolid* solid_THGEM1_diel_tmp = new G4SubtractionSolid("solid_THGEM1_diel_tmp", solid_THGEM1_diel_box, solid_THGEM1_diel_hole1, 0, hole_1_pos);
  G4SubtractionSolid* solid_THGEM1_diel = new G4SubtractionSolid("solid_THGEM1_diel", solid_THGEM1_diel_tmp, solid_THGEM1_diel_hole2, 0, hole_2_pos);
  logic_THGEM1_cell_FR4 = new G4LogicalVolume(solid_THGEM1_diel, G4Material::GetMaterial("FR4"), "logic_THGEM1_cell_FR4", 0, 0, 0);
  G4VPhysicalVolume* phys_THGEM1_cell_FR4 = new G4PVPlacement(0, zero, logic_THGEM1_cell_FR4, "phys_THGEM1_cell_FR4",
      logic_THGEM1_cell_LAr, false, 0, fCheckOverlaps);

  G4Box* solid_THGEM1_cu_box = new G4Box("solid_THGEM1_cu_box", cell_size_x / 2.0, cell_size_y / 2.0, cu_size_z / 2.0);
  G4Tubs* solid_THGEM1_cu_hole1 = new G4Tubs("solid_THGEM1_cu_hole", 0, radius_cu, cu_size_z / 1.9, 179.*deg, 271.*deg);
  G4Tubs* solid_THGEM1_cu_hole2 = new G4Tubs("solid_THGEM1_cu_hole", 0, radius_cu, cu_size_z / 1.9, -1.*deg, 91.*deg);
  G4SubtractionSolid* solid_THGEM1_cu_tmp = new G4SubtractionSolid("solid_THGEM1_cu_tmp", solid_THGEM1_cu_box, solid_THGEM1_cu_hole1, 0, hole_1_pos);
  G4SubtractionSolid* solid_THGEM1_cu = new G4SubtractionSolid("solid_THGEM1_cu", solid_THGEM1_cu_tmp, solid_THGEM1_cu_hole2, 0, hole_2_pos);
  logic_THGEM1_cell_copper = new G4LogicalVolume(solid_THGEM1_cu, G4Material::GetMaterial("FR4"), "logic_THGEM1_cell_copper", 0, 0, 0);
  G4VPhysicalVolume* phys_THGEM1_cell_copper_top = new G4PVPlacement(0, cu_top_pos, logic_THGEM1_cell_copper, "phys_THGEM1_cell_copper_top",
      logic_THGEM1_cell_LAr, false, 0, fCheckOverlaps);
  G4VPhysicalVolume* phys_THGEM1_cell_copper_bot = new G4PVPlacement(0, cu_bot_pos, logic_THGEM1_cell_copper, "phys_THGEM1_cell_copper_bot",
      logic_THGEM1_cell_LAr, false, 0, fCheckOverlaps);
}

void DetectorConstructionTHGEM1::defineSurfaces()
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
	FR4_MaterialProperty->AddProperty("REFLECTIVITY", ener, gPars::general.no_reflections ? zero : FR4_Materialrefl, 2);
	FR4_MaterialProperty->AddProperty("EFFICIENCY", ener, FR4_Materialeff, 2);
	FR4_unified->SetMaterialPropertiesTable(FR4_MaterialProperty);
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
	Cu_THGEM_MaterialProperty->AddProperty("REFLECTIVITY", ener, gPars::general.no_reflections ? zero : Cu_THGEM_Materialrefl, 2);
	Cu_THGEM_MaterialProperty->AddProperty("EFFICIENCY", ener, Cu_THGEM_Materialeff, 2);
	Cu_THGEM->SetMaterialPropertiesTable(Cu_THGEM_MaterialProperty);
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

void DetectorConstructionTHGEM1::SetSizeAndPosition()
{
	HalfWorldLength = 2 * cm / 2.0;

	x_size_THGEM1_container = 3 * mm;
  y_size_THGEM1_container = 5.5 * mm;
  z_size_THGEM1_container = gPars::det_dims.THGEM1_container_width * mm;

  position_SingleTHGEMCell = G4ThreeVector(8 * mm, 8 * mm, 8 * mm);
  position_THGEM1_container = G4ThreeVector(0 * mm, 0 * mm, 0 * mm);

  gPars::det_dims.THGEM1_active_area_size = x_size_THGEM1_container;
  gPars::det_dims.z_bottom_THGEM1 = position_THGEM1_container.x()/ mm - gPars::det_dims.THGEM1_width_total / 2.0;
  gPars::det_dims.THGEM1_single_cell_position = position_SingleTHGEMCell / mm;
  gPars::det_dims.n_PMTs = 0;
  gPars::det_dims.n_SiPMs_rows = 0;
  gPars::general.THGEM1_hole_center = G4ThreeVector(0, 0, 0);
  gPars::general.EL_gap_center = G4ThreeVector(0, 0, -3.2);
  gPars::source.z_center = gPars::det_dims.z_bottom_THGEM1 - 2.9;
}

void DetectorConstructionTHGEM1::defineMaterials()
{
	G4NistManager* man = G4NistManager::Instance();

	G4Element *C = man->FindOrBuildElement("C");
	G4Element *H = man->FindOrBuildElement("H");
	G4Element *O = man->FindOrBuildElement("O");
	G4Element *Ar = man->FindOrBuildElement("Ar");
	G4Element *Cu = man->FindOrBuildElement("Cu");

	G4double a;  // atomic mass
	G4double z;  // atomic number
	G4double density;

	const G4int numentries = 2;
	G4double energies[numentries] = { 0.1*eV, 10.0*eV };
	//------------------------------
	// LAr
	matLAr = new G4Material("LAr", 1.400*g / cm3, 1);
	matLAr->AddElement(Ar, 1);
	const G4int LAr_numentries = 12;
	G4double LAr_energies[LAr_numentries] = {1*eV, 2.95*eV, 3.2*eV, 3.4*eV, 3.8*eV, 4*eV, 5.63*eV, 6.89*eV, 7.75*eV, 8.86*eV, 9.69*eV, 10.33*eV};
	G4double LAr_rindices[LAr_numentries] = {1.23, 1.23, 1.23, 1.23, 1.23, 1.23, 1.26, 1.29, 1.31, 1.34, 1.36, 1.45}; //doi:10.1088/1748-0221/15/09/P09009
	//G4double LAr_absorpti[LAr_numentries] = { 2 * m, 2 * m }; // avoid infinite light-paths
	G4double LAr_absorpti[LAr_numentries] = { 20000 * m, 20000 * m };
	G4MaterialPropertiesTable* LAr_prop = new G4MaterialPropertiesTable();
	LAr_prop->AddProperty("ABSLENGTH", energies, LAr_absorpti, numentries);
	LAr_prop->AddProperty("RINDEX", LAr_energies, LAr_rindices, LAr_numentries);
	matLAr->SetMaterialPropertiesTable(LAr_prop);
	//------------------------------
	//FR4
	//I don't know a real chemical composition. So it's a dummy
	matFr4 = new G4Material("FR4", 1.850*g / cm3 /* * (1 - 0.28)*/ /*THGEM transparence*/, 3);
	matFr4->AddElement(C, 5); matFr4->AddElement(O, 2); matFr4->AddElement(H, 8);
	const G4int numentries_FR4 = 2;
	G4double energies_FR4[numentries_FR4] = { 0.1*eV, 10.0*eV };
	G4double rindices_FR4[numentries_FR4] = { 1.5, 1.5 };
	G4double absorpti_FR4[numentries_FR4] = { 1 * um, 1 * um };
	G4MaterialPropertiesTable* prop_FR4 = new G4MaterialPropertiesTable();
	prop_FR4->AddProperty("ABSLENGTH", energies_FR4, absorpti_FR4, numentries_FR4);
	prop_FR4->AddProperty("RINDEX", energies_FR4, rindices_FR4, numentries_FR4);
	matFr4->SetMaterialPropertiesTable(prop_FR4);
}
