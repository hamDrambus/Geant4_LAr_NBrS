#include <geant4/detector/VDetectorConstruction.hh>

VDetectorConstruction::VDetectorConstruction() : fCheckOverlaps(gPars::general.check_geometry_overlap),
    phys_THGEM1_cell_LAr(nullptr), phys_THGEM1_container(nullptr)
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

VDetectorConstruction::~VDetectorConstruction()
{
    delete rotX_90;
    delete rotY_90;
    delete rotZ_90;
    delete rotZ_180;
    delete rotZ_270;
}

DataVector VDetectorConstruction::convert_lambda_to_E(DataVector &in) const
{
    DataVector out = in; // To copy interpolation order, etc.
    out.clear();
    for (std::size_t i = 0, i_end_ = in.size(); i!=i_end_; ++i) {
        double lambda = hc / in[i].first;
        out.insert(lambda, in[i].second);
    }
    return out;
}

std::vector<G4PhysicalVolumesSearchScene::Findings> VDetectorConstruction::LocatePV(G4VPhysicalVolume* volume)
{
    std::vector<G4PhysicalVolumesSearchScene::Findings> findingsVector;
    if (physiWorld == nullptr) {
        std::cerr<<"VDetectorConstruction::LocatePV:Error:\n"
                "\tNULL world volume. This function must be called after or during detector construction in G4VPhysicalVolume::Construct()." << std::endl;
        return findingsVector;
    }
    if (volume == nullptr)
        return findingsVector;
    G4PhysicalVolumeModel searchModel(physiWorld); // Unlimited depth.
    G4ModelingParameters mp; // Default - no culling.
    searchModel.SetModelingParameters(&mp);
    G4PhysicalVolumesSearchScene searchScene(&searchModel, volume->GetName());
    searchModel.DescribeYourselfTo(searchScene); // Initiate search.
    for (const auto& findings: searchScene.GetFindings())
        findingsVector.push_back(findings);
    return findingsVector;
}

void VDetectorConstruction::SetupTHGEMsMapping()
{
    std::vector<G4PhysicalVolumesSearchScene::Findings> cells = LocatePV(phys_THGEM1_cell_LAr);
    if (cells.size() == 0) {
        std::cout<<"********************************"<<std::endl;
        std::cout<<"VDetectorConstruction::SetupTHGEM1Mapping:Warning:\n"
                "\tTHGEM1 cell is not found. Geometry without THGEM1 mapping is used."<<std::endl;
        return;
    }
    if (cells.size() > 1) {
        std::cout<<"********************************"<<std::endl;
        std::cout<<"VDetectorConstruction::SetupTHGEM1Mapping:Warning:\n"
                "\tfound several cell physical volumes (name \""<<phys_THGEM1_cell_LAr->GetName()<<")\"!"<<std::endl;
        std::cout<<"\tUsing only the first one."<<std::endl;
    }
    std::vector<G4PhysicalVolumesSearchScene::Findings> containers = LocatePV(phys_THGEM1_container);
    if (containers.size() == 0) {
        std::cout<<"********************************"<<std::endl;
        std::cout<<"VDetectorConstruction::SetupTHGEM1Mapping:Warning:\n"
                "\tTHGEM1 container is not found. Geometry without THGEM1 mapping is used."<<std::endl;
        return;
    }
    if (containers.size() > 1) {
        std::cout<<"********************************"<<std::endl;
        std::cout<<"VDetectorConstruction::SetupTHGEM1Mapping:Warning:\n"
                "\tfound several THGEM containers' physical volumes (name \""<<phys_THGEM1_container->GetName()<<")\"!"<<std::endl;
        std::cout<<"\tUsing only the first one."<<std::endl;
    }
    G4ThreeVector cell_pos = cells[0].fFoundObjectTransformation.getTranslation();
    G4ThreeVector thgem1_pos = containers[0].fFoundObjectTransformation.getTranslation();
    G4VSolid* cell_box = phys_THGEM1_cell_LAr->GetLogicalVolume()->GetSolid();
    G4VSolid* thgem1_box = phys_THGEM1_container->GetLogicalVolume()->GetSolid();
    G4ThreeVector bMin, bMax;
    cell_box->BoundingLimits(bMin, bMax);
    G4ThreeVector cell_sizes = bMax - bMin;
    thgem1_box->BoundingLimits(bMin, bMax);
    G4ThreeVector thgem1_sizes = bMax - bMin;
    HexagonalMapping thgem1_map("THGEM1", thgem1_pos, cell_pos, thgem1_sizes, cell_sizes, true);
    thgem1_map.AddTrigger(MappingTrigger(phys_THGEM1_cell_LAr, false, "THGEM1_leaving_cell"));
    thgem1_map.AddTrigger(MappingTrigger(phys_THGEM1_container, true, "THGEM1_entering_container"));
    gData.mapping_manager.AddMapping(thgem1_map);
}

void VDetectorConstruction::CreateTHGEM1Cell(G4Material* medium) //Must be the same as in gmsh-elmer simulation
{
    double cell_size_x = gPars::det_dims->THGEM1_hole_pitch / 2.0;
    double cell_size_y = cell_size_x * std::sqrt(3.0);
    double cell_size_z = gPars::det_dims->THGEM1_container_width;
    double diel_size_z = gPars::det_dims->THGEM1_dielectric_thickness;
    double radius = gPars::det_dims->THGEM1_hole_radius;
    double radius_center = gPars::det_dims->THGEM1_dielectric_radius;
    double radius_cu = radius + gPars::det_dims->THGEM1_hole_rim;
    double cu_size_z = gPars::det_dims->THGEM1_copper_thickness;

    G4ThreeVector zero(0.0, 0.0, 0.0);
    G4ThreeVector hole_1_pos(cell_size_x / 2.0, cell_size_y / 2.0, 0.0);
    G4ThreeVector hole_2_pos(-cell_size_x / 2.0, -cell_size_y / 2.0, 0.0);
    G4ThreeVector cu_top_pos(0.0, 0.0, diel_size_z / 2.0 + cu_size_z / 2.0);
    G4ThreeVector cu_bot_pos(0.0, 0.0, -diel_size_z / 2.0 - cu_size_z / 2.0);

    G4Box* solid_THGEM1_cell_isolation = new G4Box("solid_THGEM1_cell_isolation", cell_size_x, cell_size_y, cell_size_z); //so that cell is in the same material
    logic_THGEM1_cell = new G4LogicalVolume(solid_THGEM1_cell_isolation, medium, "logic_THGEM1_cell_isolation", 0, 0, 0);

    G4Box* solid_THGEM1_cell_LAr = new G4Box("solid_THGEM1_cell_LAr", cell_size_x / 2.0, cell_size_y / 2.0, cell_size_z / 2.0);
    logic_THGEM1_cell_LAr = new G4LogicalVolume(solid_THGEM1_cell_LAr, medium, "logic_THGEM1_cell_LAr", 0, 0, 0);
    phys_THGEM1_cell_LAr = new G4PVPlacement(0, zero, logic_THGEM1_cell_LAr, "phys_THGEM1_cell",
        logic_THGEM1_cell, false, 0, fCheckOverlaps);

    G4Box* solid_THGEM1_diel_box = new G4Box("solid_THGEM1_diel_box", cell_size_x / 2.0, cell_size_y / 2.0, diel_size_z / 2.0);
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
    double rs[] = {0, 0, 0};
    G4Polycone* solid_diel_hole1 = new G4Polycone("solid_diel_hole1", 177.*deg, 273.*deg, 3, Zs, rs, Rs);
    G4Polycone* solid_diel_hole2 = new G4Polycone("solid_diel_hole2", -3.*deg, 93.*deg, 3, Zs, rs, Rs);
    G4SubtractionSolid* solid_THGEM1_diel_tmp = new G4SubtractionSolid("solid_THGEM1_diel_tmp", solid_THGEM1_diel_box, solid_diel_hole1, 0, hole_1_pos);
    G4SubtractionSolid* solid_THGEM1_diel = new G4SubtractionSolid("solid_THGEM1_diel", solid_THGEM1_diel_tmp, solid_diel_hole2, 0, hole_2_pos);
    logic_THGEM1_cell_FR4 = new G4LogicalVolume(solid_THGEM1_diel, matFR4, "logic_THGEM1_cell_FR4", 0, 0, 0);
    G4VPhysicalVolume* phys_THGEM1_cell_FR4 = new G4PVPlacement(0, zero, logic_THGEM1_cell_FR4, "phys_THGEM1_cell_FR4",
        logic_THGEM1_cell_LAr, false, 0, fCheckOverlaps);

    G4Box* solid_THGEM1_cu_box = new G4Box("solid_THGEM1_cu_box", cell_size_x / 2.0, cell_size_y / 2.0, cu_size_z / 2.0);
    G4Tubs* solid_THGEM1_cu_hole1 = new G4Tubs("solid_THGEM1_cu_hole", 0, radius_cu, cu_size_z / 1.9, 177.*deg, 273.*deg);
    G4Tubs* solid_THGEM1_cu_hole2 = new G4Tubs("solid_THGEM1_cu_hole", 0, radius_cu, cu_size_z / 1.9, -3.*deg, 93.*deg);
    G4SubtractionSolid* solid_THGEM1_cu_tmp = new G4SubtractionSolid("solid_THGEM1_cu_tmp", solid_THGEM1_cu_box, solid_THGEM1_cu_hole1, 0, hole_1_pos);
    G4SubtractionSolid* solid_THGEM1_cu = new G4SubtractionSolid("solid_THGEM1_cu", solid_THGEM1_cu_tmp, solid_THGEM1_cu_hole2, 0, hole_2_pos);
    logic_THGEM1_cell_copper = new G4LogicalVolume(solid_THGEM1_cu, matFR4, "logic_THGEM1_cell_copper", 0, 0, 0);
    G4VPhysicalVolume* phys_THGEM1_cell_copper_top = new G4PVPlacement(0, cu_top_pos, logic_THGEM1_cell_copper, "phys_THGEM1_cell_copper_top",
        logic_THGEM1_cell_LAr, false, 0, fCheckOverlaps);
    G4VPhysicalVolume* phys_THGEM1_cell_copper_bot = new G4PVPlacement(0, cu_bot_pos, logic_THGEM1_cell_copper, "phys_THGEM1_cell_copper_bot",
        logic_THGEM1_cell_LAr, false, 0, fCheckOverlaps);

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
    logic_THGEM1_cell_copper->SetVisAttributes(Cu_VisAtt);
    logic_THGEM1_cell_LAr->SetVisAttributes(Invisible);
    logic_THGEM1_cell_FR4->SetVisAttributes(FR4_VisAtt);
    logic_THGEM1_cell->SetVisAttributes(AlmostInvisible);
}

void VDetectorConstruction::defineSurfaces()
{
    //*********************************************************************************
    //unified model
    G4double ener[2] = {.1*eV, 10.*eV};
    G4double zero[2] = {0.0, 0.0};
    //-------------------------------------------------------------------------------
    //Liquid Ar - gaseous Ar surface
    LAr_OpticalSurface = new G4OpticalSurface("LAr_OpticalSurface");
    LAr_OpticalSurface->SetType(dielectric_dielectric);
    LAr_OpticalSurface->SetModel(unified);
    LAr_OpticalSurface->SetFinish(ground);
    LAr_OpticalSurface->SetSigmaAlpha(gPars::det_opt.LAr_SigmaAlpha);//alpha in degrees, from 0 to 90.
    G4MaterialPropertiesTable *LAr_MaterialProperty = new G4MaterialPropertiesTable();
    G4double LAr_Materialeff[2] = { 0, 0 };
    LAr_MaterialProperty->AddProperty("EFFICIENCY", ener, LAr_Materialeff, 2);
    LAr_OpticalSurface->SetMaterialPropertiesTable(LAr_MaterialProperty);
    //-------------------------------------------------------------------------------
    //PMMA (Acrylic), both across SiPMs and across PMTs
    PMMA_OpticalSurface = new G4OpticalSurface("PMMA_OpticalSurface");
    PMMA_OpticalSurface->SetType(dielectric_dielectric);
    PMMA_OpticalSurface->SetModel(unified);
    PMMA_OpticalSurface->SetFinish(ground);
    PMMA_OpticalSurface->SetSigmaAlpha(gPars::det_opt.PMMA_SigmaAlpha);//alpha in degrees, from 0 to 90.
    G4MaterialPropertiesTable *PMMA_MaterialProperty = new G4MaterialPropertiesTable();
    G4double PMMA_Materialeff[2] = { 0, 0 };
    PMMA_MaterialProperty->AddProperty("EFFICIENCY", ener, PMMA_Materialeff, 2);
    PMMA_OpticalSurface->SetMaterialPropertiesTable(PMMA_MaterialProperty);
    //-------------------------------------------------------------------------------
    //FR4_unified
    FR4_unified = new G4OpticalSurface("FR4_unified");
    FR4_unified->SetType(dielectric_metal);
    FR4_unified->SetModel(unified);
    FR4_unified->SetFinish(ground);
    FR4_unified->SetSigmaAlpha(gPars::det_opt.FR4_SigmaAlpha);//alpha in degrees, from 0 to 90.
    G4MaterialPropertiesTable *FR4_MaterialProperty = new G4MaterialPropertiesTable();
    G4double FR4_Materialrefl[2] = {gPars::det_opt.FR4_reflectivity, gPars::det_opt.FR4_reflectivity};
    G4double FR4_Materialeff[2] = { 0, 0 };
    FR4_MaterialProperty->AddProperty("REFLECTIVITY", ener, gPars::general.no_reflections ? zero : FR4_Materialrefl, 2);
    FR4_MaterialProperty->AddProperty("EFFICIENCY", ener, FR4_Materialeff, 2);
    FR4_unified->SetMaterialPropertiesTable(FR4_MaterialProperty);
    //-------------------------------------------------------------------------------
    //Anode_Wire (Bronze-Berillium)
    Anode_wire_unified = new G4OpticalSurface("Anode_wire_unified");
    Anode_wire_unified->SetType(dielectric_metal);
    Anode_wire_unified->SetModel(unified);
    Anode_wire_unified->SetFinish(ground);
    Anode_wire_unified->SetSigmaAlpha(gPars::det_opt.Wire_SigmaAlpha * degree);//alpha in degrees, from 0 to 90.
    G4MaterialPropertiesTable *Anode_wire_MaterialProperty = new G4MaterialPropertiesTable();
    G4double Anode_wire_Materialrefl[2] = { gPars::det_opt.Wire_reflectivity, gPars::det_opt.Wire_reflectivity};
    G4double Anode_wire_Materialeff[2] = { 0, 0 };
    Anode_wire_MaterialProperty->AddProperty("REFLECTIVITY", ener, gPars::general.no_reflections ? zero : Anode_wire_Materialrefl, 2);
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
    G4double Cu_THGEM_Materialrefl[2] = {gPars::det_opt.Cu_reflectivity, gPars::det_opt.Cu_reflectivity};
    G4double Cu_THGEM_Materialeff[2] = {0, 0};
    Cu_THGEM_MaterialProperty->AddProperty("REFLECTIVITY", ener, gPars::general.no_reflections ? zero : Cu_THGEM_Materialrefl, 2);
    Cu_THGEM_MaterialProperty->AddProperty("EFFICIENCY", ener, Cu_THGEM_Materialeff, 2);
    Cu_THGEM->SetMaterialPropertiesTable(Cu_THGEM_MaterialProperty);
    //-------------------------------------------------------------------------------
    //Cu_Cathode
    const double r_factor = 1.0 - 0.28; // 1 - cathode THGEM transparency
    Cu_cathode = new G4OpticalSurface("Cu_cathode");
    Cu_cathode->SetType(dielectric_metal);
    Cu_cathode->SetModel(unified);
    Cu_cathode->SetFinish(ground);
    Cu_cathode->SetSigmaAlpha(gPars::det_opt.Cu_SigmaAlpha * degree);//alpha in degrees, from 0 to 90.
    G4MaterialPropertiesTable *Cu_Cathode_MaterialProperty = new G4MaterialPropertiesTable();
    G4double Cu_Cathode_Materialrefl[2] = {gPars::det_opt.Cu_reflectivity * r_factor, gPars::det_opt.Cu_reflectivity * r_factor};
    G4double Cu_Cathode_Materialeff[2] = {0, 0};
    Cu_Cathode_MaterialProperty->AddProperty("REFLECTIVITY", ener, gPars::general.no_reflections ? zero : Cu_Cathode_Materialrefl, 2);
    Cu_Cathode_MaterialProperty->AddProperty("EFFICIENCY", ener, Cu_Cathode_Materialeff, 2);
    Cu_cathode->SetMaterialPropertiesTable(Cu_Cathode_MaterialProperty);
    //-------------------------------------------------------------------------------
    //stainlessSteel
    stainlessSteel = new G4OpticalSurface("stainlessSteel");
    stainlessSteel->SetType(dielectric_metal);
    stainlessSteel->SetModel(unified);
    stainlessSteel->SetFinish(ground);
    stainlessSteel->SetSigmaAlpha(gPars::det_opt.StainlessSteel_SigmaAlpha * degree);//alpha in degrees, from 0 to 90.
    G4MaterialPropertiesTable *stainlessSteelMaterialProperty = new G4MaterialPropertiesTable();
    G4double stainlessSteelMaterialrefl[2] = {gPars::det_opt.StainlessSteel_reflectivity, gPars::det_opt.StainlessSteel_reflectivity};
    G4double stainlessSteelMaterialeff[2] = {0, 0};
    stainlessSteelMaterialProperty->AddProperty("REFLECTIVITY", ener, gPars::general.no_reflections ? zero : stainlessSteelMaterialrefl, 2);
    stainlessSteelMaterialProperty->AddProperty("EFFICIENCY", ener, stainlessSteelMaterialeff, 2);
    stainlessSteel->SetMaterialPropertiesTable(stainlessSteelMaterialProperty);
    //-----------------------------------------------------------------------------
    // PMT cathode
    PMT_OpticalSurface = new G4OpticalSurface("PMT_cathode", unified);
    PMT_OpticalSurface->SetType(dielectric_metal);
    PMT_OpticalSurface->SetModel(unified);
    PMT_OpticalSurface->SetFinish(ground);
    PMT_OpticalSurface->SetSigmaAlpha(0.);
    G4MaterialPropertiesTable* PMT_cathodeMaterialProperty = new G4MaterialPropertiesTable();
    G4double cathoderefl[2] = {0, 0};
    G4double cathodeeff[2] = {1, 1};
    PMT_cathodeMaterialProperty->AddProperty("EFFICIENCY", ener, cathodeeff, 2);//dummy
    PMT_cathodeMaterialProperty->AddProperty("REFLECTIVITY", ener, cathoderefl, 2);
    PMT_OpticalSurface->SetMaterialPropertiesTable(PMT_cathodeMaterialProperty);
    //-----------------------------------------------------------------------------
    // SiPM_OpticalSurface
    SiPM_OpticalSurface = new G4OpticalSurface("SiPM_OpticalSurface", unified);
    SiPM_OpticalSurface->SetType(dielectric_metal);
    SiPM_OpticalSurface->SetModel(unified);
    SiPM_OpticalSurface->SetFinish(ground);
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

void VDetectorConstruction::defineMaterials()
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

    const G4int numentries_metals = 2;
    G4double energies_metals[numentries_metals] = { 0.1*eV, 10.0*eV };
    G4double absorpti_metals[numentries_metals] = { 0.01 * um, 0.01 * um };
    G4MaterialPropertiesTable* prop_metals = new G4MaterialPropertiesTable();
    prop_metals->AddProperty("ABSLENGTH", energies_metals, absorpti_metals, numentries_metals);
    matAl->SetMaterialPropertiesTable(prop_metals);
    matFe->SetMaterialPropertiesTable(prop_metals);

    // Air
    matGas = new G4Material("gasAr", 0.005777*g / cm3, 1, G4State::kStateGas, 87.3 * kelvin, 1.0 * atmosphere); //TODO: use density from settings.
    matGas->AddElement(Ar, 1);
    const G4int numentries = 2;
    G4double energies[numentries] = { 0.1*eV, 10.0*eV };
    G4double vacrindices[numentries] = { 1., 1. };
    //G4double airabsorpti[numentries] = { 10*m, 10*m }; // avoid infinite light-paths
    G4double airabsorpti[numentries] = { 10000 * m, 10000 * m };
    G4MaterialPropertiesTable* airprop = new G4MaterialPropertiesTable();
    airprop->AddProperty("ABSLENGTH", energies, airabsorpti, numentries);
    airprop->AddProperty("RINDEX", energies, vacrindices, numentries);
    matGas->SetMaterialPropertiesTable(airprop);
    //------------------------------
    // LAr
    matLAr = new G4Material("LAr", 1.397*g / cm3, 1, G4State::kStateLiquid, 87.3 * kelvin, 1.0 * atmosphere); //TODO: use density from settings.
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
    matFR4 = new G4Material("FR4", 1.850*g / cm3, 3);
    matFR4->AddElement(C, 5); matFR4->AddElement(O, 2); matFR4->AddElement(H, 8);
    const G4int numentries_FR4 = 2;
    G4double energies_FR4[numentries_FR4] = { 0.1*eV, 10.0*eV };
    G4double rindices_FR4[numentries_FR4] = { 1.5, 1.5 };
    G4double absorpti_FR4[numentries_FR4] = { 1 * um, 1 * um };
    G4MaterialPropertiesTable* prop_FR4 = new G4MaterialPropertiesTable();
    prop_FR4->AddProperty("ABSLENGTH", energies_FR4, absorpti_FR4, numentries_FR4);
    prop_FR4->AddProperty("RINDEX", energies_FR4, rindices_FR4, numentries_FR4);
    matFR4->SetMaterialPropertiesTable(prop_FR4);
    //------------------------------
    // PMMA
    matPMMA = new G4Material("PMMA", 1.18*g / cm3, 3);
    matPMMA->AddElement(C, 5);	matPMMA->AddElement(O, 2); matPMMA->AddElement(H, 8);
    DataVector PMMA_ABSLENGTH(gPars::det_opt.pmma_absorption_length_filename);
    DataVector PMMA_RINDEX(gPars::det_opt.pmma_rindex_filename);
    G4MaterialPropertiesTable* prop_PMMA = new G4MaterialPropertiesTable();
    prop_PMMA->AddProperty("ABSLENGTH", PMMA_ABSLENGTH.get_Xs(eV), PMMA_ABSLENGTH.get_Ys(mm));
    prop_PMMA->AddProperty("RINDEX", PMMA_RINDEX.get_Xs(eV), PMMA_RINDEX.get_Ys());
    matPMMA->SetMaterialPropertiesTable(prop_PMMA);
    //prop_PMMA->DumpTable();
    //------------------------------
    // TPB
    matTPB = new G4Material("TPB", 1.18*g / cm3, 3);
    matTPB->AddElement(C, 5); matTPB->AddElement(O, 2); matTPB->AddElement(H, 8);
    DataVector TPB_rindex_data(gPars::det_opt.TPB_rindex_filename);
    TPB_rindex_data.scaleX(nm); TPB_rindex_data = convert_lambda_to_E(TPB_rindex_data);
    DataVector TPB_absorption_data(gPars::det_opt.TPB_abs_length_filename);
    TPB_absorption_data.scaleX(nm); TPB_absorption_data.scaleY(mm);
    TPB_absorption_data = convert_lambda_to_E(TPB_absorption_data);
    DataVector TPB_emission_data(gPars::det_opt.TPB_emission_spectrum_filename);
    TPB_emission_data.scaleX(nm); TPB_emission_data = convert_lambda_to_E(TPB_emission_data);
    DataVector TPB_efficiency_data(gPars::det_opt.TPB_efficiency_filename);
    TPB_efficiency_data.scaleX(nm); TPB_efficiency_data = convert_lambda_to_E(TPB_efficiency_data);

    G4MaterialPropertiesTable* tableTPB = new G4MaterialPropertiesTable();
    tableTPB->AddProperty("RINDEX", TPB_rindex_data.get_Xs(), TPB_rindex_data.get_Ys());
    tableTPB->AddProperty("WLSABSLENGTH", TPB_absorption_data.get_Xs(), TPB_absorption_data.get_Ys());
    tableTPB->AddProperty("WLSCOMPONENT", TPB_emission_data.get_Xs(), TPB_emission_data.get_Ys());
    //->AddProperty("WLSMEANNUMBERPHOTONS", TPB_efficiency_data.get_Xs(), TPB_efficiency_data.get_Ys());
    tableTPB->AddConstProperty("WLSMEANNUMBERPHOTONS", 0.5); // TODO: implement WLS process with dependence of WLS efficiency on incident photon's energy
    tableTPB->AddConstProperty("WLSTIMECONSTANT", 0.1*ns);
    matTPB->SetMaterialPropertiesTable(tableTPB);
    //------------------------------
    // PMMA_UV
    matPMMA_UV = new G4Material("PMMA_UV", 1.18*g / cm3, 3);
    matPMMA_UV->AddElement(C, 5);  matPMMA_UV->AddElement(O, 2);  matPMMA_UV->AddElement(H, 8);
    DataVector PMMA_UV_ABSLENGTH(gPars::det_opt.pmma_uv_absorption_length_filename);
    const G4int numentries_PMMA_UV = 2;
    G4double energies_PMMA_UV[numentries_PMMA_UV] = { 0.1*eV, 10.0*eV };
    G4double rindices_PMMA_UV[numentries_PMMA_UV] = { 1.5, 1.5 };
    G4MaterialPropertiesTable* prop_PMMA_UV = new G4MaterialPropertiesTable();
    prop_PMMA_UV->AddProperty("ABSLENGTH", PMMA_UV_ABSLENGTH.get_Xs(eV), PMMA_UV_ABSLENGTH.get_Ys(mm));
    prop_PMMA_UV->AddProperty("RINDEX", energies_PMMA_UV, rindices_PMMA_UV, numentries_PMMA_UV);
    matPMMA_UV->SetMaterialPropertiesTable(prop_PMMA_UV);
    //------------------------------
    // FSQ_RG715 - infrared filter (λ > 715 nm) with size of about 50x50x3 mm
    matFSQ_RG715 = new G4Material("FSQ_RG715", 2.20*g / cm3, 2);
    matFSQ_RG715->AddElement(Si, 1);  matFSQ_RG715->AddElement(O, 2);
    DataVector FSQ_RG715_absorption_data(gPars::det_opt.FSQ_RG715_absorption_length_filename);
    DataVector FSQ_RG715_rindex_data(gPars::det_opt.FSQ_RG715_rindex_filename);
    G4MaterialPropertiesTable* prop_FSQ_RG715 = new G4MaterialPropertiesTable();
    prop_FSQ_RG715->AddProperty("ABSLENGTH", FSQ_RG715_absorption_data.get_Xs(eV), FSQ_RG715_absorption_data.get_Ys(mm));
    prop_FSQ_RG715->AddProperty("RINDEX", FSQ_RG715_rindex_data.get_Xs(eV), FSQ_RG715_rindex_data.get_Ys());
    matFSQ_RG715->SetMaterialPropertiesTable(prop_FSQ_RG715);
}
