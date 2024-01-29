#include <geant4/detector/Detector_full_y2024_NIR.hh>

Detector_full_y2024_NIR::Detector_full_y2024_NIR() :
  Detector_full_y2022()
{}

Detector_full_y2024_NIR::~Detector_full_y2024_NIR()
{}

G4VPhysicalVolume * Detector_full_y2024_NIR::Construct()
{
    physiWorld = Detector_full_y2022::Construct();
    //--------------------------------------------------------------------------------
    // Remove PMMA plate
    logicWorld->RemoveDaughter(phys_PMMA_plate);

    // Create optical filter with its support framer.
    // Create filter
    auto solid_optical_filter = new G4Box("solid_optical_filter", filter_size_xy / 2.0, filter_size_xy / 2.0, filter_size_z / 2.0);
    auto logic_optical_filter = new G4LogicalVolume(solid_optical_filter, matFSQ_RG715, "logic_optical_filter", 0, 0, 0);
    auto phys_optical_filter = new G4PVPlacement(0, position_filter, logic_optical_filter, "phys_optical_filter",
        logicWorld, false, 0, fCheckOverlaps);
    // Create filter support frame
    auto solid_optical_filter_substrate = new G4Box("optical_filter_substrate", filter_frame_size_xy / 2.0, filter_frame_size_xy / 2.0, filter_frame_thickness / 2.0);
    auto solid_optical_filter_hole = new G4Box("optical_filter_hole", filter_size_xy / 2.0, filter_size_xy / 2.0, filter_frame_thickness / 2.0);
    auto solid_optical_filter_frame = new G4SubtractionSolid("solid_optical_filter_frame", solid_optical_filter_substrate, solid_optical_filter_hole);
    auto logic_optical_filter_frame = new G4LogicalVolume(solid_optical_filter_frame, matFR4, "l_optical_filter_frame", 0, 0, 0);
    auto phys_anode_grid_frame = new G4PVPlacement(0, position_filter_frame, logic_optical_filter_frame, "p_optical_filter_frame",
        logicWorld, false, 0);

    if (std::fabs(SiPMs_z_offset) > 0) {
        // Shift SiPM matrix
        G4ThreeVector pos(phys_SiPM_container->GetTranslation());
        phys_SiPM_container->SetTranslation(pos + G4ThreeVector(0, 0, SiPMs_z_offset));
        pos = phys_SiPMFR4->GetTranslation();
        phys_SiPMFR4->SetTranslation(pos + G4ThreeVector(0, 0, SiPMs_z_offset));
        if (fCheckOverlaps) {
            phys_SiPM_container->CheckOverlaps();
            phys_SiPMFR4->CheckOverlaps();
        }
    }
    
    //--------------------------------------------------------------------------------
    // Setting visualization
    G4VisAttributes FSQ_RG715_VisAtt(G4Colour(0.7, 0.2, 0.3, 0.1));
    G4VisAttributes FR4_VisAtt(G4Colour(1.0, 0.5, 0.2, 0.4));
    logic_optical_filter_frame->SetVisAttributes(FR4_VisAtt);
    logic_optical_filter->SetVisAttributes(FSQ_RG715_VisAtt);

    SetupTHGEMsMapping();
    return physiWorld;
}

void Detector_full_y2024_NIR::SetSizeAndPosition()
{
    // x = 0, y = 0 is on the detector central axis.
    // z = 0 (global) is the cathode's top.
    Detector_full_y2022::SetSizeAndPosition();
    // Optical filter, placed before SiPMs instead of PMMA plate
    // FSQ RG715 filter plate (frame, support)
	filter_frame_thickness = 1.0 * mm;
    // If filter_frame_thickness is not equal to PMMA_plate_size_z, moving SiPMs is required.
    SiPMs_z_offset = filter_frame_thickness - PMMA_plate_size_z;
	filter_frame_size_xy = THGEM1_size_xy;
	filter_frame_z_bottom = position_PMMA_plate.z() - PMMA_plate_size_z / 2.0;
	// FSQ RG715 filter itself
	filter_size_xy = 50.0 * mm;
	filter_size_z = 3.0 * mm;
	position_filter_frame = G4ThreeVector(0, 0, filter_frame_z_bottom + filter_frame_thickness / 2.0);
    // Filter is positioned so that its top surface is coincident with support grid's top surface.
	position_filter = G4ThreeVector(0, 0, filter_frame_z_bottom + filter_frame_thickness - filter_size_z / 2.0);

    if (position_filter.z() - filter_size_z / 2.0 < position_anode_grid_container.z() + anode_grid_cont_size_z / 2.0) {
        G4Exception("DetectorConstruction::SetSizeAndPosition: ",
            "InvalidSetup", FatalException, "Optical filter is intersecting with anode wires.");
    }
    if (filter_size_xy > anode_wire_length
        && (position_filter.z() - filter_size_z / 2.0) < position_anode_grid_frame.z() + anode_grid_thickness / 2.0) {
        G4Exception("DetectorConstruction::SetSizeAndPosition: ",
            "InvalidSetup", FatalException, "Optical filter does not fit into anode hole for wires."
            "As a result it is also intersecting with anode support frame.");
    }
}
