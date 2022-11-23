#ifndef DetectorSettings_h
#define DetectorSettings_h

#include <string>
#include <sstream>
#include <limits>
#include <globals.hh>

#include <G4SystemOfUnits.hh>
#include <G4ThreeVector.hh>

#include <utilities/GlobalUtilities.hh>

// Class containing input and resulting detector settings. Some values
// (such as cell position, THGEM1 hole center, EL gap center) are
// initialized after concrete detector construction is finished.
// All parameters are in Geant4 units
class VDetectorDimensions
{
protected:
  VDetectorDimensions() {}
public:
  enum DetectorType {
    Full_detector,
    THGEM1_detailed,
    THGEM1_SiPM_shading,
		Full_detector_y2022
  };
  DetectorType detector_type;

  double THGEM1_hole_pitch;
	double THGEM1_hole_radius; // at the ends of dielectric hole.
	double THGEM1_hole_rim; // Copper hole' r = THGEM1_hole_radius + THGEM1_hole_rim
  double THGEM1_dielectric_radius; // at hole center. By default is the same as THGEM1_hole_radius-THGEM1_hole_rim
  double THGEM1_dielectric_thickness;
  double THGEM1_copper_thickness;
  double THGEM1_width_total;
  double THGEM1_container_width; //area which triggers mapping to cell when hit. +-0.01 mm from both sides of THGEM1

  // For convenience and/or testing purposes. Set in VDetectorConstruction::SetSizeAndPosition()
  G4ThreeVector THGEM1_single_cell_position;
  G4ThreeVector THGEM_hole_center; // may be THGEM1 or THGEM0, depending on where NBrS is expected from.
  G4ThreeVector THGEM1_center;
  G4ThreeVector EL_gap_center;
  G4ThreeVector Cathode_top_center;
  unsigned int n_SiPMs;
  unsigned int n_PMTs;
  // This one is set after field map is loaded
  G4ThreeVector drift_start_center; // Should be x = y = 0, and z is where electric fields starts (depends on loaded field map)

  static const std::string SiPM_device_name;
  static const std::string PMT_device_name;
};

VDetectorDimensions* CreateDetectorSettings(std::string filename);

class DetectorDimsFull : public VDetectorDimensions
{
public:
  DetectorDimsFull() { detector_type = Full_detector; }
};

class DetectorDimsTGHEM1Detailed : public VDetectorDimensions
{
public:
  DetectorDimsTGHEM1Detailed() {detector_type = THGEM1_detailed; }
};

class DetectorDimsTGHEM1Shading : public VDetectorDimensions
{
public:
  DetectorDimsTGHEM1Shading() { detector_type = THGEM1_SiPM_shading; }
};

class DetectorDimsFullY2022 : public VDetectorDimensions
{
public:
	DetectorDimsFullY2022() { detector_type = Full_detector_y2022; }
	double THGEM0_hole_pitch;
	double THGEM0_hole_radius; // at the ends of dielectric hole.
	double THGEM0_hole_rim; // Copper hole' r = THGEM0_hole_radius + THGEM0_hole_rim
	double THGEM0_dielectric_radius; // at hole center. By default is the same as THGEM0_hole_radius-THGEM0_hole_rim
	double THGEM0_dielectric_thickness;
	double THGEM0_copper_thickness;
	double THGEM0_width_total;
	double THGEM0_container_width; //area which triggers mapping to cell when hit. +-0.01 mm from both sides of THGEM1
	bool is_NBrS_in_THGEM0;
	bool has_WLS;

	// For convenience and/or testing purposes. Set in VDetectorConstruction::SetSizeAndPosition()
	G4ThreeVector THGEM0_single_cell_position;
	G4ThreeVector THGEM0_center;
};

#endif // DetectorSettings_h
