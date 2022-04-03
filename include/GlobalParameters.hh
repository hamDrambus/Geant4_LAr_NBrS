/*  Global parameters such as filenames, detector dimensions, etc., shared between all classes.
 *  Also stores results of simulation. All parameters except results must be read only during
 *  simulation after their initialization.
 *  TODO: Make results thread-safe.
 *  TODO?: move parameters which are not initialized in gPars::InitGlobals() elsewhere
 */
#ifndef GlobalParameters_h
#define GlobalParameters_h
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <deque>

#include <G4SystemOfUnits.hh>
#include "G4ThreeVector.hh"

#include "PolynomialFit.hh"
#include "HexagonalMapping.hh"
#include "PhotonHit.hh"

//#define TEMP_CODE_

struct DriftElectron {
  int index;
  G4ThreeVector position;
  std::string seed_info;
};

struct GeneratedData {
  DriftElectron electron;
  std::deque<PhotonHit> photons;
};

void AddToFile (std::deque<GeneratedData> data, std::string filename);

namespace gPars
{
  struct ProgramSetups {
    bool check_geometry_overlap;
    std::string track_mapping_info_class;
    int teleportation_verbosity; // 0 is quiet
    bool no_reflections;
    bool no_diffused_reflections;
    double photon_max_time; // to forcefully kill photons stuck in full internal reflections
    G4ThreeVector THGEM1_hole_center; // For convenience and testing purposes
    G4ThreeVector EL_gap_center; // For convenience and testing purposes
  };

	struct Source {
		double energy_line;
		PDF_routine energy_spectrum;
		std::string energy_spectrum_filename;
		double x_center;
		double y_center;
		double z_center;
		double xy_radius;
		double z_width;
		int N_events;
	};

	struct DetectorDimensions {
		double external_collimator_diameter;

		double THGEM1_hole_pitch;
		double THGEM1_hole_radius;
		double THGEM1_hole_rim;
		double THGEM1_dielectric_thickness;
		double THGEM1_copper_thickness;
		double THGEM1_width_total;
		double THGEM1_container_width; //area which triggers mapping to cell when hit. +-0.01 mm from both sides of THGEM1
		double THGEM1_active_area_size;
		double width_interface_grid_support; //effective width which lengthens whole inner structure
		double width_interface_grid_frame; //real frame is much thicker, but has holes for support pillars
		G4ThreeVector THGEM1_single_cell_position;
		double EL_gap_thickness;
		double z_top_interface_grid;
		double z_bottom_THGEM1;
		double PMMA_width;
		double LAr_dead_width;
		double THGEM_cathode_width;
		double Al_window_width;

		unsigned int n_PMTs; // fixed to 4
		unsigned int n_SiPMs_rows; //total number = n_SiPMs_rows^2
		std::string SiPM_device_name;
		std::string PMT_device_name;
		std::string THGEM1_cell_name;
		std::string THGEM1_cell_container_name;
	};

	struct DetectorOptics {
		double FR4_SigmaAlpha;
		double StainlessSteel_SigmaAlpha;
		double Cu_SigmaAlpha;
		double Wire_SigmaAlpha;
		std::string pmma_absorption_length_filename;
		std::string pmma_rindex_filename;
		std::string pmma_uv_absorption_length_filename;
	};

	class Results {
	public:
		std::deque<unsigned int> SiPM_photon_n; //per each device
		std::deque<unsigned int> PMT_photon_n;
		std::deque<G4ThreeVector> SiPM_positions; //Set automatically in RunAction
		std::deque<G4ThreeVector> PMT_positions;
		std::deque<GeneratedData> generated_photons; //TODO:per thread. + Merge
		std::deque<GeneratedData> recorded_photons;
		unsigned int n_reflections; // TODO: can move to some user action and add n_reflection to hit info
		unsigned int GetNGeneratedPhotons(void) const;
		unsigned int GetNRecordedPhotons(void) const;
		std::string generated_filename;
		std::string recorded_filename;
	};

	void InitGlobals(void);

	extern bool doView;

	extern std::string this_path;
	extern std::string data_path;

	extern ProgramSetups debugging;
	extern Source source;
	extern DetectorDimensions det_dims;
	extern DetectorOptics det_opt;
	extern Results results;

	extern HexagonalMapping* THGEM1_mapping; // Set after geometry construction in RunAction. Thread safe.
}

#endif //GlobalParameters_h
