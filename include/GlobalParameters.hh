/*  Global parameters such as filenames, detector dimensions, etc., shared between all classes.
 *  The difference from GlobalData is that this class only stores simple parameters which all
 *  may be easily changed without silent breaking of simulation.
 *  All parameters must be read only during simulation after their initialization.
 */
#ifndef GlobalParameters_h
#define GlobalParameters_h
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <deque>

#include <G4SystemOfUnits.hh>
#include <G4PhysicalConstants.hh>
#include <G4ThreeVector.hh>

#include "PolynomialFit.hh"
#include "HexagonalMapping.hh"
#include "PhotonHit.hh"

//#define TEMP_CODE_

static constexpr double Td = 1.e-17 * volt * cm * cm;
static constexpr double e_mass_SI = 9.10938370e-31; // electron mass in kg
static constexpr double hc = hbarc * twopi; // for recalculating photon energy to wavelength and back. E*lambda = h*c = const.

namespace gPars
{
  struct ProgramSetups {
    std::string this_path;
    std::string data_path;
    std::string output_folder;
    std::string gnuplot_bin;
    int thread_number;
    bool doView;
    bool doViewElectronDrift;
    bool check_geometry_overlap;
    std::string track_mapping_info_class;
    int teleportation_verbosity; // 0 is quiet
    bool no_reflections;
    bool no_diffused_reflections;
    double photon_max_time; // to forcefully kill photons stuck in full internal reflections
    double electron_max_time; // to forcefully kill stuck drifting electrons. For debugging only, must be DBL_MAX in real simulation!
    bool enable_e_diffusion;
    G4ThreeVector THGEM1_hole_center; // For convenience and testing purposes
    G4ThreeVector EL_gap_center; // For convenience and testing purposes
    double surface_tolerance;
    bool record_electrons; // if false, only photons are recorded to files and kept in memory
    bool record_detailed; // if false, only number of hits per channel is kept track of
    bool print_drift_track; // works only for Detector_THGEM1_detailed and GenElectronsPatterns testing classes
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
		double NBrS_yield_factor; // To avoid drifting electrons w/o photon emission
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

		double SiPM_size;
		unsigned int n_SiPMs_rows; //total number = n_SiPMs_rows^2
		std::string SiPM_device_name;

		unsigned int n_PMTs; // fixed to 4
		std::string PMT_device_name;
		std::string THGEM1_cell_name;
		std::string THGEM1_cell_container_name;
	};

	struct ElmerFieldMap {
	  std::string elmer_mesh_folder;
    std::string elmer_solution_filename;
    G4ThreeVector elmer_mesh_center;
    std::string LAr_drift_velocity; // Used only when there is no theoretically calculated one
    std::string LAr_diffusion_longitudinal; // Used only when there is no theoretically calculated one
    std::string LAr_diffusion_transversal; // Used only when there is no theoretically calculated one
    double drift_step_size;
    double mesh_tolerance; // relative, necessary to get field near mesh boundary
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

	struct ArgonProperties {
	  bool print_calculations;
	  bool pedantic_calculations; // if true then incorrect computations (such as division by 0,
	  //infinities, DBL_MAX) throw error. Otherwise default values are used.
	  std::string XS_energy_transfer_filename;
	  std::string XS_momentum_transfer_filename;
	  std::string cache_folder;
	  double energy_max; // In Geant4 units. Maximal energy until which LAr properties are known and to be calculated to.
	  double field_max; // In Geant4 units. Maximal electric field for which argon properties should be calculated.
	  double atomic_density; // in Geant4 units [mm^-3]
	  double m_to_M; // m_electron / M_Ar_atom ratio
	};

	class Results {
	public:
		std::string generated_filename;
		std::string recorded_filename;
	};

	void InitGlobals(void);

	extern ProgramSetups general;
	extern Source source;
	extern DetectorDimensions det_dims;
	extern ElmerFieldMap field_map;
	extern DetectorOptics det_opt;
	extern ArgonProperties Ar_props;
	extern Results results;
}

#endif //GlobalParameters_h
