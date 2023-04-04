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

#include <boost/foreach.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/algorithm/string.hpp>

#include <G4SystemOfUnits.hh>
#include <G4PhysicalConstants.hh>
#include <G4ThreeVector.hh>

#include "utilities/PolynomialFit.hh"
#include "geant4/detector/DetectorSettings.hh"
#include "geant4/generator/SourceSettings.hh"

//#define TEMP_CODE_

static constexpr double Td = 1.e-17 * volt * cm * cm;
static constexpr double e_mass_SI = 9.10938370e-31; // electron mass in kg
static constexpr double hc = hbarc * twopi; // for recalculating photon energy to wavelength and back. E*lambda = h*c = const.

namespace gPars
{
  struct ProgramSetups {
    std::string this_path;
    std::string data_path;
    std::string settings_filename;
    std::string output_folder;
    std::string gnuplot_bin;
    int thread_number;
    bool doView;
    bool doViewElectronDrift;
    bool check_geometry_overlap;
    int teleportation_verbosity; // 0 is quiet
    bool no_reflections;
    bool no_diffused_reflections;
    double photon_max_time; // to forcefully kill photons stuck in full internal reflections
    double electron_max_time; // to forcefully kill stuck drifting electrons. For debugging only, must be DBL_MAX in real simulation!
    bool enable_e_diffusion;
    bool record_electrons; // if false, only photons are recorded to files and kept in memory
    bool record_detailed; // if false, only number of hits per channel is kept track of
    bool print_drift_track; // works only for Detector_THGEM1_detailed and GenElectronsPatterns testing classes
    G4long initial_seed; // if not specified in setups random time(NULL) is used as starting seed for the simulation
  };

	struct ElmerFieldMap {
	  std::string elmer_mesh_folder;
    std::string elmer_solution_filename;
    G4ThreeVector elmer_mesh_center;
    double drift_step_size;
    double mesh_tolerance; // relative, necessary to get field near mesh boundary
    double max_rel_field_change;
	};

	struct DetectorOptics {
		double FR4_SigmaAlpha;
		double LAr_SigmaAlpha;
		double PMMA_SigmaAlpha;
		double StainlessSteel_SigmaAlpha;
		double Cu_SigmaAlpha;
		double Wire_SigmaAlpha;
		double FR4_reflectivity;
		double Cu_reflectivity;
		double StainlessSteel_reflectivity;
		double Wire_reflectivity;
		std::string pmma_absorption_length_filename;
		std::string pmma_rindex_filename;
		std::string pmma_uv_absorption_length_filename;
		std::string TPB_rindex_filename;
		std::string TPB_abs_length_filename;
		std::string TPB_efficiency_filename;
		std::string TPB_emission_spectrum_filename;
	};

	enum NBrSFormula {
	  ElasticXS,
	  TransferXS
	};

	enum Medium {
		LiquidAr,
		LiquidKr,
		LiquidXe
	};

	struct MediumProperties {
		Medium medium_type;
	  bool print_calculations;
	  bool pedantic_calculations; // if true then incorrect computations (such as division by 0,
	  //infinities, DBL_MAX) throw error. Otherwise default values are used.
	  std::string XS_energy_transfer_filename;
	  std::string XS_momentum_transfer_filename;
	  std::string cache_folder;
	  std::string exp_drift_velocity; // Used only when there is no theoretically calculated one
    std::string exp_diffusion_longitudinal; // Used only when there is no theoretically calculated one
    std::string exp_diffusion_transversal; // Used only when there is no theoretically calculated one
	  double atomic_density; // in Geant4 units [mm^-3]
    double m_to_M; // m_electron / M_Ar_atom ratio (eV to eV)
    NBrSFormula NBrS_formula;
    bool force_recalculation;
    double distributions_energy_step_factor;
    double field_step_factor;
    double spectra_step_factor;
	};

	class Results {
	public:
		std::string generated_filename;
		std::string recorded_filename;
	};

	bool InitGlobals(std::string filename);
	bool LoadSettings(std::string filename);

	extern ProgramSetups general;
	extern VDetectorDimensions *det_dims;
	extern VSourceSettings *source;
	extern ElmerFieldMap field_map;
	extern DetectorOptics det_opt;
	extern MediumProperties medium_props;
	extern Results results;

	std::string MediumName(void);
}

#endif //GlobalParameters_h
