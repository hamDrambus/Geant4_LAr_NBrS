#include <GlobalParameters.hh>
#include <geant4/detector/DetectorSettings.hh>
#include <geant4/generator/SourceSettings.hh>

namespace gPars
{
	ProgramSetups general;
	VDetectorDimensions *det_dims;
	VSourceSettings *source;
	ElmerFieldMap field_map;
	DetectorOptics det_opt;
	ArgonProperties Ar_props;
	Results results;

	bool LoadSettings(std::string fname)
	{
	  std::cout << "Loading settings \"" << fname << "\"..." << std::endl;
    // Create an empty property tree object
    using boost::property_tree::ptree;
    using boost::property_tree::ptree_bad_data;
    using boost::property_tree::ptree_bad_path;
    using boost::property_tree::ptree_error;
    ptree pt;
    try {
      read_xml(fname, pt);
      {
      ptree gen = pt.get_child("Settings.General");
      general.data_path = gen.get<std::string>("data_path");
      general.output_folder = gen.get<std::string>("output_folder");
      general.doView = gen.get<bool>("doView", false);
      general.doViewElectronDrift = gen.get<bool>("doViewElectronDrift", false);
      general.record_electrons = gen.get<bool>("record_electrons", true);
      general.record_detailed = gen.get<bool>("record_detailed", true);
      results.generated_filename = gen.get<std::string>("generated_filename", "generated.dat");
      results.recorded_filename = gen.get<std::string>("recorded_filename", "recorded.dat");
      general.enable_e_diffusion = gen.get<bool>("enable_e_diffusion", true);
      general.thread_number = gen.get<int>("thread_number", 8);
      general.photon_max_time = gen.get<double>("photon_max_time_us", 3.0) * us;
      general.initial_seed = gen.get<long>("initial_seed", time(NULL));
      }
      {
      ptree deb = pt.get_child("Settings.Debug");
      general.check_geometry_overlap = deb.get<bool>("check_geometry_overlap", true);
      general.electron_max_time = deb.get<double>("electron_max_time_ns", DBL_MAX);
      if (general.electron_max_time != DBL_MAX)
        general.electron_max_time *= ns;
      general.print_drift_track = deb.get<bool>("print_drift_track", false);
      general.teleportation_verbosity = deb.get<int>("teleportation_verbosity", 0);
      Ar_props.print_calculations = deb.get<bool>("Ar_print_calculations", true);
      Ar_props.pedantic_calculations = deb.get<bool>("Ar_pedantic_calculations", true);
      }
      {
      ptree field = pt.get_child("Settings.FieldMap");
      field_map.elmer_mesh_folder = field.get<std::string>("elmer_mesh_folder");
      field_map.elmer_solution_filename = field.get<std::string>("elmer_solution_filename");
      field_map.drift_step_size = std::fabs(field.get<double>("drift_step_size_um") * um);

      field_map.mesh_tolerance = field.get<double>("mesh_tolerance", 1e-10);
      field_map.max_rel_field_change = std::fabs(field.get<double>("max_rel_field_change", 0.05));
      field_map.elmer_mesh_center.setX(field.get<double>("mesh_center_x_mm", 0) * mm);
      field_map.elmer_mesh_center.setY(field.get<double>("mesh_center_y_mm", 0) * mm);
      field_map.elmer_mesh_center.setZ(field.get<double>("mesh_center_z_mm", 0) * mm);
      }
      {
      ptree optics = pt.get_child("Settings.DetectorOptics");
      det_opt.FR4_SigmaAlpha = optics.get<double>("FR4_SigmaAlpha_deg") * deg;
      det_opt.LAr_SigmaAlpha = optics.get<double>("LAr_SigmaAlpha_deg", 0) * deg;
      det_opt.PMMA_SigmaAlpha = optics.get<double>("PMMA_SigmaAlpha_deg", 0) * deg;
      det_opt.StainlessSteel_SigmaAlpha = optics.get<double>("StainlessSteel_SigmaAlpha_deg") * deg;
      det_opt.Cu_SigmaAlpha = optics.get<double>("Cu_SigmaAlpha_deg") * deg;
      det_opt.Wire_SigmaAlpha = optics.get<double>("Wire_SigmaAlpha_deg") * deg;
      det_opt.FR4_reflectivity = optics.get<double>("FR4_reflectivity");
      det_opt.Cu_reflectivity = optics.get<double>("Cu_reflectivity");
      det_opt.StainlessSteel_reflectivity = optics.get<double>("StainlessSteel_reflectivity", 0.5); // doi:10.1063/1.2202915
      det_opt.Wire_reflectivity =  optics.get<double>("Wire_reflectivity", 0.5); //approximately https://nvlpubs.nist.gov/nistpubs/bulletin/07/nbsbulletinv7n2p197_A2b.pdf
      det_opt.pmma_absorption_length_filename = optics.get<std::string>("pmma_absorption_length_filename");
      det_opt.pmma_rindex_filename = optics.get<std::string>("pmma_rindex_filename");
      det_opt.pmma_uv_absorption_length_filename = optics.get<std::string>("pmma_uv_absorption_length_filename");
      det_opt.TPB_rindex_filename = optics.get<std::string>("TPB_rindex_filename");
      det_opt.TPB_abs_length_filename = optics.get<std::string>("TPB_abs_length_filename");
      det_opt.TPB_efficiency_filename = optics.get<std::string>("TPB_efficiency_filename");
      det_opt.TPB_emission_spectrum_filename = optics.get<std::string>("TPB_emission_spectrum_filename");

      general.no_reflections = optics.get<bool>("no_reflections", false);
      general.no_diffused_reflections = optics.get<bool>("no_diffused_reflections", false);
      }
      {
      ptree ar_props = pt.get_child("Settings.ArgonProperties");
      Ar_props.XS_energy_transfer_filename = ar_props.get<std::string>("XS_energy_transfer_filename");
      Ar_props.XS_momentum_transfer_filename = ar_props.get<std::string>("XS_momentum_transfer_filename");

      // Not changeable parameters:
      Ar_props.atomic_density = ar_props.get<double>("atomic_density_cm3", 2.10e22) / (cm3);
      Ar_props.cache_folder = ar_props.get<std::string>("cache_folder", "");
      Ar_props.LAr_drift_velocity = ar_props.get<std::string>("LAr_drift_velocity", "");
      Ar_props.LAr_diffusion_transversal = ar_props.get<std::string>("LAr_diffusion_transversal", "");
      Ar_props.LAr_diffusion_longitudinal = ar_props.get<std::string>("LAr_diffusion_longitudinal", "");
      std::string NBrS_formula_str = ar_props.get<std::string>("NBrS_formula", "ElasticXS");
      boost::optional<NBrSFormula> NBrS_formula;
      if (NBrS_formula_str == "ElasticXS")
        NBrS_formula = NBrSFormula::ElasticXS;
      if (NBrS_formula_str == "TransferXS")
        NBrS_formula = NBrSFormula::TransferXS;
      if (boost::none == NBrS_formula) {
        std::cerr<<"LoadSettings:Warning: NBrS_formula=\""<<NBrS_formula_str<<"\" is not supported! Default \"ElasticXS\" is used."<<std::endl;
        Ar_props.NBrS_formula = NBrSFormula::ElasticXS;
      } else
        Ar_props.NBrS_formula = *NBrS_formula;
      Ar_props.force_recalculation = ar_props.get<bool>("force_recalculation", false);
      Ar_props.distributions_energy_step_factor = ar_props.get<double>("distributions_energy_step_factor", 1.0);
      Ar_props.field_step_factor = ar_props.get<double>("field_step_factor", 1.0);
      Ar_props.spectra_step_factor = ar_props.get<double>("spectra_step_factor", 1.0);
      }

      det_dims = CreateDetectorSettings(fname);
      if (nullptr == det_dims)
        return false;
      source = CreateSourceSettings(fname);
      if (nullptr == source)
        return false;

      field_map.elmer_mesh_folder = general.data_path + field_map.elmer_mesh_folder;
      field_map.elmer_solution_filename = general.data_path + field_map.elmer_solution_filename;

      det_opt.pmma_absorption_length_filename = general.data_path+det_opt.pmma_absorption_length_filename;
      det_opt.pmma_rindex_filename = general.data_path+det_opt.pmma_rindex_filename;
      det_opt.pmma_uv_absorption_length_filename = general.data_path+det_opt.pmma_uv_absorption_length_filename;
      det_opt.TPB_abs_length_filename = general.data_path+det_opt.TPB_abs_length_filename;
      det_opt.TPB_efficiency_filename = general.data_path+det_opt.TPB_efficiency_filename;
      det_opt.TPB_emission_spectrum_filename = general.data_path+det_opt.TPB_emission_spectrum_filename;
      det_opt.TPB_rindex_filename = general.data_path+det_opt.TPB_rindex_filename;

      Ar_props.XS_energy_transfer_filename = general.data_path + Ar_props.XS_energy_transfer_filename;
      Ar_props.XS_momentum_transfer_filename = general.data_path + Ar_props.XS_momentum_transfer_filename;
      Ar_props.LAr_drift_velocity = general.data_path + Ar_props.LAr_drift_velocity;
      Ar_props.LAr_diffusion_transversal = general.data_path + Ar_props.LAr_diffusion_transversal;
      Ar_props.LAr_diffusion_longitudinal = general.data_path + Ar_props.LAr_diffusion_longitudinal;

    } catch (ptree_bad_path& e) {
      std::cout << "LoadSettings: ptree_bad_path exception:" << std::endl;
      std::cout << e.what() << std::endl;
      goto fail_load;
    } catch (ptree_bad_data& e) {
      std::cout << "LoadSettings: ptree_bad_data exception:" << std::endl;
      std::cout << e.what() << std::endl;
      goto fail_load;
    } catch (ptree_error& e) {
      std::cout << "LoadSettings: ptree_error exception:" << std::endl;
      std::cout << e.what() << std::endl;
      goto fail_load;
    } catch (std::exception& e) {
      std::cout << "LoadSettings: std::exception:" << std::endl;
      std::cout << e.what() << std::endl;
      goto fail_load;
    }
    std::cout << "LoadSettings: Successfully loaded \"" << fname << "\"" << std::endl;
    return true;
  fail_load:
    std::cout << "LoadSettings: Failed to load settings \"" << fname << "\"" << std::endl;
    return false;
	}

	bool InitGlobals(std::string filename)
	{
	  // Possible defaults are set in LoadSettings
	  // They may be overwritten in LoadSettings().
	  // If default is not possible and is absent in settings,
	  // initialization fails.
	  det_dims = nullptr;
	  source = nullptr;
		char path[FILENAME_MAX];
#if defined(__WIN32__)
		general.this_path = _getcwd(path, FILENAME_MAX);
#else
		general.this_path = getcwd(path, FILENAME_MAX);
#endif
		if (!general.this_path.empty())
			if (general.this_path.back()!='/')
			  general.this_path.push_back('/');

#if defined(_WIN32)||defined(_WIN64)
		general.gnuplot_bin = "\"%GNUPLOT%\\gnuplot.exe\"";
#else
		general.gnuplot_bin = "gnuplot";
#endif
		general.settings_filename = filename;
		if (!LoadSettings(filename))
		  return false;

		//=====================================================================================================================
		// Consistency checks below

    std::cout<<"This path: \""<<general.this_path<<"\""<<std::endl;
    std::cout<<"Data path: \""<<general.this_path+general.data_path<<"\""<<std::endl;

    if (gPars::general.thread_number < 1) {
      std::cerr<<"gPars::InitGlobals(): Warning, thread number < 1, setting to 1."<<std::endl;
      gPars::general.thread_number = 1;
    }

    if (!general.doView)
      general.doViewElectronDrift = false;

    if (general.doViewElectronDrift && !general.record_electrons) {
      G4Exception("gPars::InitGlobals(): ",
            "InvalidSetup", FatalException, "To view electron drift electrons must be recorded");
      return false;
    }

    if (general.no_diffused_reflections) {
      det_opt.FR4_SigmaAlpha = 0.0;
      det_opt.LAr_SigmaAlpha = 0.0;
      det_opt.PMMA_SigmaAlpha = 0.0;
      det_opt.StainlessSteel_SigmaAlpha = 0.0;
      det_opt.Cu_SigmaAlpha = 0.0;
      det_opt.Wire_SigmaAlpha = 0.0;
    }
    if (source->N_events > 200 && (general.doView || general.doViewElectronDrift)) {
      G4Exception("gPars::InitGlobals(): ",
            "InvalidSetup", FatalException, "Viewing large number of events is disabled.");
      return false;
    }
    return true;
	}
}
