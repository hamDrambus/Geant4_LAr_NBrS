#include "GlobalParameters.hh"

namespace gPars
{
	ProgramSetups general;
	Source source;
	DetectorDimensions det_dims;
	ElmerFieldMap field_map;
	DetectorOptics det_opt;
	ArgonProperties Ar_props;
	Results results;

	void InitGlobals(void)
	{
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

		general.doView = false;
		general.doViewElectronDrift = false;
		general.thread_number = 8;
		general.data_path = "../../NBrS_THGEM_LAr_v0/data/";
		general.output_folder = "../../NBrS_THGEM_LAr_v0/tests/test_11_fieldmap_speed/v01_tetra/";
		general.check_geometry_overlap = false;
		general.no_reflections = false;
		general.no_diffused_reflections = false;
		general.enable_e_diffusion = true;
		general.teleportation_verbosity = 0;
		general.photon_max_time = 3.0 * ns; // about 1 meter full path.
		general.electron_max_time = DBL_MAX;//3 * 3.5e-4 * ns; // 3 times normal drift time
		general.surface_tolerance = 1e-8 * mm;
		general.record_electrons = false; // if false, only photons are recorded to files and kept in memory
		general.record_detailed = false; // if false, only number of hits per channel is kept track of
		general.print_drift_track = false; // works only for Detector_THGEM1_detailed and GenElectronsPatterns testing classes

		// THGEM CERN 28%, in [mm]
		det_dims.THGEM1_copper_thickness = 0.03;
		det_dims.THGEM1_hole_radius = 0.25;
		det_dims.THGEM1_dielectric_thickness = 0.4;
		det_dims.THGEM1_hole_pitch = 0.9;
		det_dims.THGEM1_hole_rim = 0.1;
		//THGEM Electroconnect 75%, in [mm]
    //det_dims.THGEM1_copper_thickness = 0.03;
    //det_dims.THGEM1_hole_radius = 0.5;
    //det_dims.THGEM1_dielectric_thickness = 0.96;
    //det_dims.THGEM1_hole_pitch = 1.1;
    //det_dims.THGEM1_hole_rim = 0.0;

		det_dims.THGEM1_active_area_size = 100;
		det_dims.THGEM1_width_total = 2*det_dims.THGEM1_copper_thickness + det_dims.THGEM1_dielectric_thickness;
		det_dims.THGEM1_container_width = det_dims.THGEM1_width_total + 0.01 * 2; // 0.01 mm from each side. Different value from Elmer simulation

		det_dims.SiPM_size = 6.0;
		det_dims.n_SiPMs_rows = 5;
		det_dims.SiPM_device_name = "phys_SiPM";

		det_dims.width_interface_grid_support = 1.4;
		det_dims.width_interface_grid_frame = 5;
		det_dims.THGEM1_single_cell_position = G4ThreeVector(150 * mm, 150 * mm, 150 * mm);
		det_dims.external_collimator_diameter = 50; //No collimator if >= diameter_size_Al_window
		//det_dims.EL_gap_thickness = 13;//double phase
		det_dims.EL_gap_thickness = -4;//single phase. From THGEM1 real bottom.
		det_dims.z_top_interface_grid = 48.0 + det_dims.width_interface_grid_support;
		det_dims.z_bottom_THGEM1 = det_dims.z_top_interface_grid + 22; //=71.4
		det_dims.PMMA_width = 3;
		det_dims.LAr_dead_width = 2;
		det_dims.THGEM_cathode_width = 0.5;
		det_dims.Al_window_width = 23;
		det_dims.n_PMTs = 4;
		det_dims.PMT_device_name = "phys_PMT";
		det_dims.THGEM1_cell_name = "phys_THGEM1_cell";
		det_dims.THGEM1_cell_container_name = "phys_THGEM1_cell_container";

		source.energy_line = -1;
    source.energy_spectrum_filename = "energy_spectrum/YAP_Ce_energies_eV_1.dat";
    source.energy_spectrum.read(general.data_path+source.energy_spectrum_filename);
    source.x_center = 0.00;
    source.y_center = 0.00;
    source.z_center = det_dims.z_bottom_THGEM1 - 1.5;
    source.xy_radius = 7;
    source.z_width = 0;
    source.N_events = 500000;
    source.NBrS_yield_factor = 1;

		field_map.elmer_mesh_folder = general.data_path + "../singleTHGEM28_LAr/v00.02_THGEM1/";
		field_map.elmer_solution_filename = general.data_path + "../singleTHGEM28_LAr/Elmer_v00.02/case_6180v.result";
		field_map.elmer_mesh_center = G4ThreeVector(0,0,0);
		field_map.LAr_drift_velocity = general.data_path + "LAr_drift/LAr_drift_velocity_data.txt";
		field_map.LAr_diffusion_transversal = general.data_path + "LAr_drift/LAr_diffusionT_data.txt";
		field_map.LAr_diffusion_longitudinal = general.data_path + "LAr_drift/LAr_diffusionL_data.txt";
		field_map.mesh_tolerance = 1e-10;
		field_map.drift_step_size = 2 * um;
		field_map.max_rel_field_change = 0.05; // 5%

		det_opt.FR4_SigmaAlpha = 50;
		det_opt.StainlessSteel_SigmaAlpha = 10;
		det_opt.Cu_SigmaAlpha = 30;
		det_opt.Wire_SigmaAlpha = 5;
		det_opt.FR4_reflectivity = 0.20; //model_25b
		//det_opt.FR4_reflectivity = 0.05; //https://www.cetem.gov.br/images/congressos/2008/CAC00560008.pdf
		det_opt.Cu_reflectivity = 0.36; // Bass M. Handbook of optics, Vol.4 Edition3
		det_opt.pmma_absorption_length_filename = general.data_path+"absorption_length/PMMA_absorption_length_eV_mm.dat";
		det_opt.pmma_rindex_filename = general.data_path+"refractive_index/PMMA_rindex_eV_1.dat";
		det_opt.pmma_uv_absorption_length_filename = general.data_path+"absorption_length/PMMA_VU_absorption_length_eV_mm.dat";

		Ar_props.print_calculations = true;
		Ar_props.pedantic_calculations = false;
		Ar_props.XS_energy_transfer_filename = general.data_path + "LAr_XS/XS_energy_transfer.txt";
		Ar_props.XS_momentum_transfer_filename = general.data_path + "LAr_XS/XS_momentum_transfer.txt";
		Ar_props.cache_folder = "../../NBrS_THGEM_LAr_v0/data_cache/";
		Ar_props.energy_max = 6 * eV; //XSs are known up to 11.5 eV, but with field under 100kV/cm electron energy is < 4 eV.
		Ar_props.field_max = 100 * kilovolt / cm; // Maximal electric field for which argon properties should be calculated.
		Ar_props.atomic_density = 2.10e22 / (cm * cm * cm);
		Ar_props.m_to_M = 5.109989461e5 / 3.726e10; // m_electron / M_Ar_atom ratio (eV to eV)

		results.generated_filename = "generated.dat";
		results.recorded_filename = "recorded.dat";

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
      return;
    }

		if ((-det_dims.EL_gap_thickness < (det_dims.THGEM1_width_total + (det_dims.THGEM1_container_width - det_dims.THGEM1_width_total)*0.5)) &&
        (-det_dims.EL_gap_thickness > (-(det_dims.THGEM1_container_width - det_dims.THGEM1_width_total)*0.5))) {
      G4Exception("gPars::InitGlobals(): ",
            "InvalidSetup", FatalException, "LAr level intersects with THGEM1 (or its container)! This case is not supported.");
      return;
    }
    if (det_dims.n_PMTs != 4) {
      G4Exception("gPars::InitGlobals(): ",
            "InvalidSetup", FatalException, "Number of PMTs must be 4.");
      return;
    }

    if (general.no_diffused_reflections) {
      det_opt.FR4_SigmaAlpha = 0.0;
      det_opt.StainlessSteel_SigmaAlpha = 0.0;
      det_opt.Cu_SigmaAlpha = 0.0;
      det_opt.Wire_SigmaAlpha = 0.0;
    }
    if (source.N_events > 200 && (general.doView || general.doViewElectronDrift)) {
      G4Exception("gPars::InitGlobals(): ",
            "InvalidSetup", FatalException, "Viewing large number of events is disabled.");
      return;
    }
	}
}
