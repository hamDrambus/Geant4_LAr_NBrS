#include "GlobalParameters.hh"

void AddToFile (std::deque<GeneratedData> data, std::string filename)
{
  std::ofstream str;
  open_output_file(filename, str, std::ios_base::ate|std::ios_base::out);
  if (!str.is_open()) {
    std::cerr<<"AddToFile:Error: Failed to open file"<<std::endl;
    std::cerr<<"\t\""<<filename<<"\""<<std::endl;
    return;
  }
  for (std::size_t e = 0, e_end_ = data.size(); e!=e_end_; ++e) {
    str<<data[e].electron.index<<"\t"<<data[e].photons.size()<<"\t"<<data[e].electron.position.x() / mm<<"\t"
        <<data[e].electron.position.y() / mm<<"\t"<<data[e].electron.position.z() / mm<<"\t"
        <<data[e].electron.seed_info<<std::endl;
    for (std::size_t p = 0, p_end_ = data[e].photons.size(); p!=p_end_; ++p) {
      str<<data[e].photons[p]._energy / eV<<"\t"<<data[e].photons[p]._pos.x() / mm<<"\t"
          <<data[e].photons[p]._pos.y() / mm<<"\t"<<data[e].photons[p]._pos.z() / mm<<"\t"
          <<data[e].photons[p]._time / us<<"\t"<<data[e].photons[p]._momentum.x()<<"\t"
          <<data[e].photons[p]._momentum.y()<<"\t"<<data[e].photons[p]._momentum.z()<<std::endl;
    }
  }
}

namespace gPars
{
  unsigned int Results::GetNGeneratedPhotons(void) const
  {
    unsigned int out = 0;
    for (std::size_t e = 0, e_end_ = generated_photons.size(); e!=e_end_; ++e) {
      out += generated_photons[e].photons.size();
    }
    return out;
  }
  unsigned int Results::GetNRecordedPhotons(void) const
  {
    unsigned int out = 0;
    for (std::size_t e = 0, e_end_ = recorded_photons.size(); e!=e_end_; ++e) {
      out += recorded_photons[e].photons.size();
    }
    return out;
  }

	std::string this_path;
	std::string data_path;
	bool doView;

	ProgramSetups debugging;
	Source source;
	DetectorDimensions det_dims;
	DetectorOptics det_opt;
	Results results;

	HexagonalMapping* THGEM1_mapping;

	void InitGlobals(void)
	{
		char path[FILENAME_MAX];
#if defined(__WIN32__)
		this_path = _getcwd(path, FILENAME_MAX);
#else
		this_path = getcwd(path, FILENAME_MAX);
#endif //__WIN32__
		if (!this_path.empty())
			if (this_path.back()!='/')
				this_path.push_back('/');

		data_path = "../../NBrS_THGEM_LAr_v0/data/";
		doView = true;

		std::cout<<"This path: \""<<this_path<<"\""<<std::endl;
		std::cout<<"Data path: \""<<this_path+data_path<<"\""<<std::endl;

		debugging.check_geometry_overlap = false;
		debugging.track_mapping_info_class = "TrackMappingInfo";
		debugging.teleportation_verbosity = 0;
		debugging.no_reflections = false;
		debugging.no_diffused_reflections = false;
		debugging.photon_max_time = 3 * ns; // about 1 meter full path.

		source.energy_line = -1;
		source.energy_spectrum_filename = "energy_spectrum/YAP_Ce_energies_eV_1.dat";
		source.energy_spectrum.read(data_path+source.energy_spectrum_filename);
		source.x_center = 0.0; //1 * 0.9;
		source.y_center = 0.0;
		source.z_center = 0.0;
		source.xy_radius = 3;//*sqrt(2);
		source.z_width = 0;
		source.N_events = 100;

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

		det_dims.width_interface_grid_support = 1.4;
		det_dims.width_interface_grid_frame = 5;
		det_dims.xyz_position_SingleTHGEMHole = 150;
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
		det_dims.n_SiPMs_rows = 5;
		det_dims.SiPM_device_name = "phys_SiPM";
		det_dims.PMT_device_name = "phys_PMT";
		det_dims.THGEM1_cell_name = "phys_THGEM1_cell";
		det_dims.THGEM1_cell_container_name = "phys_THGEM1_cell_container";

		det_opt.FR4_SigmaAlpha = 50;
		det_opt.StainlessSteel_SigmaAlpha = 10;
		det_opt.Cu_SigmaAlpha = 30;
		det_opt.Wire_SigmaAlpha = 5;
		det_opt.pmma_absorption_length_filename = data_path+"absorption_length/PMMA_absorption_length_eV_mm.dat";
		det_opt.pmma_rindex_filename = data_path+"refractive_index/PMMA_rindex_eV_1.dat";
		det_opt.pmma_uv_absorption_length_filename = data_path+"absorption_length/PMMA_VU_absorption_length_eV_mm.dat";

		results.generated_filename = "generated.dat";
		results.recorded_filename = "recorded.dat";;

		THGEM1_mapping = nullptr;

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

    if (debugging.no_diffused_reflections) {
      det_opt.FR4_SigmaAlpha = 0.0;
      det_opt.StainlessSteel_SigmaAlpha = 0.0;
      det_opt.Cu_SigmaAlpha = 0.0;
      det_opt.Wire_SigmaAlpha = 0.0;
    }
	}
}
