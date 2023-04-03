#include <boost/foreach.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/algorithm/string.hpp>
#include <geant4/generator/SourceSettings.hh>
#include <geant4/generator/GenElectronsPatterns.hh>
#include <geant4/generator/GenElectronsPatternsNBrS.hh>
#include <geant4/generator/GenNBrS_InTHGEM.hh>
#include <geant4/generator/GenPhotonsDirectly.hh>


VSourceSettings* CreateSourceSettings(std::string filename)
{
  VSourceSettings *out = nullptr;
  // Create an empty property tree object
  using boost::property_tree::ptree;
  using boost::property_tree::ptree_bad_data;
  using boost::property_tree::ptree_bad_path;
  using boost::property_tree::ptree_error;
  ptree pt;
  try {
    read_xml(filename, pt);
    boost::optional<ptree&> src = pt.get_child_optional("Settings.Source.NBrS_Generator");
    if (src) {
      SettingsNBrSGenerator *settings = new SettingsNBrSGenerator();
      out = settings;
      settings->NBrS_yield_factor = src->get<double>("NBrS_yield_factor", 1.0);
      settings->xy_radius_smearing = std::fabs(src->get<double>("xy_radius_smearing_mm", 0.0));
      goto read_common;
    }
    src = pt.get_child_optional("Settings.Source.PhotonsDirectly");
    if (src) {
      SettingsDirectPhotons *settings = new SettingsDirectPhotons();
      out = settings;
      settings->energy_spectrum_filename = src->get<std::string>("energy_spectrum_filename", "");
      if (settings->energy_spectrum_filename.empty()) { // Spectrum not specified, use fixed energy
        settings->energy = src->get<double>("energy_line_eV") * eV;
      } else {
        settings->energy_spectrum.read(gPars::general.data_path+settings->energy_spectrum_filename);
        if (!settings->energy_spectrum.isValid()) {
          std::cerr<<"CreateSourceSettings: Error: energy spectrum is invalid. Attempting to load fixed energy." <<std::endl;
          settings->energy = src->get<double>("energy_line_eV") * eV;
        }
      }
      std::string pattern = src->get<std::string>("Pattern");
      if (pattern == "THGEM1_hole_center") {
        settings->pattern = GenPhotonsDirectly::THGEM1_hole_center;
      } else if (pattern == "EL_gap_center") {
        settings->pattern = GenPhotonsDirectly::EL_gap_center;
      } else if (pattern == "Cathode_center") {
        settings->pattern = GenPhotonsDirectly::Cathode_center;
      } else if (pattern == "Cathode_14mm_coll") {
        settings->pattern = GenPhotonsDirectly::Cathode_14mm_coll;
      } else if (pattern == "SiPM_shading") {
        settings->pattern = GenPhotonsDirectly::SiPM_shading;
      } else if (pattern == "By_source") {
        settings->pattern = GenPhotonsDirectly::By_source;
      } else {
        G4Exception("CreateSourceSettings: ",
            "InvalidSetup", FatalException, ("Invalid Settings.Source.PhotonsDirectly.Pattern '" + pattern + "' in settings.").c_str());
      }
      settings->angle = src->get<double>("angle_to_z_deg", 0) * deg;
      goto read_common;
    }
    src = pt.get_child_optional("Settings.Source.ElectronPattern");
    if (src) {
      SettingsElectronPattern *settings = new SettingsElectronPattern();
      out = settings;
      //, RandomSquare, UniformLineX, UniformLineY, UniformSquareGrid, Uniform1Ring, Uniform2Rings, Uniform3Rings
      std::string pattern = src->get<std::string>("Pattern");
      if (pattern == "RandomCircle") {
        settings->pattern = GenElectronsPatterns::RandomCircle;
      } else if (pattern == "UniformLineX") {
        settings->pattern = GenElectronsPatterns::UniformLineX;
      } else if (pattern == "UniformLineY") {
        settings->pattern = GenElectronsPatterns::UniformLineY;
      } else if (pattern == "UniformSquareGrid") {
        settings->pattern = GenElectronsPatterns::UniformSquareGrid;
      } else if (pattern == "Uniform1Ring") {
        settings->pattern = GenElectronsPatterns::Uniform1Ring;
      } else if (pattern == "Uniform2Rings") {
        settings->pattern = GenElectronsPatterns::Uniform2Rings;
      } else if (pattern == "Uniform3Rings") {
        settings->pattern = GenElectronsPatterns::Uniform3Rings;
      }  else {
        G4Exception("CreateSourceSettings: ",
            "InvalidSetup", FatalException, ("Invalid Settings.Source.ElectronPattern.Pattern '" + pattern + "' in settings.").c_str());
      }
      goto read_common;
    }
    src = pt.get_child_optional("Settings.Source.ElectronPatternNBrS");
		if (src) {
			SettingsElectronPatternNBrS *settings = new SettingsElectronPatternNBrS();
			out = settings;
			//, RandomSquare, UniformLineX, UniformLineY, UniformSquareGrid, Uniform1Ring, Uniform2Rings, Uniform3Rings
			std::string pattern = src->get<std::string>("Pattern");
			if (pattern == "RandomCircle") {
				settings->pattern = GenElectronsPatterns::RandomCircle;
			} else if (pattern == "UniformLineX") {
				settings->pattern = GenElectronsPatterns::UniformLineX;
			} else if (pattern == "UniformLineY") {
				settings->pattern = GenElectronsPatterns::UniformLineY;
			} else if (pattern == "UniformSquareGrid") {
				settings->pattern = GenElectronsPatterns::UniformSquareGrid;
			} else if (pattern == "Uniform1Ring") {
				settings->pattern = GenElectronsPatterns::Uniform1Ring;
			} else if (pattern == "Uniform2Rings") {
				settings->pattern = GenElectronsPatterns::Uniform2Rings;
			} else if (pattern == "Uniform3Rings") {
				settings->pattern = GenElectronsPatterns::Uniform3Rings;
			}  else {
				G4Exception("CreateSourceSettings: ",
						"InvalidSetup", FatalException, ("Invalid Settings.Source.ElectronPattern.Pattern '" + pattern + "' in settings.").c_str());
			}
			settings->NBrS_yield_factor = src->get<double>("NBrS_yield_factor", 1.0);
			goto read_common;
		}
    std::cerr<<"CreateSourceSettings: no source type (NBrS_Generator, PhotonsDirectly or ElectronPattern) was found."<<std::endl;
    goto fail_load;
   read_common:
    // Common settings
    out->N_events = src->get<unsigned int>("N_events");
    out->x_center = src->get<double>("x_center_mm") * mm;
    out->y_center = src->get<double>("y_center_mm") * mm;
    out->z_center = src->get<double>("z_center_mm") * mm;
    out->xy_radius = src->get<double>("xy_radius_mm") * mm;
    out->z_width = src->get<double>("z_width_mm") * mm;

  } catch (ptree_bad_path& e) {
    std::cout << "CreateSourceSettings: ptree_bad_path exception:" << std::endl;
    std::cout << e.what() << std::endl;
    goto fail_load;
  } catch (ptree_bad_data& e) {
    std::cout << "CreateSourceSettings: ptree_bad_data exception:" << std::endl;
    std::cout << e.what() << std::endl;
    goto fail_load;
  } catch (ptree_error& e) {
    std::cout << "CreateSourceSettings: ptree_error exception:" << std::endl;
    std::cout << e.what() << std::endl;
    goto fail_load;
  } catch (std::exception& e) {
    std::cout << "CreateSourceSettings: std::exception:" << std::endl;
    std::cout << e.what() << std::endl;
    goto fail_load;
  }
  return out;
fail_load:
  std::cerr << "CreateSourceSettings: Failed to load source settings \"" << filename << "\"" << std::endl;
  delete out;
  return nullptr;
}

