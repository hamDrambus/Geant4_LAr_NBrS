#include <boost/foreach.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/algorithm/string.hpp>
#include <geant4/detector/DetectorSettings.hh>

const std::string VDetectorDimensions::SiPM_device_name = "phys_SiPM";
const std::string VDetectorDimensions::PMT_device_name = "phys_PMT";

VDetectorDimensions* CreateDetectorSettings(std::string filename)
{
  VDetectorDimensions *out = nullptr;
  // Create an empty property tree object
  using boost::property_tree::ptree;
  using boost::property_tree::ptree_bad_data;
  using boost::property_tree::ptree_bad_path;
  using boost::property_tree::ptree_error;
  ptree pt;
  try {
    read_xml(filename, pt);
    ptree det = pt.get_child("Settings.Detector");
    std::string type = det.get<std::string>("DetectorType");
    if (type == "Full") {
      DetectorDimsFull *settings = new DetectorDimsFull();
      out = settings;
    } else if (type == "THGEM1 detailed") {
      DetectorDimsTGHEM1Detailed *settings = new DetectorDimsTGHEM1Detailed();
      out = settings;
    } else if (type == "THGEM1 shading") {
      DetectorDimsTGHEM1Shading *settings = new DetectorDimsTGHEM1Shading();
      out = settings;
    } else if (type == "Full_y2022" || type == "Full_y2024_NIR") {
    	DetectorDimsFullY2022 *settings = (type == "Full_y2022" ?  new DetectorDimsFullY2022() : new DetectorDimsFullY2024NIR());
      out = settings;
      settings->THGEM0_copper_thickness = det.get<double>("THGEM0_copper_thickness_mm") * mm;
      settings->THGEM0_hole_radius = det.get<double>("THGEM0_hole_radius_mm") * mm;
      settings->THGEM0_dielectric_thickness = det.get<double>("THGEM0_dielectric_thickness_mm") * mm;
      settings->THGEM0_hole_pitch = det.get<double>("THGEM0_hole_pitch_mm") * mm;
      settings->THGEM0_hole_rim = det.get<double>("THGEM0_hole_rim_mm") * mm;
      settings->THGEM0_dielectric_radius = det.get<double>("THGEM0_dielectric_radius_mm", (settings->THGEM0_hole_radius) / mm) * mm;
      settings->is_NBrS_in_THGEM0 = det.get<bool>("is_NBrS_in_THGEM0", true);
      settings->has_WLS = det.get<bool>("has_WLS", true);

			// Finalizing
      settings->THGEM0_width_total = 2*settings->THGEM0_copper_thickness + settings->THGEM0_dielectric_thickness;
      settings->THGEM0_container_width = settings->THGEM0_width_total + 0.01 * 2; // 0.01 mm from each side. Different value from Elmer simulation
      if (settings->THGEM0_dielectric_radius <= 0) {
      	std::cerr<<"THGEM0 has negative or zero radius of dielectric hole! r = "<<settings->THGEM0_dielectric_radius<<std::endl;
      	goto fail_load;
      }
    } else {
      std::cerr<<"Unknown detector type is specified: '"<<type<<"'!"<<std::endl;
      goto fail_load;
    }
    // Common settings
    out->THGEM1_copper_thickness = det.get<double>("THGEM1_copper_thickness_mm") * mm;
    out->THGEM1_hole_radius = det.get<double>("THGEM1_hole_radius_mm") * mm;
    out->THGEM1_dielectric_thickness = det.get<double>("THGEM1_dielectric_thickness_mm") * mm;
    out->THGEM1_hole_pitch = det.get<double>("THGEM1_hole_pitch_mm") * mm;
    out->THGEM1_hole_rim = det.get<double>("THGEM1_hole_rim_mm") * mm;
    out->THGEM1_dielectric_radius = det.get<double>("THGEM1_dielectric_radius_mm", (out->THGEM1_hole_radius) / mm) * mm;

    // Finalizing
    out->THGEM1_width_total = 2*out->THGEM1_copper_thickness + out->THGEM1_dielectric_thickness;
    out->THGEM1_container_width = out->THGEM1_width_total + 0.01 * 2; // 0.01 mm from each side. Different value from Elmer simulation
    if (out->THGEM1_dielectric_radius <= 0) {
			std::cerr<<"THGEM1 has negative or zero radius of dielectric hole! r = "<<out->THGEM1_dielectric_radius<<std::endl;
			goto fail_load;
		}
  } catch (ptree_bad_path& e) {
    std::cout << "CreateDetectorSettings: ptree_bad_path exception:" << std::endl;
    std::cout << e.what() << std::endl;
    goto fail_load;
  } catch (ptree_bad_data& e) {
    std::cout << "CreateDetectorSettings: ptree_bad_data exception:" << std::endl;
    std::cout << e.what() << std::endl;
    goto fail_load;
  } catch (ptree_error& e) {
    std::cout << "CreateDetectorSettings: ptree_error exception:" << std::endl;
    std::cout << e.what() << std::endl;
    goto fail_load;
  } catch (std::exception& e) {
    std::cout << "CreateDetectorSettings: std::exception:" << std::endl;
    std::cout << e.what() << std::endl;
    goto fail_load;
  }
  return out;
fail_load:
  std::cerr << "CreateDetectorSettings: Failed to load detector settings \"" << filename << "\"" << std::endl;
  delete out;
  return nullptr;
}

