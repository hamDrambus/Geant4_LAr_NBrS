#include <chrono>
#include <thread>

#include "KryptonPropertiesTables.hh"


void KryptonPropertiesTables::Initialize(void)
{
	pedantic = gPars::medium_props.pedantic_calculations;
	verbosity = gPars::medium_props.print_calculations ? 2 : verbosity;
	if (verbosity > 0)
		std::cout<<"Starting "<<gPars::MediumName()<<" properties initialization..."<<std::endl;

	std::ifstream str;
	str.open(gPars::medium_props.XS_energy_transfer_filename);
	if (str.is_open()) {
		XS_energy_transfer.read(str);
		XS_energy_transfer.scaleXY(eV, cm*cm); // See data files' units
		str.close();
	}
	if (!XS_energy_transfer.isValid()) {
		G4Exception((std::string("GlobalData::") + classname + "::Initialize: ").c_str(),
					"InvalidData", FatalException, "Electron-medium energy-transfer cross section is not loaded!");
		return;
	}
	str.open(gPars::medium_props.XS_momentum_transfer_filename);
	if (str.is_open()) {
		XS_momentum_transfer.read(str);
		XS_momentum_transfer.scaleXY(eV, cm*cm); // See data files' units
		str.close();
	}
	if (!XS_momentum_transfer.isValid()) {
		G4Exception((std::string("GlobalData::") + classname + "::Initialize: ").c_str(),
					"InvalidData", FatalException, "Electron-medium momentum-transfer cross section is not loaded!");
		return;
	}

	if (!gPars::medium_props.force_recalculation) {
		if (LoadCached()) {
			if (verbosity > 0)
				std::cout<<"Finished "<<gPars::MediumName()<<" properties initialization."<<std::endl;
			return;
		}
	}

	copy_file(gPars::general.settings_filename, gPars::medium_props.cache_folder + "settings_used_for_calculations.xml");
	const double kVcm = kilovolt / cm;
	double sfc = gPars::medium_props.field_step_factor; // Step factor to test results dependencies on step size.
	interval_field_NBrS = IntegrationInterval(80 * kVcm, 1050.0 * kVcm, sfc * 5 * kVcm);
	interval_field_NBrS += IntegrationInterval(10.0 * kVcm, 80.0 * kVcm, sfc * 1 * kVcm);
	interval_field_NBrS += IntegrationInterval(3.0 * kVcm, 12.0 * kVcm, sfc * 0.2 * kVcm);

	interval_field_drift = IntegrationInterval(0.01 * kVcm, 1.0 * kVcm, sfc * 0.02 * kVcm);
	interval_field_drift += IntegrationInterval(1 * kVcm, 1.5 * kVcm, sfc * 0.05 * kVcm);
	interval_field_drift += IntegrationInterval(1.5 * kVcm, 3.0 * kVcm, sfc * 0.1 * kVcm);
	interval_field_drift += interval_field_NBrS;

	low_high_threshold_field = interval_field_NBrS.min();
	max_field = interval_field_drift.max();
	min_field = interval_field_drift.min();

	sfc = gPars::medium_props.distributions_energy_step_factor;
	interval_XS = IntegrationInterval(0 * eV, 0.001 * eV, sfc * 0.0002 * eV);
	interval_XS += IntegrationInterval(0 * eV, 0.01 * eV, sfc * 0.001 * eV);
	interval_XS += IntegrationInterval(0.01 * eV, 1 * eV, sfc * 0.04 * eV);
	interval_XS += IntegrationInterval(1 * eV, 30 * eV, sfc * 0.4 * eV);

	sfc = gPars::medium_props.spectra_step_factor;
	interval_photon_En = IntegrationInterval(hc / (1000 * nm), hc / (50 * nm), sfc * 0.2 * 0.02 * eV); // from 50 to 1000 nm photons
	interval_photon_En += IntegrationInterval(hc / (1000 * nm), hc / (700 * nm), sfc * 0.2 * 0.002 * eV);
	min_f_value_probability = 1e-10 * pow(joule, -1.5); // to avoid division by 0

	if (2 == verbosity)
		PlotInputXSs();

	bool recalculate = gPars::medium_props.force_recalculation;
	str.open(gPars::medium_props.cache_folder + filename_electron_distributions, std::ios_base::binary);
	if (!str.is_open())
		recalculate = true;
	else {
		electron_distributions.read(str);
		if (electron_distributions.is_empty())
			recalculate = true;
		str.close();
	}
	if (recalculate) {
		CalcElectronDistributions();
		electron_distributions.write(gPars::medium_props.cache_folder + filename_electron_distributions);
	}
	if (electron_distributions.is_empty()) {
		G4Exception((std::string("GlobalData::") + classname + "::Initialize: ").c_str(),
					"InvalidData", FatalException, "Electron energy distributions are empty.");
		return;
	}
	recalculate = gPars::medium_props.force_recalculation;

	if (!drift_velocity.isValid()) {
		CalcDriftVelocity(); // F_distributions depend on drift velocity
		drift_velocity.write(gPars::medium_props.cache_folder + filename_drift_velocity, "Field[MV/mm]\tVelocity[mm/ns]\t(Geant4 units)");
	}

	str.open(gPars::medium_props.cache_folder + filename_F_distributions, std::ios_base::binary);
	if (!str.is_open())
		recalculate = true;
	else {
		F_distributions.read(str);
		if (F_distributions.is_empty())
			recalculate = true;
		str.close();
	}
	if (recalculate) {
		CalcFDistributions();
		F_distributions.write(gPars::medium_props.cache_folder + filename_F_distributions);
	}
	if (F_distributions.is_empty()) {
		G4Exception((std::string("GlobalData::") + classname + "::Initialize: ").c_str(),
					"InvalidData", FatalException, "F distributions are empty.");
		return;
	}
	recalculate = gPars::medium_props.force_recalculation;

	if (!drift_diffusion_T.isValid() || recalculate) {
		CalcDiffusionT();
		drift_diffusion_T.write(gPars::medium_props.cache_folder + filename_drift_diffusion_T, "Field[MV/mm]\tDiffusion transverse[mm^2/ns]\t(Geant4 units)");
	}
	if (!drift_diffusion_L.isValid() || recalculate) {
		CalcDiffusionL();
		drift_diffusion_L.write(gPars::medium_props.cache_folder + filename_drift_diffusion_L, "Field[MV/mm]\tDiffusion longitudinal[mm^2/ns]\t(Geant4 units)");
	}
	if (verbosity == 2) {
		std::cout<<"Plotting electron diffusion coefficients..."<<std::endl;
		PlotDiffusions();
	}
	if (!yield.isValid() || spectra.is_empty() || recalculate) {
		CalcYeildAndSpectra();
		spectra.write(gPars::medium_props.cache_folder + filename_spectra);
		yield.write(gPars::medium_props.cache_folder + filename_yield, "Field[MV/mm]\tYield[photons/mm]\t(Geant4 units)");
	}

	F_distributions.clear();
	electron_distributions.clear();
	if (verbosity > 0)
		std::cout<<"Finished "<<gPars::MediumName()<<" properties initialization..."<<std::endl;
}

IntegrationRange KryptonPropertiesTables::GetIntervalElectronDistributions(double field) const
{
  const double kVcm = kilovolt / cm;
  const double density = 1.73e22 / cm3; // The parameters were determined at this LKr density.
  // The values below are determined using iterative calculations of electron distributions
  // E2s here must be > E3s in GetIntervalFDistributions
  std::vector<double> fields = {0, 0.01, 0.07,  0.11, 0.5, 1,    3,    6,    10,   20,   30,   40,   70,   150,  300,  600,  1050, 1200};
  std::vector<double> E1s =    {0, 1e-3, 1e-2,  3e-2, 0.2, 0.4,  0.8,  2.0,  2.0,  2.5,  2.6,  2.8,  3.2,  4.0,  4.5,  6.0,  7.8,  8.4};
  std::vector<double> E2s =    {0, 4e-2, 0.15,  3e-1, 0.9, 1.6,  3.0,  2.6,  3.0,  3.7,  3.8,  4.2,  5.0,  6.8,  10.2, 16.0, 21.5, 23.3};
  DataVector Es_flat(fields, E1s, 2, 3); // 2nd order interpolation!
  DataVector Es_end(fields, E2s, 1, 2);
  Es_flat.scaleXY(kVcm / density, eV);
  Es_end.scaleXY(kVcm / density, eV);
  const double ref_field = 150 * kVcm / density;
  const double E1 = Es_flat(ref_field);
  const double E2 = Es_end(ref_field);
  // Integration interval, describing region of interest for electron distribution
  // at ~1 Td (150 kV/cm)
  const double sfc = gPars::medium_props.distributions_energy_step_factor; // Step factor to test results dependency on step size.
  IntegrationInterval reference_flat = IntegrationInterval(0 * eV, E1, sfc * 0.08 * eV); // flat part
  IntegrationInterval reference_decline = IntegrationInterval(E1, E2, sfc * 0.10 * eV); // fast changing (large derivative)
  const double reduced_field = field / gPars::medium_props.atomic_density;
  reference_flat.Rescale(0, Es_flat(reduced_field));
  reference_decline.Rescale(Es_flat(reduced_field), Es_end(reduced_field));
  return reference_flat + reference_decline;
}

IntegrationRange KryptonPropertiesTables::GetIntervalFDistributions(double field) const
{
  const double kVcm = kilovolt / cm;
  const double density = 1.73e22 / cm3; // The parameters were determined at this LKr density.
  // The values below are determined using iterative calculations of electron distributions
  std::vector<double> fields = {0, 0.01,    0.07,   0.11,   0.515,  1,    3,    6,    10,   20,   30,   40,   70,   150,  300,  600,  1050, 1200};
  std::vector<double> E1s =    {0, 1e-3,    7e-3,   1e-2,   6e-2,   0.1,  0.3,  0.7,  1.0,  1.5,  1.7,  2.0,  2.0,  2.4,  3.0,  4.1,  5.0,  5.1};
  std::vector<double> E2s =    {0, 5e-3,    4e-2,   5e-2,   0.25,   0.55, 1.5,  1.7,  2.1,  2.4,  2.7,  2.9,  3.3,  4.0,  5.5,  7.0,  9.0,  9.5};
  std::vector<double> E3s =    {0, 2e-2,    1e-1,   2e-1,   1.0,    1.2,  2.0,  4.0,  4.5,  5.0,  5.5,  6.0,  9.0,  11.5, 15.0, 20.0, 35.0, 37.0};
  DataVector Es_1(fields, E1s, 1, 2); Es_1.scaleXY(kVcm / density, eV);
  DataVector Es_2(fields, E2s, 1, 2); Es_2.scaleXY(kVcm / density, eV);
  DataVector Es_3(fields, E3s, 1, 2); Es_3.scaleXY(kVcm / density, eV);
  const double ref_field = 70 * kVcm / density;
  const double E1 = Es_1(ref_field);
  const double E2 = Es_2(ref_field);
  const double E3 = Es_3(ref_field);
  // Integration interval, describing region of interest for electron distribution
  // at 0.405 Td (70 kV/cm)
  const double sfc = gPars::medium_props.distributions_energy_step_factor; // Step factor to test results dependency on step size.
  IntegrationInterval reference_rise = IntegrationInterval(0 * eV, E1, sfc * 0.02 * eV); // fast changing (large derivative)
  IntegrationInterval reference_flat = IntegrationInterval(E1, E2, sfc * 0.12 * eV);
  IntegrationInterval reference_decline = IntegrationInterval(E2, E3, sfc * 0.045 * eV); // fast changing (large derivative)
  const double reduced_field = field / gPars::medium_props.atomic_density;
  reference_rise.Rescale(0, Es_1(reduced_field));
  reference_flat.Rescale(Es_1(reduced_field), Es_2(reduced_field));
  reference_decline.Rescale(Es_2(reduced_field), Es_3(reduced_field));
  return reference_rise + reference_flat + reference_decline;
}
