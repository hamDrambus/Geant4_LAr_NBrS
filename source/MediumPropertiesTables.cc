#include <chrono>
#include <thread>

#include <boost/numeric/odeint.hpp>
#include "utilities/DenseObserver.hpp"
#include "utilities/IntegrationConditional.hh"

#include "MediumPropertiesTables.hh"

namespace odeint = boost::numeric::odeint;

const std::string MediumPropertiesTables::filename_drift_velocity = "drift_velocity.dat";
const std::string MediumPropertiesTables::filename_drift_diffusion_L = "diffusion_longitudinal.dat";
const std::string MediumPropertiesTables::filename_drift_diffusion_T = "diffusion_transverse.dat";
const std::string MediumPropertiesTables::filename_yield = "NBrS_yield.dat";
const std::string MediumPropertiesTables::filename_spectra = "NBrS_spectra.dat";
const std::string MediumPropertiesTables::filename_electron_distributions = "electron_energy_distributions.dat";
const std::string MediumPropertiesTables::filename_electron_distributions_derivative = "electron_energy_distributions_deriv.dat";
const std::string MediumPropertiesTables::filename_F_distributions = "F_distributions.dat";
const std::string MediumPropertiesTables::filename_F_distributions_derivative = "F_distributions_deriv.dat";
const std::string MediumPropertiesTables::filename_XS_energy = "XS_energy_transfer_spline.dat";
const std::string MediumPropertiesTables::filename_XS_momentum = "XS_momentum_transfer_spline.dat";

bool MediumPropertiesTables::LoadCached(void)
{
  bool calculations_required = false;
  std::ifstream str;
  str.open(gPars::medium_props.cache_folder + filename_drift_velocity);
  if (str.is_open()) {
    drift_velocity.read(str);
    if (!drift_velocity.isValid())
      calculations_required = true;
    str.close();
  } else
    calculations_required = true;

  str.open(gPars::medium_props.cache_folder + filename_drift_diffusion_L);
  if (str.is_open()) {
    drift_diffusion_L.read(str);
    if (!drift_diffusion_L.isValid())
      calculations_required = true;
    str.close();
  } else
    calculations_required = true;

  str.open(gPars::medium_props.cache_folder + filename_drift_diffusion_T);
  if (str.is_open()) {
    drift_diffusion_T.read(str);
    if (!drift_diffusion_T.isValid())
      calculations_required = true;
    str.close();
  } else
    calculations_required = true;

  str.open(gPars::medium_props.cache_folder + filename_yield);
  if (str.is_open()) {
    yield.read(str);
    if (!yield.isValid())
      calculations_required = true;
    str.close();
  } else
    calculations_required = true;

  str.open(gPars::medium_props.cache_folder + filename_spectra, std::ios_base::binary);
  if (str.is_open()) {
    spectra.read(str);
    if (spectra.is_empty())
      calculations_required = true;
    str.close();
  } else
    calculations_required = true;

  return !calculations_required;
}

void MediumPropertiesTables::Initialize(void)
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

	InitializeFields();

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
	str.open(gPars::medium_props.cache_folder + filename_electron_distributions_derivative, std::ios_base::binary);
	if (!str.is_open())
		recalculate = true;
	else {
		electron_distributions_derivative.read(str);
		if (electron_distributions_derivative.is_empty())
			recalculate = true;
		str.close();
	}
	if (recalculate) {
		CalcElectronDistributions();
		electron_distributions.write(gPars::medium_props.cache_folder + filename_electron_distributions);
		electron_distributions_derivative.write(gPars::medium_props.cache_folder + filename_electron_distributions_derivative);
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
	str.open(gPars::medium_props.cache_folder + filename_F_distributions_derivative, std::ios_base::binary);
	if (!str.is_open())
		recalculate = true;
	else {
		F_distributions_derivative.read(str);
		if (F_distributions_derivative.is_empty())
			recalculate = true;
		str.close();
	}
	if (recalculate) {
		CalcFDistributions();
		F_distributions.write(gPars::medium_props.cache_folder + filename_F_distributions);
		F_distributions_derivative.write(gPars::medium_props.cache_folder + filename_F_distributions_derivative);
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
	F_distributions_derivative.clear();
	electron_distributions.clear();
	electron_distributions_derivative.clear();
	if (verbosity > 0)
		std::cout<<"Finished "<<gPars::MediumName()<<" properties initialization..."<<std::endl;
}

void MediumPropertiesTables::CalcElectronDistributions(void)
{
  if (verbosity > 0)
    std::cout<<"Calculating electron energy distributions..."<<std::endl;
  for (std::size_t i = 0, i_end_ = interval_field_drift.NumOfIndices(); i!=i_end_; ++i) {
    double field = interval_field_drift.Value(i);
    std::pair<DataVector, DataVector> e_distr = CalcElectronDistribution(field);
    electron_distributions.push(field, e_distr.first);
    electron_distributions_derivative.push(field, e_distr.second);
  }
  if (verbosity == 2) {
    std::cout<<"Plotting electron energy distributions..."<<std::endl;
    PlotElectronDistributions();
  }
}

void MediumPropertiesTables::CalcDriftVelocity(void)
{
  if (verbosity > 0)
    std::cout<<"Calculating electron drift velocity..."<<std::endl;
  drift_velocity.clear();
  drift_velocity.set_leftmost(0.0); // Velocity is 0 at negative absolute field.
  drift_velocity.use_leftmost(false); drift_velocity.use_rightmost(true);
  drift_velocity.setOrder(1); drift_velocity.setNused(2); // Linear interpolation (energy points are quite dense).
  for (std::size_t i = 0, i_end_ = interval_field_drift.NumOfIndices(); i!=i_end_; ++i) {
    double field = interval_field_drift.Value(i);
    double vel = CalcDriftVelocity(field);
    drift_velocity.insert(field, vel);
  }
  if (verbosity == 2) {
    std::cout<<"Plotting electron drift velocity..."<<std::endl;
    PlotDriftVelocity();
  }
}

void MediumPropertiesTables::CalcFDistributions(void)
{
  if (verbosity > 0)
    std::cout<<"Calculating F distributions..."<<std::endl;
  for (std::size_t i = 0, i_end_ = interval_field_drift.NumOfIndices(); i!=i_end_; ++i) {
    double field = interval_field_drift.Value(i);
    std::pair<DataVector, DataVector> F_distr = CalcFDistribution(field);
    F_distributions.push(field, F_distr.first);
    F_distributions_derivative.push(field, F_distr.second);
  }
  if (verbosity == 2) {
    std::cout<<"Plotting F distributions..."<<std::endl;
    PlotFDistributions();
  }
}

void MediumPropertiesTables::CalcDiffusionT(void)
{
  if (verbosity > 0)
    std::cout<<"Calculating electron transverse diffusion..."<<std::endl;
  drift_diffusion_T.clear();
  drift_diffusion_T.set_leftmost(0.0); // Velocity is 0 at negative absolute field.
  drift_diffusion_T.use_leftmost(false); drift_diffusion_T.use_rightmost(true);
  drift_diffusion_T.setOrder(1); drift_diffusion_T.setNused(2); // Linear interpolation (energy points are quite dense).
  for (std::size_t i = 0, i_end_ = interval_field_drift.NumOfIndices(); i!=i_end_; ++i) {
    double field = interval_field_drift.Value(i);
    double diff = CalcDiffusionT(field);
    if (diff == DBL_MAX || diff == -DBL_MAX)
      continue;
    drift_diffusion_T.insert(field, diff);
  }
}

void MediumPropertiesTables::CalcDiffusionL(void)
{
  if (verbosity > 0)
    std::cout<<"Calculating electron longitudinal diffusion..."<<std::endl;
  drift_diffusion_L.clear();
  drift_diffusion_L.set_leftmost(0.0); // Velocity is 0 at negative absolute field.
  drift_diffusion_L.use_leftmost(false); drift_diffusion_L.use_rightmost(true);
  drift_diffusion_L.setOrder(1); drift_diffusion_L.setNused(2); // Linear interpolation (energy points are quite dense).
  for (std::size_t i = 0, i_end_ = interval_field_drift.NumOfIndices(); i!=i_end_; ++i) {
    double field = interval_field_drift.Value(i);
    double diff = CalcDiffusionL(field);
    if (diff == DBL_MAX || diff == -DBL_MAX)
      continue;
    drift_diffusion_L.insert(field, diff);
  }
}

DataVector MediumPropertiesTables::CalcDiffXS(double field) const
{
  switch (gPars::medium_props.NBrS_formula) {
  case(gPars::NBrSFormula::ElasticXS):
    return CalcDiffXS_ElasticFormula(field);
  case(gPars::NBrSFormula::TransferXS):
    return CalcDiffXS_TransferFormula(field);
  default:
    G4Exception((classname + "::CalcDiffXS: ").c_str(),
          "InvalidValue", FatalException, (std::string("gPars::medium_props.NBrS_formula type ")
          + std::to_string(gPars::medium_props.NBrS_formula) + " is not implemented!").c_str()); // TODO: boost can reflect enums to text
    return DataVector();
  }
}

void MediumPropertiesTables::CalcYeildAndSpectra(void)
{
  if (verbosity > 0)
    std::cout<<"Calculating NBrS yield and spectra..."<<std::endl;
  for (std::size_t i = 0, i_end_ = interval_field_NBrS.NumOfIndices(); i!=i_end_; ++i) {
    double field = interval_field_NBrS.Value(i);
    DataVector diff_XS = CalcDiffXS(field);
    if (!diff_XS.isValid())
      continue;
    double yeild_v = diff_XS.getY(diff_XS.size() - 1);
    yield.push_back(field, yeild_v);
    if (yeild_v > 0) {
      diff_XS.scaleY(1.0 / yeild_v);
      spectra.push(field, diff_XS);
    }
  }
  if (verbosity == 2) {
    std::cout<<"Plotting NBrS yield and spectra..."<<std::endl;
    PlotYield();
    PlotSpectra();
    //PlotSpectraRaw();
  }
}

std::pair<DataVector, DataVector> MediumPropertiesTables::CalcElectronDistribution(double field) const // e distribution by their energy, f'.
{
  // Calculate eq. 10 in Borisova2021 (doi:doi:10.1209/0295-5075/ac4c03)
	// The equation was rewritten in a such way that integrands, energies and results are close to 1.0, so
	// that numerical algorithms are stable. To this end, cross-section(e) = XS(e) * cross-section(0).
  DataVector out;
  out.set_out_value(0.0); // Distribution returns 0 at energies outside specified energy range.
  out.use_leftmost(false); out.use_rightmost(false);
  out.setOrder(1); out.setNused(2); // Linear interpolation (energy points are quite dense).
  DataVector f_tmp = out; // To calculate normalization.
  DataVector out_deriv = out;
  DataVector integral;
  integral.use_leftmost(true); integral.use_rightmost(true);
  integral.setOrder(1); integral.setNused(2); // Linear interpolation (energy points should be quite dense).
  DataVector integral_dense = integral;

  const double delta = 2.0 * gPars::medium_props.m_to_M;
	const double atomic_density = gPars::medium_props.atomic_density * m*m*m;
	const double field_SI = field / volt * m;
	const double eq_10_coeff = 3.0 / field_SI / field_SI * delta / (m*m*m*m) * atomic_density * atomic_density;

  double reltol = gPars::medium_props.distributions_relative_tolerance, min_f = gPars::medium_props.distributions_relative_minimum;
  // min_f = f(En_max)/f(0); defines maximum energy. reltol defines precision of integration.
  // reltol for f is absolute tolerance for the integral!

  min_f = std::min(min_f, 1e-1);
  double eV_max = -2.0 * std::log(min_f) /
  		(3.0 * delta * atomic_density * atomic_density * (XS_energy_transfer(0) / (m*m)) * (XS_momentum_transfer(0) / (m*m)));
  eV_max = sqrt(eV_max) * field_SI; // First approximation of maximum energy. If XSs increase with energy, then real E max < E max approx.
  const double max_integral = -std::log(min_f);

  auto observer = [&] (const double &x, double energy) {
  	integral.insert(energy, x);
  };

  auto ode = [&] (const double &, double &dxdt, double energy) { // integrand function ( d{x}/d{energy} = dxdt(x, energy) )
  	//energy is in eV NOT in Geant4 units
		double Y = eq_10_coeff * energy * XS_energy_transfer(energy * eV) * XS_momentum_transfer(energy * eV);
		dxdt = Y;
	};
  auto stop_condition = [&] (const double &x, double , double ) -> bool {
  	return x >= max_integral;
  };
	double x = 0.0;
	odeint::integrate_adaptive_conditional_dense(odeint::runge_kutta_dopri5<double>(), reltol, 0.0,
				ode, x, 0.0, stop_condition, eV_max/100.0, observer);
	if (!integral.size())
		return std::pair<DataVector, DataVector> (out, out_deriv);

	// Stepping for integral x(e) may be too sparse for exponent of it, hence step size is adjusted using target reltol
	x = integral.getY(0);
	double integrand = 0.0;
	double e = integral.getX(0), e_next, x_next;
	out.push_back(e, std::exp(-x));
	ode(x, integrand, e);
	out_deriv.push_back(e, -std::exp(-x) * integrand);
	f_tmp.push_back(e, std::exp(-x)*sqrt(e));
	for (std::size_t i = 1, i_end_ = integral.size(); i != i_end_; ++i) {
		e_next = integral.getX(i);
		x_next = integral.getY(i);
		double e0 = e, e1 = e_next;
		double x0 = x, x1 = x_next;
		do { // TODO: optimize by saving evaluations of exp
			bool subdivided = false;
			double e_extr = (-(e1 - e0) * std::log((std::exp(-x1) - std::exp(-x0))/(x0 - x1)) + x1 * e0 - x0 * e1 ) / (x1 - x0);
			// the point at which the error of interpolating f(e) = exp(-x(e)) = ~exp(-x0 - (x1-x0) (e-e0)/(e1+e0)) linearly is maximal.
			double f_exact = std::exp(-x0 - (x1 - x0) * (e_extr - e0) / (e1 - e0));
			double f_approx = std::exp(-x0) + (std::exp(-x1) - std::exp(-x0)) * (e_extr - e0) / (e1 - e0);
			double error = std::fabs(f_exact - f_approx);
			//double ref = 1.0 / std::exp(-x0);
			while (error / f_exact > reltol) {
				subdivided = true;
				e1 = 0.5*(e0 + e1);
				x1 = 0.5*(x0 + x1);
				e_extr = (-(e1 - e0) * std::log((std::exp(-x1) - std::exp(-x0))/(x0 - x1)) + x1 * e0 - x0 * e1 ) / (x1 - x0);
				f_exact = std::exp(-x0 - (x1 - x0) * (e_extr - e0) / (e1 - e0));
				f_approx = std::exp(-x0) + (std::exp(-x1) - std::exp(-x0)) * (e_extr - e0) / (e1 - e0);
				error = std::fabs(f_exact - f_approx);
			}
			out.push_back(e1, std::exp(-x1));
			ode(x1, integrand, e1);
			out_deriv.push_back(e1, -std::exp(-x1) * integrand);
			f_tmp.push_back(e1, std::exp(-x1)*sqrt(e1));
			if (!subdivided)
				break;
			e0 = e1;
			x0 = x1;
			e1 = e_next;
			x1 = x_next;
		} while (true);
		e = e_next;
		x = x_next;
	}

  f_tmp.integrate();
  double normalization = f_tmp.getY(f_tmp.size() - 1);
  out.scaleY(1.0 / normalization);
  out_deriv.scaleY(1.0 / normalization);
  return std::pair<DataVector, DataVector> (out, out_deriv);
}

double MediumPropertiesTables::CalcDriftVelocity(double field) const
{ //eq. 11 in Borisova2021 (doi:doi:10.1209/0295-5075/ac4c03)
	// The equation was rewritten in a such way that integrands, energies and results are close to 1.0, so
	// that numerical algorithms are stable. To this end, cross-section(e) = XS(e) * cross-section(0).
  double reltol = gPars::medium_props.drift_velocity_relative_tolerance;
  auto indexes = electron_distributions_derivative.getX_indices(field);
  if (!indexes || indexes->first != indexes->second) {
  	if (pedantic)
			G4Exception((classname +"::CalcDriftVelocity: ").c_str(),
				"InvalidValue", FatalException, "Could not determine integration range (empty derivative of electron distribution f).");
		else
			return 0.0;
  }
  const DataVector & f_distr_deriv = electron_distributions_derivative.getY_data(indexes->first);

  double E_min = f_distr_deriv.minX();
	double E_max = f_distr_deriv.maxX();
	std::size_t En_range_size = f_distr_deriv.size();

  double dfde_max = std::max(std::fabs(f_distr_deriv.minY()), std::fabs(f_distr_deriv.maxY()));
  const double inv_dfde_max = 1.0 / dfde_max;

  const double atomic_density = gPars::medium_props.atomic_density * m*m*m;
  const double field_SI = field / volt * m;
  const double XS0 = XS_momentum_transfer(0);
  const double eq_11_coeff = -2.0 / 3.0 * e_SI * field_SI / e_mass_SI
  		* sqrt(e_mass_SI / 2.0) * dfde_max / atomic_density / (XS0 / (m*m)) / sqrt(e_SI); // [m/s * eV^1/2]


  auto ode = [&] (const double &x, double &dxdt, double energy) { // integrand function
  	//energy is eV NOT in Geant4 units
  	if (energy < 0) {
  		dxdt = 0;
  		return;
  	}
		double dfde = f_distr_deriv(energy);
		dxdt = energy * dfde * inv_dfde_max * XS0 / XS_momentum_transfer(energy * electronvolt);
  };

  double dE = std::max((E_max - E_min)/En_range_size, 5.0 * std::numeric_limits<double>::epsilon() * std::max(E_min, E_max) );
  double integral = 0.0;
	odeint::integrate_adaptive(odeint::make_controlled(0.0, reltol, dE, odeint::runge_kutta_dopri5<double>()),
			ode, integral, E_min, E_max, dE);
	double integral_rev = 0.0;
	odeint::integrate_adaptive(odeint::make_controlled(0.0, reltol, -dE, odeint::runge_kutta_dopri5<double>()),
			ode, integral_rev, E_max, E_min, -dE);

  double val = 0.5 * (integral - integral_rev) * eq_11_coeff * m / s; // return in Geant4 units
  if (isnan(val) || isinf(val) || val < 0) {
    if (pedantic) {
    	std::string fld = "E = " + dbl_to_str(field * cm / kilovolt, 3)
							+ " kV/cm, E/N = " + dbl_to_str(field / gPars::medium_props.atomic_density / Td, 3)
							+ " Td";
			std::string message = "Invalid value obtained (negative, nan, infinity or DBL_MAX) at " + fld;
    	G4Exception((classname +"::CalcDriftVelocity: ").c_str(),
        "InvalidValue", FatalException, message.c_str());
    }
    else
      val = 0;
  }
  return val;
}

std::pair<DataVector, DataVector> MediumPropertiesTables::CalcFDistribution(double field) const
{ // Calculate eq. 5 in Atrazhev1985 (doi:10.1088/0022-3719/18/6/015)
	// The equation was rewritten in a such way that integrands, energies and results are close to 1.0, so
	// that numerical algorithms are stable. To this end, cross-section(e) = XS(e) * cross-section(0).
  DataVector out;
  out.set_out_value(0.0); // Distribution returns 0 at energies outside specified energy range.
  out.use_leftmost(false); out.use_rightmost(false);
  out.setOrder(1); out.setNused(2); // Linear interpolation (energy points are quite dense).
  DataVector integrand = out; // for inner integral over dy, function of x. Has 1/sqrt(x) singularity when x->0
  DataVector integrand_reverse = out; // for inner integral over dy, not from 0 to x but from x to E max. Required due to limited precision.
  // otherwise, integrand(Emax) is not 0 as it must be in exact solution.
  DataVector out_deriv = out;

  auto indexes = electron_distributions.getX_indices(field);
	if (!indexes || indexes->first != indexes->second) {
		if (pedantic)
				G4Exception((classname +"::CalcFDistribution: ").c_str(),
					"InvalidValue", FatalException, "Could not determine integration range (empty electron distribution f).");
			else
				return std::pair<DataVector, DataVector>(out, out_deriv);
	}
	const DataVector & f_distr = electron_distributions.getY_data(indexes->first);
	indexes = electron_distributions_derivative.getX_indices(field);
	if (!indexes || indexes->first != indexes->second) {
		if (pedantic)
				G4Exception((classname +"::CalcFDistribution: ").c_str(),
					"InvalidValue", FatalException, "Could not determine integration range (empty electron distribution f derivative).");
			else
				return std::pair<DataVector, DataVector>(out, out_deriv);
	}
	const DataVector & f_distr_deriv = electron_distributions_derivative.getY_data(indexes->first);

  // It is always good to reduce number of multiplications inside loops. So constant factors are cached.
	const double atomic_density = gPars::medium_props.atomic_density * m*m*m;
	const double field_SI = field / volt * m;
	const double XS0 = XS_momentum_transfer(0);
	const double inv_XS0 = 1.0 / XS_momentum_transfer(0);
	const double eq_5_coeff = 1.5 * e_mass_SI / e_SI / field_SI * (drift_velocity(field) * s / m)
			* (XS0 / (m*m)) * atomic_density * sqrt(2.0 / e_mass_SI) * sqrt(e_SI); // in eV-1/2 (NOT geant4 unit)
	double reltol = gPars::medium_props.F_distributions_relative_tolerance;

  auto ode1 = [&] (const double &I, double &dIdy, double y) { // integrand function
  	// y is energy in eV NOT in Geant4 units
  	if (y < 0.0) {
  		dIdy = 0.0;
  		return;
  	}
  	double dfdy = f_distr_deriv(y);
  	double f = f_distr(y);
  	double V1 = y * XS0 * dfdy / XS_momentum_transfer(y * electronvolt);
		double V2 = sqrt(y) * f;
		dIdy = (V1 + V2 * eq_5_coeff); // in eV-3/2 (NOT Geant4 units)
  };

  auto observer1 = [&] (const double &I, double x) {
  	// x is energy in eV NOT in Geant4 units
		integrand.insert(x, I);
	};

  auto observer1_reverse = [&] (const double &I, double x) {
		// x is energy in eV NOT in Geant4 units
  	integrand_reverse.insert(x, I);
	};

  std::size_t En_range_sz = f_distr.size();
	double E_min = f_distr.minX();
	double E_max = f_distr.maxX();

	//estimating scale of integral for absolute error
	double En = E_min + (E_max-E_min)/100.0;
	double f0 = f_distr(0);
	double fmin = f0 * gPars::medium_props.distributions_relative_minimum;
	double V = std::min( E_max * (fmin * XS0 / XS_momentum_transfer(E_max * electronvolt)), En * f0);

	double dE = std::max((E_max - E_min)/(En_range_sz), 5 * std::numeric_limits<double>::epsilon() * std::max(E_max, E_min));

  double integral = 0.0;
  //Inner integral over dy from 0 to x
	odeint::integrate_adaptive_dense(odeint::runge_kutta_dopri5<double>(), reltol * V, reltol,
			ode1, integral, E_min, E_max, dE, observer1);
	integral = 0.0;
	//Inner integral over dy from Emax to x
	odeint::integrate_adaptive_dense(odeint::runge_kutta_dopri5<double>(), reltol * V, reltol,
			ode1, integral, E_max, E_min, -dE, observer1_reverse);
	integral = 0.0;
	double integrand_max = 0.0, E_max_integrand = 0.5*(E_max + E_min), E_max_integrand_rev = 0.5*(E_max + E_min);
	for (std::size_t i = 0, i_end_ = integrand.size(); i!=i_end_; ++i) {
		if (std::fabs(integrand[i].second) > integrand_max) {
			integrand_max = std::fabs(integrand[i].second);
			E_max_integrand = integrand[i].first;
		}
	}
	integrand_max = 0.0;
	for (std::size_t i = 0, i_end_ = integrand_reverse.size(); i!=i_end_; ++i) {
		if (std::fabs(integrand_reverse[i].second) > integrand_max) {
			integrand_max = std::fabs(integrand_reverse[i].second);
			E_max_integrand_rev = integrand_reverse[i].first;
		}
	}
	double E_crit = 0.5 * (E_max_integrand + E_max_integrand_rev);

	auto ode2 = [&] (const double &F, double &dFdy, double x) { // integrand function
		// x is energy in eV NOT in Geant4 units
		double V;
		if (x <= 0) {
			V = 1.0;
		} else {
			double f = f_distr(x);
			V = 1.0  + XS_momentum_transfer(x * electronvolt) * inv_XS0 / x / f * (x > E_crit ? integrand_reverse(x) : integrand(x));
			// x > E_max_integrand ... is making integrand behave correctly both near 0 (Emin) and Emax
			if (isnan(V) || isinf(V)) {
				V = 1.0;
			}
		}
		dFdy = V; // unitless
	};
	auto observer2 = [&] (const double &F, double x) {
		// x is energy in eV NOT in Geant4 units
		double f = f_distr(x);
		double dfde = f_distr_deriv(x);
		out.insert(x, f * F);
		double integrand;
		ode2(F, integrand, x);
		out_deriv.insert(x, f * integrand + dfde * F);
	};

	//Outer integral over dx
	odeint::integrate_adaptive_dense(odeint::runge_kutta_dopri5<double>(), fmin * std::pow(e_SI, 1.5) * reltol, reltol,
				ode2, integral, E_min, E_max, dE, observer2);

  return std::pair<DataVector, DataVector>(out, out_deriv);
}

double MediumPropertiesTables::CalcDiffusionT(double field) const
{
  //eq. 4 in Atrazhev1985 (doi:10.1088/0022-3719/18/6/015), rewritten without nu_m
	// The equation was rewritten in a such way that integrands, energies and results are close to 1.0, so
	// that numerical algorithms are stable. To this end, cross-section(e) = XS(e) * cross-section(0).

	auto indexes = electron_distributions.getX_indices(field);
	if (!indexes || indexes->first != indexes->second) {
		if (pedantic)
				G4Exception((classname +"::CalcDiffusionT: ").c_str(),
					"InvalidValue", FatalException, "Could not determine integration range (empty electron distribution f).");
			else
				return 0.0;
	}
	const DataVector & f_distr = electron_distributions.getY_data(indexes->first);

  const double atomic_density = gPars::medium_props.atomic_density * m*m*m;
  const double XS0 = XS_momentum_transfer(0);
  const double eq_4_coeff = sqrt(2.0) / 3.0 / sqrt(e_mass_SI) / atomic_density * std::pow(e_SI, 0.5) / (XS0 / (m*m));

  double reltol = gPars::medium_props.diffusion_relative_tolerance;

	std::size_t En_range_sz = f_distr.size();
	double E_min = f_distr.minX();
	double E_max = f_distr.maxX();

	auto ode = [&] (const double &x, double &dxdt, double energy) { // integrand function
		//energy is in eV NOT in Geant4 units
		if (energy < 0.0) {
			dxdt = 0.0;
			return;
		}
		double f = electron_distributions(field, energy);
		dxdt = energy * f * XS0 / XS_momentum_transfer(energy * electronvolt);
	};

	double dE = std::max((E_max - E_min)/(En_range_sz), 5 * std::numeric_limits<double>::epsilon() * std::max(E_max, E_min));
  double integral = 0.0;
  odeint::integrate_adaptive(odeint::make_controlled(0.0, reltol, dE, odeint::runge_kutta_dopri5<double>()),
  			ode, integral, E_min, E_max, dE);

  double val = integral * eq_4_coeff * m * m / s; // return in Geant4 units
  if (isnan(val) || val < 0 || isinf(val)) {
    if (pedantic)
      G4Exception((classname + "::CalcDiffusionT: ").c_str(),
        "InvalidValue", FatalException, "Invalid value obtained (negative, nan, infinity or DBL_MAX).");
    else
      val = 0;
  }
  return val;
}

double MediumPropertiesTables::CalcDiffusionL(double field) const
{
  // Eq. 5 in Atrazhev1985 (doi:10.1088/0022-3719/18/6/015)
  // !!! Formula for longitudial diffusion there is incorrect!
  // Second integral term has incorrect units. The source from which
  // the formula was derived (Parker and Lowke 1969 doi:10.1103/PhysRev.181.290)
  // seems correct. Fixed formula is in data/Atrazhev1985_corrected_equation_5.png
  // Or in latex format:
  // D_{||} = D_{T} - \frac{2}{3m} \int_{0}^{\infty}\left( \frac{\epsilon^{3/2}}{\nu_{m}(\epsilon)} \frac{d F (\epsilon)}{d \epsilon} + \frac{3mW}{2e\mathscr{E}} \sqrt{\epsilon} F(\epsilon) \right) d \epsilon
	// The equation was rewritten in a such way that integrands, energies and results are close to 1.0, so
	// that numerical algorithms are stable. To this end, cross-section(e) = XS(e) * cross-section(0).

	auto indexes = F_distributions.getX_indices(field);
	if (!indexes || indexes->first != indexes->second) {
		if (pedantic)
				G4Exception((classname +"::CalcDiffusionL: ").c_str(),
					"InvalidValue", FatalException, "Could not determine integration range (empty electron distribution F).");
			else
				return drift_diffusion_T(field);
	}
	const DataVector & F_distr = F_distributions.getY_data(indexes->first);
	indexes = F_distributions_derivative.getX_indices(field);
	if (!indexes || indexes->first != indexes->second) {
		if (pedantic)
				G4Exception((classname +"::CalcDiffusionL: ").c_str(),
					"InvalidValue", FatalException, "Could not determine integration range (empty electron distribution F derivative).");
			else
				return drift_diffusion_T(field);
	}
	const DataVector & F_distr_deriv = F_distributions_derivative.getY_data(indexes->first);

	std::size_t En_range_sz = F_distr.size();
	double E_min = F_distr.minX();
	double E_max = F_distr.maxX();

  const double atomic_density = gPars::medium_props.atomic_density * m*m*m;
  const double field_SI = field / volt * m;
  const double XS0 = XS_momentum_transfer(0);
  const double int1_coeff = sqrt(2.0) / 3.0 / atomic_density / (XS0 / (m*m) ) / sqrt(e_mass_SI) * std::pow(e_SI, 0.5);
  const double int2_coeff = (drift_velocity(field) * s / m) / e_SI / field_SI * e_SI;
  const double int2_to_int1_coeff = int2_coeff / int1_coeff;

  double reltol = gPars::medium_props.diffusion_relative_tolerance;

	auto ode = [&] (const double &x, double &dxdt, double energy) { // integrand function
		//energy is in eV NOT in Geant4 units
		if (energy < 0.0) {
			dxdt = 0.0;
			return;
		}
		double dFde = F_distr_deriv(energy);
		double F = F_distr(energy);
		dxdt = energy * dFde * XS0 / XS_momentum_transfer(energy * electronvolt) + sqrt(energy) * F * int2_to_int1_coeff;
	};

//	auto ode1 = [&] (const double &x, double &dxdt, double energy) { // integrand function
//		//energy is in eV NOT in Geant4 units
//		if (energy < 0.0) {
//			dxdt = 0.0;
//			return;
//		}
//		double dFde = F_distr_deriv(energy);
//		dxdt = energy * dFde * XS0 / XS_momentum_transfer(energy * electronvolt);
//	};
//
//	auto ode2 = [&] (const double &x, double &dxdt, double energy) { // integrand function
//		//energy is in eV NOT in Geant4 units
//		if (energy < 0.0) {
//			dxdt = 0.0;
//			return;
//		}
//		double F = F_distr(energy);
//		dxdt = sqrt(energy) * F;
//	};

	double dE = std::max((E_max - E_min)/(En_range_sz), 5 * std::numeric_limits<double>::epsilon() * std::max(E_max, E_min));
	double integral1 = 0.0, integral2 = 0.0;
	//odeint::integrate_adaptive(odeint::make_controlled(0.0, reltol, dE, odeint::runge_kutta_dopri5<double>()),
	//			ode1, integral1, E_min, E_max, dE);
	//odeint::integrate_adaptive(odeint::make_controlled(0.0, reltol, dE, odeint::runge_kutta_dopri5<double>()),
	//			ode2, integral2, E_min, E_max, dE);
  //double val = drift_diffusion_T(field) - (integral1 * int1_coeff + integral2 * int2_coeff) * m * m / s ; // return in Geant4 units;
	odeint::integrate_adaptive(odeint::make_controlled(0.0, reltol, dE, odeint::runge_kutta_dopri5<double>()),
				ode, integral1, E_min, E_max, dE);
	double val = drift_diffusion_T(field) - int1_coeff * (integral1) * m * m / s ; // return in Geant4 units;

  if (isnan(val) || val < 0 || isinf(val)) {
    if (pedantic) {
    	std::string fld = "E = " + dbl_to_str(field * cm / kilovolt, 3)
    	        + " kV/cm, E/N = " + dbl_to_str(field / gPars::medium_props.atomic_density / Td, 3)
    	        + " Td";
    	std::string message = "Invalid value obtained (negative, nan, infinity or DBL_MAX) at " + fld;
      G4Exception((classname + "::CalcDiffusionL: ").c_str(),
        "InvalidValue", FatalException, fld.c_str());
    }
    else
      val = drift_diffusion_T(field);
  }
  return val;
}

// Default: uses elastic XS for NBrS XS calculation
DataVector MediumPropertiesTables::CalcDiffXS_ElasticFormula(double field) const
{
	return CalcDiffXS_XSFormula(field, XS_energy_transfer);
}

// Uses momentum transfer XS for NBrS XS calculation (correct, but only for hv << E of electron https://doi.org/10.48550/arXiv.2206.01388)
DataVector MediumPropertiesTables::CalcDiffXS_TransferFormula(double field) const
{
  return CalcDiffXS_XSFormula(field, XS_momentum_transfer);
}

DataVector MediumPropertiesTables::CalcDiffXS_XSFormula(double field, const DataVector & cross_section) const
{
	// Calculate eq. 6 but with d(lambda) replaced by d(h*nu) in Borisova2021 (doi:doi:10.1209/0295-5075/ac4c03)
	// but using momentum transfer XS as discussed in https://doi.org/10.48550/arXiv.2206.01388
	// The equation was rewritten in a such way that integrands, energies and results are close to 1.0, so
	// that numerical algorithms are stable. To this end, cross-section(e) = XS(e) * cross-section(0)
	// and f(E) = f(0) * f'(E)
	DataVector out;
	out.set_out_value(0.0); // Distribution returns 0 at energies outside specified energy range.
	out.use_leftmost(false); out.use_rightmost(false);
	out.setOrder(1); out.setNused(2); // Linear interpolation (energy points are quite dense).
	DataVector Y_cdf1 = out, Y_cdf2 = out;

	auto indexes = electron_distributions.getX_indices(field);
	if (!indexes || indexes->first!=indexes->second) {
		if (pedantic)
			G4Exception((classname +"::CalcDiffXS_XSFormula: ").c_str(),
				"InvalidValue", FatalException, "Could not determine integration range (empty electron distribution f).");
		else
			return out;
	}
	const DataVector & e_distr = electron_distributions.getY_data(indexes->first);
	if (!e_distr.size()) {
		if (pedantic)
			G4Exception((classname +"::CalcDiffXS_XSFormula: ").c_str(),
				"InvalidValue", FatalException, "Could not determine integration range (empty electron distribution f).");
		else
			return out;
	}

	const double lambda_max = gPars::medium_props.maximum_lambda;
	const double rel_tol = gPars::medium_props.yield_relative_tolerance;
	const double atomic_density = gPars::medium_props.atomic_density * m*m*m;
	const double XS0 = cross_section(0);
	const double invXS0 = 1.0 / XS0;
	const double f0 = e_distr(0) / rel_tol; // To make sure the integrands are large numbers
	const double invf0 = 1.0/f0;
	const double invf0XS0 = invXS0 * invf0;
	const double eq_6_coeff = 8.0 / 3.0 * atomic_density * (classic_electr_radius/m) / (hc/(m*joule)) * sqrt(2.0 / e_mass_SI)
				/ (drift_velocity(field) * s / m) * std::pow(e_SI, 1.5) * XS0 / (m*m) * f0;


	double E_min = e_distr.minX();
	double E_max = e_distr.maxX();
	double dE = (E_max - E_min)/e_distr.size();
	E_min = hc / lambda_max / eV; // in eV NOT in Geant4 units
	if (E_min/E_max > 1.0 - 2.0 * std::numeric_limits<double>::epsilon())
		return out;
	dE = std::max(std::min(dE, (E_max - E_min)/100.0), 5 * std::numeric_limits<double>::epsilon() * std::max(E_min, E_max));

	auto ode2 = [&] (const double &Y, double &dYdw, double energy) { // integrand function
		//energy is in eV NOT in Geant4 units
		double photon_energy = energy;
		if (energy < E_min || (photon_energy / E_max > 1.0 - 2.0 * std::numeric_limits<double>::epsilon())) {
			dYdw = 0.0;
			return;
		}

		auto ode1 = [&] (const double &dYdw, double &dYdwde, double energy) { // integrand function
			//energy is in eV NOT in Geant4 units
			if (photon_energy / E_max > 1.0 - 2.0 * std::numeric_limits<double>::epsilon() || energy < photon_energy) {
				dYdwde = 0.0;
				return;
			}
			const double energy_diff = energy - photon_energy;
			const double XS = sqrt(energy_diff/energy) / photon_energy
					* (energy_diff * cross_section(energy * electronvolt) +
						energy * cross_section(energy_diff * electronvolt));
			const double f = e_distr(energy);
			dYdwde = energy * f * XS * invf0XS0;
		};

		double dE_loc = std::max(std::min(dE, (E_max - photon_energy)/100.0), 5 * std::numeric_limits<double>::epsilon() * E_max );
		double integral = 0.0;
		double val1, val2, val3;
		ode1(integral, val1, (0.8*photon_energy + 0.2*E_max));
		ode1(integral, val2, (0.5*photon_energy + 0.5*E_max));
		ode1(integral, val3, (0.0*photon_energy + 1.0*E_max));
		double abs_tol = rel_tol * std::max({val1, val2, val3}) * (E_max - photon_energy) * 0.1;
		odeint::integrate_adaptive(
					odeint::make_controlled(abs_tol, rel_tol, dE_loc, odeint::runge_kutta_dopri5<double>()),
					ode1, integral, photon_energy, E_max, dE_loc);
		dYdw = integral;
	};

	double prev_val = 0;
	auto observer2 = [&] (const double &Y, double energy) {
		// energy is in eV NOT in Geant4 units
		double val = Y;
		if (isnan(val) || val < 0 || isinf(val)) {
			if (pedantic)
				G4Exception((classname + "::CalcDiffXS_XSFormula: ").c_str(),
					"InvalidValue", FatalException, "Invalid value obtained (negative, nan, infinity or DBL_MAX).");
			else
				val = 0;
		}
		val = std::max(val, prev_val); // required only due to limited precision. On paper, val must always be greater than prev_val
		prev_val = val;
		Y_cdf1.insert(energy * eV, val * eq_6_coeff / m); // store in Geant4 units
	};

	auto observer2_rev = [&] (const double &Y, double energy) {
		// energy is in eV NOT in Geant4 units
		double val = Y;
		if (isnan(val) || val < 0 || isinf(val)) {
			if (pedantic)
				G4Exception((classname + "::CalcDiffXS_XSFormula: ").c_str(),
					"InvalidValue", FatalException, "Invalid value obtained (negative, nan, infinity or DBL_MAX).");
			else
				val = 0;
		}
		val = std::max(val, prev_val); // required only due to limited precision. On paper, val must always be greater than prev_val
		prev_val = val;
		Y_cdf2.insert(energy * eV, val * eq_6_coeff / m); // store in Geant4 units
	};

	double integral = 0.0;
	// Approximate value of the integral near Emax threshold (lambda min) which behaves as dE^(5/2) where dE = Emax-W
	odeint::integrate_adaptive_dense(odeint::runge_kutta_dopri5<double>(), 0.0, rel_tol,
			ode2, integral, E_min, E_max, dE, observer2);
	if (!Y_cdf1.size())
		return Y_cdf1;
	const double fraction = 0.85;
	double half_Y = Y_cdf1.getY(Y_cdf1.size()-1) * fraction;
	indexes = Y_cdf1.getY_indices(half_Y);
	if (!indexes || indexes->first == 0)
		return Y_cdf1;
	half_Y = Y_cdf1.getY(indexes->first);
	double half_E = Y_cdf1.getX(indexes->first) / eV;
	integral = 0.0;
	prev_val = 0.0;
	// This duplicated integration is required so that spectrum is smooth at small lambdas (large energies)
	odeint::integrate_adaptive_dense(odeint::runge_kutta_dopri5<double>(), 0.0, rel_tol,
			ode2, integral, half_E, E_max, dE, observer2_rev);
	if (!Y_cdf2.size()) {
		return Y_cdf1;
	}
	// Merging solutions
	out.reserve(indexes->first + Y_cdf2.size());
	for(std::size_t i = 0, i_end_ = indexes->first + 1; i!=i_end_; ++i)
		out.push_back(Y_cdf1[i]);
	for(std::size_t i = 1, i_end_ = Y_cdf2.size(); i!=i_end_; ++i)
		out.push_back(Y_cdf2.getX(i), half_Y + Y_cdf2.getY(i));
	return out;
}

DataVector MediumPropertiesTables::CalcDiffXS_ExactFormula(double field) const
{
  return DataVector();
}

bool MediumPropertiesTables::IsReady(void) const
{
  return (drift_velocity.isValid() && drift_diffusion_L.isValid()
      && drift_diffusion_T.isValid() && yield.isValid()
      && !spectra.is_empty());
}

double MediumPropertiesTables::GenPhotonEnergy(double field, double rnd) const
{
  if (spectra.is_empty())
    return 0.0;
  double f_min = spectra.getX(0);
  if (field < f_min)
    return 0.0;
  return spectra.getY(field, rnd);
}

void MediumPropertiesTables::PlotInputXSs(void) const
{
  // Note: png terminal cannot plot several graphs with replot.
  // So all data must be plotted in single plot command with comma.
  std::ofstream str;
  FILE* gnuplotPipe = popen(gPars::general.gnuplot_bin.c_str(), "w");
  if (nullptr == gnuplotPipe) {
    if (verbosity > 0) {
      std::cerr<<classname<<"::PlotInputXSs:Error: could not open gnuplot."<<std::endl;
      std::cerr<<"\tNot plotting."<<std::endl;
    }
    return;
  }
  std::stringstream gnuplotCommands;
  gnuplotCommands
      <<"set terminal png size 1000,800 enhanced"<<std::endl
      <<"set output '"+gPars::medium_props.cache_folder+"electron-atom_XSs.png'"<<std::endl
      <<"set grid"<<std::endl
      <<"set logscale xy"<<std::endl
      <<"set xrange[0.001:50]"<<std::endl
      <<"set xlabel \"Electron energy [eV]\""<<std::endl
      <<"set yrange [1e-18:1e-14]"<<std::endl
      <<"set ylabel \"Cross-section [cm^{2}]"<<std::endl <<std::endl;
  fputs(gnuplotCommands.str().c_str(), gnuplotPipe);
  fflush(gnuplotPipe);
  gnuplotCommands.str(std::string());
  gnuplotCommands
    << "undefine $Mydata1"<<std::endl
    << "$Mydata1 << EOD"<<std::endl;
  for(std::size_t j = 0, j_end_ = XS_energy_transfer.size(); j!=j_end_; ++j) {
    double X = XS_energy_transfer.getX(j) / eV;
    double Y = XS_energy_transfer.getY(j) / cm / cm;
    gnuplotCommands // Passing data through pipe together with commands
    << X << "\t" << Y << ""<<std::endl;
  }
  gnuplotCommands << "EOD"<<std::endl;
  fputs(gnuplotCommands.str().c_str(), gnuplotPipe);
  fflush(gnuplotPipe);

  open_output_file(gPars::medium_props.cache_folder + filename_XS_energy, str, std::ios_base::trunc);
  if (str.is_open())
    str << "//Energy[eV]\tElastic XS[cm^2]" <<std::endl;
  gnuplotCommands.str(std::string());
  gnuplotCommands
    << "undefine $Mydata2"<<std::endl
    << "$Mydata2 << EOD"<<std::endl;
  IntegrationRange interval_XS = IntegrationInterval(0, 0.01, 0.002);
  interval_XS += IntegrationInterval(0.01, 0.1, 0.01);
  interval_XS += IntegrationInterval(0.1, 1, 0.1);
  interval_XS += IntegrationInterval(1, 20, 0.01);
  interval_XS += IntegrationInterval(20, 100, 0.5);
  interval_XS.Trim(0.0,
  		std::max(XS_energy_transfer.getX(XS_energy_transfer.size()-1) / eV,
  				XS_momentum_transfer.getX(XS_momentum_transfer.size()-1) / eV));
  for(std::size_t j = 0, j_end_ = interval_XS.NumOfIndices(); j!=j_end_; ++j) {
    double X = interval_XS.Value(j);
    double Y = XS_energy_transfer(interval_XS.Value(j) * eV) / cm / cm;
    gnuplotCommands // Passing data through pipe together with commands
    << X << "\t" << Y << ""<<std::endl;
    if (str.is_open())
      str << X << "\t" << Y <<std::endl;
  }
  gnuplotCommands << "EOD"<<std::endl;
  fputs(gnuplotCommands.str().c_str(), gnuplotPipe);
  fflush(gnuplotPipe);
  str.close();

  gnuplotCommands.str(std::string());
  gnuplotCommands
    << "undefine $Mydata3"<<std::endl
    << "$Mydata3 << EOD"<<std::endl;
  for(std::size_t j = 0, j_end_ = XS_momentum_transfer.size(); j!=j_end_; ++j) {
    gnuplotCommands // Passing data through pipe together with commands
    << XS_momentum_transfer.getX(j) / eV
    << "\t" << XS_momentum_transfer.getY(j) / cm / cm
    << ""<<std::endl;
  }
  gnuplotCommands << "EOD"<<std::endl;
  fputs(gnuplotCommands.str().c_str(), gnuplotPipe);
  fflush(gnuplotPipe);
  gnuplotCommands.str(std::string());

  open_output_file(gPars::medium_props.cache_folder + filename_XS_momentum, str, std::ios_base::trunc);
  if (str.is_open())
    str << "//Energy[eV]\tMomentun transfer XS[cm^2]" <<std::endl;
  gnuplotCommands
    << "undefine $Mydata4"<<std::endl
    << "$Mydata4 << EOD"<<std::endl;
  for(std::size_t j = 0, j_end_ = interval_XS.NumOfIndices(); j!=j_end_; ++j) {
    double X = interval_XS.Value(j);
    double Y = XS_momentum_transfer(interval_XS.Value(j) * eV) / cm / cm;
    gnuplotCommands // Passing data through pipe together with commands
    << X << "\t" << Y << ""<<std::endl;
    if (str.is_open())
      str << X << "\t" << Y <<std::endl;
  }
  gnuplotCommands << "EOD"<<std::endl;
  fputs(gnuplotCommands.str().c_str(), gnuplotPipe);
  fflush(gnuplotPipe);
  str.close();
  gnuplotCommands.str(std::string());

  std::string title1 = "XS energy transfer input";
  std::string title2 = "XS energy transfer spline";
  std::string title3 = "XS momentum transfer input";
  std::string title4 = "XS momentum transfer spline";
  gnuplotCommands << "plot $Mydata1 using 1:2 w p pt 5 ps 2 lc rgb 'black' title \"" << title1 <<"\",";
  gnuplotCommands << " $Mydata2 using 1:2 w l lt 1 lc rgb 'black' title \"" << title2 <<"\",";
  gnuplotCommands << " $Mydata3 using 1:2 w p pt 7 ps 2 lc rgb 'black' title \"" << title3 <<"\",";
  gnuplotCommands << " $Mydata4 using 1:2 w l lt 0 lc rgb 'black' title \"" << title4 <<"\""<<std::endl;
  fputs(gnuplotCommands.str().c_str(), gnuplotPipe);
  fflush(gnuplotPipe);

  pclose(gnuplotPipe);
}

void MediumPropertiesTables::PlotElectronDistributions(void) const
{
  FILE* gnuplotPipe = popen(gPars::general.gnuplot_bin.c_str(), "w");
  if (nullptr == gnuplotPipe) {
    if (verbosity > 0) {
      std::cerr<<classname<<"::PlotElectronDistributions:Error: could not open gnuplot."<<std::endl;
      std::cerr<<"\tNot plotting."<<std::endl;
    }
    return;
  }
  delete_file(gPars::medium_props.cache_folder+"electron_distributions_low_E.gif");
  std::stringstream gnuplotCommands;
  auto En_range = electron_distributions.getRangeAtX(low_high_threshold_field);
	if (!En_range) {
		return;
	}
  gnuplotCommands
      << "set terminal gif size 800,800 animate delay 100 enhanced"<<std::endl
      <<"set output '"+gPars::medium_props.cache_folder+"electron_distributions_low_E.gif'"<<std::endl
      <<"set grid"<<std::endl
      <<"set xrange [0:"<< En_range->second<<"]"<<std::endl
      <<"set xlabel \"Energy [eV]\""<<std::endl
      <<"set logscale y"<<std::endl
      <<"set yrange [1e-4:1e3]"<<std::endl
      <<"set ylabel \"f\\' [eV^{-3/2}]\""<<std::endl <<std::endl;
  fputs(gnuplotCommands.str().c_str(), gnuplotPipe);
  fflush(gnuplotPipe);
  for(std::size_t i = 0, i_end_ = electron_distributions.size(); i != i_end_; ++i) {
    if (electron_distributions.getX(i) > interval_field_NBrS.min())
      continue;
    gnuplotCommands.str(std::string());
    gnuplotCommands
      << "undefine $Mydata"<<std::endl
      << "$Mydata << EOD"<<std::endl;
    for(std::size_t j = 0, j_end_ = electron_distributions.getY_data(i).size(); j!=j_end_; ++j) {
      double X = electron_distributions.getY_data(i).getX(j);
      double Y = electron_distributions.getY_data(i).getY(j);
      gnuplotCommands // Passing data through pipe together with commands
      << X << "\t" << Y << std::endl;
    }
    gnuplotCommands << "EOD"<<std::endl;
    std::string title = "E = " + dbl_to_str(electron_distributions.getX(i) * cm / kilovolt, 3)
        + " kV/cm, E/N = " + dbl_to_str(electron_distributions.getX(i) / gPars::medium_props.atomic_density / Td, 3)
        + " Td";
    gnuplotCommands << "plot $Mydata using 1:2 w l lw 2 lc rgb 'black' title \"" << title <<"\""<<std::endl;
    fputs(gnuplotCommands.str().c_str(), gnuplotPipe);
    fflush(gnuplotPipe);
  }

  auto data = electron_distributions.getY_data(electron_distributions.size()-1);
  double En_max = data.getX(data.size()-1);
  gnuplotCommands.str(std::string());
  delete_file(gPars::medium_props.cache_folder+"electron_distributions_high_E.gif");
  gnuplotCommands
      << "set terminal gif size 800,800 animate delay 100 enhanced"<<std::endl
      <<"set output '"+gPars::medium_props.cache_folder+"electron_distributions_high_E.gif'"<<std::endl
      <<"set xrange [0:"<< En_max<<"]"<<std::endl
      <<"set yrange [1e-5:10]"<<std::endl<<std::endl;
  fputs(gnuplotCommands.str().c_str(), gnuplotPipe);
  fflush(gnuplotPipe);
  for(std::size_t i = 0, i_end_ = electron_distributions.size(); i != i_end_; ++i) {
    if (electron_distributions.getX(i) < interval_field_NBrS.min())
      continue;
    gnuplotCommands.str(std::string());
    gnuplotCommands
      << "undefine $Mydata"<<std::endl
      << "$Mydata << EOD"<<std::endl;
    for(std::size_t j = 0, j_end_ = electron_distributions.getY_data(i).size(); j!=j_end_; ++j) {
      double X = electron_distributions.getY_data(i).getX(j);
      double Y = electron_distributions.getY_data(i).getY(j);
      gnuplotCommands // Passing data through pipe together with commands
      << X << "\t" << Y << std::endl;
    }
    gnuplotCommands << "EOD"<<std::endl;
    std::string title = "E = " + dbl_to_str(electron_distributions.getX(i) * cm / kilovolt, 3)
        + " kV/cm, E/N = " + dbl_to_str(electron_distributions.getX(i) / gPars::medium_props.atomic_density / Td, 3)
        + " Td";
    gnuplotCommands << "plot $Mydata using 1:2 w l lw 2 lc rgb 'black' title \"" << title <<"\""<<std::endl;
    fputs(gnuplotCommands.str().c_str(), gnuplotPipe);
    fflush(gnuplotPipe);
  }

  gnuplotCommands.str(std::string());
  delete_file(gPars::medium_props.cache_folder+"electron_distributions_dynamic.gif");
  gnuplotCommands
      << "set terminal gif size 800,800 animate delay 100 enhanced"<<std::endl
      <<"set output '"+gPars::medium_props.cache_folder+"electron_distributions_dynamic.gif'"<<std::endl
      <<"set xrange [0:*]"<<std::endl
      <<"set yrange [ 1e-7 < * :*]"<<std::endl<<std::endl;
  fputs(gnuplotCommands.str().c_str(), gnuplotPipe);
  fflush(gnuplotPipe);
  for(std::size_t i = 0, i_end_ = electron_distributions.size(); i != i_end_; ++i) {
    gnuplotCommands.str(std::string());
    gnuplotCommands
      << "undefine $Mydata"<<std::endl
      << "$Mydata << EOD"<<std::endl;
    for(std::size_t j = 0, j_end_ = electron_distributions.getY_data(i).size(); j!=j_end_; ++j) {
      double X = electron_distributions.getY_data(i).getX(j);
      double Y = electron_distributions.getY_data(i).getY(j);
      gnuplotCommands // Passing data through pipe together with commands
      << X << "\t" << Y << std::endl;
    }
    gnuplotCommands << "EOD"<<std::endl;
    std::string title = "E = " + dbl_to_str(electron_distributions.getX(i) * cm / kilovolt, 3)
        + " kV/cm, E/N = " + dbl_to_str(electron_distributions.getX(i) / gPars::medium_props.atomic_density / Td, 3)
        + " Td";
    gnuplotCommands << "plot $Mydata using 1:2 w l lw 2 lc rgb 'black' title \"" << title <<"\""<<std::endl;
    fputs(gnuplotCommands.str().c_str(), gnuplotPipe);
    fflush(gnuplotPipe);
  }
  pclose(gnuplotPipe);
}

void MediumPropertiesTables::PlotDriftVelocity(void) const
{
  FILE* gnuplotPipe = popen(gPars::general.gnuplot_bin.c_str(), "w");
  if (nullptr == gnuplotPipe) {
    if (verbosity > 0) {
      std::cerr<<classname<<"::PlotDriftVelocity:Error: could not open gnuplot."<<std::endl;
      std::cerr<<"\tNot plotting."<<std::endl;
    }
    return;
  }
  std::stringstream gnuplotCommands;
  gnuplotCommands
      <<"set terminal png size 1000,800 enhanced"<<std::endl
      <<"set output '"+gPars::medium_props.cache_folder+"electron_drift_velocity.png'"<<std::endl
      <<"set grid"<<std::endl
      <<"set xrange["<<interval_field_drift.min() / gPars::medium_props.atomic_density / Td
      <<":"<<interval_field_drift.max() / gPars::medium_props.atomic_density / Td << "]"<<std::endl
      <<"set xlabel \"Field [Td]\""<<std::endl
      <<"set logscale x"<<std::endl
      <<"set yrange [0:1e6]"<<std::endl
      <<"set ylabel \"Drift velocity [cm/s]"<<std::endl <<std::endl;
  fputs(gnuplotCommands.str().c_str(), gnuplotPipe);
  fflush(gnuplotPipe);
  gnuplotCommands.str(std::string());
  gnuplotCommands
    << "undefine $Mydata"<<std::endl
    << "$Mydata << EOD"<<std::endl;
  for(std::size_t j = 0, j_end_ = drift_velocity.size(); j!=j_end_; ++j) {
    double X = drift_velocity.getX(j) / gPars::medium_props.atomic_density / Td; // Geant4 field to Townsends
    double Y = drift_velocity.getY(j) * s / cm;
    gnuplotCommands // Passing data through pipe together with commands
    << X << "\t" << Y << std::endl;
  }
  gnuplotCommands << "EOD"<<std::endl;
  std::string title = gPars::MediumName() + " drift velocity";
  gnuplotCommands << "plot $Mydata using 1:2 w l lw 2 lc rgb 'black' title '" << title <<"'"<<std::endl;
  fputs(gnuplotCommands.str().c_str(), gnuplotPipe);
  fflush(gnuplotPipe);
  pclose(gnuplotPipe);
}

void MediumPropertiesTables::PlotFDistributions(void) const
{
  FILE* gnuplotPipe = popen(gPars::general.gnuplot_bin.c_str(), "w");
  if (nullptr == gnuplotPipe) {
    if (verbosity > 0) {
      std::cerr<<classname<<"::PlotFDistributions:Error: could not open gnuplot."<<std::endl;
      std::cerr<<"\tNot plotting."<<std::endl;
    }
    return;
  }
  auto En_range = F_distributions.getRangeAtX(low_high_threshold_field);
	if (!En_range) {
		return;
	}
  std::stringstream gnuplotCommands;
  delete_file(gPars::medium_props.cache_folder+"electron_F_distributions_low_E.gif");
  gnuplotCommands
      <<"set terminal gif size 800,800 animate delay 100 enhanced"<<std::endl
      <<"set output '"+gPars::medium_props.cache_folder+"electron_F_distributions_low_E.gif'"<<std::endl
      <<"set grid"<<std::endl
      <<"set xrange[0:" << En_range->second << "]"<<std::endl
      <<"set xlabel \"Energy [eV]\""<<std::endl
      <<"set logscale y"<<std::endl
      <<"set yrange [1e-4:1e3]"<<std::endl
      <<"set ylabel \"F\\' [eV^{-1/2}]\""<<std::endl <<std::endl;
  fputs(gnuplotCommands.str().c_str(), gnuplotPipe);
  fflush(gnuplotPipe);
  for(std::size_t i = 0, i_end_ = F_distributions.size(); i != i_end_; ++i) {
    if (F_distributions.getX(i) > interval_field_NBrS.min())
      continue;
    gnuplotCommands.str(std::string());
    gnuplotCommands
      << "undefine $Mydata"<<std::endl
      << "$Mydata << EOD"<<std::endl;
    for(std::size_t j = 0, j_end_ = F_distributions.getY_data(i).size(); j!=j_end_; ++j) {
      gnuplotCommands // Passing data through pipe together with commands
      << F_distributions.getY_data(i).getX(j)
      << "\t" << F_distributions.getY_data(i).getY(j)
      << std::endl;
    }
    gnuplotCommands << "EOD"<<std::endl;
    std::string title = "E = " + dbl_to_str(F_distributions.getX(i) * cm / kilovolt, 3)
        + " kV/cm, E/N = " + dbl_to_str(F_distributions.getX(i) / gPars::medium_props.atomic_density / Td, 3)
        + " Td";
    gnuplotCommands << "plot $Mydata using 1:2 w l lw 2 lc rgb 'black' title \"" << title <<"\""<<std::endl;
    fputs(gnuplotCommands.str().c_str(), gnuplotPipe);
    fflush(gnuplotPipe);
  }
  auto data = F_distributions.getY_data(F_distributions.size()-1);
	double En_max = data.getX(data.size()-1);
  gnuplotCommands.str(std::string());
  delete_file(gPars::medium_props.cache_folder+"electron_F_distributions_high_E.gif");
  gnuplotCommands
      <<"set terminal gif size 800,800 animate delay 100 enhanced"<<std::endl
      <<"set output '"+gPars::medium_props.cache_folder+"electron_F_distributions_high_E.gif'"<<std::endl
      <<"set xrange[0:" << En_max << "]"<<std::endl
      <<"set yrange [1e-4:10]"<<std::endl
      <<"set ylabel \"F\\' [eV^{-1/2}]\""<<std::endl;
  fputs(gnuplotCommands.str().c_str(), gnuplotPipe);
  fflush(gnuplotPipe);
  for(std::size_t i = 0, i_end_ = F_distributions.size(); i != i_end_; ++i) {
    if (F_distributions.getX(i) < interval_field_NBrS.min())
      continue;
    gnuplotCommands.str(std::string());
    gnuplotCommands
      << "undefine $Mydata"<<std::endl
      << "$Mydata << EOD"<<std::endl;
    for(std::size_t j = 0, j_end_ = F_distributions.getY_data(i).size(); j!=j_end_; ++j) {
      double X = F_distributions.getY_data(i).getX(j);
      double Y = F_distributions.getY_data(i).getY(j);
      gnuplotCommands // Passing data through pipe together with commands
      << X <<"\t" << Y << std::endl;
    }
    gnuplotCommands << "EOD"<<std::endl;
    std::string title = "E = " + dbl_to_str(F_distributions.getX(i) * cm / kilovolt, 3)
        + " kV/cm, E/N = " + dbl_to_str(F_distributions.getX(i) / gPars::medium_props.atomic_density / Td, 3)
        + " Td";
    gnuplotCommands << "plot $Mydata using 1:2 w l lw 2 lc rgb 'black' title \"" << title <<"\""<<std::endl;
    fputs(gnuplotCommands.str().c_str(), gnuplotPipe);
    fflush(gnuplotPipe);
  }

  gnuplotCommands.str(std::string());
  delete_file(gPars::medium_props.cache_folder+"electron_F_distributions_dynamic.gif");
  gnuplotCommands
      <<"set terminal gif size 800,800 animate delay 100 enhanced"<<std::endl
      <<"set output '"+gPars::medium_props.cache_folder+"electron_F_distributions_dynamic.gif'"<<std::endl
      <<"set xrange [*:*]"<<std::endl
      <<"set logscale x"<<std::endl
      <<"set yrange [ 1e-7< * :*]"<<std::endl<<std::endl;
  fputs(gnuplotCommands.str().c_str(), gnuplotPipe);
  fflush(gnuplotPipe);
  for(std::size_t i = 0, i_end_ = F_distributions.size(); i != i_end_; ++i) {
    gnuplotCommands.str(std::string());
    gnuplotCommands
      << "undefine $Mydata"<<std::endl
      << "$Mydata << EOD"<<std::endl;
    for(std::size_t j = 0, j_end_ = F_distributions.getY_data(i).size(); j!=j_end_; ++j) {
      double X = F_distributions.getY_data(i).getX(j);
      double Y = F_distributions.getY_data(i).getY(j);
      gnuplotCommands // Passing data through pipe together with commands
      << X <<"\t" << Y << std::endl;
    }
    gnuplotCommands << "EOD"<<std::endl;
    std::string title = "E = " + dbl_to_str(F_distributions.getX(i) * cm / kilovolt, 3)
        + " kV/cm, E/N = " + dbl_to_str(F_distributions.getX(i) / gPars::medium_props.atomic_density / Td, 3)
        + " Td";
    gnuplotCommands << "plot $Mydata using 1:2 w l lw 2 lc rgb 'black' title \"" << title <<"\""<<std::endl;
    fputs(gnuplotCommands.str().c_str(), gnuplotPipe);
    fflush(gnuplotPipe);
  }
  pclose(gnuplotPipe);
}

void MediumPropertiesTables::PlotYield(void) const
{
  FILE* gnuplotPipe = popen(gPars::general.gnuplot_bin.c_str(), "w");
  if (nullptr == gnuplotPipe) {
    if (verbosity > 0) {
      std::cerr<<classname<<"::PlotYield:Error: could not open gnuplot."<<std::endl;
      std::cerr<<"\tNot plotting."<<std::endl;
    }
    return;
  }
  std::stringstream gnuplotCommands;
  gnuplotCommands
      <<"set terminal png size 1200,800 enhanced"<<std::endl
      <<"set output '"+gPars::medium_props.cache_folder+"NBrS_yield.png'"<<std::endl
      <<"set grid"<<std::endl
      <<"set xrange[0:"<<interval_field_NBrS.max() / kilovolt * cm << "]"<<std::endl
      <<"set xlabel \"Electric field [kV/cm]\""<<std::endl
      <<"set yrange [0:*]"<<std::endl
      <<"set ylabel \"Absolute yield [photons/cm]"<<std::endl <<std::endl;
  fputs(gnuplotCommands.str().c_str(), gnuplotPipe);
  fflush(gnuplotPipe);
  gnuplotCommands.str(std::string());
  gnuplotCommands
    << "undefine $Mydata"<<std::endl
    << "$Mydata << EOD"<<std::endl;
  for(std::size_t j = 0, j_end_ = yield.size(); j!=j_end_; ++j) {
    gnuplotCommands // Passing data through pipe together with commands
    << yield.getX(j) / kilovolt * cm // Geant4 field to kV/cm
    << "\t" << yield.getY(j) * cm
    << std::endl;
  }
  gnuplotCommands << "EOD"<<std::endl;
  std::string title = "NBrS yield in " + gPars::MediumName();
  gnuplotCommands << "plot $Mydata using 1:2 w l lw 2 lc rgb 'black' title '" << title <<"'"<<std::endl;
  fputs(gnuplotCommands.str().c_str(), gnuplotPipe);
  fflush(gnuplotPipe);
  pclose(gnuplotPipe);
}

void MediumPropertiesTables::PlotSpectra(void) const
{
  FILE* gnuplotPipe = popen(gPars::general.gnuplot_bin.c_str(), "w");
  if (nullptr == gnuplotPipe) {
    if (verbosity > 0) {
      std::cerr<<classname<<"::PlotElectronDistributions:Error: could not open gnuplot."<<std::endl;
      std::cerr<<"\tNot plotting."<<std::endl;
    }
    return;
  }
  std::stringstream gnuplotCommands;
  gnuplotCommands
      << "set terminal gif size 1200,800 animate delay 100 enhanced"<<std::endl
      <<"set output '"+gPars::medium_props.cache_folder+"NBrS_spectra.gif'"<<std::endl
      <<"set grid"<<std::endl
      <<"set xrange[0:"<<gPars::medium_props.maximum_lambda / nm <<"]"<<std::endl
      <<"set xlabel \"Wavelength [nm]\""<<std::endl
      <<"set logscale y"<<std::endl
      <<"set yrange [1e-8:5e-2]"<<std::endl
      <<"set ylabel \"Normalized to 1 spectrum [nm^{-1}]"<<std::endl <<std::endl;
  fputs(gnuplotCommands.str().c_str(), gnuplotPipe);
  fflush(gnuplotPipe);
  double rel_tol = gPars::medium_props.yield_relative_tolerance;
  for(std::size_t i = 0, i_end_ = spectra.size(); i != i_end_; ++i) {
    gnuplotCommands.str(std::string());
    gnuplotCommands
      << "undefine $Mydata"<<std::endl
      << "$Mydata << EOD"<<std::endl;
    const DataVector & spectrum = spectra.getY_data(i);
    for(std::size_t j = 0, j_end_ = spectra.getY_data(i).size(); j!=j_end_; ++j) {
      // Spectra are stored as CDF (cumulative distribution function) so it needs to be differentiated.
      // Then spectra are plotted as a function of wavelength, not energy.
    	double E = spectrum.getX(j);
    	double dE = E * rel_tol;
      double pdf = j == 0 ? (spectrum(E+dE) - spectrum(E))/dE :
      		(j == (j_end_ - 1) ? (spectrum(E) - spectrum(E - dE))/dE :
      		(spectrum(E+dE) - spectrum(E-dE))/(2*dE));
      double lambda = hc / E;
      double spec = pdf * hc / lambda / lambda; // pdf transformation when variable is changed: E = h mu = h c / lambda, |dE/dl| = hc / l / l;
      gnuplotCommands // Passing data through pipe together with commands
      << lambda / nm // Geant4 unit to nm
      << "\t" << spec * nm << std::endl;
    }
    gnuplotCommands << "EOD"<<std::endl;
    std::string title = "E = " + dbl_to_str(spectra.getX(i) * cm / kilovolt, 3)
        + " kV/cm, E/N = " + dbl_to_str(spectra.getX(i) / gPars::medium_props.atomic_density / Td, 3)
        + " Td";
    gnuplotCommands << "plot $Mydata using 1:2 w l lw 2 lc rgb 'black' title '" << title <<"'"<<std::endl;
    fputs(gnuplotCommands.str().c_str(), gnuplotPipe);
    fflush(gnuplotPipe);
  }
  pclose(gnuplotPipe);
}

void MediumPropertiesTables::PlotSpectraRaw(void) const
{
	FILE* gnuplotPipe = popen(gPars::general.gnuplot_bin.c_str(), "w");
	if (nullptr == gnuplotPipe) {
		if (verbosity > 0) {
			std::cerr<<classname<<"::PlotElectronDistributions:Error: could not open gnuplot."<<std::endl;
			std::cerr<<"\tNot plotting."<<std::endl;
		}
		return;
	}
	std::stringstream gnuplotCommands;
	auto En_range = electron_distributions.getRangeAtX(interval_field_NBrS.max());
	if (!En_range) {
		return;
	}
	gnuplotCommands
			<< "set terminal gif size 1200,800 animate delay 100 enhanced"<<std::endl
			<<"set output '"+gPars::medium_props.cache_folder+"NBrS_spectra_raw.gif'"<<std::endl
			<<"set grid"<<std::endl
			<<"set logscale x"<<std::endl
			<<"set xrange[1:"<<8.0<<"]"<<std::endl
			<<"set xlabel \"Energy [eV]\""<<std::endl
			<<"set logscale y"<<std::endl
			<<"set yrange [1e-2:1.2]"<<std::endl
			<<"set ylabel \"CDF [unitless]"<<std::endl <<std::endl;
	fputs(gnuplotCommands.str().c_str(), gnuplotPipe);
	fflush(gnuplotPipe);
	for(std::size_t i = 0, i_end_ = spectra.size(); i != i_end_; ++i) {
		gnuplotCommands.str(std::string());
		gnuplotCommands
			<< "undefine $Mydata"<<std::endl
			<< "$Mydata << EOD"<<std::endl;
		for(std::size_t j = 0, j_end_ = spectra.getY_data(i).size(); (j+1)<j_end_; ++j) {
			double E = spectra.getY_data(i).getX(j) / electronvolt;
			double cdf = spectra.getY_data(i).getY(j);
			gnuplotCommands // Passing data through pipe together with commands
			<< E << "\t" << cdf << std::endl;
		}
		gnuplotCommands << "EOD"<<std::endl;
		std::string title = "E = " + dbl_to_str(spectra.getX(i) * cm / kilovolt, 3)
				+ " kV/cm, E/N = " + dbl_to_str(spectra.getX(i) / gPars::medium_props.atomic_density / Td, 3)
				+ " Td";
		gnuplotCommands << "plot $Mydata using 1:2 w l lw 2 lc rgb 'black' title '" << title <<"'"<<std::endl;
		fputs(gnuplotCommands.str().c_str(), gnuplotPipe);
		fflush(gnuplotPipe);
	}
	pclose(gnuplotPipe);
}

void MediumPropertiesTables::PlotDiffusions(void) const
{
  FILE* gnuplotPipe = popen(gPars::general.gnuplot_bin.c_str(), "w");
  if (nullptr == gnuplotPipe) {
    if (verbosity > 0) {
      std::cerr<<classname<<"::PlotDiffusions:Error: could not open gnuplot."<<std::endl;
      std::cerr<<"\tNot plotting."<<std::endl;
    }
    return;
  }
  std::stringstream gnuplotCommands;
  gnuplotCommands
      <<"set terminal png size 1200,800 enhanced"<<std::endl
      <<"set output '"+gPars::medium_props.cache_folder+"electron_diffusion.png'"<<std::endl
      <<"set grid"<<std::endl
      <<"set xrange["<<interval_field_drift.min() / gPars::medium_props.atomic_density / Td
      <<":"<<interval_field_drift.max() / gPars::medium_props.atomic_density / Td << "]"<<std::endl
      <<"set logscale x"<<std::endl
      <<"set xlabel \"Field (Td)\""<<std::endl
      <<"set yrange [0:60]"<<std::endl
      <<"set ylabel \"Diffusion coefficient (cm^{2}/s)"<<std::endl <<std::endl;
  fputs(gnuplotCommands.str().c_str(), gnuplotPipe);
  fflush(gnuplotPipe);
  gnuplotCommands.str(std::string());
  gnuplotCommands
    << "undefine $Mydata1"<<std::endl
    << "$Mydata1 << EOD"<<std::endl;
  for(std::size_t j = 0, j_end_ = drift_diffusion_T.size(); j!=j_end_; ++j) {
    gnuplotCommands // Passing data through pipe together with commands
    << drift_diffusion_T.getX(j) / gPars::medium_props.atomic_density / Td // Geant4 field to Townsends
    << "\t" << drift_diffusion_T.getY(j) * s / cm / cm
    << ""<<std::endl;
  }
  gnuplotCommands << "EOD"<<std::endl;
  fputs(gnuplotCommands.str().c_str(), gnuplotPipe);
  fflush(gnuplotPipe);
  gnuplotCommands.str(std::string());

  gnuplotCommands
    << "undefine $Mydata2"<<std::endl
    << "$Mydata2 << EOD"<<std::endl;
  for(std::size_t j = 0, j_end_ = drift_diffusion_L.size(); j!=j_end_; ++j) {
    gnuplotCommands // Passing data through pipe together with commands
    << drift_diffusion_L.getX(j) / gPars::medium_props.atomic_density / Td // Geant4 field to Townsends
    << "\t" << drift_diffusion_L.getY(j) * s / cm / cm
    << ""<<std::endl;
  }
  gnuplotCommands << "EOD"<<std::endl;
  fputs(gnuplotCommands.str().c_str(), gnuplotPipe);
  fflush(gnuplotPipe);
  gnuplotCommands.str(std::string());

  std::string title1 = gPars::MediumName() + " transversal diffusion";
  std::string title2 = gPars::MediumName() + " longitudinal diffusion";
  gnuplotCommands << "plot $Mydata1 using 1:2 w l lw 2 lt 1 lc rgb 'black' title '" << title1 <<"',";
  gnuplotCommands << " $Mydata2 using 1:2 w l lw 2 lt 4 lc rgb 'medium-blue' title '" << title2 <<"',";
  fputs(gnuplotCommands.str().c_str(), gnuplotPipe);
  fflush(gnuplotPipe);
  pclose(gnuplotPipe);
}
