#include <chrono>
#include <thread>

#include "MediumPropertiesTables.hh"

const std::string MediumPropertiesTables::filename_drift_velocity = "drift_velocity.dat";
const std::string MediumPropertiesTables::filename_drift_diffusion_L = "diffusion_longitudinal.dat";
const std::string MediumPropertiesTables::filename_drift_diffusion_T = "diffusion_transverse.dat";
const std::string MediumPropertiesTables::filename_yield = "NBrS_yield.dat";
const std::string MediumPropertiesTables::filename_spectra = "NBrS_spectra.dat";
const std::string MediumPropertiesTables::filename_electron_distributions = "electron_energy_distributions.dat";
const std::string MediumPropertiesTables::filename_F_distributions = "F_distributions.dat";
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

void MediumPropertiesTables::CalcElectronDistributions(void)
{
  if (verbosity > 0)
    std::cout<<"Calculating electron energy distributions..."<<std::endl;
  for (std::size_t i = 0, i_end_ = interval_field_drift.NumOfIndices(); i!=i_end_; ++i) {
    double field = interval_field_drift.Value(i);
    DataVector e_distr = CalcElectronDistribution(field);
    electron_distributions.push(field, e_distr);
  }
  if (verbosity == 2) {
    std::cout<<"Plotting electron energy distributions..."<<std::endl;
    PlotElectronDistributions();
  }
}

DataVector MediumPropertiesTables::CalcElectronDistribution(double field) const // e distribution by their energy, f'.
{
  // Calculate eq. 10 in Borisova2021 (doi:doi:10.1209/0295-5075/ac4c03)
  DataVector out;
  out.set_out_value(0.0); // Distribution returns 0 at energies outside specified energy range.
  out.use_leftmost(false); out.use_rightmost(false);
  out.setOrder(1); out.setNused(2); // Linear interpolation (energy points are quite dense).
  DataVector f_tmp = out; // To calculate normalization.
  const double delta = 2.0 * gPars::medium_props.m_to_M;
  const double eV_to_vel = 2.0 / eV * e_SI / e_mass_SI;
  const double atomic_density = gPars::medium_props.atomic_density * m*m*m;
  const double field_SI = field / volt * m;
  const double eq_10_coeff = 1.5 * e_mass_SI / e_SI / e_SI / field_SI / field_SI;

  IntegrationRange e_distr_interval = GetIntervalElectronDistributions(field);
  IntegrationRange E_interval = interval_XS + e_distr_interval;
  E_interval.Trim(0.0, e_distr_interval.max());

  double integral = 0.0;
  double En_prev = E_interval.Value(0) / joule;
  double Y_prev = 0; // Y == value under integral, not distribution value
  for (std::size_t i = 0, i_end_ = E_interval.NumOfIndices(); i!=i_end_; ++i) {
    double En = E_interval.Value(i); // in Geant4 units
    double nu_e = delta * atomic_density * (XS_energy_transfer(En) / (m*m)) * sqrt(En * eV_to_vel); // eq. 7 in Borisova2021. In SI.
    double nu_m = atomic_density * (XS_momentum_transfer(En) / (m*m)) * sqrt(En * eV_to_vel);// eq. 8 in Borisova2021. In SI.
    double Y = eq_10_coeff * nu_e * nu_m;
    double En_SI = En / joule; // Must be in SI in integral.
    integral += 0.5 * (Y + Y_prev) * (En_SI - En_prev);
    Y_prev = Y;
    En_prev = En_SI;
    if (isnan(integral) || integral < 0 || isinf(integral)) {
      if (pedantic)
        G4Exception((classname +"::CalcElectronDistribution: ").c_str(),
          "InvalidValue", FatalException, "Invalid value obtained (negative, nan, infinity or DBL_MAX).");
      else
        integral = 1e60;
    }
    out.push_back(En_SI, exp(-integral));
    f_tmp.push_back(En_SI, exp(-integral)*sqrt(En_SI));
  }
  double normalization = f_tmp.integrate(e_distr_interval / joule);
  out.scaleY(1.0 / normalization);
  double max_f = out.getY(0);
  //for (std::size_t i = 0, i_end_ = out.size(); i!=i_end_; ++i)
  //  if (out.getY(i) < max_f * min_e_rel_probability)
  //    out[i].second = 0;
  return out;
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

double MediumPropertiesTables::CalcDriftVelocity(double field) const
{ //eq. 11 in Borisova2021 (doi:doi:10.1209/0295-5075/ac4c03)
  const double delta = 2.0 * gPars::medium_props.m_to_M;
  const double eV_to_vel = 2.0 / eV * e_SI / e_mass_SI;
  const double sqrt_vel = sqrt(joule * eV_to_vel);
  const double atomic_density = gPars::medium_props.atomic_density * m*m*m;
  const double field_SI = field / volt * m;
  const double eq_11_coeff = 2.0 / 3.0 * e_SI * field_SI / e_mass_SI;
  IntegrationRange range = GetIntervalElectronDistributions(field) + interval_XS;
  range.Trim(0.0, GetIntervalElectronDistributions(field).max());

  double integral = 0.0;
  double f_prev = electron_distributions(field, range.Value(0) / joule); // Distribution energy is in SI, interval_XS is in Gean4 units
  double Y_prev = 0; // Y == value under integral
  for (std::size_t i = 0, i_end_ = range.NumOfIndices(); i!=i_end_; ++i) {
    double En = range.Value(i); // in Geant4 units
    double nu_m_sqrtE = atomic_density * (XS_momentum_transfer(En) / (m*m)) * sqrt_vel;// eq. 8 in Borisova2021 / sqrt(E). In SI.
    //nu_m is divided by sqrt(E) to avoid division by 0 below! Take note of units.
    double En_SI = En / joule; // Must be in SI in integral.
    double Y = En_SI / nu_m_sqrtE;
    double f = electron_distributions(field, En_SI);
    integral += -0.5 * (Y + Y_prev) * (f - f_prev);
    Y_prev = Y;
    f_prev = f;
  }
  double val = integral * eq_11_coeff * m / s; // return in Geant4 units
  if (isnan(val) || isinf(val) || val < 0) {
    if (pedantic)
      G4Exception((classname +"::CalcDriftVelocity: ").c_str(),
        "InvalidValue", FatalException, "Invalid value obtained (negaive, nan, infinity or DBL_MAX).");
    else
      val = 0;
  }
  return val;
}

void MediumPropertiesTables::CalcFDistributions(void)
{
  if (verbosity > 0)
    std::cout<<"Calculating F distributions..."<<std::endl;
  for (std::size_t i = 0, i_end_ = interval_field_drift.NumOfIndices(); i!=i_end_; ++i) {
    double field = interval_field_drift.Value(i);
    DataVector F_distr = CalcFDistribution(field);
    F_distributions.push(field, F_distr);
  }
  if (verbosity == 2) {
    std::cout<<"Plotting F distributions..."<<std::endl;
    PlotFDistributions();
  }
}

DataVector MediumPropertiesTables::CalcFDistribution(double field) const
{ // Calculate eq. 5 in Atrazhev1985 (doi:10.1088/0022-3719/18/6/015)
  DataVector out;
  out.set_out_value(0.0); // Distribution returns 0 at energies outside specified energy range.
  out.use_leftmost(false); out.use_rightmost(false);
  out.setOrder(1); out.setNused(2); // Linear interpolation (energy points are quite dense).
  IntegrationRange range = GetIntervalElectronDistributions(field) + interval_XS + GetIntervalFDistributions(field);
  range.Trim(0.0, GetIntervalFDistributions(field).max());

  // It is always good to reduce number of multiplications inside loops. So constant factors are cached.
  const double delta = 2.0 * gPars::medium_props.m_to_M;
  const double eV_to_vel = 2.0 / eV * e_SI / e_mass_SI;
  const double sqrt_vel = sqrt(joule * eV_to_vel);
  const double atomic_density = gPars::medium_props.atomic_density * m*m*m;
  const double field_SI = field / volt * m;
  const double eq_5_coeff = 1.5 * e_mass_SI / e_SI / field_SI * (drift_velocity(field) * s / m);

  double integral = 0.0;
  double En_prev = range.Value(0) / joule;
  double V_prev = 0; // V == value under integral
  for (std::size_t i = 0, i_end_ = range.NumOfIndices(); i!=i_end_; ++i) { // integral over dx
    double En = range.Value(i); // in Geant4 units
    double En_SI = En / joule; // Must be in SI in integral.
    double nu_m_e = atomic_density * (XS_momentum_transfer(En) / (m*m)) * sqrt(En * eV_to_vel); // eq. 8 in Borisova2021. In SI.
    double f_e = electron_distributions(field, En_SI);

    if (0 >= En || f_e < min_f_value_probability) {
      integral += 0.0;
      out.insert(En_SI, 0.0);
      continue;
    }

    IntegrationRange interval_y = range;
    interval_y.Trim(0.0, En);
    double integral1 = 0; // 2nd term inside integral over dx of eq. 5 in Atrazhev1985
    double integral2 = 0; // 3rd term inside integral over dx of eq. 5 in Atrazhev1985
    double y_prev = interval_y.Value(0) / joule;
    double f_prev = electron_distributions(field, y_prev);
    double V1_prev = 0; // value under integral1
    double V2_prev = 0; // value under integral2
    for (std::size_t j = 0, j_end_ = interval_y.NumOfIndices(); j!=j_end_; ++j) { // Take both integrals over dy
      double y = interval_y.Value(j);
      double y_SI = y / joule; // Must be in SI in integral.
      double nu_m_sqrty = atomic_density * (XS_momentum_transfer(y) / (m*m)) * sqrt_vel; // eq. 8 in Borisova2021 / sqrt(y_SI). In SI.
      //nu_m is divided by sqrt(y) to avoid division by 0! Take note of units.
      double f = electron_distributions(field, y_SI);
      double V1 = y_SI / nu_m_sqrty;
      double V2 = sqrt(y_SI) * f;
      integral1 += 0.5 * (V1 + V1_prev) * (f - f_prev); // df/de * de == df
      integral2 += 0.5 * (V2 + V2_prev) * (y_SI - y_prev);
      y_prev = y_SI;
      f_prev = f;
      V1_prev = V1;
      V2_prev = V2;
    }

    double V = 1.0  + nu_m_e / pow(En_SI, 1.5) / f_e * (integral1 + eq_5_coeff * integral2);
    integral += 0.5 * (V + V_prev) * (En_SI - En_prev);
    V_prev = V;
    En_prev = En_SI;
    double val =  f_e * integral;
    if (isnan(val) || isinf(val)) {
      if (pedantic)
        G4Exception((classname + "::CalcFDistribution: ").c_str(),
          "InvalidValue", FatalException, "Invalid value obtained (nan, infinity or DBL_MAX).");
      else
        val = 0;
    }
    out.insert(En_SI, val);
  }
  return out;
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
    double yeild_v = diff_XS.integrate(interval_photon_En);
    yield.push_back(field, yeild_v);
    diff_XS.integrate();
    yeild_v = diff_XS.getY(diff_XS.size() - 1);
    if (yeild_v > 0) {
      diff_XS.scaleY(1.0 / yeild_v);
      spectra.push(field, diff_XS);
    }
  }
  if (verbosity == 2) {
    std::cout<<"Plotting NBrS yield and spectra..."<<std::endl;
    PlotYield();
    PlotSpectra();
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

// Default: uses elastic XS for NBrS XS calculation
DataVector MediumPropertiesTables::CalcDiffXS_ElasticFormula(double field) const
{
  // Calculate eq. 6 but with d(lambda) replaced by d(h*nu) in Borisova2021 (doi:doi:10.1209/0295-5075/ac4c03)
  DataVector out;
  out.set_out_value(0.0); // Distribution returns 0 at energies outside specified energy range.
  out.use_leftmost(false); out.use_rightmost(false);
  out.setOrder(1); out.setNused(2); // Linear interpolation (energy points are quite dense).
  if (interval_photon_En.NumOfIndices() <= 0)
    return out;
  const double delta = 2.0 * gPars::medium_props.m_to_M;
  const double eV_to_vel = 2.0 / eV * e_SI / e_mass_SI;
  const double atomic_density = gPars::medium_props.atomic_density * m*m*m;
  const double field_SI = field / volt * m;
  const double eq_6_coeff = 8.0 / 3.0 * atomic_density * (classic_electr_radius/m) / (hc/(m*joule)) * sqrt(2.0 / e_mass_SI)
      / (drift_velocity(field) * s / m); // In SI

  for (std::size_t i = 0, i_end_ = interval_photon_En.NumOfIndices(); i!=i_end_; ++i) {
    double photon_En = interval_photon_En.Value(i);
    double photon_En_SI = photon_En / joule;
    double integral = 0.0;
    IntegrationRange range = GetIntervalElectronDistributions(field) + interval_XS;
    range.Trim(0.0, GetIntervalElectronDistributions(field).max());
    if (photon_En < range.max()) {
      range.Trim(photon_En, range.max());
      double En_prev = range.Value(0) / joule;
      double Y_prev = 0; // Y == value under integral
      for (std::size_t j = 0, j_end_ = range.NumOfIndices(); j!=j_end_; ++j) {
        double En = range.Value(j);
        double En_SI = En / joule;
        double XS = sqrt((En_SI - photon_En_SI)/En_SI) / photon_En_SI
            * ((En_SI - photon_En_SI) * XS_energy_transfer(En) / (m*m) + En_SI * XS_energy_transfer(En - photon_En) / (m*m));
        double Y = En_SI * electron_distributions(field, En_SI) * XS;
        integral += 0.5 * (Y + Y_prev) * (En_SI - En_prev);
        Y_prev = Y;
        En_prev = En_SI;
      }
    }
    double val = eq_6_coeff * integral / m / joule;  // store in Geant4 units
    if (isnan(val) || val < 0 || isinf(val)) {
      if (pedantic)
        G4Exception((classname + "::CalcDiffXS: ").c_str(),
          "InvalidValue", FatalException, "Invalid value obtained (negative, nan, infinity or DBL_MAX).");
      else
        val = 0;
    }
    out.push_back(photon_En, val);
  }
  return out;
}

// Uses momentum transfer XS for NBrS XS calculation (correct, but only for hv << E of electron https://doi.org/10.48550/arXiv.2206.01388)
DataVector MediumPropertiesTables::CalcDiffXS_TransferFormula(double field) const
{
  // Calculate eq. 6 but with d(lambda) replaced by d(h*nu) in Borisova2021 (doi:doi:10.1209/0295-5075/ac4c03)
  // but using momentum transfer XS as discussed in https://doi.org/10.48550/arXiv.2206.01388
  DataVector out;
  out.set_out_value(0.0); // Distribution returns 0 at energies outside specified energy range.
  out.use_leftmost(false); out.use_rightmost(false);
  out.setOrder(1); out.setNused(2); // Linear interpolation (energy points are quite dense).
  if (interval_photon_En.NumOfIndices() <= 0)
    return out;
  const double delta = 2.0 * gPars::medium_props.m_to_M;
  const double eV_to_vel = 2.0 / eV * e_SI / e_mass_SI;
  const double atomic_density = gPars::medium_props.atomic_density * m*m*m;
  const double field_SI = field / volt * m;
  const double eq_6_coeff = 8.0 / 3.0 * atomic_density * (classic_electr_radius/m) / (hc/(m*joule)) * sqrt(2.0 / e_mass_SI)
      / (drift_velocity(field) * s / m); // In SI

  for (std::size_t i = 0, i_end_ = interval_photon_En.NumOfIndices(); i!=i_end_; ++i) {
    double photon_En = interval_photon_En.Value(i);
    double photon_En_SI = photon_En / joule;
    double integral = 0.0;
    IntegrationRange range = GetIntervalElectronDistributions(field) + interval_XS;
    range.Trim(0.0, GetIntervalElectronDistributions(field).max());
    if (photon_En < range.max()) {
      range.Trim(photon_En, range.max());
      double En_prev = range.Value(0) / joule;
      double Y_prev = 0; // Y == value under integral
      for (std::size_t j = 0, j_end_ = range.NumOfIndices(); j!=j_end_; ++j) {
        double En = range.Value(j);
        double En_SI = En / joule;
        double XS = sqrt((En_SI - photon_En_SI)/En_SI) / photon_En_SI
            * ((En_SI - photon_En_SI) * XS_momentum_transfer(En) / (m*m) + En_SI * XS_momentum_transfer(En - photon_En) / (m*m));
        double Y = En_SI * electron_distributions(field, En_SI) * XS;
        integral += 0.5 * (Y + Y_prev) * (En_SI - En_prev);
        Y_prev = Y;
        En_prev = En_SI;
      }
    }
    double val = eq_6_coeff * integral / m / joule;  // store in Geant4 units
    if (isnan(val) || val < 0 || isinf(val)) {
      if (pedantic)
        G4Exception((classname + "::CalcDiffXS: ").c_str(),
          "InvalidValue", FatalException, "Invalid value obtained (negative, nan, infinity or DBL_MAX).");
      else
        val = 0;
    }
    out.push_back(photon_En, val);
  }
  return out;
}

DataVector MediumPropertiesTables::CalcDiffXS_ExactFormula(double field) const
{
  return DataVector();
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

double MediumPropertiesTables::CalcDiffusionT(double field) const
{
  //eq. 4 in Atrazhev1985 (doi:10.1088/0022-3719/18/6/015), rewritten without nu_m
  const double delta = 2.0 * gPars::medium_props.m_to_M;
  const double atomic_density = gPars::medium_props.atomic_density * m*m*m;
  const double field_SI = field / volt * m;
  const double eq_4_coeff = sqrt(2.0) / 3.0 / sqrt(e_mass_SI) / atomic_density;

  double integral = 0.0;
  IntegrationRange range = GetIntervalElectronDistributions(field) + interval_XS;
  range.Trim(0, GetIntervalElectronDistributions(field).max()); // To avoid singularity at 0
  double En_prev = range.Value(0) / joule; // Distribution energy is in SI, interval_XS is in Gean4 units
  double Y_prev = 0; // Y == value under integral
  for (std::size_t i = 0, i_end_ = range.NumOfIndices(); i!=i_end_; ++i) {
    double En = range.Value(i); // in Geant4 units
    double En_SI = En / joule; // Must be in SI in integral.
    double Y = En_SI / (XS_momentum_transfer(En) / (m*m)) * electron_distributions(field, En_SI);
    integral += 0.5 * (Y + Y_prev) * (En_SI - En_prev);
    En_prev = En_SI;
    Y_prev = Y;
  }
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

double MediumPropertiesTables::CalcDiffusionL(double field) const
{
  // Eq. 5 in Atrazhev1985 (doi:10.1088/0022-3719/18/6/015)
  // !!! Formula for longitudial diffusion there is incorrect!
  // Second integral term has incorrect units. The source from which
  // the formula was derived (Parker and Lowke 1969 doi:10.1103/PhysRev.181.290)
  // seems correct. Fixed formula is in data/Atrazhev1985_corrected_equation_5.png
  // Or in latex format:
  // D_{||} = D_{0} - \frac{2}{3m} \int_{0}^{\infty}\left( \frac{\epsilon^{3/2}}{\nu_{m}(\epsilon)} \frac{d F (\epsilon)}{d \epsilon} + \frac{3mW}{2e\mathscr{E}} \sqrt{\epsilon} F(\epsilon) \right) d \epsilon
  const double delta = 2.0 * gPars::medium_props.m_to_M;
  const double eV_to_vel = 2.0 / eV * e_SI / e_mass_SI;
  const double sqrt_vel = sqrt(joule * eV_to_vel);
  const double atomic_density = gPars::medium_props.atomic_density * m*m*m;
  const double field_SI = field / volt * m;
  const double eq_5_coeff = 2.0 / 3.0 / e_mass_SI;
  const double int2_coeff = 3.0 / 2.0 * (drift_velocity(field) * s / m) * e_mass_SI / e_SI / field_SI;

  IntegrationRange range = GetIntervalFDistributions(field) + interval_XS;
  range.Trim(0.0, GetIntervalFDistributions(field).max());
  double integral1 = 0.0;
  double integral2 = 0.0;
  double En_prev = range.Value(0) / joule; // Distribution energy is in SI, interval_XS is in Gean4 units
  double F_prev = F_distributions(field, En_prev);
  double V1_prev = 0; // V == value under integral
  double V2_prev = 0; // V == value under integral
  for (std::size_t i = 0, i_end_ = range.NumOfIndices(); i!=i_end_; ++i) {
    double En = range.Value(i); // in Geant4 units
    double nu_m_sqrtE = atomic_density * (XS_momentum_transfer(En) / (m*m)) * sqrt_vel;// eq. 8 in Borisova2021 / sqrt(E). In SI.
    //nu_m is divided by sqrt(E) in order to avoid division by 0 below! Take note of units.
    double En_SI = En / joule; // Must be in SI in integral.
    double F = F_distributions(field, En_SI);
    double V1 = En_SI / nu_m_sqrtE;
    double V2 = sqrt(En_SI) * F;
    integral1 += 0.5 * (V1 + V1_prev) * (F - F_prev);
    integral2 += 0.5 * (V2 + V2_prev) * (En_SI - En_prev);
    En_prev = En_SI;
    F_prev = F;
    V1_prev = V1;
    V2_prev = V2;
  }
  double val = drift_diffusion_T(field) - (integral1 + integral2 * int2_coeff) * eq_5_coeff * m * m / s ; // return in Geant4 units;
  if (isnan(val) || val < 0 || isinf(val)) {
    if (pedantic)
      G4Exception((classname + "::CalcDiffusionL: ").c_str(),
        "InvalidValue", FatalException, "Invalid value obtained (negative, nan, infinity or DBL_MAX).");
    else
      val = drift_diffusion_T(field);
  }
  return val;
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
      <<"set output '"+gPars::medium_props.cache_folder+"electron_Ar_XSs.png'"<<std::endl
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
  for(std::size_t j = 0, j_end_ = interval_XS.NumOfIndices(); j!=j_end_; ++j) {
    double X = interval_XS.Value(j) / eV;
    double Y = XS_energy_transfer(interval_XS.Value(j)) / cm / cm;
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
    double X = interval_XS.Value(j) / eV;
    double Y = XS_momentum_transfer(interval_XS.Value(j)) / cm / cm;
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
  gnuplotCommands
      << "set terminal gif size 800,800 animate delay 100 enhanced"<<std::endl
      <<"set output '"+gPars::medium_props.cache_folder+"electron_distributions_low_E.gif'"<<std::endl
      <<"set grid"<<std::endl
      <<"set xrange [0:"<< GetIntervalElectronDistributions(low_high_threshold_field).max() / eV<<"]"<<std::endl
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
      double X = electron_distributions.getY_data(i).getX(j) / e_SI; // Joule to eV
      double Y = electron_distributions.getY_data(i).getY(j) * std::pow(e_SI, 1.5);
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
  delete_file(gPars::medium_props.cache_folder+"electron_distributions_high_E.gif");
  gnuplotCommands
      << "set terminal gif size 800,800 animate delay 100 enhanced"<<std::endl
      <<"set output '"+gPars::medium_props.cache_folder+"electron_distributions_high_E.gif'"<<std::endl
      <<"set xrange [0:"<< GetIntervalElectronDistributions(max_field).max() / eV<<"]"<<std::endl
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
      double X = electron_distributions.getY_data(i).getX(j) / e_SI; // Joule to eV
      double Y = electron_distributions.getY_data(i).getY(j) * std::pow(e_SI, 1.5);
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
      double X = electron_distributions.getY_data(i).getX(j) / e_SI; // Joule to eV
      double Y = electron_distributions.getY_data(i).getY(j) * std::pow(e_SI, 1.5);
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
    double X = drift_velocity.getX(j) / gPars::medium_props.atomic_density / Td;
    double Y = drift_velocity.getY(j) * s / cm;
    gnuplotCommands // Passing data through pipe together with commands
    << drift_velocity.getX(j) / gPars::medium_props.atomic_density / Td // Geant4 field to Townsends
    << "\t" << drift_velocity.getY(j) * s / cm
    << std::endl;
  }
  gnuplotCommands << "EOD"<<std::endl;
  std::string title = "LAr drift velocity";
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
  std::stringstream gnuplotCommands;
  delete_file(gPars::medium_props.cache_folder+"electron_F_distributions_low_E.gif");
  gnuplotCommands
      <<"set terminal gif size 800,800 animate delay 100 enhanced"<<std::endl
      <<"set output '"+gPars::medium_props.cache_folder+"electron_F_distributions_low_E.gif'"<<std::endl
      <<"set grid"<<std::endl
      <<"set xrange[0:" << GetIntervalFDistributions(low_high_threshold_field).max() / eV << "]"<<std::endl
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
      << F_distributions.getY_data(i).getX(j) / e_SI // Joule to eV
      << "\t" << F_distributions.getY_data(i).getY(j) * std::pow(e_SI, 0.5)
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
  gnuplotCommands.str(std::string());
  delete_file(gPars::medium_props.cache_folder+"electron_F_distributions_high_E.gif");
  gnuplotCommands
      <<"set terminal gif size 800,800 animate delay 100 enhanced"<<std::endl
      <<"set output '"+gPars::medium_props.cache_folder+"electron_F_distributions_high_E.gif'"<<std::endl
      <<"set xrange[0:" << GetIntervalFDistributions(max_field).max() / eV << "]"<<std::endl
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
      double X = F_distributions.getY_data(i).getX(j) / e_SI; // Joule to eV
      double Y = F_distributions.getY_data(i).getY(j) * std::pow(e_SI, 0.5);
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
      double X = F_distributions.getY_data(i).getX(j) / e_SI; // Joule to eV
      double Y = F_distributions.getY_data(i).getY(j) * std::pow(e_SI, 0.5);
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
  std::string title = "NBrS yield in LAr";
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
      <<"set xrange[0:1000]"<<std::endl
      <<"set xlabel \"Wavelength [nm]\""<<std::endl
      <<"set logscale y"<<std::endl
      <<"set yrange [1e-4:5e-2]"<<std::endl
      <<"set ylabel \"Normalized to 1 spectrum [nm^{-1}]"<<std::endl <<std::endl;
  fputs(gnuplotCommands.str().c_str(), gnuplotPipe);
  fflush(gnuplotPipe);
  for(std::size_t i = 0, i_end_ = spectra.size(); i != i_end_; ++i) {
    gnuplotCommands.str(std::string());
    gnuplotCommands
      << "undefine $Mydata"<<std::endl
      << "$Mydata << EOD"<<std::endl;
    for(std::size_t j = 0, j_end_ = spectra.getY_data(i).size(); (j+1)<j_end_; ++j) {
      // Spectra are stored as CDF (cumulative distribution function) so it needs to be differentiated.
      // Then spectra are plotted as a function of wavelength, not energy.
      double pdf = (spectra.getY_data(i).getY(j + 1) - spectra.getY_data(i).getY(j))
          / (spectra.getY_data(i).getX(j + 1) - spectra.getY_data(i).getX(j));
      double lambda = 2.0 * hc / (spectra.getY_data(i).getX(j)+spectra.getY_data(i).getX(j + 1));
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
      <<"set xrange["<<min_field / gPars::medium_props.atomic_density / Td
      <<":"<<max_field / gPars::medium_props.atomic_density / Td << "]"<<std::endl
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

  std::string title1 = "LAr transversal diffusion";
  std::string title2 = "LAr longitudinal diffusion";
  gnuplotCommands << "plot $Mydata1 using 1:2 w l lw 2 lt 1 lc rgb 'black' title '" << title1 <<"',";
  gnuplotCommands << " $Mydata2 using 1:2 w l lw 2 lt 4 lc rgb 'medium-blue' title '" << title2 <<"',";
  fputs(gnuplotCommands.str().c_str(), gnuplotPipe);
  fflush(gnuplotPipe);
  pclose(gnuplotPipe);
}
