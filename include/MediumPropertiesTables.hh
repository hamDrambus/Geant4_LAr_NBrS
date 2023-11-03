#ifndef MediumPropertiesTables_h
#define MediumPropertiesTables_h

#include <string>

#include <G4SystemOfUnits.hh>
#include <G4PhysicalConstants.hh>

#include "utilities/FunctionTable.hh"
#include "utilities/GlobalUtilities.hh"
#include "utilities/IntegrationInterval.hh"
#include "utilities/PolynomialFit.hh"
#include "GlobalParameters.hh"

// Abstract class for medium.
//TODO: add multithreading to calculations. (CalcYeildAndSpectra may be quite long for high precision).
//TODO: optimize Calc...(field) methods. E.g. get const DataVector& distributions at specified field instead of calling distr(field, energy).
class MediumPropertiesTables {
  // All public properties are in Geant4 units.
public:
  // Properties below are calculated from XSs as in Atrazhev1985 (doi:10.1088/0022-3719/18/6/015)
  // and Borisova2021 (doi:10.1209/0295-5075/ac4c03). They are cached on disk and are not re-calculated
  // for each program execution.
	// The properties are in SI in Geant4 units.
  DataVector drift_velocity; // f(electric field)
  DataVector drift_diffusion_L; // f(electric field)
  DataVector drift_diffusion_T; // f(electric field)
  DataVector yield; // f(electric field). Total photon yield for 1 electron in given spectrum per unit length.
  FunctionTable spectra; // spectrum CDF as a function of electric field. Each CDF is a function of photon energy.

  MediumPropertiesTables() : classname("MediumPropertiesTables") {}
  virtual ~MediumPropertiesTables() = default;
  virtual void Initialize(void);
  virtual void InitializeFields(void) = 0; // TODO: move to settings
  bool IsReady(void) const;
  virtual double GenPhotonEnergy(double field, double rnd) const;
  inline void SetVerboseLevel(int verboseLevel)
  { verbosity = verboseLevel; }
protected:
  std::string classname;
  // XS = cross section, loaded from input files. In Geant4 units
  DataVector XS_energy_transfer; // f(electron energy),
  DataVector XS_momentum_transfer; // f(electron energy)
  FunctionTable XS_NBrS_diff_exact; // f(electron energy, photon energy) in Geant4 units [mm^2/eV]

  virtual bool LoadCached(void);
  virtual void CalcElectronDistributions(void);
  virtual void CalcFDistributions(void);
  virtual void CalcYeildAndSpectra(void);
  virtual void CalcDriftVelocity(void);
  virtual void CalcDiffusionT(void);
  virtual void CalcDiffusionL(void); // dependent on transverse diffusion.

  // Plotting functions have some hard-coded parameters for convenience.
  virtual void PlotInputXSs(void) const;
  virtual void PlotElectronDistributions(void) const;
  virtual void PlotDriftVelocity(void) const;
  virtual void PlotFDistributions(void) const;
  virtual void PlotYield(void) const;
  virtual void PlotSpectra(void) const;
  virtual void PlotSpectraRaw(void) const;
  virtual void PlotDiffusions(void) const;

  virtual std::pair<DataVector, DataVector> CalcElectronDistribution(double field) const; // e distribution by their energy, f'. f' is denoted as f in Atrazhev1985.
  //Normalized as integral{0, +inf} {e^(1/2)*f'(e) de} = 1. Usual distribution f(e) = e^(1/2)*f'(e).
  virtual std::pair<DataVector, DataVector> CalcFDistribution(double field) const; // eq. 5 in Atrazhev1985
  virtual DataVector CalcDiffXS(double field) const; // Spectrum normalized to yield. Integrate to get yield, normalize to unit area to get spectrum.
  
  /**
  * Calculates probability of NBrS photon production by hot electrons at given electric
  * field in the medium according to approximate theory (https://doi.org/10.1209/0295-5075/ac4c03)
  * using energy-transfer electron-atom cross section.
  * The approximation is correct only for hv << E of electron (https://doi.org/10.1016/j.nimb.2022.09.012).
  * Used by default if formula to use is not specifield in xml settings.
  *
  * @param field electric field in Geant4 units.
  * @returns Integrated cross section of photon production vs photon energy (spectrum)
  * which is equivalent to cumulative distribution function normalized to the total yeild
  * in number of photons per unit drift length. Photon energy is in eV in Geant4 units.
  */
  virtual DataVector CalcDiffXS_ElasticFormula(double field) const;
  
  /**
  * Calculates probability of NBrS photon production by hot electrons at given electric
  * field in the medium according to approximate theory (https://doi.org/10.1209/0295-5075/ac4c03)
  * using momentum-transfer electron-atom cross section.
  * The approximation is correct only for hv << E of electron (https://doi.org/10.1016/j.nimb.2022.09.012).
  *
  * @param field electric field in Geant4 units.
  * @returns Integrated cross section of photon production vs photon energy (spectrum)
  * which is equivalent to cumulative distribution function normalized to the total yeild
  * in number of photons per unit drift length. Photon energy is in eV in Geant4 units.
  */
  virtual DataVector CalcDiffXS_TransferFormula(double field) const;
  
  /**
  * Calculates probability of NBrS photon production by hot electrons at given electric
  * field in the medium according to approximate theory (https://doi.org/10.1209/0295-5075/ac4c03).
  * The approximation is correct only for hv << E of electron (https://doi.org/10.1016/j.nimb.2022.09.012).
  *
  * @param field electric field in Geant4 units.
  * @param XS tabulated electron-atom cross section vs electron energy whic is to be used
  * for approximation of NBrS cross section (see Eq. 6 in Borisova2021
  * https://doi.org/10.1209/0295-5075/ac4c03).
  * Energy is in eV and both it and cross section are in Geant4 units. 
  * @returns Integrated cross section of photon production vs photon energy (spectrum)
  * which is equivalent to cumulative distribution function normalized to the total yeild
  * in number of photons per unit drift length. Photon energy is in eV in Geant4 units.
  */
  virtual DataVector CalcDiffXS_XSFormula(double field, const DataVector & XS) const;
  
  /**
  * Calculates probability of NBrS photon production by hot electrons at given electric
  * field in the medium according to exact theory (https://doi.org/10.1016/j.nimb.2022.09.012).
  * Uses calculations conducted in Wolfram Mathematica v13.0 and tabulated in XS_NBrS_diff_exact.
  *
  * @param field electric field in Geant4 units.
  * @returns Integrated cross section of photon production vs photon energy (spectrum)
  * which is equivalent to cumulative distribution function normalized to the total yeild
  * in number of photons per unit drift length. Photon energy is in eV in Geant4 units.
  */
  virtual DataVector CalcDiffXS_ExactFormula(double field) const;

  virtual double CalcDriftVelocity(double field) const;
  virtual double CalcDiffusionT(double field) const;
  virtual double CalcDiffusionL(double field) const;

  // These distributions are deleted after initialization as they are required only for calculations of the public properties.
  // Field is in Geant4 unit, distributions themselves are in eVs (NOT in Geant4 eVs!) .
  FunctionTable electron_distributions; // Electron f' PDFs as a function of electric field. Each PDF is a function of photon energy. [PDF] = [energy^(-3/2)].
  FunctionTable electron_distributions_derivative; // Derivative over the energy of the above.
  FunctionTable F_distributions; // F from eq. 5 in Atrazhev1985 as a function of electric field. [F] = [energy^(-1/2)].
  FunctionTable F_distributions_derivative; // Derivative over the energy of the above.

  // Fields for NBrS and drift parameters are selected manually in individual daughter classes.
  // In case of another materials, input XSs, etc. these may need to be changed.
  IntegrationRange interval_field_NBrS; // only E > ~3kV/cm. Lower fields have no NBrS and quite different energy distributions.
  IntegrationRange interval_field_drift; // wide field range. Drift parameters are calculated for all fields.
  double low_high_threshold_field; // no NBrS below this, also distributions for display are separated by this value.

  int verbosity;
  bool pedantic; // If true then unexpected issues with calculations throw error. Otherwise they are somehow fixed but may lead to wrong results.

  static const std::string filename_drift_velocity;
  static const std::string filename_drift_diffusion_L;
  static const std::string filename_drift_diffusion_T;
  static const std::string filename_yield;
  static const std::string filename_spectra;
  static const std::string filename_electron_distributions;
  static const std::string filename_electron_distributions_derivative;
  static const std::string filename_F_distributions;
  static const std::string filename_F_distributions_derivative;
  static const std::string filename_XS_energy;
  static const std::string filename_XS_momentum;
};

#endif // MediumPropertiesTables_h
