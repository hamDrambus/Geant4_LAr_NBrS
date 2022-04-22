#ifndef ArgonPropertiesTables_h
#define ArgonPropertiesTables_h

#include <string>

#include <G4SystemOfUnits.hh>
#include <G4PhysicalConstants.hh>

#include "utilities/FunctionTable.hh"
#include "utilities/GlobalUtilities.hh"
#include "utilities/IntegrationInterval.hh"
#include "utilities/PolynomialFit.hh"
#include "GlobalParameters.hh"

//TODO: add multithreading to calculations. (CalcYeildAndSpectra is quite long for small step sizes).
class ArgonPropertiesTables {
  // All public properties are in Geant4 units.
public:
  // Properties below are calculated from XSs as in Atrazhev1985 (doi:10.1088/0022-3719/18/6/015)
  // and Borisova2021 (doi:10.1209/0295-5075/ac4c03). They are cached on disk and are not re-calculated
  // for each program execution.
  DataVector drift_velocity; // f(electric field)
  DataVector drift_diffusion_L; // f(electric field)
  DataVector drift_diffusion_T; // f(electric field)
  DataVector yield; // f(electric field). Total photon yield for 1 electron in given spectrum per unit length.
  FunctionTable spectra; // spectrum CDF as a function of electric field. Each CDF is a function of photon energy.

  ArgonPropertiesTables() : verbosity(0), pedantic(true) {}
  void Initialize(void);
  bool IsReady(void) const;
  double GenPhotonEnergy(double field, double rnd) const;
  inline void SetVerboseLevel(int verboseLevel)
  { verbosity = verboseLevel; }
protected:
  // XS = cross-section, loaded from input files. In Geant4 units
  DataVector XS_energy_transfer; // f(electron energy),
  DataVector XS_momentum_transfer; // f(electron energy)

  bool LoadCached(void);
  void CalcElectronDistributions(void);
  void CalcFDistributions(void);
  void CalcYeildAndSpectra(void);
  void CalcDriftVelocity(void);
  void CalcDiffusionT(void);
  void CalcDiffusionL(void); // dependent on transverse diffusion.

  void PlotInputXSs(void) const;
  void PlotElectronDistributions(void) const;
  void PlotDriftVelocity(void) const;
  void PlotFDistributions(void) const;
  void PlotYield(void) const;
  void PlotSpectra(void) const;
  void PlotDiffusions(void) const;

  DataVector CalcElectronDistribution(double field) const; // e distribution by their energy, f'. f' is denoted as f in Atrazhev1985.
  //Normalized as integral{0, +inf} {e^(1/2)*f'(e) de} = 1. Usual distribution f(e) = e^(1/2)*f'(e).
  DataVector CalcFDistribution(double field) const; // eq. 5 in Atrazhev1985
  DataVector CalcDiffXS(double field) const; // Spectrum normalized to yield. Integrate to get yield, normalize to unit area to get spectrum.
  double CalcDriftVelocity(double field) const;
  double CalcDiffusionT(double field) const;
  double CalcDiffusionL(double field) const;

  // These distributions are cleared after initialization if calculation of public parameters is required.
  // Field is in Geant4 unit, distributions themselves are in SI.
  FunctionTable electron_distributions; // Electron f' PDFs as a function of electric field. Each PDF is a function of photon energy. [PDF] = [energy^(-3/2)].
  FunctionTable F_distributions; // F from eq. 5 in Atrazhev1985 as a function of electric field. [F] = [energy^(-1/2)].

  // Integration intervals are selected manually after looking at integrated functions behavior.
  // In case another material or input XSs are used these may need to be changed.
  // TODO: intervals depending on electric field. Two regions are required at least E: <~ 4kV/cm and E > ~4kV/cm.
  IntegrationRange interval_XS; // for calculating electron distribution and F(e).
  IntegrationRange interval_photon_En;
  IntegrationRange interval_e_distr_low_E; // for calculating velocities, F(e) and diffusion.
  IntegrationRange interval_e_distr_high_E;
  IntegrationRange interval_F_low_E; // for calculating diffusion.
  IntegrationRange interval_F_high_E;
  IntegrationRange interval_field_NBrS; // only E > ~3kV/cm. Lower fields have no NBrS and quite different energy distributions.
  IntegrationRange interval_field_drift; // = E < ~3kV/cm + interval_field_NBrS. Drift parameters are calculated for all fields.
  double min_e_rel_probability;

  int verbosity;
  bool pedantic;

  static const std::string filename_drift_velocity;
  static const std::string filename_drift_diffusion_L;
  static const std::string filename_drift_diffusion_T;
  static const std::string filename_yield;
  static const std::string filename_spectra;
  static const std::string filename_electron_distributions;
  static const std::string filename_F_distributions;
};

#endif // ArgonPropertiesTables_h
