/*  Analog of Garfield++ Medium.hh
 *  This class is used to store empirical drift parameters in LAr
 *  and is used during electron drift.
 */

#ifndef DRIFT_MEDIUM_H_
#define DRIFT_MEDIUM_H_

#include <G4SystemOfUnits.hh>

#include <utilities/PolynomialFit.hh>


class DriftMedium {
 public:
  DriftMedium(std::string name);
  DriftMedium(std::string name, DataVector &drift_velocity, DataVector &diffusion_L, DataVector &diffusion_T);
  ~DriftMedium();

  void Initialize(DataVector &drift_velocity, DataVector &diffusion_L, DataVector &diffusion_T);
  std::string GetName() const { return m_name; }
  void SetName(std::string name) { m_name = name; }

  void SetDriftable(bool driftable) { m_driftable = driftable; }
  bool IsDriftable() const { return m_driftable; }

  // Transport parameters for electrons in Geant4 units! Not in Garfield++ units!
  // Drift velocity [mm / ns]
  bool ElectronVelocity(const double ex, const double ey, const double ez,
                        double& vx, double& vy, double& vz);
  // Longitudinal and transverse diffusion coefficients [mm1/2]
  bool ElectronDiffusion(const double ex, const double ey, const double ez,
                        double& dl, double& dt);

  DataVector& GetVelocityData(void) { return m_drift_velocity; };
  DataVector& GetLongDiffutionData(void) { return m_diffusion_longitudinal; };
  DataVector& GetTransDiffutionData(void) { return m_diffusion_transverse; };

  // Switch on/off debugging  messages
  void EnableDebugging() { m_debug = true; }
  void DisableDebugging() { m_debug = false; }

 protected:
  std::string m_className;
  std::string m_name;
  bool m_debug;

  DataVector m_drift_velocity;
  DataVector m_diffusion_longitudinal;
  DataVector m_diffusion_transverse;

  bool m_driftable;
};

#endif // DRIFT_MEDIUM_H_
