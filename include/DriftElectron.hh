/*  Heavily simplified copy of Garfield++ AvalancheMC
 */

#ifndef DRIFT_ELECTRON_H_
#define DRIFT_ELECTRON_H_

#include <vector>
#include <string>
#include <cmath>
#include <iostream>
#include <fstream>

#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4RandomTools.hh"
#include "G4Polyline.hh"
#include "G4VisAttributes.hh"

#include "PolynomialFit.hh"
#include "DriftMedium.hh"
#include "FieldElmerMap.hh"
#include "GlobalData.hh"

// Status codes for drift lines
static const int StatusAlive = 0;
static const int StatusLeftDriftArea = -1;
static const int StatusTooManySteps = -2;
static const int StatusCalculationAbandoned = -3;
static const int StatusLeftDriftMedium = -5;
static const int StatusAttached = -7;
static const int StatusBelowTransportCut = -16;
static const int StatusOutsideTimeWindow = -17;

static const double BoundaryDistance = 1.e-8;

class DriftElectron {

 public:
  DriftElectron();
  ~DriftElectron();

  // Switch on/off diffusion (default: enabled)
  void EnableDiffusion() { m_useDiffusion = true; }
  void DisableDiffusion() { m_useDiffusion = false; }

  void SetTimeWindow(const double t0, const double t1) {
    if (fabs(t1 - t0) < m_Small) {
      std::cerr << m_className << "::SetTimeWindow:\n";
      std::cerr << "    Time interval must be greater than zero.\n";
      return;
    }
    m_hasTimeWindow = true;
    m_tMin = std::min(t0, t1);
    m_tMax = std::max(t0, t1);
  }
  void UnsetTimeWindow(void) { m_hasTimeWindow = false; }

  // Stepping model
  // Fixed time step (default 20 ps)
  void SetTimeSteps(const double d = 20.0 * ps);
  // Fixed distance step (default 10 um)
  void SetDistanceSteps(const double d = 1 * um);
  // Add limitation on field change along step |(E1 - E2)/E1| < relative_delta
  void SetFieldChangeLimit(const double relative_delta = 0.01);

  const DriftTrack& GetDriftTrack(void) const { return m_drift; }

  void GetElectronEndpoint(double& x1, double& y1, double& z1,
                           double& t1, int& status) const;
  bool DoDriftElectron(const double x0, const double y0, const double z0, const double t0);
  // Switch on/off debugging messages
  void EnableDebugging() { m_debug = true; }
  void DisableDebugging() { m_debug = false; }

  void Draw(void) const;
  void WriteDriftTrack(std::string filename) const;

 protected:
  std::string m_className;
  double m_Small;
  bool m_useDiffusion;

  unsigned int m_nDrift;
  DriftTrack m_drift;

  // Step size model
  int m_stepModel;
  // Fixed time step
  double m_tMc;
  // Fixed distance step
  double m_dMc;
  double m_field_rel_delta;

  // Time window
  bool m_hasTimeWindow;
  double m_tMin, m_tMax;

  bool m_debug;

  // Compute a drift line with starting point (x0, y0, z0, t0)
  bool DriftLine(const double x0, const double y0, const double z0,
                 const double t0, const int type);
};

#endif // DRIFT_ELECTRON_H_
