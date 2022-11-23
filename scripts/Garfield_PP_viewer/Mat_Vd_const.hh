#ifndef MAT_VD_CONST_HH
#define MAT_VD_CONST_HH

#include <string>
#include <iostream>
#include <fstream>
#include <cmath>
#include "Medium.hh"

namespace Garfield {
class Mat_Vd_const : public Medium
{
protected:
  double V_drift;
public:
  Mat_Vd_const(double V_dr);
  ~Mat_Vd_const();
  // Transport parameters for electrons
  // Drift velocity [cm / ns] = const
  virtual bool ElectronVelocity(const double ex, const double ey,
                                const double ez, const double bx,
                                const double by, const double bz, double& vx,
                                double& vy, double& vz);
  // Longitudinal and transverse diffusion coefficients [cm1/2] ==0
  virtual bool ElectronDiffusion(const double ex, const double ey,
                                 const double ez, const double bx,
                                 const double by, const double bz, double& dl,
                                 double& dt);
};
}
#endif
