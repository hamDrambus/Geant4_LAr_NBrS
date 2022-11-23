#include "Mat_Vd_const.hh"

namespace Garfield {
Mat_Vd_const::Mat_Vd_const(double V_d): Medium(), V_drift(V_d)
{
  m_className="Material_Vd_const";
  m_driftable = true;
}

Mat_Vd_const::~Mat_Vd_const()
{}

bool Mat_Vd_const::ElectronVelocity(const double ex, const double ey,
                                const double ez, const double bx,
                                const double by, const double bz, double& vx,
                                double& vy, double& vz)
{
  vx=vy=vz=0.0;
  if (!((0==bx)&&(0==by)&&(0==bz))) {
     std::cerr<<"Magnetic field is not supported"<<std::endl;
     return false;
  }
  double E = std::sqrt(ex*ex+ey*ey+ez*ez);
  double V = V_drift*1e-9; //velocity in Garfield is in cm/ns, not in cm/s
  vx = -ex*V/E;
  vy = -ey*V/E;
  vz = -ez*V/E;
  return true;
}

bool Mat_Vd_const::ElectronDiffusion(const double ex, const double ey,
                                 const double ez, const double bx,
                                 const double by, const double bz, double& dl,
                                 double& dt)
{
  dt = dl = 0.0;
  if (!((0==bx)&&(0==by)&&(0==bz))) {
     std::cerr<<"Magnetic field is not supported"<<std::endl;
     return false;
  }
  return true;
}
}
