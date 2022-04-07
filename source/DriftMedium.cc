#include "DriftMedium.hh"

DriftMedium::DriftMedium(std::string name) :
    m_className("LAr_drift_medium"), m_name(name), m_driftable(false), m_debug(false)
{}

DriftMedium::DriftMedium(std::string name, DataVector &drift_velocity, DataVector &diffusion_L, DataVector &diffusion_T) :
   m_className("LAr_drift_medium"), m_name(name), m_driftable(false), m_debug(false)
{
  Initialize(drift_velocity, diffusion_L, diffusion_T);
}

DriftMedium::~DriftMedium()
{}

void DriftMedium::Initialize(DataVector &drift_velocity, DataVector &diffusion_L, DataVector &diffusion_T)
{
  m_drift_velocity = drift_velocity;
  m_diffusion_longitudinal = diffusion_L;
  m_diffusion_transverse = diffusion_T;
  if (!m_drift_velocity.isValid()) {
    std::cerr << "DriftMedium::Initialize:Warning: Drift velocity data is invalid." << std::endl;
  }
  if (!m_diffusion_longitudinal.isValid()) {
    std::cerr << "DriftMedium::Initialize:Warning: Longitudinal diffusion data is invalid." << std::endl;
  }
  if (!m_diffusion_transverse.isValid()) {
    std::cerr << "DriftMedium::Initialize:Warning: Transverse diffusion data is invalid." << std::endl;
  }
}

// Transport parameters for electrons in Geant4 units! Not in Garfield++ units!
// Drift velocity [mm / ns]
bool DriftMedium::ElectronVelocity(const double ex, const double ey, const double ez,
                      double& vx, double& vy, double& vz)
{
  vx = vy = vz = 0.0;
  if (!m_drift_velocity.isValid())
    return false;
  double E = std::sqrt(ex*ex + ey*ey + ez*ez);
  double V = m_drift_velocity(E); // is responsible for handling too high and too low fields
  if (V <= 0 && E != 0.0) {
    std::cerr<<"DriftMedium::ElectronVelocity:Error: material " << m_name <<std::endl;
    std::cerr<< "\tvelocity "<<V * us / mm <<"[mm/us] is <= 0 for electric field "<<E * cm / volt <<"[V/cm]"<<std::endl;
    return false;
  }
  const double mu = -1 * V / E;
  vx = mu * ex;
  vy = mu * ey;
  vz = mu * ez;
  return true;
}

// Longitudinal and transverse diffusion coefficients [mm1/2]
bool DriftMedium::ElectronDiffusion(const double ex, const double ey, const double ez,
                      double& dl, double& dt)
{
  dl = dt = 0.;
  if (!m_diffusion_longitudinal.isValid() && !m_diffusion_transverse.isValid())
      return true;
  double E = std::sqrt(ex*ex + ey*ey + ez*ez);
  double Ld = m_diffusion_longitudinal.isValid() ? m_diffusion_longitudinal(E) : 0; // is responsible for handling too high and too low fields
  double Td = m_diffusion_transverse.isValid() ? m_diffusion_transverse(E) : 0; // is responsible for handling too high and too low fields
  if (Ld < 0) {
    std::cerr<<"DriftMedium::ElectronDiffusion:Error: material " << m_name <<std::endl;
    std::cerr<< "\tlongitudinal diffusion "<<Ld / std::sqrt(mm) <<"[mm^1/2] is < 0 for electric field "<<E * cm / volt <<"[V/cm]"<<std::endl;
    return false;
  }
  if (Td < 0) {
    std::cerr<<"DriftMedium::ElectronDiffusion:Error: material " << m_name <<std::endl;
    std::cerr<< "\ttransverse diffusion "<<Td / std::sqrt(mm) <<"[mm^1/2] is < 0 for electric field "<<E * cm / volt <<"[V/cm]"<<std::endl;
    return false;
  }
  if (!m_diffusion_longitudinal.isValid())
    Ld = Td;
  if (!m_diffusion_transverse.isValid())
    Td = Ld;
  dl = Ld;
  dt = Td;
  return true;
}
