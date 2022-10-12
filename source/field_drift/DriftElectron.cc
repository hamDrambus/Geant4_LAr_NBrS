#include <field_drift/DriftElectron.hh>

DriftElectron::DriftElectron()
    : m_Small(1e-25),
      m_nDrift(0),
      m_stepModel(2),
      m_tMc(20.0 * ps),
      m_dMc(1 * um),
      m_useDiffusion(true),
      m_debug(false),
      m_field_rel_delta(DBL_MAX),
      m_hasTimeWindow(false),
      m_tMin(-DBL_MAX),
      m_tMax(DBL_MAX),
      stuck_check_counter(0)
{
  m_className = "DriftElectron";
  m_drift.track.reserve(10000);
}

DriftElectron::~DriftElectron()
{}

void DriftElectron::SetTimeSteps(const double d)
{
  m_stepModel = 0;
  if (d < m_Small) {
    std::cerr << m_className << "::SetTimeSteps:\n";
    std::cerr << "    Specified step size is too small.\n";
    std::cerr << "    Using default (20 ps) instead.\n";
    m_tMc = 20.0 * ps;
  } else {
    if (m_debug) {
      std::cout << m_className << "::SetTimeSteps:\n";
      std::cout << "    Step size set to " << d / ns << " [ns].\n";
    }
    m_tMc = d;
  }
}

void DriftElectron::SetDistanceSteps(const double d)
{
  m_stepModel = 1;
  if (d < m_Small) {
    std::cerr << m_className << "::SetDistanceSteps:\n";
    std::cerr << "    Specified step size is too small.\n";
    std::cerr << "    Using default (1 um) instead.\n";
    m_dMc = 1 * um;
  } else {
    if (m_debug) {
      std::cout << m_className << "::SetDistanceSteps:\n";
      std::cout << "    Step size set to " << d / um << " um.\n";
    }
    m_dMc = d;
  }
}

void DriftElectron::SetFieldChangeLimit(const double relative_delta)
{
  if (relative_delta < m_Small) {
    std::cerr << m_className << "::SetFieldChangeLimit:\n";
    std::cerr << "    Specified relative change is too small.\n";
    std::cerr << "    Using default (1 %) instead.\n";
    m_field_rel_delta = 0.01;
  } else {
    if (m_debug) {
      std::cout << m_className << "::SetFieldChangeLimit:\n";
      std::cout << "    Step field maximum relative change to " << relative_delta << ".\n";
    }
    m_field_rel_delta = relative_delta;
  }
}

void DriftElectron::GetElectronEndpoint(double& x1, double& y1, double& z1,
                                        double& t1, int& status) const
{
  if (m_drift.track.empty()) {
    std::cerr << m_className << "::GetElectronEndpoint:\n";
    std::cerr << "    Endpoint does not exist.\n";
    status = -1;
    return;
  }
  status = 0;
  x1 = m_drift.track.back().pos.x();
  y1 = m_drift.track.back().pos.y();
  z1 = m_drift.track.back().pos.z();
  t1 = m_drift.track.back().time;
}

bool DriftElectron::DoDriftElectron(const double x0, const double y0,
                                const double z0, const double t0)
{
  if (!gData.mapping_manager.HasFieldMapping()) {
    std::cerr << m_className << "::DriftElectron:\n";
    std::cerr << "    THGEM1 field mapping is not available.\n";
    return false;
  }
  if (!gData.field_map || !gData.field_map->IsReady()) {
    std::cerr << m_className << "::DriftElectron:\n";
    std::cerr << "    Field map is not available.\n";
    return false;
  } //1.8004977683123558  -3.8561206537818222 68.410000000000011
  if (!DriftLine(x0, y0, z0, t0, -1)) {
    if (m_debug)
      std::cout<<"Electron drift failed."<<std::endl;
    return false;
  }
  if (m_debug)
    std::cout<<"Electron drifted."<<std::endl;
  return true;
}

bool DriftElectron::DriftLine(const double x0, const double y0,
                            const double z0, const double t0,
                            const int type)
{
  // Current position
  double x = x0, y = y0, z = z0;
  // Time step
  double delta;
  // Medium
  DriftMedium* medium = 0;
  // Electric field
  int status = 0;
  double ex = 0., ey = 0., ez = 0.;
  // Drift velocity
  double vx = 0., vy = 0., vz = 0., v = 0., vt = 0.;
  // Longitudinal and transverse diffusion coefficients
  double dl = 0., dt = 0.;
  // Diffusion vector
  double dx = 0., dy = 0., dz = 0., d = 0.;
  // Rotation angles
  double phi, cphi, sphi;
  double theta, ctheta, stheta;
  // Collision time
  double tau = 0.;

  // Reset the drift line.
  m_drift.track.clear();
  stuck_check_counter = 0;
  // Add the starting point to the drift line.
  DriftTrack::driftPoint point;
  point.pos = G4Point3D(x0, y0, z0);
  point.time = t0;

  bool ok = true;
  bool trapped = false;
  bool validAlphaEta = false;
  int abortReason = 0;

  if (m_hasTimeWindow && (t0 < m_tMin || t0 > m_tMax)) {
    std::cerr << m_className << "::DriftLine:\n";
    std::cerr << "    Starting time " << t0 / us << " is outside the specified\n";
    std::cerr << "    time window (" << m_tMin / us << ", " << m_tMax / us << ").\n";
    ok = false;
    abortReason = StatusOutsideTimeWindow;
  }
  // Get the electric field at the starting point.
  G4ThreeVector E = gData.GetFieldAtGlobal(G4ThreeVector(x, y, z), medium, &status);
  ex = E.x(); ey = E.y(); ez = E.z();
  // Make sure the starting point is inside a drift medium.
  if (status != 0) {
    std::cerr << m_className << "::DriftLine:\n";
    std::cerr << "    No drift medium at initial position (" << x / mm << ", " << y / mm
              << ", " << z / mm << ").\n";
    ok = false;
    abortReason = StatusLeftDriftMedium;
  }

  double e = sqrt(ex * ex + ey * ey + ez * ez);
  if (e < m_Small) {
    std::cerr << m_className << "::DriftLine:\n";
    std::cerr << "    Electric field at initial position is too small:\n";
    std::cerr << "      ex = " << ex * cm / volt << " V/cm\n";
    std::cerr << "      ey = " << ey * cm / volt << " V/cm\n";
    std::cerr << "      ez = " << ez * cm / volt << " V/cm\n";
    ok = false;
    abortReason = StatusCalculationAbandoned;
  }
  point.field = G4ThreeVector(ex, ey, ez);
  m_drift.track.push_back(point);
  m_nDrift = 1;
  while (ok) {
    if (type ==-1) {
      if (!medium->ElectronVelocity(ex, ey, ez, vx, vy, vz) ||
          !medium->ElectronDiffusion(ex, ey, ez, dl, dt)) {
        std::cerr << m_className << "::DriftLine:\n";
        std::cerr << "    Error calculating electron"
                  << " velocity or diffusion\n";
        std::cerr << "    at (" << x / mm << ", " << y / mm << ", " << z / mm << ")\n";
        ok = false;
        abortReason = StatusCalculationAbandoned;
        break;
      } else {
        if ((vx*ex+vy*ey+vz*ez)>0) {//interpolation error for LAr, electron cannot move along E
          vx=-vx;
          vy=-vy;
          vz=-vz;
          std::cerr << m_className << "::DriftLine:\n";
          std::cerr << "    Warning! Electron moves along the field\n";
          std::cerr << "    at (" << x / mm << ", " << y / mm << ", " << z / mm << ")\n";
          std::cerr << "    at E "<<std::sqrt(ex*ex+ey*ey+ez*ez) * cm / volt
                    <<"=|(" << ex * cm / volt << ", " << ey * cm / volt << ", " << ez * cm / volt << ")|\n";
        }
      }
    } else {
      std::cerr << m_className << "::DriftLine:\n";
      std::cerr << "    Unknown drift line type (" << type << ").\n";
      std::cerr << "    Program bug!\n";
      ok = false;
      abortReason = StatusCalculationAbandoned;
      return false;
    }

    if (m_debug) {
      std::cout << m_className << "::DriftLine:\n";
      std::cout << "    Drift velocity [mm/us] at " << x / mm << ", " << y / mm << ", " << z / mm
                << ": " << vx * us / mm << ", " << vy * us / mm << ", " << vz * us / mm << "\n";
    }
    v = sqrt(vx * vx + vy * vy + vz * vz);
    if (v < m_Small) {
      std::cerr << m_className << "::DriftLine:\n";
      std::cerr << "    Drift velocity at (" << x / mm << ", " << y / mm << ", " << z / mm
                << ") is too small:\n";
      std::cerr << "      vx = " << vx * us / mm << " [mm/us]\n";
      std::cerr << "      vy = " << vy * us / mm << " [mm/us]\n";
      std::cerr << "      vz = " << vz * us / mm << " [mm/us]\n";
      ok = false;
      abortReason = StatusCalculationAbandoned;
      break;
    }

    // Determine the time step.
    switch (m_stepModel) {
      case 0:
        // Fixed time steps
        delta = m_tMc;
        break;
      case 1:
        // Fixed distance steps
        delta = m_dMc / v;
        break;
      default:
        std::cerr << m_className << "::DriftLine:\n";
        std::cerr << "    Unknown stepping model.\n";
        return false;
    }

    if (DBL_MAX != m_field_rel_delta) { // limit step by electric field change along it
      G4ThreeVector new_pos (x + delta * vx, y + delta * vy, z + delta * vz);
      E = gData.GetFieldAtGlobal(new_pos, medium, &status);
      if (status == 0) {
        double new_e = std::sqrt(E*E);
        double rel_change = std::fabs(new_e - e) / e;
        if (rel_change > m_field_rel_delta) {
          delta *= m_field_rel_delta / rel_change;
        }
      }
    }

    // Draw a random diffusion direction in the particle frame.
    if (m_useDiffusion) {
      d = sqrt(v * delta);
      dx = d * CLHEP::RandGaussQ::shoot(0., dl);
      dy = d * CLHEP::RandGaussQ::shoot(0., dt);
      dz = d * CLHEP::RandGaussQ::shoot(0., dt);
    }
    if (m_debug) {
      std::cout << m_className << "::DriftLine:\n";
      std::cout << "    Adding diffusion step "
                << dx / mm << ", " << dy / mm << ", " << dz / mm << "\n";
    }
    // Compute the rotation angles to align the diffusion
    // and drift velocity vectors
    vt = sqrt(vx * vx + vy * vy);
    if (vt < m_Small) {
      phi = 0.;
      theta = halfpi;
      if (vz < 0.) theta = -theta;
    } else {
      phi = atan2(vy, vx);
      theta = atan2(vz, vt);
    }
    cphi = cos(phi);
    sphi = sin(phi);
    ctheta = cos(theta);
    stheta = sin(theta);

    // Compute the proposed end-point of this step.
    x += delta * vx + cphi * ctheta * dx - sphi * dy - cphi * stheta * dz;
    y += delta * vy + sphi * ctheta * dx + cphi * dy - sphi * stheta * dz;
    z += delta * vz + stheta * dx + ctheta * dz;

    if (m_debug) {
      std::cout << m_className << "::DriftLine:\n";
      std::cout << "    New point: "
                << x / mm << ", " << y / mm << ", " << z / mm << "\n";
      std::cout << "    Velocity: "
                << vx * us / mm << ", " << vy * us / mm << ", " << vz * us / mm << "\n";
    }
    // Compute the electric field at the new point.
    G4ThreeVector E = gData.GetFieldAtGlobal(G4ThreeVector(x, y, z), medium, &status);
    ex = E.x(); ey = E.y(), ez = E.z();
    // Check if the new position is inside a drift medium.
    if (status != 0) {
      // Try to terminate the drift line
      // close to the boundary.
      dx = x - point.pos.x();
      dy = y - point.pos.y();
      dz = z - point.pos.z();
      d = sqrt(dx * dx + dy * dy + dz * dz);
      if (d > 0.) {
        dx /= d;
        dy /= d;
        dz /= d;
      }
      while (d > BoundaryDistance) {
        delta *= 0.5;
        d *= 0.5;
        x = point.pos.x() + dx * d;
        y = point.pos.y() + dy * d;
        z = point.pos.z() + dz * d;
        // Check if the mid-point is inside the drift medium.
        E = gData.GetFieldAtGlobal(G4ThreeVector(x, y, z), medium, &status);
        if (status == 0) {
          point.pos = G4Point3D(x, y, z);
          point.time += delta;
        }
      }
      // Place the particle OUTSIDE the drift medium.
      point.pos = point.pos + G4Point3D(dx, dy, dz) * d;
      point.field = G4ThreeVector(0, 0, 0);
      m_drift.track.push_back(point);
      ++m_nDrift;
      abortReason = StatusLeftDriftMedium;
      if (m_debug) {
        std::cout << m_className << "::DriftLine:\n";
        std::cout << "    Particle left the drift medium.\n";
        std::cout << "    At " << point.pos.x() / mm << ", " << point.pos.y() / mm << ", " << point.pos.z() / mm
                  << "\n";
      }
      break;
    }

    e = sqrt(ex * ex + ey * ey + ez * ez);
    if (e < m_Small) {
      std::cerr << m_className << "::DriftLine:\n";
      std::cerr << "    Electric field at (" << x / mm << ", " << y / mm << ", " << z / mm
                << ") is too small:\n";
      std::cerr << "      ex = " << ex * cm / volt << " V/cm\n";
      std::cerr << "      ey = " << ey * cm / volt << " V/cm\n";
      std::cerr << "      ez = " << ez * cm / volt << " V/cm\n";
      ok = false;
      abortReason = StatusCalculationAbandoned;
      break;
    }
    // Add the new point to drift line.
    point.pos = G4Point3D(x, y, z);
    point.time += delta;
    point.field = G4ThreeVector(ex, ey, ez);
    m_drift.track.push_back(point);
    ++m_nDrift;

    // Check if the time is still within the specified interval.
    if (m_hasTimeWindow && point.time > m_tMax) {
      abortReason = StatusOutsideTimeWindow;
      ok = false;
      break;
    }
    if (IsStuck()) {
      std::cerr << m_className <<"::DriftLine(\n";
      std::cerr << "    Drifting particle is stuck! Aborting drift.\n"
                << "    Starting position: (" << m_drift.track.begin()->pos.x() / mm << ", " << m_drift.track.begin()->pos.y() / mm << ", "
                << m_drift.track.begin()->pos.z() / mm << ")\n"
                << "    Current position: (" << m_drift.track.back().pos.x() / mm << ", " << m_drift.track.back().pos.y() / mm << ", "
                << m_drift.track.back().pos.z() / mm << ").\n"
                << "    Diffusion:" << (m_useDiffusion ? "true" : "false") << "\n";
      abortReason = StatusStuck;
      ok = false;
      break;
    }
  }

  if (m_debug) {
    std::cout << m_className << "::DriftLine:\n";
    std::cout << "    Finished drifting electron from\n"
              << "      (" << m_drift.track.begin()->pos.x() / mm << ", " << m_drift.track.begin()->pos.y() / mm << ", "
              << m_drift.track.begin()->pos.z() / mm << ") to \n"
              << "      (" << m_drift.track.back().pos.x() / mm << ", " << m_drift.track.back().pos.y() / mm << ", "
              << m_drift.track.back().pos.z() / mm << ").\n";
  }

  if (!ok) return false;
  return true;
}

bool DriftElectron::IsStuck(void)
{
  if (stuck_check_counter < StuckCheckFrequency) {
    ++stuck_check_counter;
    return false;
  }
  stuck_check_counter = 0;
  if (m_drift.track.size() < StuckCheckIterations || StuckCheckIterations < 3)
    return false;
  double length = 0;
  double min_x = m_drift.track.back().pos.x(), max_x = min_x;
  double min_y = m_drift.track.back().pos.y(), max_y = min_y;
  double min_z = m_drift.track.back().pos.z(), max_z = min_z;
  for (std::size_t i = m_drift.track.size() - 2, i_end_ = m_drift.track.size() - StuckCheckIterations; i!=i_end_; --i) {
    length += (m_drift.track[i].pos - m_drift.track[i+1].pos).mag();
    min_x = std::min(min_x, m_drift.track[i].pos.x());
    max_x = std::max(max_x, m_drift.track[i].pos.x());
    min_y = std::min(min_y, m_drift.track[i].pos.y());
    max_y = std::max(max_y, m_drift.track[i].pos.y());
    min_z = std::min(min_z, m_drift.track[i].pos.z());
    max_z = std::max(max_z, m_drift.track[i].pos.z());
  }
  double bound = (max_x - min_x + max_y - min_y + max_z - min_z) / 3;
  if (length > 5 * bound * sqrt(StuckCheckIterations)) // In random walk (length ~=~ bound * sqrt(StuckCheckIterations))
    return true;
  return false;
}

void DriftElectron::Draw(void) const
{
  m_drift.Draw();
}

void DriftElectron::WriteDriftTrack(std::string filename) const
{
  m_drift.Write(filename);
}
