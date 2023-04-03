#include <geant4/generator/GenElectronsPatterns.hh>

GenElectronsPatterns::GenElectronsPatterns() :
  VGeneratePrimaries(), mPattern(PatternElectron::RandomCircle)
{}

GenElectronsPatterns::GenElectronsPatterns(PatternElectron pattern) :
  VGeneratePrimaries(), mPattern(pattern)
{}

GenElectronsPatterns::~GenElectronsPatterns()
{}

void GenElectronsPatterns::GeneratePrimaries(G4Event* anEvent)
{
  SetupElectronDrift();
  ClearRecords();
  int ID = anEvent->GetEventID();
  do {
    G4long seed = GetAndFixSeed(); // Results may be reproduced by fixating seed to this value and copying the rest of parameters.
    G4ThreeVector e_pos = GenPosition(ID);
    bool ok = mElectronDrift->DoDriftElectron(e_pos.x(), e_pos.y(), e_pos.z(), 0);
    if (ok) {
      // Turns out, as of v10 Geant4, G4VVisManager::Draw only works after worker threads have finished
      // (doi:10.1088/1742-6596/513/2/022005 page 7). So electron track has to be saved in Run results.
      //if (gPars::general.doViewElectronDrift)
      //  mElectronDrift->Draw();
      if (gPars::general.doViewElectronDrift)
        RecordElectron(e_pos, ID, seed, mElectronDrift->GetDriftTrack());
      else
        RecordElectron(e_pos, ID, seed);
      if (gPars::general.print_drift_track) {
        std::string filename = gPars::general.output_folder + "e_track_T";
        filename += std::to_string(G4Threading::G4GetThreadId());
        filename += "_E" + std::to_string(ID) + ".txt";
        mElectronDrift->WriteDriftTrack(filename);
      }
    }
    ++ID;
  } while ((anEvent->GetEventID() == gPars::source->N_events - 1) && ID != (gPars::source->N_events + ExtraEventsN()));
  // Loop triggers only at the end of beamOn and when extra events are required.
}

G4ThreeVector GenElectronsPatterns::GenPosition(int event_number) const
{
  switch(mPattern) {
  case RandomCircle:
    return GenPosition_RandomCircle(event_number);
  case RandomSquare:
    return GenPosition_RandomSquare(event_number);
  case UniformLineX:
    return GenPosition_UniformLineX(event_number);
  case UniformLineY:
    return GenPosition_UniformLineY(event_number);
  case UniformSquareGrid:
    return GenPosition_UniformSquareGrid(event_number);
  case Uniform1Ring:
    return GenPosition_Uniform1Ring(event_number);
  case Uniform2Rings:
    return GenPosition_Uniform2Rings(event_number);
  case Uniform3Rings:
    return GenPosition_Uniform3Rings(event_number);
  default: {
    G4Exception("GenElectronsPatterns::GenPosition: ",
          "InvalidSetup", FatalException, "Chosen electron pattern is not implemented!");
  }
  }
  return G4ThreeVector(0, 0, 0);
}

G4ThreeVector GenElectronsPatterns::GenPosition_RandomCircle(int event_number) const
{
  double x, y, z;
  double phi = 2 * M_PI * G4UniformRand();
  double R = std::sqrt(G4UniformRand()) * gPars::source->xy_radius;
  x = gPars::source->x_center + R*std::cos(phi);
  y = gPars::source->y_center + R*std::sin(phi);
  z = gPars::source->z_center;
  return gPars::det_dims->drift_start_center + G4ThreeVector(x, y, z);
}

G4ThreeVector GenElectronsPatterns::GenPosition_RandomSquare(int event_number) const
{
  double x, y, z;
  x = gPars::source->x_center + (G4UniformRand() - 0.5) * gPars::source->xy_radius;
  y = gPars::source->y_center + (G4UniformRand() - 0.5) * gPars::source->xy_radius;
  z = gPars::source->z_center;
  return gPars::det_dims->drift_start_center + G4ThreeVector(x, y, z);
}

G4ThreeVector GenElectronsPatterns::GenPosition_UniformLineX(int event_number) const
{
  double x, y, z;
  double step = 2.0 * gPars::source->xy_radius / (gPars::source->N_events - 1);
  x = gPars::source->x_center + event_number * step - (gPars::source->N_events - 1) / 2.0 * step;
  y = gPars::source->y_center;
  z = gPars::source->z_center;
  return gPars::det_dims->drift_start_center + G4ThreeVector(x, y, z);
}

G4ThreeVector GenElectronsPatterns::GenPosition_UniformLineY(int event_number) const
{
  double x, y, z;
  double step = 2.0 * gPars::source->xy_radius / (gPars::source->N_events - 1);
  x = gPars::source->x_center;
  y = gPars::source->y_center + event_number * step - (gPars::source->N_events - 1) / 2.0 * step;
  z = gPars::source->z_center;
  return gPars::det_dims->drift_start_center + G4ThreeVector(x, y, z);
}

G4ThreeVector GenElectronsPatterns::GenPosition_UniformSquareGrid(int event_number) const
{
  double x, y, z;
  int N = gPars::source->N_events;
  int real_n = std::pow(std::sqrt(N), 2) < N ?
      std::pow(std::sqrt(N) + 1, 2) : N;
  N = std::sqrt(real_n);
  double stepX = 2.0 * gPars::source->xy_radius / (N - 1);
  double stepY = stepX;
  x = gPars::source->x_center + (event_number % N) * stepX - (N - 1) / 2.0 * stepX;
  y = gPars::source->y_center + (event_number / N) * stepY - (N - 1) / 2.0 * stepY;
  z = gPars::source->z_center;
  return gPars::det_dims->drift_start_center + G4ThreeVector(x, y, z);
}

G4ThreeVector GenElectronsPatterns::GenPosition_Uniform1Ring(int event_number) const
{
  double x, y, z;
  double step_phi = 2.0 * pi / gPars::source->N_events;
  double R = gPars::source->xy_radius;
  double phi = event_number * step_phi;
  x = gPars::source->x_center + R*std::cos(phi);
  y = gPars::source->y_center + R*std::sin(phi);
  z = gPars::source->z_center;
  return gPars::det_dims->drift_start_center + G4ThreeVector(x, y, z);
}

G4ThreeVector GenElectronsPatterns::GenPosition_Uniform2Rings(int event_number) const
{
  double x, y, z;
  int N = gPars::source->N_events;
  int real_n = N + N % 2;
  N = real_n / 2;
  double step_phi = 2.0 * pi / N;
  double R = gPars::source->xy_radius / (1 + event_number / N);
  double phi = (event_number % N) * step_phi;
  x = gPars::source->x_center + R*std::cos(phi);
  y = gPars::source->y_center + R*std::sin(phi);
  z = gPars::source->z_center;
  return gPars::det_dims->drift_start_center + G4ThreeVector(x, y, z);
}

G4ThreeVector GenElectronsPatterns::GenPosition_Uniform3Rings(int event_number) const
{
  double x, y, z;
  int N = gPars::source->N_events;
  int real_n = N + (N % 3 ? 0 : (3 - N % 3));
  N = real_n / 3;
  double step_phi = 2.0 * pi / N;
  double R = gPars::source->xy_radius / (1 + event_number / N);
  double phi = (event_number % N) * step_phi;
  x = gPars::source->x_center + R*std::cos(phi);
  y = gPars::source->y_center + R*std::sin(phi);
  z = gPars::source->z_center;
  return gPars::det_dims->drift_start_center + G4ThreeVector(x, y, z);
}

int GenElectronsPatterns::ExtraEventsN() const
{
  switch(mPattern) {
  case UniformSquareGrid: {
    int N = gPars::source->N_events;
    int real_n = std::pow(std::sqrt(N), 2) < N ?
        std::pow(std::sqrt(N) + 1, 2) : N;
    return real_n - N;
  }
  case Uniform2Rings: {
    int N = gPars::source->N_events;
    int real_n = N + N % 2;
    return real_n - N;
  }
  case Uniform3Rings: {
    int N = gPars::source->N_events;
    int real_n = N + (N % 3 ? 0 : (3 - N % 3));
    return real_n - N;
  }
  default: {
    return 0;
  }
  }
}

