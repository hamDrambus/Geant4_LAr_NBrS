#include "VGeneratePrimaries.hh"

VGeneratePrimaries::VGeneratePrimaries(double energy)
{
  SetParticleEnergySpectrum(energy);
  mNavigator = nullptr;
  mElectronDrift = nullptr;
}
VGeneratePrimaries::VGeneratePrimaries(PDF_routine& energy_spectrum)
{
  SetParticleEnergySpectrum(energy_spectrum);
  mNavigator = nullptr;
  mElectronDrift = nullptr;
}

VGeneratePrimaries::~VGeneratePrimaries()
{
  if (nullptr!=mNavigator)
    delete mNavigator;
  if (nullptr!=mElectronDrift)
    delete mElectronDrift;
}

void VGeneratePrimaries::SetParticleEnergySpectrum(double energy)
{
  PDF_routine spectrum;
  spectrum.insert(energy, 1.0);
  spectrum.pdf_to_cdf();
  SetParticleEnergySpectrum(spectrum);
}

void VGeneratePrimaries::SetParticleEnergySpectrum(PDF_routine energy_spectrum)
{
  mEnergySpectrum = energy_spectrum;
  if (!mEnergySpectrum.isValid()) {
    G4Exception("VGeneratePrimaries::SetParticleEnergySpectrum()", "Event0101",
      FatalException, "Invalid energy spectrum is given.");
  }
}

void VGeneratePrimaries::SetupNavigator(void)
{
  G4Navigator* theNavigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
  G4VPhysicalVolume* world = theNavigator->GetWorldVolume();
  if (nullptr != mNavigator)
    delete mNavigator;
  mNavigator = new G4Navigator();
  mNavigator->SetWorldVolume(world);
}

void VGeneratePrimaries::SetupElectronDrift(void)
{
  if (nullptr == gData.LAr_medium || !gData.LAr_medium->GetVelocityData().isValid()) {
    G4Exception("PrimaryGeneratorAction::SetupElectronDrift: ",
      "InvalidSetup", FatalException, "Invalid drift medium");
    return;
  }
  if (nullptr != mElectronDrift)
      delete mElectronDrift;
  mElectronDrift = new DriftElectron();
  mElectronDrift->SetDistanceSteps(gPars::field_map.drift_step_size);
  if ((!gData.LAr_medium->GetLongDiffutionData().isValid() && !gData.LAr_medium->GetTransDiffutionData().isValid())
      || !gPars::general.enable_e_diffusion)
    mElectronDrift->DisableDiffusion();
  else
    mElectronDrift->EnableDiffusion();
  if (gPars::general.electron_max_time != DBL_MAX)
    mElectronDrift->SetTimeWindow(0.0, gPars::general.electron_max_time);
}

G4long VGeneratePrimaries::GetAndFixSeed(void)
{
  G4long seed = G4Random::getTheEngine()->operator unsigned int();
  G4Random::setTheSeed(seed);
  return seed;
}

void VGeneratePrimaries::RecordElectron(G4ThreeVector position, int index, G4long seed)
{
  mGeneratedInfo.push_back(GeneratedData());
  mGeneratedInfo.back().electron.index = index;
  mGeneratedInfo.back().electron.position = position;
  mGeneratedInfo.back().electron.seed_info = "\"" + std::to_string(seed) + "\"";
}

void VGeneratePrimaries::RecordPhoton(PhotonHit photon)
{
  if (mGeneratedInfo.empty())
    mGeneratedInfo.push_back(GeneratedData());
  mGeneratedInfo.back().photons.push_back(photon);
}

void VGeneratePrimaries::ClearRecords(void)
{
  mGeneratedInfo.clear();
}