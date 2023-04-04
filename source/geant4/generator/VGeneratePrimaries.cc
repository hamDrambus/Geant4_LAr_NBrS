#include <geant4/generator/VGeneratePrimaries.hh>

VGeneratePrimaries::VGeneratePrimaries()
{
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
  if (nullptr == gData.drift_medium || !gData.drift_medium->GetVelocityData().isValid()) {
    G4Exception("PrimaryGeneratorAction::SetupElectronDrift: ",
      "InvalidSetup", FatalException, "Invalid drift medium");
    return;
  }
  if (nullptr != mElectronDrift)
      delete mElectronDrift;
  mElectronDrift = new DriftElectron();
  mElectronDrift->SetDistanceSteps(gPars::field_map.drift_step_size);
  if ((!gData.drift_medium->GetLongDiffutionData().isValid() && !gData.drift_medium->GetTransDiffutionData().isValid())
      || !gPars::general.enable_e_diffusion)
    mElectronDrift->DisableDiffusion();
  else
    mElectronDrift->EnableDiffusion();
  if (gPars::general.electron_max_time != DBL_MAX)
    mElectronDrift->SetTimeWindow(0.0, gPars::general.electron_max_time);
  if (gPars::field_map.max_rel_field_change > 0)
    mElectronDrift->SetFieldChangeLimit(gPars::field_map.max_rel_field_change);
}

G4long VGeneratePrimaries::GetAndFixSeed(void)
{
  G4long seed = G4Random::getTheEngine()->operator unsigned int();
  G4Random::setTheSeed(seed);
  return seed;
}

void VGeneratePrimaries::RecordElectron(G4ThreeVector position, int index, G4long seed, const DriftTrack &track)
{
  mGeneratedInfo.push_back(GeneratedData());
  mGeneratedInfo.back().electron.index = index;
  mGeneratedInfo.back().electron.position = position;
  mGeneratedInfo.back().electron.seed_info = "\"" + std::to_string(seed) + "\"";
  mGeneratedInfo.back().electron.track = track;
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
