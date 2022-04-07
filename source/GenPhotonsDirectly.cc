#include "GenPhotonsDirectly.hh"

GenPhotonsDirectly::GenPhotonsDirectly(double energy, PatternPhoton pattern) :
  mPattern(pattern)
{
  photon_spectrum = PDF_routine();
  photon_spectrum.insert(energy, 1.0);
  photon_spectrum.pdf_to_cdf();
  SetParticleEnergySpectrum(photon_spectrum);
  navigator = nullptr;
}

GenPhotonsDirectly::GenPhotonsDirectly(PDF_routine& pdf, PatternPhoton pattern) :
  mPattern(pattern)
{
  SetParticleEnergySpectrum(pdf);
  navigator = nullptr;
}


GenPhotonsDirectly::~GenPhotonsDirectly()
{
  if (nullptr!=navigator)
    delete navigator;
}

void GenPhotonsDirectly::SetParticleEnergySpectrum(PDF_routine energySpectrum)
{
  photon_spectrum = energySpectrum;
  if (!photon_spectrum.isValid()) {
    G4Exception("PrimaryGeneratorAction::SetParticleEnergySpectrum()", "Event0101",
      FatalException, "Invalid energy spectrum is given.");
  }
}

void GenPhotonsDirectly::SetupNavigator(void)
{
  G4Navigator* theNavigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
  G4VPhysicalVolume* world = theNavigator->GetWorldVolume();
  if (nullptr != navigator)
    delete navigator;
  navigator = new G4Navigator();
  navigator->SetWorldVolume(world);
}

void GenPhotonsDirectly::GeneratePrimaries(G4Event* anEvent)
{
  SetupNavigator();
	G4ParticleDefinition* particleDef;
	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	G4String particleName = "opticalphoton";
	particleDef = particleTable->FindParticle(particleName);
	if (NULL == particleDef) {
		G4Exception("GenPhotonsDirectly::GeneratePrimaries: ",
			"InvalidParticle", FatalException, "Could not find 'opticalphoton' particle");
		return;
	}

	if (gData.results.generated_photons.empty()) { // Set only 1 dummy electron
	  gData.results.generated_photons.push_back(GlobalData::GeneratedData());
    gData.results.generated_photons.back().electron.index = 0;
    gData.results.generated_photons.back().electron.position = G4ThreeVector(0, 0, 0);
    G4long seed = CLHEP::HepRandom::getTheEngine()->operator unsigned int();
    CLHEP::HepRandom::setTheSeed(seed);
    gData.results.generated_photons.back().electron.seed_info = "\"" + std::to_string(seed) + "\"";
    gData.results.recorded_photons.push_back(gData.results.generated_photons.back());
	}

  // Set particle energy
  PhotonHit photon;
  PrimaryGenerator single_particle(navigator);
  single_particle.SetParticleDefinition(particleDef);
  double energy = 0;
  if (photon_spectrum.isValid()) {
    double rnd = G4UniformRand();
    energy = photon_spectrum(rnd);
  } else {
    G4Exception("GenPhotonsDirectly::GeneratePrimaries: ",
      "InvalidParticle", FatalException, "Invalid energy spectrum is given.");
    return;
  }
  energy *= eV;
  single_particle.SetParticleEnergy(energy);
  photon._energy = energy;

  // Set particle direction
  G4ThreeVector direction;
  if (mPattern == PatternPhoton::SiPM_shading)
    direction = G4ThreeVector(0, 0, 1);
  else {
    double phi = 2 * pi * G4UniformRand();
    double cosTheta = (G4UniformRand() - 0.5) * 2;
    double sinTheta = sqrt(1 - cosTheta*cosTheta);
    direction = G4ThreeVector(sinTheta*std::cos(phi), sinTheta*std::sin(phi), cosTheta);
  }
  single_particle.SetParticleMomentumDirection(direction);
  photon._momentum = direction;

  // Set particle position
  double x, y, z;
  if (mPattern == PatternPhoton::SiPM_shading) { // Square-shaped distribution
    x = gPars::det_dims.SiPM_size*(G4UniformRand() - 0.5);
    y = gPars::det_dims.SiPM_size*(G4UniformRand() - 0.5);
    z = gPars::general.EL_gap_center.z();
  } else if (mPattern == PatternPhoton::THGEM1_hole_center) {
    x = gPars::general.THGEM1_hole_center.x();
    y = gPars::general.THGEM1_hole_center.y();
    z = gPars::general.THGEM1_hole_center.z();
  } else if (mPattern == PatternPhoton::EL_gap_center) {
    x = gPars::general.EL_gap_center.x();
    y = gPars::general.EL_gap_center.y();
    z = gPars::general.EL_gap_center.z();
  } else if (mPattern == PatternPhoton::Cathode_center) {
    x = 0;
    y = 0;
    z = gPars::det_dims.THGEM_cathode_width + 0.1;
  } else if (mPattern == PatternPhoton::Cathode_14mm_coll) {
    double phi = 2 * pi * G4UniformRand();
    double R = std::sqrt(G4UniformRand()) * 14;
    x = R * std::cos(phi);
    y = R * std::sin(phi);
    z = gPars::det_dims.THGEM_cathode_width + 0.1;
  } else if (mPattern == PatternPhoton::By_source) {
    double phi = 2 * pi * G4UniformRand();
    double R = std::sqrt(G4UniformRand()) * gPars::source.xy_radius;
    x = gPars::source.x_center + R * std::cos(phi);
    y = gPars::source.y_center + R * std::sin(phi);
    z = gPars::source.z_center + gPars::source.z_width*(G4UniformRand() - 0.5);
  } else {
    G4Exception("GenPhotonsDirectly::GeneratePrimaries: ",
        "InvalidCode", FatalException, "Used PatternPhoton is not implemented!.");
    return;
  }
  G4ThreeVector position(x * mm, y * mm, z * mm);
  single_particle.SetParticlePosition(position);
  photon._pos = G4ThreeVector(position);
  //------------------------------------
  G4ThreeVector polarization = G4PlaneVectorRand(direction);
  single_particle.SetParticlePolarization(polarization);
  photon._polarisation = polarization;

  double time = 0.0;
  single_particle.SetParticleTime(time);
  photon._time = time;

  single_particle.GeneratePrimaryVertex(anEvent);
  gData.results.generated_photons.back().photons.push_back(photon);
}
