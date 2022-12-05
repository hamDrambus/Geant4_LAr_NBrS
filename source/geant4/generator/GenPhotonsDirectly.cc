#include <geant4/generator/GenPhotonsDirectly.hh>

GenPhotonsDirectly::GenPhotonsDirectly(double energy, PatternPhoton pattern, double angle) :
  VGeneratePrimaries(), mPattern(pattern), mAngle(angle)
{
	SetParticleEnergySpectrum(energy);
}

GenPhotonsDirectly::GenPhotonsDirectly(PDF_routine& energy_spectrum, PatternPhoton pattern, double angle) :
  VGeneratePrimaries(), mPattern(pattern), mAngle(angle)
{
	SetParticleEnergySpectrum(energy_spectrum);
}

GenPhotonsDirectly::~GenPhotonsDirectly()
{}

void GenPhotonsDirectly::SetParticleEnergySpectrum(double energy)
{
  PDF_routine spectrum;
  spectrum.insert(energy, 1.0);
  spectrum.pdf_to_cdf();
  SetParticleEnergySpectrum(spectrum);
}

void GenPhotonsDirectly::SetParticleEnergySpectrum(PDF_routine energy_spectrum)
{
  mEnergySpectrum = energy_spectrum;
  if (!mEnergySpectrum.isValid()) {
    G4Exception("VGeneratePrimaries::SetParticleEnergySpectrum()", "Event0101",
      FatalException, "Invalid energy spectrum is given.");
  }
}

void GenPhotonsDirectly::GeneratePrimaries(G4Event* anEvent)
{
  SetupNavigator();
  ClearRecords();
	G4ParticleDefinition* particleDef;
	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	G4String particleName = "opticalphoton";
	particleDef = particleTable->FindParticle(particleName);
	if (NULL == particleDef) {
		G4Exception("GenPhotonsDirectly::GeneratePrimaries: ",
			"InvalidParticle", FatalException, "Could not find 'opticalphoton' particle");
		return;
	}

	if (mGeneratedInfo.empty()) {
	  G4long seed = GetAndFixSeed();
	  RecordElectron(G4ThreeVector(0, 0, 0), 0, seed);
	}

  // Set particle energy
  PhotonHit photon;
  PrimaryGenerator single_particle(mNavigator);
  single_particle.SetParticleDefinition(particleDef);
  double energy = 0;
  if (mEnergySpectrum.isValid()) {
    double rnd = G4UniformRand();
    energy = mEnergySpectrum(rnd);
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
  if (mPattern == PatternPhoton::SiPM_shading) {
    double angle = mAngle;
    direction = G4ThreeVector(std::sin(angle), 0, std::cos(angle));
  } else {
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
    x = gPars::det_dims->THGEM_hole_center.x() + gPars::source->xy_radius*(G4UniformRand() - 0.5);
    y = gPars::det_dims->THGEM_hole_center.y() + gPars::source->xy_radius*(G4UniformRand() - 0.5);
    z = gPars::det_dims->THGEM_hole_center.z() + gPars::source->z_center;
  } else if (mPattern == PatternPhoton::THGEM1_hole_center) {
    x = gPars::det_dims->THGEM_hole_center.x() + gPars::source->x_center;
    y = gPars::det_dims->THGEM_hole_center.y() + gPars::source->y_center;
    z = gPars::det_dims->THGEM_hole_center.z() + gPars::source->z_center;
  } else if (mPattern == PatternPhoton::EL_gap_center) {
    x = gPars::det_dims->EL_gap_center.x();
    y = gPars::det_dims->EL_gap_center.y();
    z = gPars::det_dims->EL_gap_center.z();
  } else if (mPattern == PatternPhoton::Cathode_center) {
    x = gPars::det_dims->Cathode_top_center.x();
    y = gPars::det_dims->Cathode_top_center.y();
    z = gPars::det_dims->Cathode_top_center.z() + 0.1;
  } else if (mPattern == PatternPhoton::Cathode_14mm_coll) {
    double phi = 2 * pi * G4UniformRand();
    double R = std::sqrt(G4UniformRand()) * 14;
    x = gPars::det_dims->Cathode_top_center.x() + R * std::cos(phi);
    y = gPars::det_dims->Cathode_top_center.x() + R * std::sin(phi);
    z = gPars::det_dims->Cathode_top_center.z() + 0.1;
  } else if (mPattern == PatternPhoton::By_source) {
    double phi = 2 * pi * G4UniformRand();
    double R = std::sqrt(G4UniformRand()) * gPars::source->xy_radius;
    x = gPars::source->x_center + R * std::cos(phi);
    y = gPars::source->y_center + R * std::sin(phi);
    z = gPars::source->z_center + gPars::source->z_width*(G4UniformRand() - 0.5);
  } else {
    G4Exception("GenPhotonsDirectly::GeneratePrimaries: ",
        "InvalidCode", FatalException, "Used PatternPhoton is not implemented!.");
    return;
  }
  G4ThreeVector position(x, y, z);
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
  RecordPhoton(photon);
}
