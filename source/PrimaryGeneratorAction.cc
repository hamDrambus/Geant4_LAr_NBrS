#include <iostream>
using namespace std;

#include "PrimaryGeneratorAction.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction(double energy)
{
  photon_spectrum = PDF_routine();
  photon_spectrum.insert(energy, 1.0);
  photon_spectrum.pdf_to_cdf();
  SetParticleEnergySpectrum(photon_spectrum);
  navigator = nullptr;
}

PrimaryGeneratorAction::PrimaryGeneratorAction(PDF_routine& pdf)
{
  SetParticleEnergySpectrum(pdf);
  navigator = nullptr;
}


PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
	delete navigator;
}

void PrimaryGeneratorAction::SetParticleEnergySpectrum(PDF_routine energySpectrum)
{
  photon_spectrum = energySpectrum;
  if (!photon_spectrum.isValid()) {
    G4Exception("PrimaryGeneratorAction::SetParticleEnergySpectrum()", "Event0101",
      FatalException, "Invalid energy spectrum is given.");
  }
}

void PrimaryGeneratorAction::SetupNavigator(void)
{
  G4Navigator* theNavigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
  G4VPhysicalVolume* world = theNavigator->GetWorldVolume();
  if (nullptr != navigator)
    delete navigator;
  navigator = new G4Navigator();
  navigator->SetWorldVolume(world);
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  SetupNavigator();
	G4ParticleDefinition* particle;
	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	G4String particleName = "opticalphoton";
	particle = particleTable->FindParticle(particleName);
	if (NULL == particle) {
		G4Exception("PrimaryGeneratorAction::GeneratePrimaries: ",
			"InvalidParticle", FatalException, "Could not find 'opticalphoton' particle");
		return;
	}

	gPars::results.generated_photons.push_back(GeneratedData());
	gPars::results.generated_photons.back().electron.index = anEvent->GetEventID();
	gPars::results.generated_photons.back().electron.position = G4ThreeVector(0, 0, 0);
	G4long seed = CLHEP::HepRandom::getTheEngine()->operator unsigned int();
	CLHEP::HepRandom::setTheSeed(seed);
	gPars::results.generated_photons.back().electron.seed_info = "\"" + std::to_string(seed) + "\"";
	gPars::results.recorded_photons.push_back(gPars::results.generated_photons.back());
	{ // Will be 'for' loop in the future, describing electron drift and emission of NBrS photons along it

	  // Set particle energy
	  PhotonHit photon;
    PrimaryGenerator single_particle(navigator);
    single_particle.SetParticleDefinition(particle);
    double energy = 0;
    if (photon_spectrum.isValid()) {
      double rnd = G4UniformRand();
      energy = photon_spectrum(rnd);
    } else {
      G4Exception("PrimaryGeneratorAction::GeneratePrimaries: ",
        "InvalidParticle", FatalException, "Invalid energy spectrum is given.");
      return;
    }
    energy *= eV;
    single_particle.SetParticleEnergy(energy);
    photon._energy = energy;

    // Set particle direction
    double phi = 2 * M_PI * G4UniformRand();
    double cosTheta = (G4UniformRand() - 0.5) * 2;
    double sinTheta = sqrt(1 - cosTheta*cosTheta);
    G4ThreeVector direction(sinTheta*std::cos(phi), sinTheta*std::sin(phi), cosTheta);
    //G4ThreeVector direction(0, 0, 1.0);
    single_particle.SetParticleMomentumDirection(direction);
    photon._momentum = direction;

    // Set particle position
    double x, y, z;
    phi = 2 * M_PI * G4UniformRand();
    double R = std::sqrt(G4UniformRand()) * gPars::source.xy_radius;
    x = gPars::source.x_center + R*std::cos(phi) * mm;
    y = gPars::source.y_center + R*std::sin(phi) * mm;
    z = gPars::source.z_center + gPars::source.z_width*(G4UniformRand() - 0.5) * mm;
    //G4ThreeVector position(x, y, z);
    //G4ThreeVector position = gPars::debugging.EL_gap_center;
    G4ThreeVector position = gPars::debugging.THGEM1_hole_center;
    single_particle.SetParticlePosition(position);
    photon._pos = G4ThreeVector(position);
    //------------------------------------
    G4ThreeVector polarization = G4PlaneVectorRand(direction);
    single_particle.SetParticlePolarization(polarization);
    photon._polarisation = polarization;

    double time = 1.0*us;
    single_particle.SetParticleTime(time);
    photon._time = time;

    single_particle.GeneratePrimaryVertex(anEvent);
    gPars::results.generated_photons.back().photons.push_back(photon);
  }
}
