#include "GenNBrS_InTHGEM.hh"

GenNBrS_InTHGEM::GenNBrS_InTHGEM(double energy) :
  VGeneratePrimaries(energy)
{}

GenNBrS_InTHGEM::GenNBrS_InTHGEM(PDF_routine& energy_spectrum) :
  VGeneratePrimaries(energy_spectrum)
{}

GenNBrS_InTHGEM::~GenNBrS_InTHGEM()
{}

void GenNBrS_InTHGEM::GeneratePrimaries(G4Event* anEvent)
{
  SetupNavigator();
  SetupElectronDrift();
  ClearRecords();
	G4ParticleDefinition* particle;
	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	G4String particleName = "opticalphoton";
	particle = particleTable->FindParticle(particleName);
	if (NULL == particle) {
		G4Exception("PrimaryGeneratorAction::GeneratePrimaries: ",
			"InvalidParticle", FatalException, "Could not find 'opticalphoton' particle");
		return;
	}
	G4long seed = GetAndFixSeed(); // Results may be reproduced by fixating seed to this value and copying the rest of parameters.

	// Set electron position
  double x, y, z;
  double phi = 2 * M_PI * G4UniformRand();
  double R = std::sqrt(G4UniformRand()) * gPars::source.xy_radius;
  x = gPars::source.x_center + R*std::cos(phi);
  y = gPars::source.y_center + R*std::sin(phi);
  z = gPars::source.z_center + gPars::source.z_width*(G4UniformRand() - 0.5);
  x += anEvent->GetEventID() * 0.01;

  G4ThreeVector electron_position = G4ThreeVector(x * mm, y * mm, z * mm);
  RecordElectron(electron_position, anEvent->GetEventID(), seed);

	mElectronDrift->DoDriftElectron(x, y, z, 0);
	if (gPars::general.doViewElectronDrift)
	  mElectronDrift->Draw();

	for (std::size_t i=1, i_end_ = mElectronDrift->GetDriftTrack().size(); i!=i_end_ && i_end_!=0; ++i) {
	  //TODO: NBrS calculations are here
	  if (false) {
      // Set particle energy
      PhotonHit photon;
      PrimaryGenerator single_particle(mNavigator);
      single_particle.SetParticleDefinition(particle);
      double energy = 0;
      if (mEnergySpectrum.isValid()) {
        double rnd = G4UniformRand();
        energy = mEnergySpectrum(rnd);
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
      phi = 2 * pi * G4UniformRand();
      double R = std::sqrt(G4UniformRand()) * gPars::source.xy_radius;
      x = gPars::source.x_center + R*std::cos(phi) * mm;
      y = gPars::source.y_center + R*std::sin(phi) * mm;
      z = gPars::source.z_center + gPars::source.z_width*(G4UniformRand() - 0.5) * mm;
      //G4ThreeVector position(x, y, z);
      G4ThreeVector position = gPars::general.EL_gap_center;
      //G4ThreeVector position = gPars::debugging.THGEM1_hole_center;
      //position.setX(0 * 0.9); position.setY(0 * 0.9 * std::sqrt(3.0));
      single_particle.SetParticlePosition(position);
      photon._pos = G4ThreeVector(position);
      //------------------------------------
      G4ThreeVector polarization = G4PlaneVectorRand(direction);
      single_particle.SetParticlePolarization(polarization);
      photon._polarisation = polarization;

      double time = 0.0*us;
      single_particle.SetParticleTime(time);
      photon._time = time;

      single_particle.GeneratePrimaryVertex(anEvent);
      RecordPhoton(photon);
	  }
  }
}
