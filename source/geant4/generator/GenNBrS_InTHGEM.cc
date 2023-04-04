#include <geant4/generator/GenNBrS_InTHGEM.hh>

GenNBrS_InTHGEM::GenNBrS_InTHGEM() :
  VGeneratePrimaries()
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
	if (!gData.medium_props->IsReady()) {
	  G4Exception("PrimaryGeneratorAction::GeneratePrimaries: ",
      "InvalidData", FatalException, "Argon NBrS and/or drift parameters are not available.");
    return;
	}

	G4ThreeVector e_pos;
	G4long seed;
	bool ok = false;
	int j = 0;
	do {
    // Generate electron
    seed = GetAndFixSeed(); // Results may be reproduced by fixating seed to this value and copying the rest of parameters.
    // Set electron position
    double x, y, z;
    double phi1 = 2 * M_PI * G4UniformRand();
    double R;
    if (!mRadiusPDF.isValid()) {
    	R = std::sqrt(G4UniformRand()) * mRadius;
    } else {
    	R = mRadiusPDF(G4UniformRand());
    }
    x = gPars::source->x_center + R*std::cos(phi1);
    y = gPars::source->y_center + R*std::sin(phi1);
    z = gPars::source->z_center + gPars::source->z_width*(G4UniformRand() - 0.5);

    e_pos = G4ThreeVector(x, y, z) + gPars::det_dims->drift_start_center;
    ok = mElectronDrift->DoDriftElectron(e_pos.x(), e_pos.y(), e_pos.z(), 0);
	  ++j;
	} while (!ok && j < 50); // j is to avoid deadlock when there are issues with settings.
  if (!ok)
    return;

	if (gPars::general.doViewElectronDrift)
    RecordElectron(e_pos, anEvent->GetEventID(), seed, mElectronDrift->GetDriftTrack());
  else
    RecordElectron(e_pos, anEvent->GetEventID(), seed);

  const std::vector<DriftTrack::driftPoint> &TR = mElectronDrift->GetDriftTrack().track;
	for (std::size_t i=1, i_end_ = TR.size(); i!=i_end_ && i_end_!=0; ++i) {
	  double dl = (TR[i].pos - TR[i-1].pos).mag();
	  double E_avg = 0.5 * (TR[i].field.mag() + TR[i-1].field.mag());
	  double yield = mNBrS_yield_factor * gData.medium_props->yield(E_avg);
	  double N_ph_avg = yield * dl; // All variables are in Geant4 units.
	  G4long N_ph_actual = G4Poisson(N_ph_avg);
	  for (G4long ph = 0; ph!=N_ph_actual; ++ph) {
      // Set particle energy
      PhotonHit photon;
      PrimaryGenerator single_particle(mNavigator);
      single_particle.SetParticleDefinition(particle);
      double energy = gData.medium_props->GenPhotonEnergy(E_avg, G4UniformRand());
      if (energy <= 0) {
        G4Exception("PrimaryGeneratorAction::GeneratePrimaries: ",
            "InvalidParticle", JustWarning, "Argon NBrS photon with incorrect energy aborted.");
        continue;
      }
      single_particle.SetParticleEnergy(energy);
      photon._energy = energy;

      // Set particle direction
      double phi = 2 * M_PI * G4UniformRand();
      double cosTheta = (G4UniformRand() - 0.5) * 2;
      double sinTheta = sqrt(1 - cosTheta*cosTheta);
      G4ThreeVector direction(sinTheta*std::cos(phi), sinTheta*std::sin(phi), cosTheta);
      single_particle.SetParticleMomentumDirection(direction);
      photon._momentum = direction;

      // Set particle position
      double t = G4UniformRand();
      G4Point3D position = TR[i-1].pos + t * (TR[i].pos - TR[i-1].pos);
      single_particle.SetParticlePosition(position);
      photon._pos = G4ThreeVector(position);
      //------------------------------------
      G4ThreeVector polarization = G4PlaneVectorRand(direction);
      single_particle.SetParticlePolarization(polarization);
      photon._polarisation = polarization;

      double time = TR[i-1].time + t * (TR[i].time - TR[i-1].time);
      single_particle.SetParticleTime(time);
      photon._time = time;

      single_particle.GeneratePrimaryVertex(anEvent);
      RecordPhoton(photon);
	  }
  }
}

void GenNBrS_InTHGEM::SetSourceXYprofile(double radius, double radius_smearing)
{
	mRadius = radius;
	mRadiusSmearing = radius_smearing;
	if (mRadiusSmearing > 0) {
		double r_max = mRadius + mRadiusSmearing * 7;
		IntegrationRange rs = IntegrationInterval(0, std::max(mRadius - mRadiusSmearing*7, 0.0), std::max(mRadius - mRadiusSmearing*7, 0.0)/200.0);
		rs += IntegrationInterval(std::max(mRadius - mRadiusSmearing*7, 0.0), r_max, std::min(mRadiusSmearing/50.0, mRadius/200));
		std::vector<double> xs(rs.NumOfIndices()), ys(rs.NumOfIndices());
		for (long int i = 0, i_end_ = rs.NumOfIndices(); i!=i_end_; ++i) {
			xs[i] = rs.Value(i);
			ys[i] = xs[i] * 1.0 / (exp((xs[i] - mRadius)/mRadiusSmearing) + 1);
		}
		mRadiusPDF = PDF_routine(xs, ys);
		if (!mRadiusPDF.isValid()) {
			G4Exception("GenNBrS_InTHGEM::SetSourceXYprofile: ",
					"InvalidData", FatalException, "Failed to set source radius distribution function.");
				return;
		}
	}
}
