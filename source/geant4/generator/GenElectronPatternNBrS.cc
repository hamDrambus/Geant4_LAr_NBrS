#include <geant4/generator/GenElectronsPatternsNBrS.hh>

GenElectronsPatternsNBrS::GenElectronsPatternsNBrS() :
	GenElectronsPatterns(), mNBrS_yield_factor(1.0)
{}

GenElectronsPatternsNBrS::~GenElectronsPatternsNBrS()
{}

void GenElectronsPatternsNBrS::GeneratePrimaries(G4Event* anEvent)
{
	SetupNavigator();
  SetupElectronDrift();
  ClearRecords();
  G4ParticleDefinition* particle;
	G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
	G4String particleName = "opticalphoton";
	particle = particleTable->FindParticle(particleName);
	if (NULL == particle) {
		G4Exception("GenElectronsPatternsNBrS::GeneratePrimaries: ",
			"InvalidParticle", FatalException, "Could not find 'opticalphoton' particle");
		return;
	}
	if (!gData.medium_props->IsReady()) {
		G4Exception("GenElectronsPatternsNBrS::GeneratePrimaries: ",
			"InvalidData", FatalException, "Argon NBrS and/or drift parameters are not available.");
		return;
	}

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

					single_particle.SetDisableMapping(gPars::det_dims->detector_type == gPars::det_dims->DetectorType::THGEM1_detailed ? true : false);
					single_particle.GeneratePrimaryVertex(anEvent);
					RecordPhoton(photon);
				}
			}
    }
    ++ID;
  } while ((anEvent->GetEventID() == gPars::source->N_events - 1) && ID != (gPars::source->N_events + ExtraEventsN()));
  // Loop triggers only at the end of beamOn and when extra events are required.
}


