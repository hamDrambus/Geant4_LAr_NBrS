#include <geant4/generator/PrimaryGenerator.hh>

PrimaryGenerator::PrimaryGenerator(G4Navigator* navig)
{
  particle_definition = nullptr;
  G4ThreeVector zero;
  particle_momentum_direction = (G4ParticleMomentum)zero;
  particle_energy = 0.0;
  particle_momentum = 0.0;
  particle_position = zero;
  particle_time = 0.0;
  particle_polarization = zero;
  particle_charge = 0.0;
  navigator = navig;
  if(navigator == nullptr) {
    G4Exception("PrimaryGenerator::PrimaryGenerator()", "Event0109",
          FatalException, "NULL navigator is given");
    return;
  }
}

PrimaryGenerator::~PrimaryGenerator()
{}

void PrimaryGenerator::SetParticleDefinition(G4ParticleDefinition* aParticleDefinition)
{
  if(aParticleDefinition == nullptr) {
    G4Exception("G4ParticleGun::SetParticleDefinition()", "Event0101",
                FatalException, "Null pointer is given.");
  }
  if(aParticleDefinition->IsShortLived()) {
    if(aParticleDefinition->GetDecayTable() == nullptr) {
      G4ExceptionDescription ED;
      ED << "PrimaryGenerator does not support shooting a short-lived "
         << "particle without a valid decay table." << G4endl;
      ED << "PrimaryGenerator::SetParticleDefinition for "
         << aParticleDefinition->GetParticleName() << " is ignored." << G4endl;
      G4Exception("G4ParticleGun::SetParticleDefinition()", "Event0102",
                  JustWarning, ED);
      return;
    }
  }
  particle_definition = aParticleDefinition;
  particle_charge = particle_definition->GetPDGCharge();
  if(particle_momentum>0.0) {
    G4double mass = particle_definition->GetPDGMass();
    particle_energy =
      std::sqrt(particle_momentum*particle_momentum+mass*mass)-mass;
  }
}

void PrimaryGenerator::SetParticleEnergy(G4double aKineticEnergy)
{
  particle_energy = aKineticEnergy;
  if(particle_momentum>0.0) {
    if(particle_definition != nullptr) {
      G4cout << "PrimaryGenerator::" << particle_definition->GetParticleName()
             << G4endl;
    } else {
      G4cout << "PrimaryGenerator::" << " " << G4endl;
    }
    G4cout << " was defined in terms of Momentum: "
           << particle_momentum/GeV << "GeV/c" << G4endl;
    G4cout << " is now defined in terms of KineticEnergy: "
           << particle_energy/GeV   << "GeV"   << G4endl;
    particle_momentum = 0.0;
  }
}

void PrimaryGenerator::SetParticleMomentum(G4double aMomentum)
{
  if(particle_energy>0.0) {
    if(particle_definition != nullptr) {
      G4cout << "PrimaryGenerator::" << particle_definition->GetParticleName()
             << G4endl;
    } else {
      G4cout << "PrimaryGenerator::" << " " << G4endl;
    }
    G4cout << " was defined in terms of KineticEnergy: "
           << particle_energy/GeV << "GeV"   << G4endl;
    G4cout << " is now defined in terms Momentum: "
           << aMomentum/GeV       << "GeV/c" << G4endl;
  }
  if(particle_definition == nullptr) {
    G4cout << "Particle Definition not defined yet for G4ParticleGun"
           << G4endl;
    G4cout << "Zero Mass is assumed" << G4endl;
    particle_momentum = aMomentum;
    particle_energy = aMomentum;
  } else {
    G4double mass = particle_definition->GetPDGMass();
    particle_momentum = aMomentum;
    particle_energy =
      std::sqrt(particle_momentum*particle_momentum+mass*mass)-mass;
  }
}

void PrimaryGenerator::SetParticleMomentum(G4ParticleMomentum aMomentum)
{
  if(particle_energy>0.0) {
    if(particle_definition != nullptr) {
      G4cout << "PrimaryGenerator::" << particle_definition->GetParticleName()
             << G4endl;
    } else {
      G4cout << "PrimaryGenerator::" << " " << G4endl;
    }
    G4cout << " was defined in terms of KineticEnergy: "
           << particle_energy/GeV << "GeV"   << G4endl;
    G4cout << " is now defined in terms Momentum: "
           << aMomentum.mag()/GeV << "GeV/c" << G4endl;
  }
  if(particle_definition == nullptr) {
    G4cout << "Particle Definition not defined yet for PrimaryGenerator"
           << G4endl;
    G4cout << "Zero Mass is assumed" << G4endl;
    particle_momentum_direction = aMomentum.unit();
    particle_momentum = aMomentum.mag();
    particle_energy = aMomentum.mag();
  } else {
    G4double mass =  particle_definition->GetPDGMass();
    particle_momentum = aMomentum.mag();
    particle_momentum_direction = aMomentum.unit();
    particle_energy =
      std::sqrt(particle_momentum*particle_momentum+mass*mass)-mass;
  }
}

void PrimaryGenerator::GeneratePrimaryVertex(G4Event* event)
{
	if(particle_definition == nullptr) {
		G4ExceptionDescription ED;
		ED << "Particle definition is not defined." << G4endl;
		ED << "PrimaryGenerator::SetParticleDefinition() has to be invoked beforehand."
		   << G4endl;
		G4Exception("PrimaryGenerator::GeneratePrimaryVertex()", "Event0109",
					FatalException, ED);
		return;
	}

	HexagonalMappingData mapping_data;
  mapping_data.position = particle_position;
  mapping_data.momentum = particle_momentum_direction;
  mapping_data.polarization = particle_polarization;
  if (gData.mapping_manager.HasMapping() && !disable_mapping) {
    G4VPhysicalVolume* volume = navigator->LocateGlobalPointAndSetup(mapping_data.position);
    mapping_data = gData.mapping_manager.GetNewState(volume, mapping_data);
  }
  G4PrimaryVertex* vertex = new G4PrimaryVertex(mapping_data.position, particle_time);
  G4double mass =  particle_definition->GetPDGMass();
  G4PrimaryParticle* particle = new G4PrimaryParticle(particle_definition);
  particle->SetKineticEnergy(particle_energy);
  particle->SetMass(mass);
  particle->SetMomentumDirection(mapping_data.momentum);
  particle->SetCharge(particle_charge);
  particle->SetPolarization(mapping_data.polarization.x(),
      mapping_data.polarization.y(), mapping_data.polarization.z());
  PrimaryParticleUserInfo *map_info = new PrimaryParticleUserInfo(mapping_data);
  particle->SetUserInformation(map_info);
  vertex->SetPrimary(particle);
  event->AddPrimaryVertex(vertex);
}
