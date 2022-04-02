#ifndef PRIMARY_GENERATOR_H_
#define PRIMARY_GENERATOR_H_

#include "globals.hh"
#include "G4VPrimaryGenerator.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleDefinition.hh"
#include "G4PrimaryVertex.hh"
#include "G4ParticleMomentum.hh"
#include "G4SystemOfUnits.hh"
#include "G4PrimaryParticle.hh"
#include "G4ParticleGunMessenger.hh"
#include "G4Navigator.hh"
#include "G4Event.hh"
#include "G4ios.hh"
#include "globals.hh"
#include "Randomize.hh"
#include <G4SystemOfUnits.hh>
#include "G4RandomTools.hh"

#include "GlobalParameters.hh"
#include "GlobalUtilities.hh"
#include "HexagonalMapping.hh"
#include "PolynomialFit.hh"
#include "PrimaryParticleUserInfo.hh"

/* This class is responsible for generating photons emitted by electron
 * during its drift in electric field. Produced photons have different
 * spectra (depending on position and field) and origin positions.
 * This class also sets additional mapping information in G4PrimaryParticle via
 * G4VUserPrimaryParticleInformation. This allows for thread-safe and logical
 * mapping parameters storage.
 * !!! In case secondary particles are produced, corresponding processes
 * must be modified to set mapping parameters into G4PrimaryParticle there as well.
 */

class PrimaryGenerator : public G4VPrimaryGenerator
{
  public:
  PrimaryGenerator() = delete;
	PrimaryGenerator(G4Navigator* navig);
  virtual ~PrimaryGenerator();

  const PrimaryGenerator& operator=(const PrimaryGenerator&) = delete;
  G4bool operator==(const PrimaryGenerator&) const = delete;
  G4bool operator!=(const PrimaryGenerator&) const = delete;

  virtual void GeneratePrimaryVertex(G4Event* evt);
    // Creates a primary vertex at the given point
    // and put primary particles to it.

  // Followings are the Set methods for the particle properties.
  // SetParticleDefinition() should be called first.
  // By using SetParticleMomentum(), both particle_momentum_direction and
  // particle_energy(Kinetic Energy) are set.
  //
  void SetParticleDefinition(G4ParticleDefinition* aParticleDefinition);
  void SetParticleEnergy(G4double aKineticEnergy);
  void SetParticleMomentum(G4double aMomentum);
  void SetParticleMomentum(G4ParticleMomentum aMomentum);
  inline void SetParticleMomentumDirection(G4ParticleMomentum aMomDirection)
  { particle_momentum_direction =  aMomDirection.unit(); }
  inline void SetParticlePolarization(G4ThreeVector aVal)
  { particle_polarization = aVal; }

  inline G4ParticleDefinition* GetParticleDefinition() const
  { return particle_definition; }
  inline PDF_routine GetParticleEnergySpectrum() const
  { return particle_energy_spectrum; }
  inline G4double GetParticleEnergy() const
  { return particle_energy; }
  inline G4double GetParticleMomentum() const
  { return particle_momentum; }
  inline G4ParticleMomentum GetParticleMomentumDirection() const
  { return particle_momentum_direction; }
  inline G4ThreeVector GetParticlePolarization() const
  { return particle_polarization; }

  protected:
  G4Navigator* navigator;

  G4ParticleDefinition* 	particle_definition = nullptr;
  PDF_routine			        particle_energy_spectrum;
  G4ParticleMomentum    	particle_momentum_direction;
  G4double              	particle_energy = 0.0;
  G4double              	particle_momentum = 0.0;
  G4double              	particle_charge = 0.0;
  G4ThreeVector         	particle_polarization;
};

#endif //PRIMARY_GENERATOR_H_
