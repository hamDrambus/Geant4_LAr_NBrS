#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include <vector>

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"
#include "Randomize.hh"
#include <G4SystemOfUnits.hh>
#include "G4RandomTools.hh"
#include "G4Navigator.hh"
#include "G4TransportationManager.hh"
#include "G4VUserPrimaryGeneratorAction.hh"

#include "GlobalParameters.hh"
#include "PolynomialFit.hh"
#include "PrimaryGenerator.hh"

//TODO: place electron drift and multiple photon emission here.
//Use primary generator to just queue particle to generate and set their correct
//parameters.

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
	public:
		PrimaryGeneratorAction(double energy); //monoenergetic case
		PrimaryGeneratorAction(PDF_routine& pdf); //continious spectrum case
		~PrimaryGeneratorAction();

	public:
		void GeneratePrimaries(G4Event* anEvent);
		void SetParticleEnergySpectrum(PDF_routine energySpectrum);
		void SetupNavigator(void);

	private:
		PDF_routine photon_spectrum;
		G4Navigator* navigator; // Separate from tracking navigator is used to locate points before particle generation.
};

#endif
