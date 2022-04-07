#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include <vector>

#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"
#include <G4SystemOfUnits.hh>
#include "Randomize.hh"
#include "G4RandomTools.hh"
#include "G4Navigator.hh"
#include "G4TransportationManager.hh"
#include "G4VUserPrimaryGeneratorAction.hh"

#include "GlobalParameters.hh"
#include "PolynomialFit.hh"
#include "PrimaryGenerator.hh"
#include "GlobalData.hh"
#include "DriftElectron.hh"

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
	public:
		PrimaryGeneratorAction(double energy); //monoenergetic case
		PrimaryGeneratorAction(PDF_routine& pdf); //continious spectrum case
		~PrimaryGeneratorAction();

	public:
		void GeneratePrimaries(G4Event* anEvent);
		void SetParticleEnergySpectrum(PDF_routine energySpectrum);
		//void PlotField(std::string filename, G4ThreeVector line_start, G4ThreeVector line_finish, int Num, std::string name="", double L_fine=0, int Num_fine=0);

	protected:
		void SetupNavigator(void);
		void SetupElectronDrift(void);

		//void SetupFieldMap(void);
		//G4ThreeVector GetFieldAtGlobal(G4ThreeVector position);

	private:
		PDF_routine photon_spectrum;
		G4Navigator* navigator; // Separate from tracking navigator is used to locate points before particle generation.
		DriftElectron* electron_drift;
};

#endif
