/*  PrimaryGeneratorAction inheritor.
 *
 *  Generates electrons in specific patterns to display
 *  electron drift (with and w/o diffusion) in THGEM cells
 *  and NBrS generation.
 *  This class is supposed to be used with cut version of the detector
 *  Detector_THGEM1_detailed which defines THGEM active area as
 *  cell assembly.
 *
 *  Number of electrons is defined by gPars::source.N_events but may be
 *  slightly larger in order to get symmetric pattern.
 *
 *  Generates NBrS optical photons.
 */

#ifndef GenElectronPatternsNBrS_h
#define GenElectronPatternsNBrS_h

#include <vector>

#include <G4Event.hh>
#include <G4ParticleTable.hh>
#include <G4ParticleDefinition.hh>
#include <globals.hh>
#include <G4SystemOfUnits.hh>
#include <Randomize.hh>
#include <G4RandomTools.hh>
#include <G4Navigator.hh>
#include <G4Poisson.hh>
#include <G4TransportationManager.hh>
#include <G4VUserPrimaryGeneratorAction.hh>

#include "GlobalParameters.hh"
#include "VGeneratePrimaries.hh"
#include "GenElectronsPatterns.hh"
#include "SourceSettings.hh"

class GenElectronsPatternsNBrS : public GenElectronsPatterns
{
public:
	GenElectronsPatternsNBrS();
  ~GenElectronsPatternsNBrS();

  void GeneratePrimaries(G4Event* anEvent) override;
	void SetNBrSYieldFactor(double yield_factor) { mNBrS_yield_factor = yield_factor; }
  double GetNBrSYieldFactor(void) const {return mNBrS_yield_factor; }

protected:
	double mNBrS_yield_factor;
};

class SettingsElectronPatternNBrS : public SettingsElectronPattern
{
public:
	double NBrS_yield_factor;
	SettingsElectronPatternNBrS() : SettingsElectronPattern(), NBrS_yield_factor(1.0)
	{ generator_type = ElectronPatternsNBrS; }
};

#endif // GenElectronPatterns_h
