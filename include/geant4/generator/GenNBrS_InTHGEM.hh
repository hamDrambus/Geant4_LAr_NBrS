#ifndef GenNBrS_InTHGEM_h
#define GenNBrS_InTHGEM_h

#include <vector>

#include <G4Event.hh>
#include <G4ParticleTable.hh>
#include <G4ParticleDefinition.hh>
#include <globals.hh>
#include <G4SystemOfUnits.hh>
#include <Randomize.hh>
#include <G4RandomTools.hh>
#include <G4Poisson.hh>
#include <G4Navigator.hh>
#include <G4TransportationManager.hh>
#include <G4VUserPrimaryGeneratorAction.hh>

#include "GlobalParameters.hh"
#include "VGeneratePrimaries.hh"
#include "SourceSettings.hh"

class GenNBrS_InTHGEM : public VGeneratePrimaries
{
public:
  GenNBrS_InTHGEM();
  ~GenNBrS_InTHGEM();

  void GeneratePrimaries(G4Event* anEvent);
  double GetNBrSYieldFactor(void) const {return mNBrS_yield_factor; }
  void SetNBrSYieldFactor(double yield_factor) { mNBrS_yield_factor = yield_factor; }
  void SetSourceXYprofile(double radius, double radius_smearing);
protected:
  double mNBrS_yield_factor;
  double mRadius;
  double mRadiusSmearing;
  PDF_routine mRadiusPDF;
};

class SettingsNBrSGenerator: public VSourceSettings
{
public:
  double NBrS_yield_factor;
  double xy_radius_smearing;
  SettingsNBrSGenerator() : NBrS_yield_factor(1), xy_radius_smearing(0.0) { generator_type = NBrS; }
};

#endif // GenNBrS_InTHGEM_h
