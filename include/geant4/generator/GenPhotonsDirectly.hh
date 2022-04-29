/*  PrimaryGeneratorAction inheritor.
 *
 *  Generates optical photons in several configurations to test
 *  detector geometry, optical parameters and cell mapping using light collection efficiency.
 *  This class should be used with DetectorConstruction geometry.
 *
 *  Does not have anything to do with electron drift of field map.
 */

#ifndef GenPhotonsDirectly_h
#define GenPhotonsDirectly_h

#include <vector>

#include <G4Event.hh>
#include <G4ParticleTable.hh>
#include <G4ParticleDefinition.hh>
#include <globals.hh>
#include <G4SystemOfUnits.hh>
#include <Randomize.hh>
#include <G4RandomTools.hh>
#include <G4Navigator.hh>
#include <G4TransportationManager.hh>
#include <G4VUserPrimaryGeneratorAction.hh>

#include "GlobalParameters.hh"
#include "VGeneratePrimaries.hh"
#include "SourceSettings.hh"

class GenPhotonsDirectly : public VGeneratePrimaries
{
public:
  enum PatternPhoton {
      THGEM1_hole_center,
      EL_gap_center,
      Cathode_center,
      Cathode_14mm_coll,
      SiPM_shading,
      By_source
  };

  GenPhotonsDirectly(double energy, PatternPhoton pattern, double angle = 0.0); //monoenergetic case
  GenPhotonsDirectly(PDF_routine& energy_spectrum, PatternPhoton pattern, double angle = 0.0); //continious spectrum case
  ~GenPhotonsDirectly();

  void SetParticleEnergySpectrum(PDF_routine energy_spectrum);
  void SetParticleEnergySpectrum(double energy);
  void SetParticleAngle(double angle) { mAngle = angle; }

  void SetPattern(PatternPhoton pattern) { mPattern = pattern; }
  void GeneratePrimaries(G4Event* anEvent);

protected:
  PDF_routine mEnergySpectrum;
  PatternPhoton mPattern;
  double mAngle;
};

class SettingsDirectPhotons: public VSourceSettings
{
public:
  GenPhotonsDirectly::PatternPhoton pattern;
  double energy;
  double angle; // Photon angle to z axis in Geant4 units. Works only for SiPM_shading pattern
  PDF_routine energy_spectrum;
  std::string energy_spectrum_filename;
  SettingsDirectPhotons() : pattern(GenPhotonsDirectly::By_source), energy(-1), angle(0)
  { generator_type = PhotonsDirectly; }
};

#endif // GenPhotonsDirectly_h
