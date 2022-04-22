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

class GenPhotonsDirectly : public VGeneratePrimaries
{
public:
  enum PatternPhoton {
      THGEM1_hole_center,
      EL_gap_center,
      Cathode_center,
      Cathode_14mm_coll,
      SiPM_shading,
      SiPM_shading_test,
      By_source
  };

  // parameter meaning depends on pattern selected.
  GenPhotonsDirectly(double energy, PatternPhoton pattern, double parameter = 0.0); //monoenergetic case
  GenPhotonsDirectly(PDF_routine& energy_spectrum, PatternPhoton pattern, double parameter = 0.0); //continious spectrum case
  ~GenPhotonsDirectly();

  void SetPattern(PatternPhoton pattern) { mPattern = pattern; }
  void GeneratePrimaries(G4Event* anEvent);

protected:
  PatternPhoton mPattern;
  double mParameter;
};

#endif // GenPhotonsDirectly_h
