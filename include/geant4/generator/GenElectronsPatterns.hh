/*  PrimaryGeneratorAction inheritor.
 *
 *  Generates electrons in specific patterns to test
 *  electron drift (with and w/o diffusion), field map and THGEM to cell mapping
 *  (DriftElectron, FieldElmerMap and HexagonalMapping classes)
 *  which were added to otherwise standard geant4 simulation.
 *  This class is best used with cut version of the detector
 *  DetectorConstructionTHGEM1 which defines THGEM active area as
 *  cell assembly.
 *
 *  Number of electrons is defined by gPars::source.N_events but may be
 *  slightly larger in order to get symmetric pattern.
 *
 *  Does not generate actual Geant4 primary particles.
 */

#ifndef GenElectronPatterns_h
#define GenElectronPatterns_h

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

class GenElectronsPatterns : public VGeneratePrimaries
{
public:
  enum PatternElectron {
    RandomCircle,
    RandomSquare,
    UniformLineX,
    UniformLineY,
    UniformSquareGrid,
    Uniform1Ring,
    Uniform2Rings,
    Uniform3Rings,
  };
  GenElectronsPatterns();
  GenElectronsPatterns(PatternElectron pattern);
  ~GenElectronsPatterns();

public:
  virtual void GeneratePrimaries(G4Event* anEvent);
  void SetPattern(PatternElectron pattern) { mPattern = pattern; }

protected:
  int ExtraEventsN() const;
  G4ThreeVector GenPosition(int event_number) const;
  G4ThreeVector GenPosition_RandomCircle(int event_number) const;
  G4ThreeVector GenPosition_RandomSquare(int event_number) const;
  G4ThreeVector GenPosition_UniformLineX(int event_number) const;
  G4ThreeVector GenPosition_UniformLineY(int event_number) const;
  G4ThreeVector GenPosition_UniformSquareGrid(int event_number) const;
  G4ThreeVector GenPosition_Uniform1Ring(int event_number) const;
  G4ThreeVector GenPosition_Uniform2Rings(int event_number) const;
  G4ThreeVector GenPosition_Uniform3Rings(int event_number) const;
  PatternElectron mPattern;
};

class SettingsElectronPattern: public VSourceSettings
{
public:
  GenElectronsPatterns::PatternElectron pattern;
  SettingsElectronPattern() : pattern(GenElectronsPatterns::RandomCircle) { generator_type= ElectronPatterns; }
};

#endif // GenElectronPatterns_h
