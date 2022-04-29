#ifndef VGeneratrePrimaries_h
#define VGeneratrePrimaries_h

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

#include <utilities/PolynomialFit.hh>
#include <field_drift/DriftElectron.hh>
#include <geant4/PhotonHit.hh>
#include "GlobalParameters.hh"
#include "PrimaryGenerator.hh"

/*
 *  Helper class. Defining useful functions which are reused in concrete generators
 */
class VGeneratePrimaries : public G4VUserPrimaryGeneratorAction
{
public:
  VGeneratePrimaries();
  ~VGeneratePrimaries();

public:
  virtual void GeneratePrimaries(G4Event* anEvent) = 0;
  inline const std::deque<GeneratedData>& GetGeneratedData(void) const
  { return mGeneratedInfo; }

protected:
  void SetupNavigator(void);
  virtual void SetupElectronDrift(void);
  G4long GetAndFixSeed(void);
  void RecordElectron(G4ThreeVector position, int index, G4long seed, const DriftTrack &track = DriftTrack());
  void RecordPhoton(PhotonHit photon);
  void ClearRecords(void);

  G4Navigator* mNavigator; // Separate from tracking navigator is used to locate points before particle generation.
  DriftElectron* mElectronDrift;

  std::deque<GeneratedData> mGeneratedInfo;
};

#endif
