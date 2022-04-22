#ifndef DetectorSensor_h
#define DetectorSensor_h

#include <vector>
#include <iostream>
#include <fstream>
#include <G4VSensitiveDetector.hh>
#include <G4HCofThisEvent.hh>
#include <G4Step.hh>
#include <G4ThreeVector.hh>
#include <G4Track.hh>
#include <G4DynamicParticle.hh>
#include <G4SDManager.hh>
#include <G4ios.hh>
#include "G4ParticleTypes.hh"
#include "G4ParticleDefinition.hh"

#include <GlobalParameters.hh>
#include "PhotonHit.hh"
#include <GlobalData.hh>

class DetectorSensor : public G4VSensitiveDetector
{
  public:
    DetectorSensor(G4String name);
    ~DetectorSensor() override;

    void Initialize(G4HCofThisEvent*);
    G4bool ProcessHits(G4Step*, G4TouchableHistory*);
    G4bool ProcessHits_Optical(const G4Step* aStep, G4TouchableHistory* );
    void EndOfEvent(G4HCofThisEvent*);

    static const std::string PhotonCollectionName;
  protected:
    PhotonHitCollection *hitCollection;
};

#endif //DetectorSensor_h
