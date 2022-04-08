#ifndef Run_h
#define Run_h

#include <G4Run.hh>
#include <globals.hh>
#include <G4StatAnalysis.hh>
#include <G4RunManager.hh>
#include <G4Event.hh>

#include <G4SDManager.hh>
#include <G4HCofThisEvent.hh>
#include <G4THitsMap.hh>
#include <G4SystemOfUnits.hh>

#include "GlobalParameters.hh"
#include "GlobalData.hh"
#include "PhotonHit.hh"
#include "VGeneratePrimaries.hh"
#include "DetectorSensor.hh"

class Run : public G4Run
{
public:
  static void AddToFile(std::deque<GeneratedData> data, std::string filename);

  class Results {
  public:
    std::deque<unsigned int> SiPM_photon_n; //per each device
    std::deque<unsigned int> PMT_photon_n;
    std::deque<G4ThreeVector> SiPM_positions; //Set automatically in RunAction
    std::deque<G4ThreeVector> PMT_positions;
    std::deque<GeneratedData> generated_photons;
    std::deque<GeneratedData> recorded_photons;
    unsigned int n_generated;
    unsigned int n_reflections; // TODO: can move to some user action (stepping?) and add n_reflection to PhotonHit
    unsigned int GetNGeneratedPhotons(void) const;
    unsigned int GetNRecordedPhotons(void) const;
    void FindSensorsCoordinates(void); //Sets global positions of SiPMs and PMTs after geometry is constructed in G4RunManager::Initialize()
    void Clear(void);
    void Init(void);
    void ClearAndInit(void);
    void Print(std::ostream &str) const;
    void Merge(const Results &with);
  };

public:
  Run();
  ~Run() override;
  virtual void RecordEvent(const G4Event*) override;
  virtual void Merge(const G4Run*) override; // Note that merged info is still in storage until local run destruction
  virtual void Merged();
  Results results;
protected:
  int hit_collection_ID;
};

#endif // Run_h


