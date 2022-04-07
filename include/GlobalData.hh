/*  Global data classes such as geometry mapping, field map, results, etc., shared between all classes.
 *  The difference from GlobalParameters is that this class has much more logic and dependent on other
 *  classes as well as Geant4 initialization
 *  Initialize these globals only after G4RunManager::Initialize()!
 *  TODO: Make all thread-safe.
 */
#ifndef GlobalData_h
#define GlobalData_h
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <deque>

#include <G4SystemOfUnits.hh>
#include "G4ThreeVector.hh"

#include "GlobalParameters.hh"
#include "PolynomialFit.hh"
#include "HexagonalMapping.hh"
#include "PhotonHit.hh"
#include "FieldElmerMap.hh"


class GlobalData {
public:
  struct DriftElectron {
    int index;
    G4ThreeVector position;
    std::string seed_info;
  };

  struct GeneratedData {
    DriftElectron electron;
    std::deque<PhotonHit> photons;
  };

  static void AddToFile(std::deque<GeneratedData> data, std::string filename);

  class Results { // TODO: Make it thread-local and safe
  public:
    std::deque<unsigned int> SiPM_photon_n; //per each device
    std::deque<unsigned int> PMT_photon_n;
    std::deque<G4ThreeVector> SiPM_positions; //Set automatically in RunAction
    std::deque<G4ThreeVector> PMT_positions;
    std::deque<GeneratedData> generated_photons;
    std::deque<GeneratedData> recorded_photons;
    unsigned int n_reflections; // TODO: can move to some user action and add n_reflection to PhotonHit
    unsigned int GetNGeneratedPhotons(void) const;
    unsigned int GetNRecordedPhotons(void) const;
    void FindSensorsCoordinates(); //Sets global positions of SiPMs and PMTs after geometry is constructed in G4RunManager::Initialize()
  } results;

  GlobalData();
  ~GlobalData();

  void Initialize(void);
  void SetupFieldMap(void);
  void SetupTHGEM1Mapping(void); //Finds THGEM1 volumes and initializes mapping class gPars::THGEM1_mapping
  G4ThreeVector GetFieldAtGlobal(G4ThreeVector position, DriftMedium* &medium, int *status = nullptr);
  void PlotField(std::string filename, G4ThreeVector line_start, G4ThreeVector line_finish, int Num, std::string name="", double L_fine=0, int Num_fine=0);

  HexagonalMapping* THGEM1_mapping; // Set after geometry construction in G4RunManager::Initialize(). Thread safe.
  FieldElmerMap* field_map; // Depends only on gPars::
  DriftMedium* LAr_medium;

  protected:
  // Experimental diffusion is in cm^2/sec which needs to be in mm^(1/2) in ElectronDrift
  // Similar situation was in Garfield++
  static void RecalculateDiffusion(DataVector& diffusion, DataVector& velocity);
};

extern GlobalData gData;

#endif // GlobalData_h
