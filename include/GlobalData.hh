/*  Global data classes such as geometry mapping, field map, etc., shared between all classes.
 *  The difference from GlobalParameters is that this class has much more logic and dependent on other
 *  classes as well as Geant4 initialization
 *  Initialize these globals only after G4RunManager::Initialize()!
 */
#ifndef GlobalData_h
#define GlobalData_h
#include <field_drift/FieldElmerMap.hh>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <deque>

#include <indicators/dynamic_progress.hpp>
#include <indicators/progress_bar.hpp>

#include <G4SystemOfUnits.hh>
#include <G4ThreeVector.hh>
#include <G4Polyline.hh>

#include "utilities/IntegrationInterval.hh"
#include "utilities/PolynomialFit.hh"
#include "GlobalParameters.hh"
#include "geant4/HexagonalMapping.hh"
#include "geant4/PhotonHit.hh"
#include "ArgonPropertiesTables.hh"

class DriftTrack {
public:
  struct driftPoint {
    G4Point3D pos;
    G4ThreeVector field;
    double time;
  };
  std::vector<driftPoint> track;
  DriftTrack() {}
  // Turns out, as of v10 Geant4, G4VVisManager::Draw only works after worker threads have finished
  // (doi:10.1088/1742-6596/513/2/022005 page 7). So electron track has to be saved in Run results.
  void Draw(void) const;
  void Write(std::string filename) const;
};

class DriftElectronInfo {
public:
  DriftElectronInfo() :
    index(-1), position(0, 0, 0), seed_info("?"), track() {}
  int index;
  G4ThreeVector position; // starting position
  std::string seed_info;
  DriftTrack track; // saved only for drawing if enabled
};

struct GeneratedData {
  DriftElectronInfo electron;
  std::deque<PhotonHit> photons;
};

class GlobalData {
public:
  class ProgressBarHelper {
  public:
    ProgressBarHelper();
    void tick(void); // Reduces update rate of progress_bar
    void set_N_events(std::size_t N_events);
    void reset(void);
    indicators::ProgressBar progress_bar;
    bool is_finished(void);
    void set_as_finished(void);
  protected:
    void start(void);
    void finish(void);
    bool has_started;
    bool has_finished;
    std::size_t max_N;
    std::size_t current_N;
    std::mutex mutex_;
  } progress_bar;

  GlobalData();
  ~GlobalData();

  void Initialize(void);
  void SetupFieldMap(void);
  void SetupTHGEM1Mapping(void); //Finds THGEM1 volumes and initializes mapping class gPars::THGEM1_mapping
  G4ThreeVector GetFieldAtGlobal(G4ThreeVector position, DriftMedium* &medium, int *status = nullptr); // Not const but thred-safe
  void PlotField(std::string filename, G4ThreeVector line_start, G4ThreeVector line_finish, int Num, std::string name="", double L_fine=0, int Num_fine=0);
  G4ThreeVector GetDriftStartCenter(void) const;

  HexagonalMapping* THGEM1_mapping; // Set after geometry construction in G4RunManager::Initialize(). Thread safe.
  FieldElmerMap* field_map; // Depends only on gPars::
  DriftMedium* LAr_medium;
  ArgonPropertiesTables Ar_props;

protected:
  // Diffusion coefficients are in mm^2/ns but need to be in mm^(1/2) in ElectronDrift
  // Similar situation was in Garfield++
  static void RecalculateDiffusion(DataVector& diffusion, DataVector& velocity);
};

extern GlobalData gData;

#endif // GlobalData_h
