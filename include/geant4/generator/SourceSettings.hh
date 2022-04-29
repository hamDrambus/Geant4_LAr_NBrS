#ifndef SourceSettings_h
#define SourceSettings_h

#include <string>
#include <sstream>
#include <limits>
#include <globals.hh>

#include <G4SystemOfUnits.hh>
#include <G4ThreeVector.hh>

#include <utilities/GlobalUtilities.hh>

// Class containing input source settings.
// All parameters are in Geant4 units
class VSourceSettings
{
protected:
  VSourceSettings() {}
public:
  enum GeneratorType {
    ElectronPatterns,
    NBrS,
    PhotonsDirectly
  };
  GeneratorType generator_type;

  unsigned int N_events;
  double x_center;
  double y_center;
  double z_center;
  double xy_radius;
  double z_width;
};

VSourceSettings* CreateSourceSettings(std::string filename);

#endif // SourceSettings_h
