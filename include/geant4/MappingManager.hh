#ifndef MAPPING_MANAGER_H_
#define MAPPING_MANAGER_H_

#include <G4ThreeVector.hh>
#include <G4Point3D.hh>
#include <G4Vector3D.hh>
#include <G4ParallelWorldProcess.hh>
#include <G4TransportationManager.hh>
#include "HexagonalMapping.hh"
#include "field_drift/FieldElmerMap.hh"

// This manager class is necessary because more than a single container(THGEM plate)<->THGEM cell mapping
// is necessary if NBrS is measured from THGEM0 instead of THGEM1.
// All mapping (in ParticleGenerators and TeleporationProcess) is done through calls to this class.
// Does not allow for nested mappings (did not bother to implement).
class MappingManager {
public:
	MappingManager();
	std::vector<HexagonalMapping> mappings;
	void AddMapping(HexagonalMapping& mapping);
	void RemoveMapping(std::string name_id);
  bool HasMapping(void) const {
  	for (auto i : mappings) {
  		if (i.isValid())
  			return true;
  	}
    return false;
  }
  bool HasFieldMapping(void) const {
		for (auto i : mappings) {
			if (i.isValid() && i.has_field_map)
				return true;
		}
		return false;
	}
  // These two overloaded functions have different logic!
	// Called in TeleportationProcess (during particle transport)
	HexagonalMappingData GetNewState(const G4Track& aTrack, const G4Step& aStep, const HexagonalMappingData& old_state) const;
	// Called in PrimaryGenerator (during particle creation)
	HexagonalMappingData GetNewState(const G4VPhysicalVolume* volume, const HexagonalMappingData& old_state) const;

	G4ThreeVector GetFieldAtGlobal(G4ThreeVector position, DriftMedium* &medium, int *status) const;
};

#endif //MAPPING_MANAGER_H_
