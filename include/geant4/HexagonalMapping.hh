#ifndef HEXAGONAL_MAPPING_H_
#define HEXAGONAL_MAPPING_H_

#include <G4ThreeVector.hh>
#include <G4RotationMatrix.hh>
#include <G4GeometryTolerance.hh>
#include <G4Transform3D.hh>
#include <G4Point3D.hh>
#include <G4Vector3D.hh>
#include <G4ParallelWorldProcess.hh>
#include <G4TransportationManager.hh>

// Particle(G4Track)-specific data used by global mapping class
// HexagonalMapping.
class HexagonalMappingData {
public:
  HexagonalMappingData() : cell_x_ind(-1), cell_y_ind(-1), mapping_id("")
  {}
  int cell_x_ind;
  int cell_y_ind;
  std::string mapping_id;
  G4Point3D position;
  G4Vector3D momentum;
  G4Vector3D polarization;
  inline bool isInCell() const {
    return (cell_x_ind >=0 && cell_y_ind >=0);
  }
  friend bool operator!=(const HexagonalMappingData& left, const HexagonalMappingData& right)
  {
    if (!left.isInCell() && !right.isInCell())
      return false;
    if (!left.isInCell() && right.isInCell())
      return true;
    if (left.isInCell() && !right.isInCell())
      return true;
    if (left.cell_x_ind == right.cell_x_ind && left.cell_y_ind == right.cell_y_ind
        && left.position == right.position && left.momentum == right.momentum && left.polarization == right.polarization)
      return false;
    return true;
  }
};

class MappingTrigger {
public:
	MappingTrigger();
	MappingTrigger(G4VPhysicalVolume* volume, bool activate_on_entering, std::string id);
	//These two overloaded functions have different logic!
	virtual bool IsSatistied(const G4Track& aTrack, const G4Step& aStep) const;
	virtual bool IsSatistied(const G4VPhysicalVolume* volume, const HexagonalMappingData& old_state) const;
	inline virtual bool IsValid() const {
		return is_set;
	}
protected:
	std::string _id;
	bool is_set;
	bool is_for_entering; // trigger is for leaving cell (volume) if false
	G4VPhysicalVolume* trigger_volume;
public:
	friend class HexagonalMapping;
	friend class MappingManager;
};

// Class calculating mapping between position in THGEM1 dummy volume and its single cell in global coordinates.
// This class is used to teleport particle to single THGEM1 cell to model its behavior there.
// This is done because THGEM typically has ~1e6 cells which is too much for detector parameterization.
// This class does not track geometry changes automatically.

// TODO: When particle enters container from its sides (not x-y planes) it is not mapped to cell
//because cell container is always bigger than cell covered area: (n * cell_size <= full_size) because n is an integer.
// TODO*: Add different symmetry cases: (-,-) (implemented), (+,+), (-,+), (+,-) for (x,y)
// TODO*: It is possible to add arbitrary rotation of cell and cell container
// TODO*: General virtual mapping class can be added.
class HexagonalMapping {
public:
	// At the moment only single field map in the program is available and thus only single HexagonalMapping has
	// associated_with_field_map set to true in MappingManager
  HexagonalMapping(std::string name, G4ThreeVector container_pos, G4ThreeVector cell_pos, G4ThreeVector container_sizes,
  		G4ThreeVector cell_sizes, bool associated_with_field_map = false);

  // Returns cell indices and position inside cell according to global position
  HexagonalMappingData MapToCell(const HexagonalMappingData& map_info, bool ignore_z = false) const;
  // Returns global position and momentum according to cell indices and position inside cell
  HexagonalMappingData MapFromCell(const HexagonalMappingData& map_info, bool forced = false) const;

  void AddTrigger(MappingTrigger trigger);
  void RemoveTrigger(std::string id);
  void ClearTriggers(void);

  bool isValid(void) const;
  int GetNcells(void) const;
  G4ThreeVector GetCellRelPosition(int x_ind, int y_ind) const;
  G4Transform3D GetCellGlobalPointTransform(int x_ind, int y_ind) const; // transformation of global position to global position in cell
  G4Transform3D GetCellRelativePointTransform(int x_ind, int y_ind) const; // transformation of relative to container position to relative position in cell
  std::pair<int, int> GetIndices(const G4ThreeVector &position) const;
  std::pair<int, int> GetIndices(int index) const;
protected:
  std::vector<MappingTrigger> map_triggers;
  G4Transform3D GetCellVectorTransform(int x_ind, int y_ind) const; // transformation of global vector to cell local frame
  HexagonalMappingData MoveToNeighbourCell(const HexagonalMappingData& map_info) const;
  HexagonalMappingData MoveFromCell(const HexagonalMappingData& map_info) const;
  bool isInCell(const HexagonalMappingData& map_info) const;

  std::string name_id;
  bool has_field_map;
  G4ThreeVector g_container_pos;//center of cell container volume (THGEM) in global coordinates
  G4ThreeVector g_cell_pos;//center of cell volume in global coordinates
  double container_x_size;
  double container_y_size;
  double container_z_size; // must be equal to cell z size in geometry construction
  double cell_x_size;
  double cell_y_size;
  int cell_x_num;
  int cell_y_num;
  double tolerance; // for determining when particle hit cell surface

  // For speeding up computations
  double x_min;
  double x_max;
  double y_min;
  double y_max;

  static const G4Transform3D mirror_X;
  static const G4Transform3D mirror_Y;
  static const G4Transform3D mirror_XY;
public:
  friend class MappingManager;
};

#endif //HEXAGONAL_MAPPING_H_
