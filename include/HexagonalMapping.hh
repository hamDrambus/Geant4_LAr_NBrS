#ifndef HEXAGONAL_MAPPING_H_
#define HEXAGONAL_MAPPING_H_

#include "G4ThreeVector.hh"
#include "G4GeometryTolerance.hh"

// Particle(G4Track)-specific data used by global mapping class
// HexagonalMapping.
class HexagonalMappingData {
public:
  HexagonalMappingData() : cell_x_ind(-1), cell_y_ind(-1)
  {}
  int cell_x_ind;
  int cell_y_ind;
  G4ThreeVector position;
  G4ThreeVector momentum;
  G4ThreeVector polarization;
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



// Class calculating mapping between position in THGEM1 dummy volume and its single cell in global coordinates.
// This class is used to teleport particle to single THGEM1 cell to model its behavior there.
// This is done because THGEM typically has ~1e6 cells which is too much for detector parameterization.
// This class does not track geometry changes automatically.

// TODO: When particle enters container from its sides (not x-y planes) it is not mapped to cell
//because cell container is always bigger that cell covered area: (n * cell_size <= full_size) because n is an integer.
// TODO*: Add different symmetry cases: (-,-) (implemented), (+,+), (-,+), (+,-) for (x,y)
// TODO*: It is possible to add arbitrary rotation of cell and cell container
class HexagonalMapping {
public:
  HexagonalMapping(G4ThreeVector container_pos, G4ThreeVector cell_pos, G4ThreeVector container_sizes, G4ThreeVector cell_sizes);

  // Returns cell indices and position inside cell according to global position
  HexagonalMappingData MapToCell(const HexagonalMappingData& map_info, bool ignore_z = false) const;
  // Returns global position and momentum according to cell indices and position inside cell
  HexagonalMappingData MapFromCell(const HexagonalMappingData& map_info, bool ignore_z = false) const;
  // Propagates particle to the next cell. Cell is changed when particle hit current cell sides only.
  HexagonalMappingData MapToNeighbourCell(const HexagonalMappingData& map_info) const;

  bool isValid(void) const;

protected:
  HexagonalMappingData GetIndices(const G4ThreeVector &position) const;

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
};

#endif //HEXAGONAL_MAPPING_H_
