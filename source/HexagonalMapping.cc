#include "HexagonalMapping.hh"

const G4Transform3D HexagonalMapping::mirror_X = G4ReflectX3D();
const G4Transform3D HexagonalMapping::mirror_Y = G4ReflectY3D();
const G4Transform3D HexagonalMapping::mirror_XY = G4ReflectX3D()*G4ReflectY3D();

HexagonalMapping::HexagonalMapping(G4ThreeVector container_pos, G4ThreeVector cell_pos, G4ThreeVector container_sizes, G4ThreeVector cell_sizes) :
  g_container_pos(container_pos), g_cell_pos(cell_pos)
{
  container_x_size = container_sizes.x();
  container_y_size = container_sizes.y();
  container_z_size = container_sizes.z();
  cell_x_size = cell_sizes.x();
  cell_y_size = cell_sizes.y();
  cell_x_num = (int)(container_x_size / cell_x_size);
  cell_y_num = (int)(container_y_size / cell_y_size);

  tolerance = std::fabs(G4GeometryTolerance::GetInstance()->GetSurfaceTolerance());
  // Further correction is done for convenience only
  // in order to have container center mapping to cell hole center.
  // The reduction in cell coverage area is negligible in practice
  cell_x_num -= (cell_x_num + 2) % 4;
  cell_y_num -= (cell_y_num + 2) % 4;
  cell_x_num = std::max(cell_x_num, 0);
  cell_y_num = std::max(cell_y_num, 0);

  x_min = -cell_x_size * cell_x_num / 2.0;
  x_max = cell_x_size * cell_x_num / 2.0;
  y_min = -cell_y_size * cell_y_num / 2.0;
  y_max = cell_y_size * cell_y_num / 2.0;
}

bool HexagonalMapping::isValid(void) const
{
  return (cell_x_num > 0 && cell_y_num > 0 && cell_x_size > 0 && cell_y_size && container_x_size > 0 && container_y_size > 0);
}

int HexagonalMapping::GetNcells(void) const
{
  return cell_x_num * cell_y_num;
}

std::pair<int, int> HexagonalMapping::GetIndices(const G4ThreeVector &position) const
{
  std::pair<int, int> out(-1, -1);
  G4ThreeVector rel_pos = position - g_container_pos;
  if (rel_pos.x() < x_min || rel_pos.x() > x_max)
    return out;
  if (rel_pos.y() < y_min || rel_pos.y() > y_max)
    return out;
  out.first = (int) ((rel_pos.x() - x_min) / cell_x_size);
  out.second = (int) ((rel_pos.y() - y_min) / cell_y_size);
  if (out.first == cell_x_num)
    out.first--;
  if (out.second == cell_y_num)
    out.second--;
  return out;
}

G4ThreeVector HexagonalMapping::GetCellRelPosition(int x_ind, int y_ind) const
{
  double cell_center_x = x_min + x_ind * cell_x_size + cell_x_size / 2.0;
  double cell_center_y = y_min + y_ind * cell_y_size + cell_y_size / 2.0;
  return G4ThreeVector(cell_center_x, cell_center_y, 0);
}

G4Transform3D HexagonalMapping::GetCellVectorTransform(int x_ind, int y_ind) const
{
  //neighbour cells are asymmetric along both x and y
  int symm_x = (x_ind % 2 ? 1 : -1);
  int symm_y = (y_ind % 2 ? 1 : -1);
  if (symm_x > 0 && symm_y > 0)
    return G4Transform3D::Identity;
  if (symm_x < 0 && symm_y > 0)
    return mirror_X;
  if (symm_x > 0 && symm_y < 0)
    return mirror_Y;
  return mirror_XY;
}

// transformation of global position to global position in cell
G4Transform3D HexagonalMapping::GetCellGlobalPointTransform(int x_ind, int y_ind) const
{
  //cell_point = VectorTransform[(global - g_container_pos - cell_pos)] + g_cell_pos
  G4ThreeVector translation1 = -g_container_pos-GetCellRelPosition(x_ind, y_ind);
  G4ThreeVector translation2 = g_cell_pos;
  G4Transform3D vector_transform = GetCellVectorTransform(x_ind, y_ind);
  return G4Transform3D(G4Translate3D(translation2)) * vector_transform * G4Transform3D(G4Translate3D(translation1)); //should not be very expensive operation
}

// transformation of relative to container position to relative position in cell
G4Transform3D HexagonalMapping::GetCellRelativePointTransform(int x_ind, int y_ind) const
{
  G4ThreeVector translation1 = GetCellRelPosition(x_ind, y_ind);
  G4Transform3D vector_transform = GetCellVectorTransform(x_ind, y_ind);
  return G4Transform3D(G4Translate3D(translation1)) * vector_transform; //should not be very expensive operation
}

std::pair<int, int> HexagonalMapping::GetIndices(int index) const
{
  return std::pair<int, int>(index % cell_x_num, index / cell_x_num);
}

bool HexagonalMapping::isInCell(const HexagonalMappingData& map_info) const
{
  return map_info.isInCell() && (map_info.cell_x_ind < cell_x_num) && (map_info.cell_y_ind < cell_y_num);
}

HexagonalMappingData HexagonalMapping::MapFromCell(const HexagonalMappingData& map_info, bool forced) const
{
  if (!isInCell(map_info)) {
    return MoveFromCell(map_info);
  }
  if (forced)
    return MoveFromCell(map_info);
  G4Transform3D transform = GetCellGlobalPointTransform(map_info.cell_x_ind, map_info.cell_y_ind);
  G4Point3D rel_pos = map_info.position - g_cell_pos;
  if (!forced && ((std::fabs(rel_pos.z() + container_z_size / 2.0) < tolerance && G4ThreeVector(0, 0, -1) * map_info.momentum > 0.0)
      || (std::fabs(rel_pos.z() - container_z_size / 2.0) < tolerance && G4ThreeVector(0, 0, 1) * map_info.momentum > 0.0))) {
    return MoveFromCell(map_info);
  }
  return MoveToNeighbourCell(map_info);
}

HexagonalMappingData HexagonalMapping::MoveToNeighbourCell(const HexagonalMappingData& map_info) const
{
  if (!isInCell(map_info)) {
    return MoveFromCell(map_info);
  }
  HexagonalMappingData out = map_info;
  int index_x_delta = 0, index_y_delta = 0;
  G4Transform3D prev_transform = GetCellGlobalPointTransform(map_info.cell_x_ind, map_info.cell_y_ind).inverse();
  G4Vector3D rel_cell_pos = map_info.position - g_cell_pos;
  rel_cell_pos = (prev_transform * rel_cell_pos); // Transformed as vector, need relative
  //position with respect to world (i.e. with reflection removed) to determine where particle actually moves in terms of indexes

  G4Vector3D global_momentum = (prev_transform *  map_info.momentum); // G4Vector3D and G4Point3D are transformed differently!
  if ((std::abs(rel_cell_pos.x() - cell_x_size / 2.0) < tolerance) && (G4ThreeVector(1, 0, 0)*global_momentum > 0))
    index_x_delta = 1;
  if ((std::abs(rel_cell_pos.x() + cell_x_size / 2.0) < tolerance) && (G4ThreeVector(-1, 0, 0)*global_momentum > 0))
    index_x_delta = -1;
  if ((std::abs(rel_cell_pos.y() - cell_y_size / 2.0) < tolerance) && (G4ThreeVector(0, 1, 0)*global_momentum > 0))
    index_y_delta = 1;
  if ((std::abs(rel_cell_pos.y() + cell_y_size / 2.0) < tolerance) && (G4ThreeVector(0, -1, 0)*global_momentum > 0))
    index_y_delta = -1;
  out.cell_x_ind += index_x_delta;
  out.cell_y_ind += index_y_delta;
  if (out.cell_x_ind == map_info.cell_x_ind && out.cell_y_ind == map_info.cell_y_ind)
    return out;
  if (!isInCell(out)) { //leaving cell-filled area in the container
    return MoveFromCell(map_info);
  }
  G4Transform3D new_transform = GetCellGlobalPointTransform(out.cell_x_ind, out.cell_y_ind);
  out.position = new_transform * (prev_transform * map_info.position); // re-map the same global position to new cell
  out.momentum = new_transform * global_momentum;
  out.polarization = new_transform * (prev_transform * map_info.polarization);
  return out;
}

HexagonalMappingData HexagonalMapping::MapToCell(const HexagonalMappingData& map_info, bool ignore_z) const
{
  HexagonalMappingData out = map_info;
  std::pair<int, int> inds = GetIndices(map_info.position);
  out.cell_x_ind = inds.first;
  out.cell_y_ind = inds.second;
  if (!isInCell(out)) {
    return out;
  }
  G4Transform3D transform = GetCellGlobalPointTransform(out.cell_x_ind, out.cell_y_ind);
  out.position = transform * map_info.position; // G4Vector3D and G4Point3D are transformed differently!
  G4ThreeVector rel_pos = out.position - g_cell_pos;
  if (!ignore_z && (rel_pos.z() < -container_z_size / 2.0 - tolerance || rel_pos.z() > container_z_size / 2.0 + tolerance)) {
    return map_info;
  }
  out.momentum = transform * map_info.momentum;
  out.polarization = transform * map_info.polarization;
  return MoveToNeighbourCell(out); // instead of 'return out;' to remedy edge cases
}

HexagonalMappingData HexagonalMapping::MoveFromCell(const HexagonalMappingData& map_info) const
{
  HexagonalMappingData out;
  out.momentum = map_info.momentum;
  out.polarization = map_info.polarization;
  out.position = G4ThreeVector(- container_x_size / 2.0, - container_y_size / 2.0, - container_z_size / 2.0) + g_container_pos;
  if (!isInCell(map_info)) {
    std::cerr<<"HexagonalMapping::MoveFromCell:Error:"<<std::endl;
    std::cerr<<"\tWrong input mapping info, unknown THGEM cell. Returning default position."<<std::endl;
    return out;
  }
  G4Transform3D transform = GetCellGlobalPointTransform(map_info.cell_x_ind, map_info.cell_y_ind).inverse();
  out.position = transform * map_info.position;
  out.momentum = transform * map_info.momentum;
  out.polarization = transform * map_info.polarization;
  return out;
}
