#include "HexagonalMapping.hh"

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

HexagonalMappingData HexagonalMapping::GetIndices(const G4ThreeVector &position) const
{
  HexagonalMappingData out;
  G4ThreeVector rel_pos = position - g_container_pos;
  if (rel_pos.x() < x_min || rel_pos.x() > x_max)
    return out;
  if (rel_pos.y() < y_min || rel_pos.y() > y_max)
    return out;
  out.cell_x_ind = (int) ((rel_pos.x() - x_min) / cell_x_size);
  out.cell_y_ind = (int) ((rel_pos.y() - y_min) / cell_y_size);
  if (out.cell_x_ind == cell_x_num)
    out.cell_x_ind--;
  if (out.cell_y_ind == cell_y_num)
    out.cell_y_ind--;
  return out;
}

HexagonalMappingData HexagonalMapping::MapToCell(const HexagonalMappingData& map_info, bool ignore_z) const
{
  HexagonalMappingData out = map_info;
  G4ThreeVector rel_pos = map_info.position - g_container_pos;
  if (!ignore_z && (rel_pos.z() < -container_z_size / 2.0 || rel_pos.z() > container_z_size / 2.0)) {
    return out;
  }
  out = GetIndices(map_info.position);
  out.momentum = map_info.momentum;
  out.polarization = map_info.polarization;
  out.position = map_info.position;
  if (!out.isInCell() || out.cell_x_ind >= cell_x_num || out.cell_y_ind >= cell_y_num) {
    return out;
  }

  double cell_center_x = x_min + out.cell_x_ind * cell_x_size + cell_x_size / 2.0;
  double cell_center_y = y_min + out.cell_y_ind * cell_y_size + cell_y_size / 2.0;
  double rel_to_cell_pos_x = (rel_pos.x() - cell_center_x) * (out.cell_x_ind % 2 ? 1 : -1); //neighbour cells are asymmetric along both x and y
  double rel_to_cell_pos_y = (rel_pos.y() - cell_center_y) * (out.cell_y_ind % 2 ? 1 : -1);
  rel_pos.setX(rel_to_cell_pos_x);
  rel_pos.setY(rel_to_cell_pos_y);
  out.position = rel_pos + g_cell_pos;
  out.momentum.setX(out.momentum.x() * (out.cell_x_ind % 2 ? 1 : -1));
  out.momentum.setY(out.momentum.y() * (out.cell_y_ind % 2 ? 1 : -1));
  out.polarization.setX(out.polarization.x() * (out.cell_x_ind % 2 ? 1 : -1));
  out.polarization.setY(out.polarization.y() * (out.cell_y_ind % 2 ? 1 : -1));
  return out;
}

HexagonalMappingData HexagonalMapping::MapFromCell(const HexagonalMappingData& map_info, bool ignore_z) const
{
  HexagonalMappingData out;
  out.momentum = map_info.momentum;
  out.polarization = map_info.polarization;
  out.position = G4ThreeVector(- container_x_size / 2.0, - container_y_size / 2.0, - container_z_size / 2.0);
  out.position += g_container_pos;
  if (!map_info.isInCell() || map_info.cell_x_ind >= cell_x_num || map_info.cell_y_ind >= cell_y_num) {
    std::cerr<<"HexagonalMapping::MapFromCell:Error:"<<std::endl;
    std::cerr<<"\tWrong input mapping info, unknown THGEM cell. Returning default position."<<std::endl;
    return out;
  }
  G4ThreeVector rel_pos = map_info.position - g_cell_pos;
  if (!ignore_z && (rel_pos.z() < -container_z_size / 2.0 - tolerance || rel_pos.z() > container_z_size / 2.0 + tolerance)) {
    std::cout<<"HexagonalMapping::MapFromCell:Warning:"<<std::endl;
      std::cout<<"\tZ coordinate lies outside cell bounds. ignore_z forced to true."<<std::endl;
    ignore_z = true;
  }

  // Exact inverse of HexagonalMapping::MapToCell
  // rel_pos.x() here is rel_to_cell_pos_x in HexagonalMapping::MapToCell
  // and rel_to_container_pos_x here is rel_pos.x() there.
  double cell_center_x = x_min + map_info.cell_x_ind * cell_x_size + cell_x_size / 2.0;
  double cell_center_y = y_min + map_info.cell_y_ind * cell_y_size + cell_y_size / 2.0;
  double rel_to_container_pos_x = (rel_pos.x() * (map_info.cell_x_ind % 2 ? 1 : -1)) + cell_center_x;
  double rel_to_container_pos_y = (rel_pos.y() * (map_info.cell_y_ind % 2 ? 1 : -1)) + cell_center_y;
  rel_pos.setX(rel_to_container_pos_x);
  rel_pos.setY(rel_to_container_pos_y);
  out.position = rel_pos + g_container_pos;
  out.momentum.setX(out.momentum.x() * (out.cell_x_ind % 2 ? 1 : -1));
  out.momentum.setY(out.momentum.y() * (out.cell_y_ind % 2 ? 1 : -1));
  out.polarization.setX(out.polarization.x() * (out.cell_x_ind % 2 ? 1 : -1));
  out.polarization.setY(out.polarization.y() * (out.cell_y_ind % 2 ? 1 : -1));
  return out;
}

HexagonalMappingData HexagonalMapping::MapToNeighbourCell(const HexagonalMappingData& map_info) const
{
  if (!map_info.isInCell() || map_info.cell_x_ind >= cell_x_num || map_info.cell_y_ind >= cell_y_num) {
    std::cerr<<"HexagonalMapping::MapToNeighbourCell:Error:"<<std::endl;
    std::cerr<<"\tWrong input mapping info, unknown THGEM cell. Returning default position."<<std::endl;
    HexagonalMappingData out;
    out.position = G4ThreeVector(- container_x_size / 2.0, - container_y_size / 2.0, - container_z_size / 2.0);
    out.position += g_container_pos;
    return out;
  }
  HexagonalMappingData out = map_info;
  G4ThreeVector rel_pos = map_info.position - g_cell_pos;
  int index_x_delta = 0;
  int index_y_delta = 0;

  if ((std::abs(rel_pos.x() - cell_x_size / 2.0) < tolerance) && (G4ThreeVector(1, 0, 0)*map_info.momentum > 0))
    index_x_delta = 1;
  if ((std::abs(rel_pos.x() + cell_x_size / 2.0) < tolerance) && (G4ThreeVector(-1, 0, 0)*map_info.momentum > 0))
    index_x_delta = -1;
  if ((std::abs(rel_pos.y() - cell_y_size / 2.0) < tolerance) && (G4ThreeVector(0, 1, 0)*map_info.momentum > 0))
    index_y_delta = 1;
  if ((std::abs(rel_pos.y() + cell_y_size / 2.0) < tolerance) && (G4ThreeVector(0, -1, 0)*map_info.momentum > 0))
    index_y_delta = -1;
  out.cell_x_ind += index_x_delta;
  out.cell_y_ind += index_y_delta;
  if (!out.isInCell() || out.cell_x_ind >= cell_x_num || out.cell_y_ind >= cell_y_num) { //leaving cell-filled area in the container
    out = MapFromCell(map_info, true);
    return out;
  }
  if (index_x_delta != 0) { // Coordinate in cell is not changed due to asymmetry
    out.momentum.setX(-out.momentum.x());
    out.polarization.setX(-out.polarization.x());
  }
  if (index_y_delta != 0) {
    out.momentum.setY(-out.momentum.y());
    out.polarization.setY(-out.polarization.y());
  }
  return out;
}
