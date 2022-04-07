#include "GlobalData.hh"

GlobalData gData;

GlobalData::GlobalData() :
    THGEM1_mapping(nullptr),
    field_map(nullptr),
    LAr_medium(nullptr)
{}

GlobalData::~GlobalData()
{
  if (nullptr != THGEM1_mapping)
    delete THGEM1_mapping;
  if (nullptr != field_map)
    delete field_map;
  if (nullptr != LAr_medium)
    delete LAr_medium;
}

void GlobalData::Initialize(void)
{
  results.FindSensorsCoordinates();
  SetupTHGEM1Mapping();
  SetupFieldMap();
}

unsigned int GlobalData::Results::GetNGeneratedPhotons(void) const
{
  unsigned int out = 0;
  for (std::size_t e = 0, e_end_ = generated_photons.size(); e!=e_end_; ++e) {
    out += generated_photons[e].photons.size();
  }
  return out;
}

unsigned int GlobalData::Results::GetNRecordedPhotons(void) const
{
  unsigned int out = 0;
  for (std::size_t e = 0, e_end_ = recorded_photons.size(); e!=e_end_; ++e) {
    out += recorded_photons[e].photons.size();
  }
  return out;
}

void GlobalData::Results::FindSensorsCoordinates(void)
{
  SiPM_positions.clear();
  PMT_positions.clear();
  std::vector<G4PhysicalVolumesSearchScene::Findings> findingsVectorSiPM = FindAllPVs(gPars::det_dims.SiPM_device_name);
  std::vector<G4PhysicalVolumesSearchScene::Findings> findingsVectorPMT = FindAllPVs(gPars::det_dims.PMT_device_name);
  G4int max_number = -1;
  for (const auto& findings: findingsVectorSiPM)
    max_number = std::max(findings.fFoundPVCopyNo, max_number);
  SiPM_positions = std::deque<G4ThreeVector>(max_number+1);
  for (const auto& findings: findingsVectorSiPM)
    SiPM_positions[findings.fFoundPVCopyNo] = findings.fFoundObjectTransformation.getTranslation()/mm;
  max_number = -1;
  for (const auto& findings: findingsVectorPMT)
    max_number = std::max(findings.fFoundPVCopyNo, max_number);
  PMT_positions = std::deque<G4ThreeVector>(max_number+1);
  for (const auto& findings: findingsVectorPMT)
    PMT_positions[findings.fFoundPVCopyNo] = findings.fFoundObjectTransformation.getTranslation()/mm;
}

void GlobalData::SetupTHGEM1Mapping(void)
{
  if (nullptr != THGEM1_mapping)
    delete THGEM1_mapping;
  boost::optional<G4PhysicalVolumesSearchScene::Findings> foundCell = FindSinglePV(gPars::det_dims.THGEM1_cell_name);
  boost::optional<G4PhysicalVolumesSearchScene::Findings> foundTHGEM1 = FindSinglePV(gPars::det_dims.THGEM1_cell_container_name);
  if (boost::none == foundCell) {
    std::cerr<<"RunAction::SetupTHGEM1Mapping:Error:"<<std::endl;
    std::cerr<<"\tTHGEM1 cell is not found. Geometry without THGEM1 mapping is used."<<std::endl;
    return;
  }
  if (boost::none == foundTHGEM1) {
    std::cerr<<"RunAction::SetupTHGEM1Mapping:Error:"<<std::endl;
    std::cerr<<"\tTHGEM1 cell container is not found. Geometry without THGEM1 mapping is used."<<std::endl;
    return;
  }
  G4ThreeVector cell_pos = foundCell->fFoundObjectTransformation.getTranslation();
  G4ThreeVector thgem1_pos = foundTHGEM1->fFoundObjectTransformation.getTranslation();
  G4VSolid* cell_box = foundCell->fpFoundPV->GetLogicalVolume()->GetSolid();
  G4VSolid* thgem1_box = foundTHGEM1->fpFoundPV->GetLogicalVolume()->GetSolid();
  G4ThreeVector bMin, bMax;
  cell_box->BoundingLimits(bMin, bMax);
  G4ThreeVector cell_sizes = bMax - bMin;
  thgem1_box->BoundingLimits(bMin, bMax);
  G4ThreeVector thgem1_sizes = bMax - bMin;
  THGEM1_mapping = new HexagonalMapping(thgem1_pos, cell_pos, thgem1_sizes, cell_sizes);
}

void GlobalData::SetupFieldMap(void)
{
  std::string prfx = gPars::field_map.elmer_mesh_folder;
  field_map = new FieldElmerMap(prfx+"mesh.header", prfx+"mesh.elements",
      prfx+"mesh.nodes",gPars::field_map.elmer_solution_filename, "mm");
  if (!field_map->IsReady()) {
    G4Exception("GlobalData::SetupFieldMap: ",
        "InvalidSetup", FatalException, "Failed to prepare field map.");
    return;
  }
  double tolerance = gPars::general.surface_tolerance;
  double x_min, x_max, y_min, y_max, z_min, z_max;
  field_map->GetBoundingBox(x_min, y_min, z_min, x_max, y_max, z_max);
  if ((gPars::det_dims.THGEM1_hole_pitch / 2.0 - (x_max- x_min)) > tolerance ||
      (gPars::det_dims.THGEM1_hole_pitch * std::sqrt(3) / 2.0 - (y_max- y_min)) > tolerance ||
      (gPars::det_dims.THGEM1_container_width > (z_max- z_min))) {
    G4Exception("GlobalData::SetupFieldMap: ",
        "InvalidSetup", FatalException, "Loaded field map has dimensions incomparable to THGEM1 cell.");
    return;
  }
  field_map->SetRelativeTolerance(gPars::field_map.mesh_tolerance);

  DataVector LAr_drift;
  std::ifstream str;
  str.open(gPars::field_map.LAr_drift_velocity);
  if (!str.is_open()) {
    G4Exception("GlobalData::SetupFieldMap: ",
      "InvalidSetup", FatalException, "Failed to open LAr drift velocity file");
  }
  LAr_drift.read(str);
  LAr_drift.scaleXY(1e3*volt / cm, cm / second); // Check units in data files
  LAr_drift.use_leftmost(true);
  LAr_drift.use_rightmost(true);
  str.close();

  DataVector LAr_longitudinal;
  str.open(gPars::field_map.LAr_diffusion_longitudinal);
  if (!str.is_open()) {
    G4Exception("GlobalData::SetupFieldMap: ",
      "InvalidSetup", JustWarning, "Failed to open LAr longitudinal diffusion file");
  } else {
    LAr_longitudinal.read(str);
    LAr_longitudinal.scaleXY(1e3*volt / cm, cm * cm / second); // Check units in data files
    RecalculateDiffusion(LAr_longitudinal, LAr_drift);
    LAr_longitudinal.use_leftmost(true);
    LAr_longitudinal.use_rightmost(true);
  }

  DataVector LAr_transversal;
  str.open(gPars::field_map.LAr_diffusion_transversal);
  if (!str.is_open()) {
    G4Exception("GlobalData::SetupFieldMap: ",
      "InvalidSetup", JustWarning, "Failed to open LAr transversal diffusion file");
  } else {
    LAr_transversal.read(str);
    LAr_transversal.scaleXY(1e3*volt / cm, cm * cm / second); // Check units in data files
    RecalculateDiffusion(LAr_transversal, LAr_drift);
    LAr_transversal.use_leftmost(true);
    LAr_transversal.use_rightmost(true);
  }

  LAr_medium = new DriftMedium("LAr", LAr_drift, LAr_longitudinal, LAr_transversal);
  LAr_medium->SetDriftable(true);
  field_map->SetMedium(0, LAr_medium);
}

G4ThreeVector GlobalData::GetFieldAtGlobal(G4ThreeVector position, DriftMedium* &medium, int *status)
{
  medium = nullptr;
  if (nullptr == field_map || !field_map->IsReady()) {
    std::cerr<<"GlobalData::GetFieldAtGlobal:Error: Field map is not available."<<std::endl;
    if (nullptr != status)
      *status = 1;
    return G4ThreeVector(0, 0, 0);
  }
  if (nullptr == THGEM1_mapping || !THGEM1_mapping->isValid()) {
    std::cerr<<"GlobalData::GetFieldAtGlobal:Error: THGEM1 mapping is not available."<<std::endl;
    if (nullptr != status)
      *status = 2;
    return G4ThreeVector(0, 0, 0);
  }
  HexagonalMappingData mapping_data;
  mapping_data.cell_x_ind = -1;
  mapping_data.cell_y_ind = -1;
  mapping_data.position = position;
  mapping_data.momentum = G4ThreeVector(0, 0, 0);
  mapping_data.polarization = G4ThreeVector(0, 0, 0);
  mapping_data = THGEM1_mapping->MapToCell(mapping_data, true);
  if (!mapping_data.isInCell()) {
    if (nullptr != status)
        *status = -6;
    return mapping_data.momentum;
  }
  G4ThreeVector in_mesh_pos = mapping_data.position - gPars::det_dims.THGEM1_single_cell_position + gPars::field_map.elmer_mesh_center;
  double ex,ey,ez;
  int status_;
  field_map->ElectricField(in_mesh_pos.x(),in_mesh_pos.y(),in_mesh_pos.z(),ex,ey,ez,medium, status_);
  mapping_data.momentum = G4ThreeVector(ex, ey, ez);
  mapping_data = THGEM1_mapping->MapFromCell(mapping_data, true);
  if (nullptr != status)
    *status = status_;
  return mapping_data.momentum;
}

void GlobalData::PlotField(std::string filename, G4ThreeVector line_start, G4ThreeVector line_finish, int Num, std::string name, double L_fine, int Num_fine)
{
  if (nullptr == field_map || !field_map->IsReady()) {
    std::cerr<<"GlobalData::PlotField:Error: Field map is not available. No plotting."<<std::endl;
    return;
  }
  if (nullptr == THGEM1_mapping || !THGEM1_mapping->isValid()) {
    std::cerr<<"GlobalData::PlotField:Error: THGEM1 mapping is not available."<<std::endl;
    return;
  }
  ensure_file(filename);
  double dphi=0;
  G4ThreeVector diff = line_finish - line_start;
  double L = std::sqrt(diff * diff);
  double s=0,x=line_start.x(),y=line_start.y(),z=line_start.z();
  double ex,ey,ez,ee;
  DriftMedium * medium;
  int status;
  int LN=(Num-Num_fine)/2.0;
  int RN=Num-LN;
  double ds=0;
  if ((L_fine<=0)||(Num_fine>=Num)||(L_fine>=s)||(RN<LN)||(Num_fine<=0)||(LN>Num)||(LN<0)) {
    Num_fine=0;
    L_fine=0;
    LN=Num;
    RN=-Num;
  }
  std::ofstream flll(filename,std::ios_base::out);
  //flll<<"//z\tx\ty\tL\tEabs\tEx\tEy\tEz"
  for (int hh=0; hh<=Num; hh++) {//hh<=Num in order to cover [x0;x1], not [x0;x1)
     x=(line_finish.x()-line_start.x())*s/L+line_start.x();
     y=(line_finish.y()-line_start.y())*s/L+line_start.y();
     z=(line_finish.z()-line_start.z())*s/L+line_start.z();
     G4ThreeVector E = GetFieldAtGlobal(G4ThreeVector(x, y, z), medium, &status);
     ex = E.x(); ey = E.y(); ez = E.z();
     ee=std::sqrt(ex*ex+ey*ey+ez*ez);
     if (ee>0) {//there was a bug: because of -0.0 comparison in elmer component
     //resulting in appearing of 0 field in arbitrary points - fixed by me, yet better to leave this condition
        if((line_finish.x()==line_start.x())&&(line_finish.y()==line_start.y()))
           flll<<z / mm<<'\t';
        else if ((line_finish.z()==line_start.z())&&(line_finish.y()==line_start.y()))
           flll<<x / mm<<'\t';
        else if ((line_finish.x()==line_start.x())&&(line_finish.z()==line_start.z()))
           flll<<y / mm<<'\t';
        else
           flll<<s<<'\t';
        flll<<ee * cm / volt <<'\t'<<ex * cm / volt<<'\t'<<ey * cm / volt<<'\t'<<ez * cm / volt<<std::endl;
     } else {
       std::cerr<<"PrimaryGeneratorAction::PlotField:Warning:"<<std::endl;
       std::cerr<<"0 field at "<<G4ThreeVector(x, y, z)<<std::endl;
     }
     if ((hh<LN)||(hh>RN)) {
        ds=(L-L_fine)/(Num-LN-RN-1);
     } else {
        ds=L_fine/(RN-LN+1);
     }
     dphi+=(ds/L)*(ex*(line_finish.x()-line_start.x())+ey*(line_finish.y()-line_start.y())+ez*(line_finish.z()-line_start.z()));
     s+=ds;
  }
  if (name.size())
      std::cout<<"Int{E(r)*dr} [V] on "<<name<<" = "<<dphi / volt<<std::endl;
  flll.close();
}

void GlobalData::AddToFile(std::deque<GeneratedData> data, std::string filename)
{
  std::ofstream str;
  open_output_file(filename, str, std::ios_base::ate|std::ios_base::out);
  if (!str.is_open()) {
    std::cerr<<"AddToFile:Error: Failed to open file"<<std::endl;
    std::cerr<<"\t\""<<filename<<"\""<<std::endl;
    return;
  }
  for (std::size_t e = 0, e_end_ = data.size(); e!=e_end_; ++e) {
    str<<data[e].electron.index<<"\t"<<data[e].photons.size()<<"\t"<<data[e].electron.position.x() / mm<<"\t"
        <<data[e].electron.position.y() / mm<<"\t"<<data[e].electron.position.z() / mm<<"\t"
        <<data[e].electron.seed_info<<std::endl;
    for (std::size_t p = 0, p_end_ = data[e].photons.size(); p!=p_end_; ++p) {
      str<<data[e].photons[p]._energy / eV<<"\t"<<data[e].photons[p]._pos.x() / mm<<"\t"
          <<data[e].photons[p]._pos.y() / mm<<"\t"<<data[e].photons[p]._pos.z() / mm<<"\t"
          <<data[e].photons[p]._time / us<<"\t"<<data[e].photons[p]._momentum.x()<<"\t"
          <<data[e].photons[p]._momentum.y()<<"\t"<<data[e].photons[p]._momentum.z()<<std::endl;
    }
  }
}

void GlobalData::RecalculateDiffusion(DataVector& diffusion, DataVector& velocity)
{
  if (!diffusion.isValid() || !velocity.isValid()) {
    std::cerr<<"GlobalData::RecalculateDiffusion: Invalid data."<<std::endl;
  }
  for (std::size_t i = 0, i_end_ = diffusion.size(); i!=i_end_; ++i) {
    diffusion[i].second = std::sqrt(2.0 * diffusion[i].second / velocity(diffusion[i].first));
  }
}
