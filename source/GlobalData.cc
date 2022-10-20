#include "GlobalData.hh"

GlobalData gData;

// Turns out, as of v10 Geant4, G4VVisManager::Draw only works after worker threads have finished
// (doi:10.1088/1742-6596/513/2/022005 page 7). So electron track has to be saved in Run results.
void DriftTrack::Draw(void) const
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager && !track.empty()) {
    G4Polyline track_obj;
    track_obj.reserve(track.size());
    for (std::size_t i = 0, i_end_ = track.size(); i!=i_end_; ++i)
      track_obj.push_back(track[i].pos);
    G4VisAttributes attribs;
    attribs.SetColour(1.0, 0.55, 0.0, 0.8);
    attribs.SetLineStyle(G4VisAttributes::unbroken);
    attribs.SetLineWidth(0.4);
    track_obj.SetVisAttributes(attribs);
    pVVisManager->Draw(track_obj);
  }
}

void DriftTrack::Write(std::string filename) const
{
  std::ofstream str;
  open_output_file(filename, str, std::ios_base::trunc);
  if (!str.is_open()) {
    std::cerr << "DriftTrack::WriteDriftTrack:Error: Failed to open file." << std::endl;
    return;
  }
  double ee = 0.0;
  double l = 0.0; // Full path
  str << "//L[mm]\tx[mm]\ty[mm]\tz[mm]\tex[kV/cm]\tey[kV/cm]\tez[kV/cm]\tE[kV/cm]\tt[us]" << std::endl;
  for (std::size_t i = 0, i_end_ = track.size(); i!=i_end_; ++i) {
    ee = sqrt(pow(track[i].field.x(), 2) + pow(track[i].field.y(), 2) + pow(track[i].field.z(), 2));
    str << l / mm << "\t" << track[i].pos.x() / mm << "\t"<< track[i].pos.y() / mm << "\t"
        << track[i].pos.z() / mm << "\t" << track[i].field.x() * cm / kilovolt << "\t"
        << track[i].field.y() * cm / kilovolt << "\t" << track[i].field.z() * cm / kilovolt << "\t"
        << ee * cm / kilovolt << "\t" << track[i].time / us << std::endl;
    G4ThreeVector dx = (i == i_end_ - 1) ? G4ThreeVector(0,0,0) : track[i+1].pos - track[i].pos;
    l += sqrt(dx*dx);
  }
}

GlobalData::ProgressBarHelper::ProgressBarHelper() :
    progress_bar(indicators::option::BarWidth{50},
        indicators::option::Fill{"█"},
        indicators::option::Lead{"█"},
        indicators::option::Remainder{"-"},
        indicators::option::ShowPercentage(true),
        indicators::option::ShowElapsedTime(true),
        indicators::option::ShowRemainingTime(true),
        indicators::option::PrefixText{" Completion:"}),
    max_N(0), current_N(0), has_started(false), has_finished(false)
{}

void GlobalData::ProgressBarHelper::tick(void)
{
  std::lock_guard<std::mutex> guard(mutex_);
  if (has_finished)
    return;
  if (!has_started)
    start();
  ++current_N;
  std::size_t current_rate = progress_bar.current();
  double fraction = (double) current_N / max_N;
  std::size_t new_rate = (std::size_t)(100 * fraction);
  while (new_rate > current_rate) {
    if (!progress_bar.is_completed())
      progress_bar.tick();
    current_rate = progress_bar.current();
  }
  if (current_N == max_N)
    finish();
}

void GlobalData::ProgressBarHelper::set_N_events(std::size_t N_events)
{
  std::lock_guard<std::mutex> guard(mutex_);
  if (!has_started) {
    max_N = N_events;
    current_N = 0;
  }
}

bool GlobalData::ProgressBarHelper::is_finished(void)
{
  std::lock_guard<std::mutex> guard(mutex_);
  return has_finished;
}

void GlobalData::ProgressBarHelper::set_as_finished(void)
{
  std::lock_guard<std::mutex> guard(mutex_);
  has_started = true;
  has_finished = true;
  if (!progress_bar.is_completed())
    progress_bar.mark_as_completed();
  current_N = max_N;
}

void GlobalData::ProgressBarHelper::reset(void)
{
  std::lock_guard<std::mutex> guard(mutex_);
  if (has_started && has_finished) {
    has_started = false;
    has_finished = false;
    current_N = 0;
    progress_bar.set_progress(0);
  }
}

void GlobalData::ProgressBarHelper::start(void)
{
  if (!has_started && !has_finished) {
    has_started = true;
    current_N = 0;
    progress_bar.set_progress(0);
  }
}

void GlobalData::ProgressBarHelper::finish(void)
{
  if (has_started && !has_finished) {
    has_finished = true;
    current_N = max_N;
    if (!progress_bar.is_completed())
      progress_bar.mark_as_completed();
  }
}


GlobalData::GlobalData() :
    field_map(nullptr),
    LAr_medium(nullptr),
    progress_bar()
{}

GlobalData::~GlobalData()
{
  if (nullptr != field_map)
    delete field_map;
  if (nullptr != LAr_medium)
    delete LAr_medium;
}

void GlobalData::Initialize(void)
{
  Ar_props.Initialize();
  SetupFieldMap();
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
  double tolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
  double x_min, x_max, y_min, y_max, z_min, z_max;
  // TODO: connect elmer field map and detector's mapping properly.
  // Maybe move field loading as a function in Detector construction?
  // Or use MappingManager's mapping with has_field_map=true?
  field_map->GetBoundingBox(x_min, y_min, z_min, x_max, y_max, z_max);
  if (gPars::det_dims->detector_type == VDetectorDimensions::Full_detector_y2022) {
  	DetectorDimsFullY2022 *dims = (DetectorDimsFullY2022*) gPars::det_dims;
  	if ((dims->THGEM0_hole_pitch / 2.0 - (x_max- x_min)) > tolerance ||
				(dims->THGEM0_hole_pitch * std::sqrt(3) / 2.0 - (y_max- y_min)) > tolerance ||
				(dims->THGEM0_container_width > (z_max- z_min))) {
			G4Exception("GlobalData::SetupFieldMap: ",
					"InvalidSetup", FatalException, "Loaded field map has dimensions incomparable to THGEM0 cell.");
			return;
		}
  } else {
		if ((gPars::det_dims->THGEM1_hole_pitch / 2.0 - (x_max- x_min)) > tolerance ||
				(gPars::det_dims->THGEM1_hole_pitch * std::sqrt(3) / 2.0 - (y_max- y_min)) > tolerance ||
				(gPars::det_dims->THGEM1_container_width > (z_max- z_min))) {
			G4Exception("GlobalData::SetupFieldMap: ",
					"InvalidSetup", FatalException, "Loaded field map has dimensions incomparable to THGEM1 cell.");
			return;
		}
  }
  field_map->SetRelativeTolerance(gPars::field_map.mesh_tolerance);

  DataVector* velocity = nullptr;
  if (!Ar_props.drift_velocity.isValid()) {
    velocity = new DataVector();
    std::ifstream str;
    str.open(gPars::Ar_props.LAr_drift_velocity);
    if (!str.is_open()) {
      G4Exception("GlobalData::SetupFieldMap: ",
        "InvalidSetup", FatalException, "Failed to open LAr drift velocity file");
    }
    velocity->read(str);
    velocity->scaleXY(1e3*volt / cm, cm / second); // Check units in data files
    velocity->use_leftmost(true);
    velocity->use_rightmost(true);
    str.close();
  } else {
    velocity = &Ar_props.drift_velocity;
  }

  DataVector* diff_longitudinal = nullptr;
  DataVector* diff_transversal = nullptr;
  if (!Ar_props.drift_diffusion_L.isValid() && !Ar_props.drift_diffusion_T.isValid()) {
    diff_longitudinal = new DataVector();
    diff_transversal = new DataVector();
    std::ifstream str;
    str.open(gPars::Ar_props.LAr_diffusion_longitudinal);
    if (!str.is_open()) {
      G4Exception("GlobalData::SetupFieldMap: ",
        "InvalidSetup", JustWarning, "Failed to open LAr longitudinal diffusion file");
    } else {
      diff_longitudinal->read(str);
      diff_longitudinal->scaleXY(1e3*volt / cm, cm * cm / second); // Check units in data files
      diff_longitudinal->use_leftmost(true);
      diff_longitudinal->use_rightmost(true);
    }

    str.open(gPars::Ar_props.LAr_diffusion_transversal);
    if (!str.is_open()) {
      G4Exception("GlobalData::SetupFieldMap: ",
        "InvalidSetup", JustWarning, "Failed to open LAr transversal diffusion file");
    } else {
      diff_transversal->read(str);
      diff_transversal->scaleXY(1e3*volt / cm, cm * cm / second); // Check units in data files
      diff_transversal->use_leftmost(true);
      diff_transversal->use_rightmost(true);
    }
  } else {
    if (Ar_props.drift_diffusion_L.isValid() && Ar_props.drift_diffusion_T.isValid()) {
      diff_longitudinal = &Ar_props.drift_diffusion_L;
      diff_transversal = &Ar_props.drift_diffusion_T;
    } else { // In only one type of diffusion is not available then use the same for both.
      if (Ar_props.drift_diffusion_L.isValid())
        diff_transversal = new DataVector(*(diff_longitudinal = &Ar_props.drift_diffusion_L));
      else
        diff_longitudinal = new DataVector(*(diff_transversal= &Ar_props.drift_diffusion_T));
    }
  }
  RecalculateDiffusion(*diff_longitudinal, *velocity);
  RecalculateDiffusion(*diff_transversal, *velocity);
  LAr_medium = new DriftMedium("LAr", *velocity, *diff_longitudinal, *diff_transversal);
  LAr_medium->SetDriftable(true);
  field_map->SetMedium(0, LAr_medium);
  gPars::det_dims->drift_start_center = GetDriftStartCenter();
}

G4ThreeVector GlobalData::GetDriftStartCenter(void) const
{
  G4ThreeVector c = gPars::field_map.elmer_mesh_center;
  if (!field_map)
    return c;
  double xmin, ymin, zmin, xmax, ymax, zmax;
  field_map->GetBoundingBox(xmin, ymin, zmin, xmax, ymax, zmax);
  if (gPars::det_dims->detector_type == VDetectorDimensions::Full_detector_y2022) {
    	DetectorDimsFullY2022 *dims = (DetectorDimsFullY2022*) gPars::det_dims;
    	return dims->THGEM0_center + G4ThreeVector(0, 0, c.z() + zmin);
  }
  return gPars::det_dims->THGEM1_center + G4ThreeVector(0, 0, c.z() + zmin);
}

G4ThreeVector GlobalData::GetFieldAtGlobal(G4ThreeVector position, DriftMedium* &medium, int *status)
{
  return mapping_manager.GetFieldAtGlobal(position, medium, status);
}

void GlobalData::PlotField(std::string filename, G4ThreeVector line_start, G4ThreeVector line_finish, int Num, std::string name, double L_fine, int Num_fine)
{
  if (nullptr == field_map || !field_map->IsReady()) {
    std::cerr<<"GlobalData::PlotField:Error: Field map is not available. No plotting."<<std::endl;
    return;
  }
  if (!mapping_manager.HasFieldMapping()) {
    std::cerr<<"GlobalData::PlotField:Error: THGEM1 field mapping is not available."<<std::endl;
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

void GlobalData::RecalculateDiffusion(DataVector& diffusion, DataVector& velocity)
{
  if (!diffusion.isValid() || !velocity.isValid()) {
    std::cerr<<"GlobalData::RecalculateDiffusion: Invalid data."<<std::endl;
  }
  for (std::size_t i = 0, i_end_ = diffusion.size(); i!=i_end_; ++i) {
    diffusion[i].second = std::sqrt(2.0 * diffusion[i].second / velocity(diffusion[i].first));
  }
}
