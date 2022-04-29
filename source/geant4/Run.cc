#include <geant4/Run.hh>
#include <geant4/generator/GenNBrS_InTHGEM.hh>

Run::Run() : hit_collection_ID(-1)
{
  results.ClearAndInit();
}

Run::~Run()
{}

std::size_t Run::Results::GetNGeneratedElectrons(void) const
{
  return generated_photons.size();
}

std::size_t Run::Results::GetNGeneratedPhotons(void) const
{
  std::size_t out = 0;
  for (std::size_t e = 0, e_end_ = generated_photons.size(); e!=e_end_; ++e) {
    out += generated_photons[e].photons.size();
  }
  return out;
}

std::size_t Run::Results::GetNRecordedPhotons(void) const
{
  std::size_t out = 0;
  for (std::size_t e = 0, e_end_ = recorded_photons.size(); e!=e_end_; ++e) {
    out += recorded_photons[e].photons.size();
  }
  return out;
}

std::size_t Run::Results::GetNRecoredPMTsRAW(void) const
{
  std::size_t out = 0;
  for (std::size_t e = 0, e_end_ = PMT_photon_n.size(); e!=e_end_; ++e) {
    out += PMT_photon_n[e];
  }
  return out;
}

std::size_t Run::Results::GetNRecoredSiPMsRAW(void) const
{
  std::size_t out = 0;
  for (std::size_t e = 0, e_end_ = SiPM_photon_n.size(); e!=e_end_; ++e) {
    out += SiPM_photon_n[e];
  }
  return out;
}

std::size_t Run::Results::GetNRecoredPMTavg(void) const
{
  if (PMT_photon_n.size() != 4)
    return 0;
  std::size_t PMT3 = PMT_photon_n[0] + PMT_photon_n[2] + PMT_photon_n[3];
  double grid_fraction = (double) PMT_photon_n[1] * 3 / PMT3;
  std::size_t out = std::round((PMT3 * grid_fraction + PMT_photon_n[1]) / 4.0);
  return out;
}

std::size_t Run::Results::GetNRecoredSiPMs23(void) const
{
  if (SiPM_photon_n.size() != 25)
    return 0;
  std::size_t SiPM43_avg = SiPM_photon_n[1] + SiPM_photon_n[3] + SiPM_photon_n[5] +
      SiPM_photon_n[9] + SiPM_photon_n[15] + SiPM_photon_n[19] + SiPM_photon_n[21] +
      SiPM_photon_n[23];
  std::size_t SiPM44_avg = SiPM_photon_n[0] + SiPM_photon_n[4] + SiPM_photon_n[20] +
        SiPM_photon_n[24];
  std::size_t out = GetNRecoredSiPMsRAW();
  out -= std::round(SiPM43_avg/8.0 + SiPM44_avg/4.0);
  return out;
}

void Run::Results::FindSensorsCoordinates(void) //Sets global positions of SiPMs and PMTs after geometry is constructed in G4RunManager::Initialize()
{
  SiPM_positions.clear();
  PMT_positions.clear();
  std::vector<G4PhysicalVolumesSearchScene::Findings> findingsVectorSiPM = FindAllPVs(gPars::det_dims->SiPM_device_name);
  std::vector<G4PhysicalVolumesSearchScene::Findings> findingsVectorPMT = FindAllPVs(gPars::det_dims->PMT_device_name);
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

void Run::Results::Clear(void)
{
  SiPM_positions.clear();
  PMT_positions.clear();
  PMT_photon_n.clear();
  SiPM_photon_n.clear();
  generated_photons.clear();
  recorded_photons.clear();
  n_reflections = 0;
  n_electrons = 0;
  n_photons = 0;
  n_events = 0;
}

void Run::Results::Init(void)
{
  FindSensorsCoordinates();
  SiPM_photon_n.resize(SiPM_positions.size(), 0);
  PMT_photon_n.resize(PMT_positions.size(), 0);
}

void Run::Results::ClearAndInit(void)
{
  Clear();
  Init();
}

void Run::Results::Print(std::ostream &str) const
{
  str<<"************************************************"<<std::endl;
  if (VSourceSettings::NBrS == gPars::source->generator_type) {
    SettingsNBrSGenerator *settings = static_cast<SettingsNBrSGenerator*>(gPars::source);
    str<<"NBrS real yield was multiplied by: "<<settings->NBrS_yield_factor<<std::endl;
  }
  str<<"Run information from detailed records:"<<std::endl;
  str<<"Electrons generated:"<<GetNGeneratedElectrons()<<std::endl;
  str<<"Photons generated:"<<GetNGeneratedPhotons()<<std::endl;
  str<<"Total photons detected:"<<GetNRecordedPhotons()<<std::endl;
  str<<std::endl;
  str<<"Run information from counters:"<<std::endl;
  str<<"Total number of events:"<<n_events<<std::endl;
  str<<"Electrons generated:"<<n_electrons<<std::endl;
  str<<"Photons generated:"<<n_photons<<std::endl;
  str<<"Total photons detected: "<<GetNRecoredPMTsRAW()<<" (PMTs) + "<<GetNRecoredSiPMsRAW() <<"(SiPMs) = "<<GetNRecoredPMTsRAW()+GetNRecoredSiPMsRAW()<<std::endl;
  //str<<"Total number of reflections:"<<n_reflections<<std::endl;
  str<<std::endl;

  if (!PMT_photon_n.empty()) {
    str<<"PMT detected photons:"<<std::endl;
    str<<"PMT#0: "<<PMT_photon_n[0]<<"\t"<<"PMT#1: "<<PMT_photon_n[1] \
      <<"\t"<<"PMT#2: "<<PMT_photon_n[2]<<"\t"<<"PMT#3: "<<PMT_photon_n[3]<<std::endl;
    str<<"PMT positions:"<<std::endl;
    str<<"X: PMT#0: "<<PMT_positions[0].getX()/mm<<"\t"<<"PMT#1: "<<PMT_positions[1].getX()/mm \
      <<"\t"<<"PMT#2: "<<PMT_positions[2].getX()/mm<<"\t"<<"PMT#3: "<<PMT_positions[3].getX()/mm<<std::endl;
    str<<"Y: PMT#0: "<<PMT_positions[0].getY()/mm<<"\t"<<"PMT#1: "<<PMT_positions[1].getY()/mm \
      <<"\t"<<"PMT#2: "<<PMT_positions[2].getY()/mm<<"\t"<<"PMT#3: "<<PMT_positions[3].getY()/mm<<std::endl;
  }

  str<<std::endl<<"SiPM results (ch, Nph, X[mm], Y[mm]):"<<std::endl;
  for (std::size_t i = 0, i_end_ = SiPM_photon_n.size(); i!=i_end_; ++i)
    str<<i<<"\t";
  str<<std::endl;
  for (std::size_t i = 0, i_end_ = SiPM_photon_n.size(); i!=i_end_; ++i)
    str<<SiPM_photon_n[i]<<"\t";
  str<<std::endl;
  for (std::size_t i = 0, i_end_ = SiPM_photon_n.size(); i!=i_end_; ++i)
    str<<SiPM_positions[i].getX()/mm<<"\t";
  str<<std::endl;
  for (std::size_t i = 0, i_end_ = SiPM_photon_n.size(); i!=i_end_; ++i)
    str<<SiPM_positions[i].getY()/mm<<"\t";
  str<<std::endl;

  std::cout<<std::endl<<"Average Npe per 1 PMT (grid is accounted for) : "<<GetNRecoredPMTavg()<<std::endl;
  std::cout<<"Npe for 23 SiPMs (no #43 and #44) : "<<GetNRecoredSiPMs23()<<std::endl;
}

void Run::Results::Merge(const Results &with)
{
  if(n_events == 0) {
    *this = with;
    return;
  }
  if ((PMT_positions.size() != with.PMT_positions.size()) ||
      (PMT_photon_n.size() != with.PMT_photon_n.size())) {
    G4cerr<<"Run::Results::Merge:Error: PMT number mismatch. Not merging."<<std::endl;
    return;
  }
  if ((SiPM_positions.size() != with.SiPM_positions.size()) ||
      (SiPM_photon_n.size() != with.SiPM_photon_n.size())) {
    G4cerr<<"Run::Results::Merge:Error: SiPM number mismatch. Not merging."<<std::endl;
    return;
  }
  generated_photons.insert(generated_photons.end(), with.generated_photons.begin(), with.generated_photons.end());
  recorded_photons.insert(recorded_photons.end(), with.recorded_photons.begin(), with.recorded_photons.end());
  for (std::size_t i = 0, i_end_ = SiPM_photon_n.size(); i!=i_end_; ++i) {
    SiPM_photon_n[i] += with.SiPM_photon_n[i];
  }
  for (std::size_t i = 0, i_end_ = PMT_photon_n.size(); i!=i_end_; ++i) {
    PMT_photon_n[i] += with.PMT_photon_n[i];
  }
  n_reflections += with.n_reflections;
  n_events += with.n_events;
  n_electrons += with.n_electrons;
  n_photons += with.n_photons;
}

void Run::RecordEvent(const G4Event* event)
{
  if (hit_collection_ID < 0 ) {
    G4HCtable* HCtable = G4SDManager::GetSDMpointer()->GetHCtable();
    if (HCtable) // Not calling G4SDManager::GetCollectionID because this approach is quiet when there is no collection
      hit_collection_ID = HCtable->GetCollectionID(DetectorSensor::PhotonCollectionName);
  }

  ++results.n_events;
  G4int n = GetNumberOfEvent();
  if (0 == n) {
    G4int N = GetNumberOfEventToBeProcessed();
    gData.progress_bar.set_N_events(N);
  }
  gData.progress_bar.tick();

  const VGeneratePrimaries* generator =
      static_cast<const VGeneratePrimaries*>(G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction());
  if (!generator) {
    std::cerr<<"Run::RecordEvent:Error: No primary generator found.\n";
    std::cerr<<"\tCan't gather run results and statistics."<<std::endl;
    return;
  }
  if (generator->GetGeneratedData().empty())
    return G4Run::RecordEvent(event);

  results.n_electrons += generator->GetGeneratedData().size();
  for (std::size_t e = 0, e_end_ = generator->GetGeneratedData().size(); e!=e_end_; ++e)
    results.n_photons += generator->GetGeneratedData()[e].photons.size();
  if (gPars::general.record_detailed) {
    if (gPars::general.record_electrons) {
      results.generated_photons.insert(results.generated_photons.end(),
          generator->GetGeneratedData().begin(), generator->GetGeneratedData().end());
      // Only 1 electron should be in one event, at least in case photons are generated.
      results.recorded_photons.push_back(GeneratedData());
      results.recorded_photons.back().electron = results.generated_photons.back().electron;
    } else {
      if (results.generated_photons.empty()) { // first event
        results.generated_photons = generator->GetGeneratedData(); // Store first electron for seed info
        results.recorded_photons.push_back(GeneratedData());
      } else {
        for (std::size_t i = 0, i_end_ = generator->GetGeneratedData().size(); i!=i_end_; ++i) {
          results.generated_photons.back().photons.insert(results.generated_photons.back().photons.end(),
              generator->GetGeneratedData()[i].photons.begin(), generator->GetGeneratedData()[i].photons.end());
        }
      }
    }
  } else if (gPars::general.record_electrons) { // record electron only
    for (std::size_t e = 0, e_end_ = generator->GetGeneratedData().size(); e!=e_end_; ++e) {
      GeneratedData electron_data;
      electron_data.electron = generator->GetGeneratedData()[e].electron;
      results.generated_photons.push_back(electron_data);
    }
  }

  G4HCofThisEvent* HCE = event->GetHCofThisEvent();
  if(!HCE || hit_collection_ID < 0)
     return G4Run::RecordEvent(event);
  G4VHitsCollection* coll = HCE->GetHC(hit_collection_ID);
  if (!coll)
    return G4Run::RecordEvent(event);
  PhotonHitCollection* hit_collection = static_cast<PhotonHitCollection*>(coll);
  if (!hit_collection)
    return G4Run::RecordEvent(event);
  for (std::size_t i = 0, i_end_ = hit_collection->GetSize(); i!=i_end_; ++i) {
    PhotonHit* hit = static_cast<PhotonHit*>(hit_collection->GetHit(i));
    if (hit) {
      if (gPars::general.record_detailed)
        results.recorded_photons.back().photons.push_back(*hit);
      if (hit->_isPMT) {
        if (hit->_channel < results.PMT_photon_n.size())
          ++results.PMT_photon_n[hit->_channel];
      } else {
        if (hit->_channel < results.SiPM_photon_n.size())
          ++results.SiPM_photon_n[hit->_channel];
      }
    }
  }
  G4Run::RecordEvent(event);
}

void Run::Merge(const G4Run* aRun)
{
  const Run* localRun = static_cast<const Run*>(aRun);
  results.Merge(localRun->results);
  G4Run::Merge(aRun);
}

void Run::Merged()
{
  results.Clear();
}

void Run::AddToFile(std::deque<GeneratedData> data, std::string filename)
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
