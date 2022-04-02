#include "RunAction.hh"

RunAction::RunAction()
{}

RunAction::~RunAction()
{}

void RunAction::BeginOfRunAction(const G4Run* aRun)
{
	gPars::results.n_reflections = 0;
	gPars::results.generated_photons.clear();
	gPars::results.recorded_photons.clear();
	unsigned int N = gPars::det_dims.n_SiPMs_rows;
	gPars::results.SiPM_photon_n = std::deque<unsigned int>(N*N, 0);
	gPars::results.PMT_photon_n = std::deque<unsigned int>(gPars::det_dims.n_PMTs, 0);
	FindSensorsCoordinates();
	SetupTHGEM1Mapping();
}

void RunAction::EndOfRunAction(const G4Run* aRun)
{
	std::cout<<"************************************************"<<std::endl;
	std::cout<<"Run finished"<<std::endl;
	std::cout<<"Photons generated:"<<gPars::results.GetNGeneratedPhotons()<<std::endl;
	std::cout<<"Total photons detected:"<<gPars::results.GetNRecordedPhotons()<<std::endl;
	std::cout<<"Total number of reflections:"<<gPars::results.n_reflections<<std::endl;

	std::cout<<std::endl<<"PMT detected photons:"<<std::endl;
	std::cout<<"PMT#0: "<<gPars::results.PMT_photon_n[0]<<"\t"<<"PMT#1: "<<gPars::results.PMT_photon_n[1] \
		<<"\t"<<"PMT#2: "<<gPars::results.PMT_photon_n[2]<<"\t"<<"PMT#3: "<<gPars::results.PMT_photon_n[3]<<std::endl;
	std::cout<<"PMT positions:"<<std::endl;
	std::cout<<"X: PMT#0: "<<gPars::results.PMT_positions[0].getX()/mm<<"\t"<<"PMT#1: "<<gPars::results.PMT_positions[1].getX()/mm \
		<<"\t"<<"PMT#2: "<<gPars::results.PMT_positions[2].getX()/mm<<"\t"<<"PMT#3: "<<gPars::results.PMT_positions[3].getX()/mm<<std::endl;
	std::cout<<"Y: PMT#0: "<<gPars::results.PMT_positions[0].getY()/mm<<"\t"<<"PMT#1: "<<gPars::results.PMT_positions[1].getY()/mm \
			<<"\t"<<"PMT#2: "<<gPars::results.PMT_positions[2].getY()/mm<<"\t"<<"PMT#3: "<<gPars::results.PMT_positions[3].getY()/mm<<std::endl;

	std::cout<<std::endl<<"SiPM results (ch, Nph, X[mm], Y[mm]):"<<std::endl;
	for (std::size_t i = 0, i_end_ = gPars::results.SiPM_photon_n.size(); i!=i_end_; ++i)
		std::cout<<i<<"\t";
	std::cout<<std::endl;
	for (std::size_t i = 0, i_end_ = gPars::results.SiPM_photon_n.size(); i!=i_end_; ++i)
		std::cout<<gPars::results.SiPM_photon_n[i]<<"\t";
	std::cout<<std::endl;
	for (std::size_t i = 0, i_end_ = gPars::results.SiPM_photon_n.size(); i!=i_end_; ++i)
		std::cout<<gPars::results.SiPM_positions[i].getX()/mm<<"\t";
	std::cout<<std::endl;
	for (std::size_t i = 0, i_end_ = gPars::results.SiPM_photon_n.size(); i!=i_end_; ++i)
		std::cout<<gPars::results.SiPM_positions[i].getY()/mm<<"\t";
	std::cout<<std::endl;

	AddToFile(gPars::results.generated_photons, gPars::results.generated_filename);
	AddToFile(gPars::results.recorded_photons, gPars::results.recorded_filename);
}

void RunAction::FindSensorsCoordinates()
{
	std::vector<G4PhysicalVolumesSearchScene::Findings> findingsVectorSiPM = FindAllPVs(gPars::det_dims.SiPM_device_name);
	std::vector<G4PhysicalVolumesSearchScene::Findings> findingsVectorPMT = FindAllPVs(gPars::det_dims.PMT_device_name);
	G4int max_number = -1;
	for (const auto& findings: findingsVectorSiPM)
		max_number = std::max(findings.fFoundPVCopyNo, max_number);
	gPars::results.SiPM_positions = std::deque<G4ThreeVector>(max_number+1);
	for (const auto& findings: findingsVectorSiPM)
		gPars::results.SiPM_positions[findings.fFoundPVCopyNo] = findings.fFoundObjectTransformation.getTranslation()/mm;
	max_number = -1;
	for (const auto& findings: findingsVectorPMT)
		max_number = std::max(findings.fFoundPVCopyNo, max_number);
	gPars::results.PMT_positions = std::deque<G4ThreeVector>(max_number+1);
	for (const auto& findings: findingsVectorPMT)
		gPars::results.PMT_positions[findings.fFoundPVCopyNo] = findings.fFoundObjectTransformation.getTranslation()/mm;
}

void RunAction::SetupTHGEM1Mapping()
{
  if (nullptr != gPars::THGEM1_mapping)
    delete gPars::THGEM1_mapping;
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
  gPars::THGEM1_mapping = new HexagonalMapping(thgem1_pos, cell_pos, thgem1_sizes, cell_sizes);
}
