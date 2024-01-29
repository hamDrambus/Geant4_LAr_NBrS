#include <geant4/MappingManager.hh>
#include <GlobalParameters.hh>
#include <GlobalData.hh>

MappingManager::MappingManager()
{}

void MappingManager::AddMapping(HexagonalMapping& mapping)
{
	for (std::size_t i = 0, i_end_ = mappings.size(); i!=i_end_; ++i) {
		if (mappings[i].name_id == mapping.name_id) {
			mappings[i] = mapping;
			if (mappings[i].has_field_map) {
				for (std::size_t m = 0, m_end_ = mappings.size(); m!=m_end_; ++m) {
					if (m != i)
						mappings[m].has_field_map = false;
				}
			}
			return;
		}
	}
	if (mapping.has_field_map) {
		for (std::size_t m = 0, m_end_ = mappings.size(); m!=m_end_; ++m)
			mappings[m].has_field_map = false;
	}
	mappings.push_back(mapping);
}

void MappingManager::RemoveMapping(std::string name_id)
{
	mappings.erase(std::remove_if(mappings.begin(), mappings.end(), [name_id](const HexagonalMapping& t) {return t.name_id == name_id;}), mappings.end());
}

void MappingManager::ClearMappings(void)
{
	mappings.clear();
}


HexagonalMappingData MappingManager::GetNewState(const G4Track& aTrack, const G4Step& aStep, const HexagonalMappingData& old_state) const
{
	HexagonalMappingData new_state = old_state;
	std::vector<bool> is_entering(mappings.size(), false);
	std::vector<bool> is_leaving(mappings.size(), false);
	for (std::size_t m = 0, m_end_ = mappings.size(); m!=m_end_; ++m) {
		for (std::size_t i = 0, i_end_ = mappings[m].map_triggers.size(); i!=i_end_; ++i) {
			if (mappings[m].map_triggers[i].IsValid() && mappings[m].map_triggers[i].IsSatistied(aTrack, aStep)) {
				if (mappings[m].map_triggers[i].is_for_entering)
					is_entering[m] = old_state.mapping_id.empty(); // Nested mappings are forbidden.
				else
					is_leaving[m] = (old_state.mapping_id == mappings[m].name_id); // Cannot leave mapping if state does not belong to it.
			}
		}
	}
	int n_entering = 0, n_leaving = 0;
	for (std::size_t m = 0, m_end_ = mappings.size(); m!=m_end_; ++m) {
		if (is_entering[m] && is_leaving[m]) {
			std::cerr<<"MappingManager::GetNewState:Error:"<<std::endl;
			std::cerr<<"\tConditions for both leaving and entering cell are satisfied\n"
					"at the same time for mapping \""<<mappings[m].name_id<<"\".\n"
					"Most likely caused by wrong geometry/mapping setup.\n"
					"This mapping is ignored."<<std::endl;
			is_entering[m] = false;
			is_leaving[m] = false;
			continue;
		}
		if (is_entering[m])
			++n_entering;
		if (is_leaving[m])
			++n_leaving;
	}
	if (n_entering > 1) {
		std::cerr<<"MappingManager::GetNewState:Error:"<<std::endl;
		std::cerr<<"\tConditions for entering several mappings are satisfied at the same time!\n"
				"Most likely caused by wrong geometry/mapping setup.\n"
				"First mapping in the array is used."<<std::endl;
	}
	if (n_leaving > 1) {
		std::cerr<<"MappingManager::GetNewState:Error:"<<std::endl;
		std::cerr<<"\tConditions for leaving several mappings are satisfied at the same time!\n"
				"Most likely caused by wrong geometry/mapping setup.\n"
				"First mapping in the array is used."<<std::endl;
	}
	if (n_leaving == 1 && n_entering == 1) {
		std::cerr<<"MappingManager::GetNewState:Error:"<<std::endl;
		std::cerr<<"\tConditions for leaving one mapping and entering another are satisfied at the same time!\n"
				"Most likely caused by wrong geometry/mapping setup.\n"
				"Mapping is ignored"<<std::endl;
		return new_state;
	}
	if (n_entering > 0) {
		for (std::size_t m = 0, m_end_ = mappings.size(); m!=m_end_; ++m) {
			if (is_entering[m])
				return mappings[m].MapToCell(old_state, true);
		}
	}
	if (n_leaving > 0) {
		for (std::size_t m = 0, m_end_ = mappings.size(); m!=m_end_; ++m) {
			if (is_leaving[m])
				return mappings[m].MapFromCell(old_state, false);
		}
	}
	return new_state;
}

HexagonalMappingData MappingManager::GetNewState(const G4VPhysicalVolume* volume, const HexagonalMappingData& old_state) const
{
	HexagonalMappingData new_state = old_state;
	std::vector<bool> is_entering(mappings.size(), false);
	std::vector<bool> is_leaving(mappings.size(), false);
	for (std::size_t m = 0, m_end_ = mappings.size(); m!=m_end_; ++m) {
		for (std::size_t i = 0, i_end_ = mappings[m].map_triggers.size(); i!=i_end_; ++i) {
			if (mappings[m].map_triggers[i].IsValid() && mappings[m].map_triggers[i].IsSatistied(volume, old_state)) {
				if (mappings[m].map_triggers[i].is_for_entering)
					is_entering[m] = old_state.mapping_id.empty(); // Nested mappings are forbidden.
				else
					is_leaving[m] = (old_state.mapping_id == mappings[m].name_id); // Cannot leave mapping if state does not belong to it.
			}
		}
	}
	int n_entering = 0, n_leaving = 0;
	for (std::size_t m = 0, m_end_ = mappings.size(); m!=m_end_; ++m) {
		if (is_entering[m] && is_leaving[m]) {
			std::cerr<<"MappingManager::GetNewState:Error:"<<std::endl;
			std::cerr<<"\tConditions for both leaving and entering cell are satisfied\n"
					"at the same time for mapping \""<<mappings[m].name_id<<"\".\n"
					"Most likely caused by wrong geometry/mapping setup.\n"
					"This mapping is ignored."<<std::endl;
			is_entering[m] = false;
			is_leaving[m] = false;
			continue;
		}
		if (is_entering[m])
			++n_entering;
		if (is_leaving[m])
			++n_leaving;
	}
	if (n_entering > 1) {
		std::cerr<<"MappingManager::GetNewState:Error:"<<std::endl;
		std::cerr<<"\tConditions for entering several mappings are satisfied at the same time!\n"
				"Most likely caused by wrong geometry/mapping setup.\n"
				"First mapping in the array is used."<<std::endl;
	}
	if (n_leaving > 1) {
		std::cerr<<"MappingManager::GetNewState:Error:"<<std::endl;
		std::cerr<<"\tConditions for leaving several mappings are satisfied at the same time!\n"
				"Most likely caused by wrong geometry/mapping setup.\n"
				"First mapping in the array is used."<<std::endl;
	}
	if (n_leaving == 1 && n_entering == 1) {
		std::cerr<<"MappingManager::GetNewState:Error:"<<std::endl;
		std::cerr<<"\tConditions for leaving one mapping and entering another are satisfied at the same time!\n"
				"Most likely caused by wrong geometry/mapping setup.\n"
				"Mapping is ignored"<<std::endl;
		return new_state;
	}
	if (n_entering > 0) {
		for (std::size_t m = 0, m_end_ = mappings.size(); m!=m_end_; ++m) {
			if (is_entering[m])
				return mappings[m].MapToCell(old_state, false);
		}
	}
	if (n_leaving > 0) {
		std::cerr<<"MappingManager::GetNewState:Error:"<<std::endl;
		std::cerr<<"\tParticle is generated inside mapping (THGEM) cell! Cell mapping is undefined. No mapping is used."<<std::endl;
	}
	return new_state;
}

G4ThreeVector MappingManager::GetFieldAtGlobal(G4ThreeVector position, DriftMedium* &medium, int *status) const
{
	medium = nullptr;
	if (nullptr == gData.field_map || !gData.field_map->IsReady()) {
		std::cerr<<"MappingManager::GetFieldAtGlobal:Error: Field map is not available."<<std::endl;
		if (nullptr != status)
			*status = 1;
		return G4ThreeVector(0, 0, 0);
	}
	if (!HasFieldMapping()) {
		std::cerr<<"MappingManager::GetFieldAtGlobal:Error: THGEM1 field mapping is not available."<<std::endl;
		if (nullptr != status)
			*status = 2;
		return G4ThreeVector(0, 0, 0);
	}
	const HexagonalMapping *map = nullptr;
	for (std::size_t i = 0, i_end_ = mappings.size(); i!=i_end_; ++i)
		if (mappings[i].has_field_map)
			map = &mappings[i];

	HexagonalMappingData mapping_data;
	mapping_data.cell_x_ind = -1;
	mapping_data.cell_y_ind = -1;
	mapping_data.position = position;
	mapping_data.momentum = G4ThreeVector(0, 0, 0);
	mapping_data.polarization = G4ThreeVector(0, 0, 0);
	mapping_data = map->MapToCell(mapping_data, true);
	if (!mapping_data.isInCell()) {
		if (nullptr != status)
				*status = -6;
		return mapping_data.momentum;
	}
	G4ThreeVector in_mesh_pos = mapping_data.position - map->g_cell_pos + gPars::field_map.elmer_mesh_center;
	double ex,ey,ez;
	int status_;
	gData.field_map->ElectricField(in_mesh_pos.x(),in_mesh_pos.y(),in_mesh_pos.z(),ex,ey,ez,medium, status_);
	mapping_data.momentum = G4ThreeVector(ex, ey, ez);
	mapping_data = map->MapFromCell(mapping_data, true);
	if (nullptr != status)
		*status = status_;
	return mapping_data.momentum;
}
