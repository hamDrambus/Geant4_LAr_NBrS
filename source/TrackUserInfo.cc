#include "TrackUserInfo.hh"

const std::string TrackUserInfo::track_info_classname = "TrackMappingInfo";

TrackUserInfo::TrackUserInfo()
{
  pType = new G4String(track_info_classname);
}

TrackUserInfo::TrackUserInfo(const G4String& infoType)
{
  pType = new G4String(infoType);
}

TrackUserInfo::TrackUserInfo(const HexagonalMappingData& mapping_info)
{
  pType = new G4String(track_info_classname);
  mapping_data = mapping_info;
}

TrackUserInfo::TrackUserInfo(const TrackUserInfo& right)
{
  if(right.pType != nullptr)
    pType = new G4String(*(right.pType));
  mapping_data = right.mapping_data;
}

TrackUserInfo& TrackUserInfo::operator=(const TrackUserInfo& right)
{
  if(this != &right)
  {
    if(pType != nullptr)
      delete pType;
    if(right.pType)
      pType = new G4String(*(right.pType));
    else
      pType = nullptr;
    mapping_data = right.mapping_data;
  }
  return *this;
}
