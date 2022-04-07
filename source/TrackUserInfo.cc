#include "TrackUserInfo.hh"

TrackUserInfo::TrackUserInfo()
{
  pType = new G4String(gPars::general.track_mapping_info_class);
}

TrackUserInfo::TrackUserInfo(const G4String& infoType)
{
  pType = new G4String(infoType);
}

TrackUserInfo::TrackUserInfo(const HexagonalMappingData& mapping_info)
{
  pType = new G4String(gPars::general.track_mapping_info_class);
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
