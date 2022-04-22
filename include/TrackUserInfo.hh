#ifndef TRACK_USER_INFO_H_
#define TRACK_USER_INFO_H_

#include "G4VUserTrackInformation.hh"

#include "GlobalParameters.hh"
#include "HexagonalMapping.hh"

class TrackUserInfo : public G4VUserTrackInformation
{
public:
  TrackUserInfo();
  TrackUserInfo(const G4String& infoType);
  TrackUserInfo(const HexagonalMappingData& mapping_info);

  inline void SetMappingData(HexagonalMappingData mapping_info)
  { mapping_data = mapping_info;}
  inline HexagonalMappingData GetMappingData(void) const
  { return mapping_data;}

  TrackUserInfo(const TrackUserInfo&);
  TrackUserInfo& operator=(const TrackUserInfo&);

protected:
  HexagonalMappingData mapping_data;
public:
  static const std::string track_info_classname;
};

#endif // TRACK_USER_INFO_H_
