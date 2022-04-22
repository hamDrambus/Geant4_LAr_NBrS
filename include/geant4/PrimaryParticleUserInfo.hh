#ifndef PRIMARY_PARTICLE_USER_INFO_H_
#define PRIMARY_PARTICLE_USER_INFO_H_

#include <G4VUserPrimaryParticleInformation.hh>
#include "HexagonalMapping.hh"

class PrimaryParticleUserInfo : public G4VUserPrimaryParticleInformation
{
public:
  PrimaryParticleUserInfo() = delete;
  PrimaryParticleUserInfo(HexagonalMappingData data)
  { mapping_data = data; }
  HexagonalMappingData mapping_data;
  void Print() const
  {}
};

#endif // PRIMARY_PARTICLE_USER_INFO_H_
