#include <geant4/TrackingAction.hh>

void TrackingAction::PreUserTrackingAction(const G4Track* track)
{
  TrackUserInfo* user_info = (TrackUserInfo*) track->GetUserInformation();
  if (nullptr != user_info) //already has mapping data
    return;
  const G4PrimaryParticle *particle = track->GetDynamicParticle()->GetPrimaryParticle();
  if (nullptr == particle) {
    std::cout<<"TrackingAction::PreUserTrackingAction:Warning:"<<std::endl;
    std::cout<<"\tTrack without primary is missing mapping data. Using default."<<std::endl;
    goto bad_finish;
  } else {
    const PrimaryParticleUserInfo* particle_data = (PrimaryParticleUserInfo*) particle->GetUserInformation();
    if (nullptr == particle_data) {
      std::cout<<"TrackingAction::PreUserTrackingAction:Warning:"<<std::endl;
      std::cout<<"\tPrimary is missing mapping data. Using default."<<std::endl;
      goto bad_finish;
    }
    TrackUserInfo* user_info = new TrackUserInfo(particle_data->mapping_data);
    track->SetUserInformation(user_info);
  }
  return G4UserTrackingAction::PreUserTrackingAction(track);
bad_finish:
  HexagonalMappingData data;
  user_info = new TrackUserInfo(data);
  track->SetUserInformation(user_info);
  return G4UserTrackingAction::PreUserTrackingAction(track);
}
