#ifndef TRACKING_ACTION_H_
#define TRACKING_ACTION_H_

#include "G4UserTrackingAction.hh"
#include "G4Track.hh"
#include "G4PrimaryParticle.hh"

#include "HexagonalMapping.hh"
#include "PrimaryParticleUserInfo.hh"
#include "TrackUserInfo.hh"

// Only required to pass mapping data from G4PrimaryParticle to G4Track
class TrackingAction : public G4UserTrackingAction
{
public:
  void PreUserTrackingAction(const G4Track* track);
};

#endif // TRACKING_ACTION_H_
