#ifndef PhotonHit_h
#define PhotonHit_h

#include <G4VHit.hh>
#include <G4THitsCollection.hh>
#include <G4Allocator.hh>
#include <G4ThreeVector.hh>
#include <G4UnitsTable.hh>
#include <G4VVisManager.hh>
#include <G4Circle.hh>
#include <G4Colour.hh>
#include <G4VisAttributes.hh>

using namespace CLHEP;

// Used to store generated photons as well
class PhotonHit : public G4VHit
{
  public:
    PhotonHit();
    virtual ~PhotonHit();
    G4int operator==(const PhotonHit&) const;
    void Draw();
    void Print(std::ostream &stream = G4cout, bool printtime = true, bool printposition = false, bool printenergy = false);

  public:
    G4double _time;
    G4double _energy;
    G4ThreeVector _pos;
    G4ThreeVector _momentum;
    G4ThreeVector _polarisation;
    unsigned int _channel;
    bool _isPMT;
};

typedef G4THitsCollection<PhotonHit> PhotonHitCollection;

#endif
