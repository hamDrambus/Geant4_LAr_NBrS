#include <geant4/PhotonHit.hh>

G4ThreadLocal G4Allocator<PhotonHit>* PhotonHitAllocator=0;

PhotonHit::PhotonHit()
  : _time(0), _energy(0), _pos(G4ThreeVector(0,0,0)), _channel(0), _isPMT(true)
{}

PhotonHit::~PhotonHit() {}

G4int PhotonHit::operator==(const PhotonHit& right) const
{
  return (this==&right) ? 1 : 0;
}

void PhotonHit::Draw()
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager)
  {
    G4Circle circle(_pos);
    circle.SetScreenSize(2.);
    circle.SetFillStyle(G4Circle::filled);
    G4Colour colour(1.,0.,0.);
    G4VisAttributes attribs(colour);
    circle.SetVisAttributes(attribs);
    circle.SetScreenDiameter(4.0);
    pVVisManager->Draw(circle);
  }
}

void PhotonHit::Print(std::ostream &stream, bool printtime, bool printposition, bool printenergy)
{
  if (printtime)
    stream << _time/ns;
  if (printposition) {
    if (printtime)
      stream << "\t";
    stream  << _pos.x()/mm << "\t" << _pos.y()/mm;
  }
  if (printenergy) {
    if (printtime || printposition)
      stream << "\t";
    stream << _energy/eV;
  }
  stream << std::endl;
}
