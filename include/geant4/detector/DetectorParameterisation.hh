#ifndef DetectorParametrisation_h
#define DetectorParametrisation_h

#include <G4ThreeVector.hh>
#include <G4Tubs.hh>
#include <G4Box.hh>
#include <G4SystemOfUnits.hh>
#include <G4VPhysicalVolume.hh>
#include <G4VPVParameterisation.hh>

///  A quite arbitrary parameterisation that describes a series of rotated elements along any of the axes.
///  The elements have equal sizes and the same rotation.
///  They are spaced an equal distance along each of the axes.
///  Does not check the validity of geometry.
class DetectorParameterisation : public G4VPVParameterisation
{
public:
	DetectorParameterisation(G4int noElemsX, G4int noElemsY, G4int noElemsZ,
		G4RotationMatrix* rot, G4ThreeVector startPos, G4ThreeVector spacing);
	void ComputeTransformation(const G4int copyNo, G4VPhysicalVolume* physVol) const;
private:
	G4int fNoElementsX;
	G4int fNoElementsY;
	G4int fNoElementsZ;
	G4ThreeVector fStartPos;
	G4ThreeVector fSpacing;
	G4RotationMatrix *fRotation;
};

#endif
