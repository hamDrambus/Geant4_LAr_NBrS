#include <geant4/detector/DetectorParameterisation.hh>

DetectorParameterisation::DetectorParameterisation(G4int noElemsX, G4int noElemsY, G4int noElemsZ,
	G4RotationMatrix* rot, G4ThreeVector startPos, G4ThreeVector spacing)
	: G4VPVParameterisation(), fNoElementsX(noElemsX), fNoElementsY(noElemsY), fNoElementsZ(noElemsZ),
	fStartPos(startPos), fSpacing(spacing), fRotation(rot)
{
	fNoElementsX = std::max(fNoElementsX, 1);
	fNoElementsY = std::max(fNoElementsY, 1);
	fNoElementsZ = std::max(fNoElementsZ, 1);
}

void DetectorParameterisation::ComputeTransformation (const G4int copyNo, G4VPhysicalVolume* physVol)  const
{
	// Note: copyNo will start with zero! copyNo <= fNoElementsX*fNoElementsY*fNoElementsZ - 1
	G4double Xposition = (-fNoElementsX / 2.0 + copyNo % fNoElementsX + 0.5) * fSpacing.getX();
	G4double Yposition = (-fNoElementsY / 2.0 + (copyNo / fNoElementsX) % fNoElementsY + 0.5) * fSpacing.getY();
	G4double Zposition = (-fNoElementsZ / 2.0 + (copyNo / (fNoElementsX * fNoElementsY)) + 0.5) * fSpacing.getZ();
	physVol->SetTranslation(G4ThreeVector(Xposition, Yposition, Zposition) + fStartPos);
	physVol->SetRotation(fRotation);
}

