#include <geant4/DetectorSensor.hh>

const std::string DetectorSensor::PhotonCollectionName = "PhotonHitCollection";

DetectorSensor::DetectorSensor(G4String name) :
	G4VSensitiveDetector(name)
{
	collectionName.insert(PhotonCollectionName);
}

DetectorSensor::~DetectorSensor()
{}

void DetectorSensor::Initialize(G4HCofThisEvent* HCE)
{
	hitCollection = new PhotonHitCollection(GetName(), collectionName[0]);
	static G4int HCID = -1;
	if (HCID < 0)
		HCID = GetCollectionID(0);
	HCE->AddHitsCollection(HCID, hitCollection);
}

G4bool DetectorSensor::ProcessHits(G4Step* aStep, G4TouchableHistory* touchHist)
{
	return ProcessHits_Optical(aStep, touchHist);
}

//Generates a hit and uses the postStepPoint's mother volume replica number
//PostStepPoint because the hit is generated manually when the photon is
//absorbed by the photocathode
G4bool DetectorSensor::ProcessHits_Optical(const G4Step* aStep, G4TouchableHistory* touchHist)
{
	//need to know if this is an optical photon
	if(aStep->GetTrack()->GetDefinition()!= G4OpticalPhoton::OpticalPhotonDefinition())
		return false;
	G4StepPoint *thePrePoint  = aStep->GetPostStepPoint();
	G4VPhysicalVolume* vol = aStep->GetPostStepPoint()->GetTouchable()->GetVolume();
	if (NULL != vol && vol->GetName() == gPars::det_dims->PMT_device_name) {
		std::size_t PMT_no = vol->GetCopyNo();
		PhotonHit *hit = new PhotonHit();
		hit->_energy = thePrePoint->GetTotalEnergy();
		hit->_time = thePrePoint->GetGlobalTime();
		hit->_pos = thePrePoint->GetPosition();
		hit->_isPMT = true;
		hit->_channel = PMT_no;
		hitCollection->insert(hit);
	}
	if (NULL != vol && vol->GetName() == gPars::det_dims->SiPM_device_name) {
		std::size_t SiPM_no = vol->GetCopyNo();
		PhotonHit *hit = new PhotonHit();
		hit->_energy = thePrePoint->GetTotalEnergy();
		hit->_time = thePrePoint->GetGlobalTime();
		hit->_pos = thePrePoint->GetPosition();
		hit->_isPMT = false;
		hit->_channel = SiPM_no;
		hitCollection->insert(hit);
	}
	return true;
}

void DetectorSensor::EndOfEvent(G4HCofThisEvent* HCE)
{}
