#include <geant4/PhysicsList.hh>

PhysicsList::PhysicsList() : G4VModularPhysicsList()
{
	G4LossTableManager::Instance();
	//// EM physics;
	defaultCutValue = 0.001*mm;
	SetVerboseLevel(1);
}


PhysicsList::~PhysicsList()
{}

void PhysicsList::ConstructParticle()
{
	G4OpticalPhoton::OpticalPhotonDefinition();
}

void PhysicsList::ConstructProcess()
{
  AddTransportation();
  ConstructOptical();
}

void PhysicsList::ConstructOptical()
{
  G4OpAbsorption* theAbsorptionProcess = new G4OpAbsorption();
  //G4OpRayleigh* theRayleighScattering=new G4OpRayleigh();

  G4OpWLS* fWLSProcess = new G4OpWLS();
  fWLSProcess->UseTimeProfile("delta");
  //fWLSProcess->UseTimeProfile("exponential");

  G4OpBoundaryProcess* theBoundaryProcess = new G4OpBoundaryProcess();
  theBoundaryProcess->SetInvokeSD(true); // Produces Hits according to EFFINIENCY of surface
  TeleportationProcess* theMappingProcess = new TeleportationProcess();
  KillProcess* theCutProcess = new KillProcess();

  G4ProcessManager * pManager = G4OpticalPhoton::OpticalPhoton()->GetProcessManager();

  pManager->AddDiscreteProcess(theAbsorptionProcess);
  //pManager->AddDiscreteProcess(theRayleighScattering);
  pManager->AddDiscreteProcess(theBoundaryProcess);
  pManager->AddDiscreteProcess(theMappingProcess);
  pManager->AddDiscreteProcess(theCutProcess);
  //pManager->AddDiscreteProcess(fWLSProcess);
}
