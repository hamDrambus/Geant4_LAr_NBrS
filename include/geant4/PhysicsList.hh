#ifndef PhysicsList_h
#define PhysicsList_h

#include <iostream>

#include <G4SystemOfUnits.hh>
#include <globals.hh>

#include <G4UnitsTable.hh>
#include <G4LossTableManager.hh>
#include <G4Region.hh>
#include <G4RegionStore.hh>

#include <G4VModularPhysicsList.hh>
#include <G4OpticalPhoton.hh>
#include <G4OpAbsorption.hh>
#include <G4OpRayleigh.hh>
#include <G4OpBoundaryProcess.hh>
#include <G4OpWLS.hh>
#include <G4ProcessManager.hh>

#include "TeleportationProcess.hh"
#include "KillProcess.hh"

class PhysicsList: public G4VModularPhysicsList
{
  public:
    PhysicsList();
    ~PhysicsList();

    void ConstructParticle();
    void ConstructOptical();
    void ConstructProcess();
};

#endif
