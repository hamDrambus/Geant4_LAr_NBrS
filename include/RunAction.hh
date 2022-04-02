#ifndef RunAction_h
#define RunAction_h

#include <iostream>

#include "G4Timer.hh"
#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "globals.hh"
#include "G4UserRunAction.hh"
#include "G4TransportationManager.hh"
#include "G4PhysicalVolumesSearchScene.hh"

#include "GlobalParameters.hh"
#include "GlobalUtilities.hh"

class RunAction : public G4UserRunAction
{
  public:
    RunAction();
    ~RunAction();

  public:
    void BeginOfRunAction(const G4Run* aRun);
    void EndOfRunAction(const G4Run* aRun);
    
  protected:
    void FindSensorsCoordinates(); //Sets global positions of SiPMs and PMTs in gPars:: after geometry is constructed
    void SetupTHGEM1Mapping();
};

#endif
