#include "TeleportationProcess.hh"

TeleportationProcess::TeleportationProcess(const G4String& processName, G4ProcessType type)
  : G4VDiscreteProcess(processName, type)
{
  Initialise();
  if(verboseLevel > 0) {
    G4cout << GetProcessName() << " is created " << G4endl;
  }
  SetProcessSubType(fTeleportTHGEM);
  fTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();
}

TeleportationProcess::~TeleportationProcess() {}

void TeleportationProcess::Initialise()
{
  SetVerboseLevel(gPars::debugging.teleportation_verbosity);
}

G4VParticleChange* TeleportationProcess::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep)
{
  aParticleChange.Initialize(aTrack);

  // Get hyperStep from  G4ParallelWorldProcess
  //  NOTE: PostStepDoIt of this process to be invoked after
  //  G4ParallelWorldProcess!
  const G4Step* pStep = &aStep;
  const G4Step* hStep = G4ParallelWorldProcess::GetHyperStep();
  if(hStep != nullptr)
    pStep = hStep;

  if(pStep->GetPostStepPoint()->GetStepStatus() != fGeomBoundary) {
    if(verboseLevel > 1)
      TeleportationProcessVerbose();
    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
  }

  G4VPhysicalVolume* prevPVolume = pStep->GetPreStepPoint()->GetPhysicalVolume();
  G4VPhysicalVolume* postPVolume = pStep->GetPostStepPoint()->GetPhysicalVolume();

  if(verboseLevel > 1) {
    G4cout << " Particle at Boundary! " << G4endl;
    if(prevPVolume != nullptr)
      G4cout << " prevPVolume:  " << prevPVolume->GetName() << G4endl;
    if(postPVolume != nullptr)
      G4cout << " postPVolume: " << postPVolume->GetName() << G4endl;
  }

  if (nullptr == gPars::THGEM1_mapping) {
     if(verboseLevel > 1) {
       G4cout << " No mapping class in global parameters, skipping process." << G4endl;
     }
     return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
   }

  TrackUserInfo *mapping_data = (TrackUserInfo*) aTrack.GetUserInformation();
  if (nullptr == mapping_data) {
    if(verboseLevel > 1) {
      G4cout << " No mapping data in track, skipping process." << G4endl;
    }
    return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
  }
  fOldState = mapping_data->GetMappingData();
  const G4DynamicParticle* aParticle = aTrack.GetDynamicParticle();
  fOldState.momentum = aParticle->GetMomentumDirection();
  fOldState.polarization = aParticle->GetPolarization();
  fOldState.position = pStep->GetPostStepPoint()->GetPosition();
  fNewState = fOldState;

  if(verboseLevel > 1) {
    G4cout << " Old Position: " << fOldState.position << G4endl
           << " Old Momentum Direction: " << fOldState.momentum << G4endl
           << " Old Polarization:       " << fOldState.polarization << G4endl;
  }

  bool valid;
  // ID of Navigator which limits step
  G4int hNavId = G4ParallelWorldProcess::GetHypNavigatorID();
  auto iNav = G4TransportationManager::GetTransportationManager()->GetActiveNavigatorsIterator();
  G4ThreeVector globalNormal = (iNav[hNavId])->GetGlobalExitNormal(fOldState.position, &valid);
  if(!valid) {
    G4ExceptionDescription ed;
    ed << " TeleportationProcess/PostStepDoIt(): "
       << " The Navigator reports that it returned an invalid normal" << G4endl;
    G4Exception(
      "TeleportationProcess::PostStepDoIt", "TeleportProc01", EventMustBeAborted, ed,
      "Invalid Surface Normal - Geometry must return valid surface normal");
  }

  if (prevPVolume->GetName() == gPars::det_dims.THGEM1_cell_name
      && globalNormal*fOldState.momentum > 0.0
      && postPVolume->GetLogicalVolume()->IsAncestor(prevPVolume)) {
    if (std::fabs(globalNormal * G4ThreeVector(0, 0, 1) - 1.0) < fTolerance
        || std::fabs(globalNormal * G4ThreeVector(0, 0, -1) - 1.0) < fTolerance) //TODO: make condition not hard-coded. Maybe move this logic to HexagonalMapping?
      fNewState = gPars::THGEM1_mapping->MapFromCell(fOldState, false);
    else
      fNewState = gPars::THGEM1_mapping->MapToNeighbourCell(fOldState);
  }
  if (postPVolume->GetName() == gPars::det_dims.THGEM1_cell_container_name
      && globalNormal*fOldState.momentum > 0.0
      && !postPVolume->GetLogicalVolume()->IsAncestor(prevPVolume)) {
    fNewState = gPars::THGEM1_mapping->MapToCell(fOldState, true);
  }

  if (fNewState != fOldState) {
    aParticleChange.ProposeTrackStatus(fStopAndKill);
    aParticleChange.ProposeLocalEnergyDeposit(0.0);
    // Generate a new particle:
    G4DynamicParticle* replacement = new G4DynamicParticle(aParticle->GetDefinition(), fNewState.momentum, aParticle->GetKineticEnergy());
    replacement->SetPolarization(fNewState.polarization);

    G4Track* replacementTrack = new G4Track(replacement, pStep->GetPostStepPoint()->GetGlobalTime(), fNewState.position);
    replacementTrack->SetTouchableHandle(nullptr); // Must treat teleported particle as primary. So navigator must relocate it from scratch.
    replacementTrack->SetParentID(aTrack.GetTrackID());
    replacementTrack->SetUserInformation(new TrackUserInfo(fNewState));
    replacementTrack->SetWeight(aTrack.GetWeight());
    replacementTrack->SetBelowThresholdFlag(aTrack.IsBelowThreshold());
    replacementTrack->SetGoodForTrackingFlag(aTrack.IsGoodForTracking());
    replacementTrack->SetLocalTime(aTrack.GetLocalTime());
    replacementTrack->SetProperTime(aTrack.GetProperTime());
    aParticleChange.SetNumberOfSecondaries(1);
    aParticleChange.AddSecondary(replacementTrack);
  }
  return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
}

void TeleportationProcess::TeleportationProcessVerbose() const
{}

G4double TeleportationProcess::GetMeanFreePath(const G4Track&, G4double, G4ForceCondition* condition)
{
  *condition = Forced;
  return DBL_MAX;
}

void TeleportationProcess::SetVerboseLevel(G4int verbose)
{
  verboseLevel = verbose;
}

inline G4bool TeleportationProcess::IsApplicable(const G4ParticleDefinition& aParticleType)
{
  return true;
}
