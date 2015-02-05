#include "Action_CreateCrd.h"
#include "CpptrajStdio.h"

Action_CreateCrd::Action_CreateCrd() {}

void Action_CreateCrd::Help() {
  mprintf("\t[<name>] [ parm <name> | parmindex <#> ]\n"
          "  Create a COORDS data set named <name> for frames associated with the\n"
          "  specified topology.\n");
}

Action::RetType Action_CreateCrd::Init(ArgList& actionArgs, TopologyList* PFL, DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  // Keywords
  Topology* parm = PFL->GetParm( actionArgs );
  if (parm == 0) {
    mprinterr("Error: createcrd: No parm files loaded.\n");
    return Action::ERR;
  }
  pindex_ = parm->Pindex();
  // DataSet
  std::string setname = actionArgs.GetStringNext();
  if (setname == "_DEFAULTCRD_") {
    // Special case: Creation of COORDS DataSet has been requested by an
    //               analysis and should already be present.
    coords_ = (DataSet_Coords*)DSL->FindSetOfType(setname, DataSet::COORDS);
  } else 
    coords_ = (DataSet_Coords*)DSL->AddSet(DataSet::COORDS, setname, "CRD");
  if (coords_ == 0) return Action::ERR;
  // Do not set topology here since it may be modified later.

  mprintf("    CREATECRD: Saving coordinates from Top %s to \"%s\"\n",
          parm->c_str(), coords_->Legend().c_str());
  return Action::OK;
}

Action::RetType Action_CreateCrd::Setup(Topology* currentParm, Topology** parmAddress) {
  // Set COORDS topology now if not already set.
  if (currentParm->Pindex() == pindex_ && coords_->Top().Natom() == 0)
    coords_->SetTopology( *currentParm );
  // If # atoms in currentParm does not match coords, warn user.
  if (currentParm->Natom() != coords_->Top().Natom())
    mprintf("Warning: # atoms in current topology (%i) != # atoms in coords set \"%s\" (%i)\n"
            "Warning:   The resulting COORDS data set may have problems.\n",
            currentParm->Natom(), coords_->Legend().c_str(), coords_->Top().Natom());
  return Action::OK;
}

Action::RetType Action_CreateCrd::DoAction(int frameNum, Frame* currentFrame, 
                                           Frame** frameAddress) 
{
  coords_->AddFrame( *currentFrame );
  return Action::OK;
}
