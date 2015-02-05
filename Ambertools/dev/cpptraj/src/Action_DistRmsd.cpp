// DISTRMSD
#include "Action_DistRmsd.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_DistRmsd::Action_DistRmsd() : drmsd_(0) {}

void Action_DistRmsd::Help() {
  mprintf("\t[<name>] [<mask>] [<refmask>] [out filename]\n"
          "\t[ first | %s |\n"
          "\t  reftraj <filename> [parm <parmname> | parmindex <#>] ]\n"
          "  Calculate distance RMSD (DME) for specified atoms.\n", DataSetList::RefArgs);
}

// Action_DistRmsd::Init()
/** Called once before traj processing. Set up reference info. */
Action::RetType Action_DistRmsd::Init(ArgList& actionArgs, TopologyList* PFL, DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  // Check for keywords
  DataFile* outfile = DFL->AddDataFile( actionArgs.GetStringKey("out"), actionArgs );
  // Reference keywords
  // TODO: Can these just be put in the InitRef call?
  bool first = actionArgs.hasKey("first");
  ReferenceFrame REF = DSL->GetReferenceFrame( actionArgs );
  std::string reftrajname = actionArgs.GetStringKey("reftraj");
  Topology* RefParm = PFL->GetParm( actionArgs );
  // Get the RMS mask string for target 
  std::string mask0 = actionArgs.GetMaskNext();
  TgtMask_.SetMaskString(mask0);
  // Get the RMS mask string for reference
  std::string mask1 = actionArgs.GetMaskNext();
  if (mask1.empty())
    mask1 = mask0;

  // Initialize reference
  if (refHolder_.InitRef(false, first, false, false, reftrajname, REF, RefParm,
                         mask1, actionArgs, "distrmsd"))
    return Action::ERR;
 
  // Set up the RMSD data set
  drmsd_ = DSL->AddSet(DataSet::DOUBLE, actionArgs.GetStringNext(),"DRMSD");
  if (drmsd_==0) return Action::ERR;
  // Add dataset to data file list
  if (outfile != 0) outfile->AddSet( drmsd_ );

  mprintf("    DISTRMSD: (%s), reference is %s\n",TgtMask_.MaskString(),
          refHolder_.RefModeString());

  return Action::OK;
}

// Action_DistRmsd::setup()
/** Called every time the trajectory changes. Set up TgtMask for the new 
  * parmtop and allocate space for selected atoms from the Frame.
  */
Action::RetType Action_DistRmsd::Setup(Topology* currentParm, Topology** parmAddress) {

  if ( currentParm->SetupIntegerMask(TgtMask_) ) return Action::ERR;
  if ( TgtMask_.None() ) {
    mprintf("    Error: DistRmsd::setup: No atoms in mask.\n");
    return Action::ERR;
  }
  // Allocate space for selected atoms in the frame. This will also put the
  // correct masses in based on the mask.
  SelectedTgt_.SetupFrameFromMask(TgtMask_, currentParm->Atoms());

  if (refHolder_.SetupRef(*currentParm, TgtMask_.Nselected(), "distrmsd"))
    return Action::ERR; 

  return Action::OK;
}

// Action_DistRmsd::action()
/** Called every time a frame is read in. Calc distance RMSD.
  * If first is true, set the first frame read in as reference.
  */
Action::RetType Action_DistRmsd::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {
  // Perform any needed reference actions
  refHolder_.ActionRef( *currentFrame, false, false );
  // Set selected frame atoms. Masses have already been set.
  SelectedTgt_.SetCoordinates(*currentFrame, TgtMask_);
  double DR = SelectedTgt_.DISTRMSD( refHolder_.SelectedRef() );
  drmsd_->Add(frameNum, &DR);
  return Action::OK;
}
