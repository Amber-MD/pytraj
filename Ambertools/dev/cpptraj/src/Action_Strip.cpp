// Action_Strip
#include "Action_Strip.h"
#include "CpptrajStdio.h"
#include "ParmFile.h"

// CONSTRUCTOR
Action_Strip::Action_Strip() :
  oldParm_(0),
  newParm_(0),
  removeBoxInfo_(false)
{
  //fprintf(stderr,"Strip Con\n");
} 

void Action_Strip::Help() {
  mprintf("\t<mask> [outprefix <name>] [nobox]\n"
          "  Strip atoms in <mask> from the system.\n");
}

void Action_Unstrip::Help() {
  mprintf("\tReturn to original topology/coordinates.\n");
}

// DESTRUCTOR
Action_Strip::~Action_Strip() {
  //fprintf(stderr,"Strip Des\n");
  // oldParm does not need dealloc because it is passed in
  if (newParm_!=0) delete newParm_;
}

// Action_Strip::Init()
Action::RetType Action_Strip::Init(ArgList& actionArgs, TopologyList* PFL, DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  // Get output stripped parm filename
  prefix_ = actionArgs.GetStringKey("outprefix");
  removeBoxInfo_ = actionArgs.hasKey("nobox");

  // Get mask of atoms to be stripped
  std::string mask1 = actionArgs.GetMaskNext();
  //mprintf("    Mask 1: %s\n",mask1);
  if (mask1.empty()) {
    mprinterr("Error: Requires atom mask.\n");
    return Action::ERR;
  }
  M1_.SetMaskString(mask1);
  // We want to strip the atoms inside the mask and keep those outside
  // the mask. Since modifyStateByMask needs to know the kept atoms,
  // invert the mask selection.
  M1_.InvertMask();

  mprintf("    STRIP: Stripping atoms in mask [%s]\n",M1_.MaskString());
  if (!prefix_.empty()) 
    mprintf("           Stripped topology will be output with prefix %s\n",prefix_.c_str());
  if (removeBoxInfo_)
    mprintf("           Any existing box information will be removed.\n");

  return Action::OK;
}

// Action_Strip::Setup()
/** Attempt to create a new stripped down version of the input parmtop
  */
Action::RetType Action_Strip::Setup(Topology* currentParm, Topology** parmAddress) {
  if (currentParm->SetupIntegerMask( M1_ )) return Action::ERR;
  if (M1_.None()) {
    mprintf("Warning: strip: Mask [%s] has no atoms.\n",M1_.MaskString());
    return Action::ERR;
  }
  int numStripped = currentParm->Natom() - M1_.Nselected();
  mprintf("\tStripping %i atoms.\n", numStripped);
  // If no atoms will be stripped, no need to use this command. SKIP
  if ( numStripped == 0 ) {
    mprintf("Warning: No atoms to strip. Skipping 'strip' for topology '%s'\n",
            currentParm->c_str());
    return Action::ERR;
  }

  // Store old parm
  oldParm_ = currentParm;

  // Attempt to create new parmtop based on mask
  if (newParm_!=0) delete newParm_;
  newParm_ = currentParm->modifyStateByMask(M1_);
  if (newParm_==0) {
    mprinterr("Error: Could not create new parmtop.\n");
    return Action::ERR;
  }
  // Remove box information if asked
  if (removeBoxInfo_)
    newParm_->SetParmBox( Box() ); 
  newParm_->Brief("Stripped parm:");
  // Allocate space for new frame
  newFrame_.SetupFrameV(newParm_->Atoms(), newParm_->ParmCoordInfo());

  // If prefix given then output stripped parm
  if (!prefix_.empty()) {
    ParmFile pfile;
    if ( pfile.WritePrefixTopology( *newParm_, prefix_, ParmFile::AMBERPARM, 0 ) )
      mprinterr("Error: Could not write out stripped parm file.\n");
  }

  // Set parm
  *parmAddress = newParm_;

  return Action::OK;  
}

// Action_Strip::DoAction()
/** Modify the coordinate frame to reflect stripped parmtop. */
Action::RetType Action_Strip::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {

  newFrame_.SetFrame(*currentFrame, M1_);

  // Set frame
  *frameAddress = &newFrame_;

  return Action::OK;
} 
