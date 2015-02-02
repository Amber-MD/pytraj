// Action_Average
#include "Action_Average.h"
#include "CpptrajStdio.h"
#include "Trajout.h"

// CONSTRUCTOR
Action_Average::Action_Average() :
  ensembleNum_(-1),
  debug_(0),
  AvgFrame_(0),
  Natom_(0),
  Nframes_(0)
{ } 

void Action_Average::Help() {
  mprintf("\t<filename> [<mask>] %s\n\t[TRAJOUT ARGS]\n"
          "  Calculate the average structure of atoms in <mask> over specified input frames.\n",
          ActionFrameCounter::HelpText);
}

// DESTRUCTOR
Action_Average::~Action_Average() {
  //fprintf(stderr,"Average Destructor.\n");
  if (AvgFrame_!=0) delete AvgFrame_;
}

// Action_Average::Init()
Action::RetType Action_Average::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  ensembleNum_ = DSL->EnsembleNum();
  debug_ = debugIn;
  // Get Keywords
  avgfilename_ = actionArgs.GetStringNext();
  if (avgfilename_.empty()) {
    mprinterr("Error: average: No filename given.\n");
    return Action::ERR;
  }
  // Get start/stop/offset args
  if (InitFrameCounter(actionArgs)) return Action::ERR;

  // Get Masks
  Mask1_.SetMaskString( actionArgs.GetMaskNext() );

  // Save all remaining arguments for setting up the trajectory at the end.
  trajArgs_ = actionArgs.RemainingArgs();

  mprintf("    AVERAGE: Averaging over coordinates in mask [%s]\n",Mask1_.MaskString());
  FrameCounterInfo();
  mprintf("\tWriting averaged coords to [%s]\n",avgfilename_.c_str());

  Nframes_ = 0;

  return Action::OK;
}

// Action_Average::Setup()
/** On first call, set up Frame according to first parmtop. This will be
  * used for coordinate output.
  * On subsequent calls, determine if the number of atoms is greater than or
  * less than the original # atoms. Never calculate more than the original
  * # atoms.
  */
Action::RetType Action_Average::Setup(Topology* currentParm, Topology** parmAddress) {

  if ( currentParm->SetupIntegerMask( Mask1_ ) ) return Action::ERR;

  if (Mask1_.None()) {
    mprinterr("Error: Cannot create average: No Atoms in mask.\n");
    return Action::ERR;
  }

  if (AvgFrame_==0) {
    mprintf("\tAveraging over %i atoms.\n",Mask1_.Nselected());
    AvgFrame_ = new Frame(Mask1_.Nselected());
    AvgFrame_->ZeroCoords();
    Natom_ = AvgFrame_->Natom(); // Initiially equal to Mask1.Nselected
    // AvgParm will be used for coordinate output
    // If the number of selected atoms is less than the current parm, strip
    // the parm for output purposes.
    if (Mask1_.Nselected()<currentParm->Natom()) {
      mprintf("Warning: Atom selection < natom, stripping parm for averaging only:\n");
      AvgParm_ = *(currentParm->modifyStateByMask(Mask1_));
      if (debug_ > 0)
        AvgParm_.Summary();
    } else 
      AvgParm_ = *currentParm; 
  } else {
    // If the frame is already set up, check to see if the current number
    // of atoms is bigger or smaller. If bigger, only average Natom coords.
    // If smaller, only average P->natom coords.
    if (Mask1_.Nselected() > AvgFrame_->Natom()) {
      Natom_ = AvgFrame_->Natom();
      mprintf("Warning: Average [%s]: Parm %s selected # atoms (%i) > original parm %s\n",
              avgfilename_.c_str(), currentParm->c_str(),
              Mask1_.Nselected(), AvgParm_.c_str());
      mprintf("Warning:   selected# atoms (%i).\n",AvgFrame_->Natom());
    } else if (Mask1_.Nselected() < AvgFrame_->Natom()) {
      Natom_ = Mask1_.Nselected();
      mprintf("Warning: Average[%s]: Parm %s selected # atoms (%i) < original parm %s\n",
              avgfilename_.c_str(), currentParm->c_str(), 
              Mask1_.Nselected(), AvgParm_.c_str());
      mprintf("Warning:   selected # atoms (%i).\n",AvgFrame_->Natom());
    } else {
      Natom_ = AvgFrame_->Natom();
    }
    mprintf("\t%i atoms will be averaged for '%s'.\n",Natom_, currentParm->c_str());
  }
        
  return Action::OK;  
}

// Action_Average::DoAction()
Action::RetType Action_Average::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) 
{
  if ( CheckFrameCounter( frameNum ) ) return Action::OK;

  if (AvgFrame_->AddByMask(*currentFrame, Mask1_)) return Action::ERR;
  ++Nframes_; 

  return Action::OK;
} 

// Action_Average::Print()
void Action_Average::Print() {
  Trajout outfile;
  double d_Nframes;

  if (Nframes_ < 1) return;
  d_Nframes = (double) Nframes_;
  AvgFrame_->Divide(d_Nframes);

  mprintf("    AVERAGE: [%s %s]\n",avgfilename_.c_str(), trajArgs_.ArgLine());

  if (outfile.InitEnsembleTrajWrite(avgfilename_, trajArgs_, &AvgParm_, 
                                    TrajectoryFile::UNKNOWN_TRAJ, ensembleNum_)) 
  {
    mprinterr("Error: AVERAGE: Could not set up %s for write.\n",avgfilename_.c_str());
    return;
  }

  outfile.PrintInfo(0);

  outfile.WriteFrame(0, &AvgParm_, *AvgFrame_);

  outfile.EndTraj();
}
