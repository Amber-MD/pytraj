// Action_Mask
#include "Action_Mask.h"
#include "CpptrajStdio.h"
#include "Trajout.h"

// CONSTRUCTOR
Action_Mask::Action_Mask() :
  CurrentParm_(0),
  debug_(0),
  trajFmt_(TrajectoryFile::PDBFILE),
  trajOpt_(0)
{ } 

void Action_Mask::Help() {
  mprintf("\t<mask1> [maskout <filename>] [maskpdb <filename> | maskmol2 <filename>]\n"
          "  Print atoms selected by <mask1> to file specified by 'maskout' and/or\n"
          "  the PDB or Mol2 file specified by 'maskpdb' or 'maskmol2'. Good for\n"
          "  distance-based masks.\n");
}

// Action_Mask::Init()
// NOTE: Could also split the arglist at maskpdb and make it so any type of 
//       file can be written out.
Action::RetType Action_Mask::Init(ArgList& actionArgs, TopologyList* PFL, DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  debug_ = debugIn;
  // Get Keywords
  std::string maskFilename = actionArgs.GetStringKey("maskout");
  maskpdb_ = actionArgs.GetStringKey("maskpdb");
  std::string maskmol2 = actionArgs.GetStringKey("maskmol2");
  if (!maskpdb_.empty()) {
    trajFmt_ = TrajectoryFile::PDBFILE;
    // Set pdb output options: multi so that 1 file per frame is written; dumpq
    // so that charges are written out.
    trajOpt_ = "multi dumpq nobox";
  } else if (!maskmol2.empty()) {
    maskpdb_ = maskmol2;
    trajFmt_ = TrajectoryFile::MOL2FILE;
    trajOpt_ = "multi nobox";
  }
  // Get Mask
  Mask1_.SetMaskString( actionArgs.GetMaskNext() );

  mprintf("    ACTIONMASK: Information on atoms in mask %s will be printed",
          Mask1_.MaskString());
  if (!maskFilename.empty())
    mprintf(" to file %s",maskFilename.c_str());
  mprintf(".\n");
  if (!maskpdb_.empty()) 
    mprintf("\t%ss of atoms in mask will be written to %s.X\n",
            TrajectoryFile::FormatString(trajFmt_), maskpdb_.c_str());

  // Open output file
  if (!maskFilename.empty()) {
    if ( outfile_.OpenEnsembleWrite( maskFilename, DSL->EnsembleNum() ) )
      return Action::ERR;
      // Header
    outfile_.Printf("%-8s %8s %4s %8s %4s %8s\n","#Frame","AtomNum","Atom",
                    "ResNum","Res", "MolNum");
  }

  return Action::OK;
}

// Action_Mask::Setup()
Action::RetType Action_Mask::Setup(Topology* currentParm, Topology** parmAddress) {
  CurrentParm_ = currentParm;
  return Action::OK;
}

// Action_Mask::DoAction()
Action::RetType Action_Mask::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {
  // Get atom selection
  if ( CurrentParm_->SetupCharMask(Mask1_, *currentFrame) ) {
    mprintf("Warning: Could not set up atom mask [%s]\n",
            Mask1_.MaskString());
    return Action::ERR;
  }
  // Print out information for every atom in the mask
  for (int atom=0; atom < CurrentParm_->Natom(); atom++) {
    if (Mask1_.AtomInCharMask(atom)) {
      int res = (*CurrentParm_)[atom].ResNum();
      if (outfile_.IsOpen())
        outfile_.Printf("%8i %8i %4s %8i %4s %8i\n", frameNum+1,
                        atom+1, (*CurrentParm_)[atom].c_str(), res+1,
                        CurrentParm_->Res(res).c_str(), (*CurrentParm_)[atom].MolNum()+1);
      /*mprintf(" Type=%4s",CurrentParm_->types[atom]);
      mprintf(" Charge=%lf",CurrentParm_->charge[atom]);
      mprintf(" Mass=%lf",CurrentParm_->mass[atom]);
      outfile.Printf("\n");*/
    }
  }

  // Optional PDB write out of selected atoms for the frame.
  if (!maskpdb_.empty()) {
    Trajout pdbout;
    // Convert Mask1 to an integer mask for use in parm/frame functions
    AtomMask Mask2 = Mask1_;
    Mask2.ConvertToIntMask();
    // Create new parm and frame based on atoms in Mask. Since we dont care
    // about advanced parm info for PDB write just do a partial modify.
    Topology* pdbParm = CurrentParm_->partialModifyStateByMask(Mask2);
    //pdbParm->Summary(); // DEBUG
    Frame pdbFrame(*currentFrame, Mask2);
    // Set up output file. 
    pdbout.SetDebug(debug_);
    // Set pdb output options: multi so that 1 file per frame is written; dumpq
    // so that charges are written out. 
    if (pdbout.InitTrajWriteWithArgs(maskpdb_,trajOpt_,pdbParm,trajFmt_)) 
    {
      mprinterr("Error: %s: Could not set up for write of frame %i.\n",
                maskpdb_.c_str(),frameNum);
    } else {
      if (debug_ > 0) pdbout.PrintInfo(0);
      pdbout.WriteFrame(frameNum,pdbParm,pdbFrame);
      pdbout.EndTraj();
    }
    delete pdbParm;
  }

  return Action::OK;
} 

// Action_Mask::Print()
/** Close the output file. */
void Action_Mask::Print() {
  outfile_.CloseFile();
}
