// Action_DihedralScan
#include <cmath>
#include <ctime> // time
#include "Action_DihedralScan.h"
#include "CpptrajStdio.h"
#include "Constants.h" // DEGRAD
#include "DistRoutines.h"
#include "TorsionRoutines.h"
// Activate DEBUG info
//#define DEBUG_DIHEDRALSCAN

// CONSTRUCTOR
Action_DihedralScan::Action_DihedralScan() :
  mode_(INTERVAL),
  outframe_(0),
  interval_(60.0),
  maxVal_(0),
  check_for_clashes_(false),
  checkAllResidues_(false),
  max_rotations_(0),
  max_factor_(2),
  cutoff_(0.64), // 0.8^2
  rescutoff_(100.0), // 10.0^2
  backtrack_(5),
  increment_(1),
  max_increment_(360),
  debug_(0),
  CurrentParm_(0),
  number_of_problems_(0)
{} 

Action_DihedralScan::~Action_DihedralScan() {
  outtraj_.EndTraj();
}

void Action_DihedralScan::Help() {
  mprintf("\tresrange <range> [{interval*|random}]");
  DihedralSearch::ListKnownTypes();
  mprintf("\t'*' denotes default.\n"
          "\tOptions for 'random': [rseed <rseed>]\n"
          "\t\t[ check [cutoff <cutoff>] [rescutoff <rescutoff>] [checkallresidues]\n"
          "\t\t  [backtrack <backtrack>] [increment <increment>] [maxfactor <max_factor>] ]\n"
          "\tOptions for 'interval': <interval deg> [outtraj <filename> [<outfmt>]]\n"
          "  Rotate specified dihedral(s) by specific intervals or by random values.\n");
}

// Action_DihedralScan::Init()
Action::RetType Action_DihedralScan::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  if (DSL->EnsembleNum() > -1) {
    mprinterr("Error: DIHEDRALSCAN currently cannot be used in ensemble mode.\n");
    return Action::ERR;
  }
  TrajectoryFile::TrajFormatType outfmt = TrajectoryFile::UNKNOWN_TRAJ;
  Topology* outtop = 0;
  int iseed = -1;

  debug_ = debugIn;
  // Get Keywords - first determine mode
  if (actionArgs.hasKey("random"))
    mode_ = RANDOM;
  else if (actionArgs.hasKey("interval"))
    mode_ = INTERVAL;
  else
    mode_ = INTERVAL;
  // Get residue range
  resRange_.SetRange(actionArgs.GetStringKey("resrange"));
  if (!resRange_.Empty())
    resRange_.ShiftBy(-1); // User res args start from 1
  // Determine which angles to search for
  dihSearch_.SearchForArgs(actionArgs);
  // If nothing is enabled, enable all 
  dihSearch_.SearchForAll();
  // For interval, get value to shift by, set max rotations and set up 
  // output trajectory.
  if ( mode_ == INTERVAL ) { 
    interval_ = actionArgs.getNextDouble(60.0);
    maxVal_ = (int) (360.0 / interval_);
    if (maxVal_ < 0) maxVal_ = -maxVal_;
    outfilename_ = actionArgs.GetStringKey("outtraj");
    if (!outfilename_.empty()) {
      outfmt = TrajectoryFile::GetFormatFromArg( actionArgs );
      outtop = PFL->GetParm( actionArgs );
      if (outtop == 0) {
        mprinterr("Error: dihedralscan: No topology for output traj.\n");
        return Action::ERR;
      }
    }
  }
  // Get 'random' args
  if (mode_ == RANDOM) {
    check_for_clashes_ = actionArgs.hasKey("check");
    checkAllResidues_ = actionArgs.hasKey("checkallresidues");
    cutoff_ = actionArgs.getKeyDouble("cutoff",0.8);
    rescutoff_ = actionArgs.getKeyDouble("rescutoff",10.0);
    backtrack_ = actionArgs.getKeyInt("backtrack",4);
    increment_ = actionArgs.getKeyInt("increment",1);
    max_factor_ = actionArgs.getKeyInt("maxfactor",2);
    // Check validity of args
    if (cutoff_ < Constants::SMALL) {
      mprinterr("Error: cutoff too small.\n");
      return Action::ERR;
    }
    if (rescutoff_ < Constants::SMALL) {
      mprinterr("Error: rescutoff too small.\n");
      return Action::ERR;
    }
    if (backtrack_ < 0) {
      mprinterr("Error: backtrack value must be >= 0\n");
      return Action::ERR;
    }
    if ( increment_<1 || (360 % increment_)!=0 ) {
      mprinterr("Error: increment must be a factor of 360.\n");
      return Action::ERR;
    }
    // Calculate max increment
    max_increment_ = 360 / increment_;
    // Seed random number gen
    iseed = actionArgs.getKeyInt("rseed",-1);
    RN_.rn_set( iseed );
  }
  // Output file for # of problems
  std::string problemFile = actionArgs.GetStringKey("out");
  // Dataset to store number of problems
  number_of_problems_ = DSL->AddSet(DataSet::INTEGER, actionArgs.GetStringNext(),"Nprob");
  if (number_of_problems_==0) return Action::ERR;
  // Add dataset to data file list
  DFL->AddSetToFile(problemFile,number_of_problems_);

  mprintf("    DIHEDRALSCAN: Dihedrals in");
  if (resRange_.Empty())
    mprintf(" all solute residues.\n");
  else
    mprintf(" residue range [%s]\n", resRange_.RangeArg());
  switch (mode_) {
    case RANDOM:
      mprintf("\tDihedrals will be rotated to random values.\n");
      if (iseed==-1)
        mprintf("\tRandom number generator will be seeded using time.\n");
      else
        mprintf("\tRandom number generator will be seeded using %i\n",iseed);
      if (check_for_clashes_) {
        mprintf("\tWill attempt to recover from bad steric clashes.\n");
        if (checkAllResidues_)
          mprintf("\tAll residues will be checked.\n");
        else
          mprintf("\tResidues up to the currenly rotating dihedral will be checked.\n");
        mprintf("\tAtom cutoff %.2f, residue cutoff %.2f, backtrack = %i\n",
                cutoff_, rescutoff_, backtrack_);
        mprintf("\tWhen clashes occur dihedral will be incremented by %i\n",increment_);
        mprintf("\tMax # attempted rotations = %i times number dihedrals.\n",
                max_factor_);
      }
      break;
    case INTERVAL:
      mprintf("\tDihedrals will be rotated at intervals of %.2f degrees.\n",
              interval_);
      if (!outfilename_.empty())
        mprintf("\tCoordinates output to %s, format %s\n",outfilename_.c_str(), 
                TrajectoryFile::FormatString(outfmt));
      break;
  }
  // Setup output trajectory
  if (!outfilename_.empty()) {
    if (outtraj_.InitTrajWrite(outfilename_, outtop, outfmt)) return Action::ERR;
    outframe_ = 0;
  } 
  // Square cutoffs to compare to dist^2 instead of dist
  cutoff_ *= cutoff_;
  rescutoff_ *= rescutoff_;
  // Increment backtrack by 1 since we need to skip over current res
  ++backtrack_;
  // Initialize CheckStructure
  ArgList cs_args("noimage nobondcheck");
  if (checkStructure_.Init( cs_args, PFL, FL, DSL, DFL, debug_) != Action::OK) {
    mprinterr("Error: Could not set up structure check for DIHEDRALSCAN.\n");
    return Action::ERR;
  }
  return Action::OK;
}

// Action_DihedralScan::Setup()
/** Determine from selected mask atoms which dihedrals will be rotated. */
Action::RetType Action_DihedralScan::Setup(Topology* currentParm, Topology** parmAddress) {
  DihedralScanType dst;
  // If range is empty (i.e. no resrange arg given) look through all 
  // solute residues.
  Range actualRange;
  if (resRange_.Empty())
    actualRange = currentParm->SoluteResidues();
  else 
    actualRange = resRange_;
  // Search for dihedrals
  if (dihSearch_.FindDihedrals(*currentParm, actualRange))
    return Action::ERR;
  // For each found dihedral, set up mask of atoms that will move upon 
  // rotation. Also set up mask of atoms in this residue that will not
  // move, including atom2.
  if (debug_>0)
    mprintf("DEBUG: Dihedrals:\n");
  AtomMask cMask;
  for (DihedralSearch::mask_it dih = dihSearch_.begin();
                               dih != dihSearch_.end(); ++dih)
  {
    dst.checkAtoms.clear();
    // Set mask of atoms that will move during dihedral rotation.
    dst.Rmask = DihedralSearch::MovingAtoms(*currentParm, (*dih).A1(), (*dih).A2());
    // If randomly rotating angles, check for atoms that are in the same
    // residue as A1 but will not move. They need to be checked for clashes
    // since further rotations will not help them.
    if (mode_ == RANDOM && check_for_clashes_) {
      cMask = dst.Rmask;
      cMask.ConvertToCharMask(); 
      int a1res = (*currentParm)[(*dih).A1()].ResNum();
      for (int maskatom = currentParm->Res(a1res).FirstAtom();
               maskatom < currentParm->Res(a1res).LastAtom(); ++maskatom)
        if (!cMask.AtomInCharMask(maskatom))
          dst.checkAtoms.push_back( maskatom );
      dst.checkAtoms.push_back((*dih).A1()); // TODO: Does this need to be added first?
      // Since only the second atom and atoms it is bonded to move during 
      // rotation, base the check on the residue of the second atom.
      dst.resnum = a1res;
    }
    dst.atom0 = (*dih).A0(); // FIXME: This duplicates info
    dst.atom1 = (*dih).A1();
    dst.atom2 = (*dih).A2();
    dst.atom3 = (*dih).A3();
    BB_dihedrals_.push_back(dst);
    // DEBUG: List dihedral info.
    if (debug_ > 0) {
      mprintf("\t%s-%s-%s-%s\n", 
              currentParm->TruncResAtomName((*dih).A0()).c_str(),
              currentParm->TruncResAtomName((*dih).A1()).c_str(),
              currentParm->TruncResAtomName((*dih).A2()).c_str(),
              currentParm->TruncResAtomName((*dih).A3()).c_str() );
      if (debug_ > 1 && mode_ == RANDOM && check_for_clashes_) {
        mprintf("\t\tCheckAtoms=");
        for (std::vector<int>::iterator ca = dst.checkAtoms.begin();
                                        ca != dst.checkAtoms.end(); ++ca)
          mprintf(" %i", *ca + 1);
        mprintf("\n");
      }
      if (debug_ > 2) {
        mprintf("\t\t");
        dst.Rmask.PrintMaskAtoms("Rmask:");
      }
    }
  }

  // Set up CheckStructure for this parm
  if (checkStructure_.Setup(currentParm, parmAddress) != Action::OK)
    return Action::ERR;

  // Set the overall max number of rotations to try
  max_rotations_ = (int) BB_dihedrals_.size();
  max_rotations_ *= max_factor_;

  // Set up simple structure check. First step is coarse; check distances 
  // between a certain atom in each residue (first, COM, CA, some other atom?)
  // to see if residues are in each others neighborhood. Second step is to 
  // check the atoms in each close residue.
  if (check_for_clashes_) {
    ResidueCheckType rct;
    int res = 0;
    for (Topology::res_iterator residue = currentParm->ResStart();
                                residue != currentParm->ResEnd(); ++residue)
    {
      rct.resnum = res++;
      rct.start = (*residue).FirstAtom();
      rct.stop = (*residue).LastAtom();
      rct.checkatom = rct.start;
      ResCheck_.push_back(rct);
    }
  }
  CurrentParm_ = currentParm;
  return Action::OK;  
}

// Action_DihedralScan::CheckResidue()
/** \return 1 if a new dihedral should be tried, 0 if no clashes
  * \return -1 if further rotations will not help.
  */
int Action_DihedralScan::CheckResidue( Frame const& FrameIn, DihedralScanType const& dih, 
                                       int nextres, double *clash ) 
{
  int resnumIn = dih.resnum;
  int rstart = ResCheck_[ resnumIn ].start;
  int rstop = ResCheck_[ resnumIn ].stop;
  int rcheck = ResCheck_[ resnumIn ].checkatom;
  // Check for clashes with self
#ifdef DEBUG_DIHEDRALSCAN
  mprintf("\tChecking residue %i\n",resnumIn+1);
  mprintf("\tATOMS %i to %i\n",rstart+1,rstop);
#endif
  for (int atom1 = rstart; atom1 < rstop - 1; atom1++) {
    for (int atom2 = atom1 + 1; atom2 < rstop; atom2++) {
      // Skip bonded atoms
      bool isBonded = false;
      for (Atom::bond_iterator bndatm = (*CurrentParm_)[atom1].bondbegin();
                               bndatm != (*CurrentParm_)[atom1].bondend(); ++bndatm)
        if (*bndatm == atom2) {
          isBonded = true;
          break;
        }
      if (!isBonded) {
        double atomD2 = DIST2_NoImage(FrameIn.XYZ(atom1), FrameIn.XYZ(atom2));
        if (atomD2 < cutoff_) {
#         ifdef DEBUG_DIHEDRALSCAN 
          mprintf("\t\tCurrent Res %i Atoms %s and %s are close (%.3lf)\n", resnumIn+1, 
                  CurrentParm_->AtomMaskName(atom1).c_str(),
                  CurrentParm_->AtomMaskName(atom2).c_str(), sqrt(atomD2));
#         endif
          *clash = atomD2;
          return 1;
        }
      }
    }
  }
  // Check for clashes with previous residues, as well as clashes up to and
  // including the next residue in which a dihedral will be rotated.
  for (int res = 0; res <= nextres; res++) {
    if (res == resnumIn) continue;
    int rstart2 = ResCheck_[ res ].start;
    int rstop2 = ResCheck_[ res ].stop;
    int rcheck2 = ResCheck_[ res ].checkatom;
    double resD2 = DIST2_NoImage(FrameIn.XYZ(rcheck), FrameIn.XYZ(rcheck2));
    // If residues are close enough check each atom
    if (resD2 < rescutoff_) { 
#ifdef DEBUG_DIHEDRALSCAN
      mprintf("\tRES %i ATOMS %i to %i\n",res+1,rstart2+2,rstop2);
#endif
      for (int atom1 = rstart; atom1 < rstop; atom1++) {
        for (int atom2 = rstart2; atom2 < rstop2; atom2++) {
          double D2 = DIST2_NoImage(FrameIn.XYZ(atom1), FrameIn.XYZ(atom2));
          if (D2 < cutoff_) {
#ifdef DEBUG_DIHEDRALSCAN
            mprintf("\t\tResCheck %i Atoms %s and %s are close (%.3lf)\n", res+1,
                    CurrentParm_->TruncResAtomName(atom1).c_str(),
                    CurrentParm_->TruncResAtomName(atom2).c_str(), sqrt(D2));
#endif
            *clash = D2;
            // If the clash involves any atom that will not be moved by further
            // rotation, indicate it is not possible to resolve clash by
            // more rotation by returning -1.
            //if (atom1 == dih.atom2 || atom1 == dih.atom1) return -1;
            for (std::vector<int>::const_iterator ca = dih.checkAtoms.begin();
                                                  ca != dih.checkAtoms.end(); ca++) 
            {
              if (atom1 == *ca) return -1;
            }
            return 1;
          }
        }
      }
    }
  }
  return 0;
} 

// Action_DihedralScan::RandomizeAngles()
void Action_DihedralScan::RandomizeAngles(Frame& currentFrame) {
  Matrix_3x3 rotationMatrix;
#ifdef DEBUG_DIHEDRALSCAN
  // DEBUG
  int debugframenum=0;
  Trajout DebugTraj;
  DebugTraj.InitTrajWrite("debugtraj.nc",CurrentParm_,TrajectoryFile::AMBERNETCDF);
  DebugTraj.WriteFrame(debugframenum++,CurrentParm_,currentFrame);
#endif
  int next_resnum;
  int bestLoop = 0;
  int number_of_rotations = 0;

  std::vector<DihedralScanType>::const_iterator next_dih = BB_dihedrals_.begin();
  next_dih++;
  for (std::vector<DihedralScanType>::const_iterator dih = BB_dihedrals_.begin();
                                                     dih != BB_dihedrals_.end(); 
                                                     ++dih, ++next_dih)
  {
    ++number_of_rotations;
    // Get the residue atom of the next dihedral. Residues up to and
    // including this residue will be checked for bad clashes 
    if (next_dih!=BB_dihedrals_.end()) 
      next_resnum = next_dih->resnum;
    else
      next_resnum = dih->resnum-1;
    // Set axis of rotation
    Vec3 axisOfRotation = currentFrame.SetAxisOfRotation(dih->atom1, dih->atom2);
    // Generate random value to rotate by in radians
    // Guaranteed to rotate by at least 1 degree.
    // NOTE: could potentially rotate 360 - prevent?
    // FIXME: Just use 2PI and rn_gen, get everything in radians
    double theta_in_degrees = ((int)(RN_.rn_gen()*100000) % 360) + 1;
    double theta_in_radians = theta_in_degrees * Constants::DEGRAD;
    // Calculate rotation matrix for random theta
    rotationMatrix.CalcRotationMatrix(axisOfRotation, theta_in_radians);
    int loop_count = 0;
    double clash = 0;
    double bestClash = 0;
    if (debug_>0) mprintf("DEBUG: Rotating dihedral %zu res %8i:\n",dih - BB_dihedrals_.begin(),
                          dih->resnum+1);
    bool rotate_dihedral = true;
    while (rotate_dihedral) {
      if (debug_>0) {
        mprintf("\t%8i %12s %12s, +%.2lf degrees (%i).\n",dih->resnum+1,
                CurrentParm_->AtomMaskName(dih->atom1).c_str(),
                CurrentParm_->AtomMaskName(dih->atom2).c_str(),
                theta_in_degrees,loop_count);
      }
      // Rotate around axis
      currentFrame.Rotate(rotationMatrix, dih->Rmask);
#ifdef DEBUG_DIHEDRALSCAN
      // DEBUG
      DebugTraj.WriteFrame(debugframenum++,CurrentParm_,currentFrame);
#endif
      // If we dont care about sterics exit here
      if (!check_for_clashes_) break;
      // Check resulting structure for issues
      int checkresidue;
      if (!checkAllResidues_)
        checkresidue = CheckResidue(currentFrame, *dih, next_resnum, &clash);
      else
        checkresidue = CheckResidue(currentFrame, *dih, CurrentParm_->Nres(), &clash);
      if (checkresidue==0)
        rotate_dihedral = false;
      else if (checkresidue==-1) {
        if (dih - BB_dihedrals_.begin() < 2) {
          mprinterr("Error: Cannot backtrack; initial structure already has clashes.\n");
          number_of_rotations = max_rotations_ + 1;
        } else {
          dih--; //  0
          dih--; // -1
          next_dih = dih;
          next_dih++;
          if (debug_>0)
            mprintf("\tCannot resolve clash with further rotations, trying previous again.\n");
        }
        break;
      }
      if (clash > bestClash) {bestClash = clash; bestLoop = loop_count;}
      //n_problems = CheckResidues( currentFrame, second_atom );
      //if (n_problems > -1) {
      //  mprintf("%i\tCheckResidues: %i problems.\n",frameNum,n_problems);
      //  rotate_dihedral = false;
      //} else if (loop_count==0) {
      if (loop_count==0 && rotate_dihedral) {
        if (debug_>0)
          mprintf("\tTrying dihedral increments of +%i\n",increment_);
        // Instead of a new random dihedral, try increments
        theta_in_degrees = (double)increment_;
        theta_in_radians = theta_in_degrees * Constants::DEGRAD;
        // Calculate rotation matrix for new theta
        rotationMatrix.CalcRotationMatrix(axisOfRotation, theta_in_radians);
      }
      ++loop_count;
      if (loop_count == max_increment_) {
        if (debug_>0)
          mprintf("%i iterations! Best clash= %.3lf at %i\n",max_increment_,
                  sqrt(bestClash),bestLoop);
        if (dih - BB_dihedrals_.begin() < backtrack_) {
          mprinterr("Error: Cannot backtrack; initial structure already has clashes.\n");
          number_of_rotations = max_rotations_ + 1;
        } else { 
          for (int bt = 0; bt < backtrack_; bt++)
            dih--;
          next_dih = dih;
          next_dih++;
          if (debug_>0)
            mprintf("\tCannot resolve clash with further rotations, trying previous %i again.\n",
                    backtrack_ - 1);
        }
        break;
        // Calculate how much to rotate back in order to get to best clash
        /*int num_back = bestLoop - 359;
        theta_in_degrees = (double) num_back;
        theta_in_radians = theta_in_degrees * Constants::DEGRAD;
        // Calculate rotation matrix for theta
        calcRotationMatrix(rotationMatrix, axisOfRotation, theta_in_radians);
        // Rotate back to best clash
        currentFrame->RotateAroundAxis(rotationMatrix, theta_in_radians, dih->Rmask);
        // DEBUG
        DebugTraj.WriteFrame(debugframenum++,currentParm,*currentFrame);
        // Sanity check
        CheckResidue(currentFrame, *dih, second_atom, &clash);
        rotate_dihedral=false;*/
        //DebugTraj.EndTraj();
        //return 1;
      }
    } // End dihedral rotation loop
    // Safety valve - number of defined dihedrals times * maxfactor
    if (number_of_rotations > max_rotations_) {
      mprinterr("Error: DihedralScan: # of rotations (%i) exceeds max rotations (%i), exiting.\n",
                number_of_rotations, max_rotations_);
//#ifdef DEBUG_DIHEDRALSCAN
//      DebugTraj.EndTraj();
//#endif
      // Return gracefully for now
      break;
      //return 1;
    }
  } // End loop over dihedrals
#ifdef DEBUG_DIHEDRALSCAN
  DebugTraj.EndTraj();
  mprintf("\tNumber of rotations %i, expected %u\n",number_of_rotations,BB_dihedrals_.size());
#endif
}

// Action_DihedralScan::IntervalAngles()
void Action_DihedralScan::IntervalAngles(Frame& currentFrame) {
  Matrix_3x3 rotationMatrix;
  double theta_in_radians = interval_ * Constants::DEGRAD;
  // Write original frame
  if (!outfilename_.empty())
    outtraj_.WriteFrame(outframe_++, CurrentParm_, currentFrame);
  for (std::vector<DihedralScanType>::iterator dih = BB_dihedrals_.begin();
                                               dih != BB_dihedrals_.end();
                                               dih++)
  {
    // Set axis of rotation
    Vec3 axisOfRotation = currentFrame.SetAxisOfRotation((*dih).atom1, (*dih).atom2);
    // Calculate rotation matrix for interval 
    rotationMatrix.CalcRotationMatrix(axisOfRotation, theta_in_radians);
    if (debug_ > 0) {
      mprintf("\tRotating Dih %s-%s by %.2f deg %i times.\n",
               CurrentParm_->TruncResAtomName( (*dih).atom1 ).c_str(), 
               CurrentParm_->TruncResAtomName( (*dih).atom2 ).c_str(), interval_, maxVal_); 
    }
    for (int rot = 0; rot < maxVal_; ++rot) {
      // Rotate around axis
      currentFrame.Rotate(rotationMatrix, (*dih).Rmask);
      // Write output trajectory
      if (outtraj_.TrajIsOpen())
        outtraj_.WriteFrame(outframe_++, CurrentParm_, currentFrame);
    }
  }
}

// Action_DihedralScan::DoAction()
Action::RetType Action_DihedralScan::DoAction(int frameNum, Frame* currentFrame, 
                                              Frame** frameAddress) 
{
  switch (mode_) {
    case RANDOM: RandomizeAngles(*currentFrame); break;
    case INTERVAL: IntervalAngles(*currentFrame); break;
  }
  // Check the resulting structure
  int n_problems = checkStructure_.CheckFrame( frameNum+1, *currentFrame );
  //mprintf("%i\tResulting structure has %i problems.\n",frameNum,n_problems);
  number_of_problems_->Add(frameNum, &n_problems);

  return Action::OK;
} 
