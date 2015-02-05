#include <cmath> // sqrt
#include "Action_Diffusion.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_Diffusion::Action_Diffusion() :
  printIndividual_(false),
  time_(1),
  hasBox_(false),
  debug_(0)
{
  boxcenter_[0] = 0;
  boxcenter_[1] = 0;
  boxcenter_[2] = 0;
}

void Action_Diffusion::Help() {
  mprintf("\t<mask> <time per frame> [average] [<prefix>]\n"
          "  Compute a mean square displacement plot for the atoms in the mask.\n"
          "  The following files are produced:\n"
          "    <prefix>_x.xmgr: Mean square displacement(s) in the X direction (in Å^2).\n"
          "    <prefix>_y.xmgr: Mean square displacement(s) in the Y direction (in Å^2).\n"
          "    <prefix>_z.xmgr: Mean square displacement(s) in the Z direction (in Å^2).\n"
          "    <prefix>_r.xmgr: Overall mean square displacement(s) (in Å^2).\n"
          "    <prefix>_a.xmgr: Total distance travelled (in Å).\n");
}

Action::RetType Action_Diffusion::Init(ArgList& actionArgs, TopologyList* PFL, DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  debug_ = debugIn;
  printIndividual_ = !(actionArgs.hasKey("average"));
  mask_.SetMaskString( actionArgs.GetMaskNext() );
  time_ = actionArgs.getNextDouble(1.0);
  if (time_ < 0) {
    mprinterr("Error: diffusion: time per frame incorrectly specified\n");
    return Action::ERR;
  }
  std::string outputNameRoot = actionArgs.GetStringNext();
  // Default filename: 'diffusion'
  if (outputNameRoot.empty()) 
    outputNameRoot.assign("diffusion");
  
  // Open output files
  std::string fname = outputNameRoot + "_x.xmgr";
  if (outputx_.OpenEnsembleWrite( fname, DSL->EnsembleNum() )) return Action::ERR;
  fname = outputNameRoot + "_y.xmgr";
  if (outputy_.OpenEnsembleWrite( fname, DSL->EnsembleNum() )) return Action::ERR;
  fname = outputNameRoot + "_z.xmgr";
  if (outputz_.OpenEnsembleWrite( fname, DSL->EnsembleNum() )) return Action::ERR;
  fname = outputNameRoot + "_r.xmgr";
  if (outputr_.OpenEnsembleWrite( fname, DSL->EnsembleNum() )) return Action::ERR;
  fname = outputNameRoot + "_a.xmgr";
  if (outputa_.OpenEnsembleWrite( fname, DSL->EnsembleNum() )) return Action::ERR;

  mprintf("    DIFFUSION:\n");
  mprintf("\tAtom Mask is [%s]\n", mask_.MaskString());
  if (printIndividual_)
    mprintf("\tThe average and individual results will be printed to:\n");
  else
    mprintf("\tOnly the average results will be printed to:\n");
  const char* onr = outputNameRoot.c_str();
  mprintf("\t  %s_x.xmgr: Mean square displacement(s) in the X direction (in Å^2).\n"
          "\t  %s_y.xmgr: Mean square displacement(s) in the Y direction (in Å^2).\n"
          "\t  %s_z.xmgr: Mean square displacement(s) in the Z direction (in Å^2).\n"
          "\t  %s_r.xmgr: Overall mean square displacement(s) (in Å^2).\n"
          "\t  %s_a.xmgr: Total distance travelled (in Å).\n",
          onr, onr, onr, onr, onr);
  mprintf("\tThe time between frames in ps is %.3f.\n", time_);
  mprintf("\tTo calculate diffusion constants from a mean squared displacement plot\n"
          "\t(i.e. {_x|_y|_z|_r}.xmgr), calculate the slope of the line and multiply\n"
          "\tby 10.0/6.0; this will give units of 1x10^-5 cm^2/s\n");

  return Action::OK;
}

Action::RetType Action_Diffusion::Setup(Topology* currentParm, Topology** parmAddress) {
  // Setup atom mask
  if (currentParm->SetupIntegerMask( mask_ )) return Action::ERR;
  if (mask_.None()) {
    mprinterr("Error: diffusion: No atoms selected.\n");
    return Action::ERR;
  }

  // Check for box
  if ( currentParm->BoxType() != Box::NOBOX ) {
    // Currently only works for orthogonal boxes
    if ( currentParm->BoxType() != Box::ORTHO ) {
      mprintf("Warning: diffusion command currently only works with orthogonal boxes.\n");
      mprintf("Warning: imaging will be disabled for this command. This may result in\n");
      mprintf("Warning: large jumps if target molecule is imaged. To counter this please\n");
      mprintf("Warning: use the 'unwrap' command prior to 'diffusion'.\n");
      hasBox_ = false;
    } else 
      hasBox_ = true;
  } else
    hasBox_ = false;

  // Allocate the distance arrays
  distancex_.resize( mask_.Nselected() );
  distancey_.resize( mask_.Nselected() );
  distancez_.resize( mask_.Nselected() );
  distance_.resize(  mask_.Nselected() );

  // Allocate the delta arrays
  deltax_.assign( mask_.Nselected(), 0 );
  deltay_.assign( mask_.Nselected(), 0 );
  deltaz_.assign( mask_.Nselected(), 0 );

  // Reserve space for the initial and previous frame arrays
  //initial_.reserve( mask_.Nselected() );
  previousx_.reserve( mask_.Nselected() );
  previousy_.reserve( mask_.Nselected() );
  previousz_.reserve( mask_.Nselected() );

  // If initial frame already set and current # atoms > # atoms in initial
  // frame this will probably cause an error.
  if (!initial_.empty() && currentParm->Natom() > initial_.Natom()) {
    mprintf("Warning: # atoms in current parm (%s, %i) > # atoms in initial frame (%i)\n",
             currentParm->c_str(), currentParm->Natom(), initial_.Natom());
    mprintf("Warning: This may lead to segmentation faults.\n");
  }

  return Action::OK;
}

Action::RetType Action_Diffusion::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {
  // Load initial frame if necessary
  if (initial_.empty()) {
    initial_ = *currentFrame;
    for (AtomMask::const_iterator atom = mask_.begin(); atom != mask_.end(); ++atom)
    {
      const double* XYZ = currentFrame->XYZ(*atom);
      //initial_.push_back( XYZ[0] );
      previousx_.push_back( XYZ[0] );
      //initial_.push_back( XYZ[1] );
      previousy_.push_back( XYZ[1] );
      //initial_.push_back( XYZ[2] );
      previousz_.push_back( XYZ[2] );
    }
  } else {
    if (hasBox_) 
      boxcenter_ = currentFrame->BoxCrd().Center();
    Vec3 boxL = currentFrame->BoxCrd().Lengths();
    // Set iterators
    std::vector<double>::iterator px = previousx_.begin();
    std::vector<double>::iterator py = previousy_.begin();
    std::vector<double>::iterator pz = previousz_.begin();
    std::vector<double>::iterator dx = deltax_.begin();
    std::vector<double>::iterator dy = deltay_.begin();
    std::vector<double>::iterator dz = deltaz_.begin();
    std::vector<double>::iterator distx = distancex_.begin();
    std::vector<double>::iterator disty = distancey_.begin();
    std::vector<double>::iterator distz = distancez_.begin();
    std::vector<double>::iterator dist  = distance_.begin();
    // For averaging over selected atoms
    double average = 0;
    double avgx = 0;
    double avgy = 0;
    double avgz = 0;
    for (AtomMask::const_iterator atom = mask_.begin(); atom != mask_.end(); ++atom)
    { // Get current and initial coords for this atom.
      const double* XYZ = currentFrame->XYZ(*atom);
      const double* iXYZ = initial_.XYZ(*atom);
      // Calculate distance to previous frames coordinates.
      double delx = XYZ[0] - *px;
      double dely = XYZ[1] - *py;
      double delz = XYZ[2] - *pz;
      // If the particle moved more than half the box, assume
      // it was imaged and adjust the distance of the total
      // movement with respect to the original frame.
      if (hasBox_) {
        if      (delx >  boxcenter_[0]) *dx -= boxL[0];
        else if (delx < -boxcenter_[0]) *dx += boxL[0];
        else if (dely >  boxcenter_[1]) *dy -= boxL[1];
        else if (dely < -boxcenter_[1]) *dy += boxL[1];
        else if (delz >  boxcenter_[2]) *dz -= boxL[2];
        else if (delz < -boxcenter_[2]) *dz += boxL[2];
      }
      // DEBUG
      if (debug_ > 2)
        mprintf("ATOM: %5i %10.3f %10.3f %10.3f",*atom,XYZ[0],delx,*dx);
      // Set the current x with reference to the un-imaged trajectory.
      double xx = XYZ[0] + *dx; 
      double yy = XYZ[1] + *dy; 
      double zz = XYZ[2] + *dz;
      // Calculate the distance between this "fixed" coordinate
      // and the reference (initial) frame.
      delx = xx - iXYZ[0];
      dely = yy - iXYZ[1];
      delz = zz - iXYZ[2];
      // DEBUG
      if (debug_ > 2)
        mprintf(" %10.3f\n", delx);
      // Store distance for this atom
      *distx = delx*delx;
      *disty = dely*dely;
      *distz = delz*delz;
      *dist = *distx + *disty + *distz;
      // Accumulate averages
      avgx += *distx;
      avgy += *disty;
      avgz += *distz;
      average += *dist;
      // Update the previous coordinate set to match the current coordinates
      *px = XYZ[0];
      *py = XYZ[1];
      *pz = XYZ[2];
      // Increment iterators
      ++px;
      ++py;
      ++pz;
      ++dx;
      ++dy;
      ++dz;
      ++distx;
      ++disty;
      ++distz;
      ++dist;
    } // END loop over selected atoms
    double dNselected = (double)mask_.Nselected();
    average /= dNselected;
    avgx /= dNselected;
    avgy /= dNselected;
    avgz /= dNselected;

    // ----- OUTPUT -----
    // Output averages
    double Time = time_ * (double)frameNum;
    outputx_.Printf("%8.3f  %8.3f", Time, avgx);
    outputy_.Printf("%8.3f  %8.3f", Time, avgy);
    outputz_.Printf("%8.3f  %8.3f", Time, avgz);
    outputr_.Printf("%8.3f  %8.3f", Time, average);
    outputa_.Printf("%8.3f  %8.3f", Time, sqrt(average));
    // Individual values
    if (printIndividual_) {
      for (int i = 0; i < mask_.Nselected(); ++i) {
        outputx_.Printf("  %8.3f", distancex_[i]);
        outputy_.Printf("  %8.3f", distancey_[i]);
        outputz_.Printf("  %8.3f", distancez_[i]);
        outputr_.Printf("  %8.3f", distance_[i]);
        outputa_.Printf("  %8.3f", sqrt(distance_[i]));
      }
    }
    // Print newlines
    outputx_.Printf("\n");
    outputy_.Printf("\n");
    outputz_.Printf("\n");
    outputr_.Printf("\n");
    outputa_.Printf("\n");
  }

  return Action::OK;
}  
