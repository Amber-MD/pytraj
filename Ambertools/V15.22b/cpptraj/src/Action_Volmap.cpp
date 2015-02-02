#include <cmath> // pow, exp, sqrt
#include <algorithm> // std::min, std::max
#include "Action_Volmap.h"
#include "Constants.h" // PI
#include "CpptrajStdio.h"

const double Action_Volmap::sqrt_8_pi_cubed = sqrt(8.0*Constants::PI*Constants::PI*Constants::PI);
const double Action_Volmap::one_over_6 = 1.0 / 6.0;
// CONSTRUCTOR
Action_Volmap::Action_Volmap() :
  ensembleNum_(-1),
  dx_(0.0), dy_(0.0), dz_(0.0),
  xmin_(0.0), ymin_(0.0), zmin_(0.0),
  Nframes_(0),
  setupGridOnMask_(false),
  grid_(0),
  peakcut_(0.05),
  buffer_(3.0),
  radscale_(1.0)
{}

void Action_Volmap::Help() {
  RawHelp();
  mprintf("    filename  : Output file name\n"
          "    dx, dy, dz: grid spacing in the x-, y-, and z-dimensions, respectively.\n"
          "  The grid size can be determined either by the size (x,y,z in Angstroms)\n"
          "  or by a rectangular prism enclosing a mask with <buffer> clearance\n"
          "  in every dimension. The density is calculated from the atoms in the\n"
          "  required <mask>. If a <buffer> is given, the grid is centered on the\n"
          "  centermask if provided, or the required mask if not.\n");
}

void Action_Volmap::RawHelp() {
  mprintf("\tfilename dx dy dz <mask> [xplor] [radscale <factor>]\n"
          "\t[ [[buffer <buffer>] [centermask <mask>]] | [center <x,y,z>] [size <x,y,z>] ]\n"
          "\t[peakcut <cutoff>] [peakfile <xyzfile>]\n");
}

// Action_Volmap::Init()
Action::RetType Action_Volmap::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  ensembleNum_ = DSL->EnsembleNum();
  // Get the required mask
  std::string reqmask = actionArgs.GetMaskNext();
  if (reqmask.empty()) {
     mprinterr("Error: Volmap: no density mask specified.\n");
     return Action::ERR;
  }
  densitymask_.SetMaskString(reqmask);
  // Get output filename
  std::string filename = actionArgs.GetStringNext();
  if (filename.empty()) {
    mprinterr("Error: Volmap: no filename specified.\n");
    return Action::ERR;
  }
  // Get grid resolutions
  dx_ = actionArgs.getNextDouble(0.0);
  dy_ = actionArgs.getNextDouble(0.0);
  dz_ = actionArgs.getNextDouble(0.0);
  // Get extra options
  peakcut_ = actionArgs.getKeyDouble("peakcut", 0.05);
  peakfilename_ = actionArgs.GetStringKey("peakfile");
  radscale_ = 1.0 / actionArgs.getKeyDouble("radscale", 1.0);
  std::string sizestr = actionArgs.GetStringKey("size");
  std::string center = actionArgs.GetStringKey("centermask");
  std::string density = actionArgs.GetStringKey("density");

  // See how we are going to be setting up our grid
  setupGridOnMask_ = false;
  std::string dsname = actionArgs.GetStringKey("data");
  if (!dsname.empty()) {
    // Get existing grid dataset
    grid_ = (DataSet_GridFlt*)DSL->FindSetOfType( dsname, DataSet::GRID_FLT );
    if (grid_ == 0) {
      mprinterr("Error: volmap: Could not find grid data set with name %s\n",
                dsname.c_str());
      return Action::ERR;
    }
  } else {
    // Create new grid.
    grid_ = (DataSet_GridFlt*)DSL->AddSet(DataSet::GRID_FLT, actionArgs.GetStringKey("name"),
                                         "VOLMAP");
    if (grid_ == 0) return Action::ERR;
    if (!sizestr.empty()) {
      // Get grid sizes from the specified arguments
      ArgList sizeargs = ArgList(sizestr, ",");
      double xsize = sizeargs.getNextDouble(0.0);
      double ysize = sizeargs.getNextDouble(0.0);
      double zsize = sizeargs.getNextDouble(0.0);
      if (xsize <= 0 || ysize <= 0 || zsize <= 0) {
        mprinterr("Error: Volmap: Illegal grid sizes [%s]\n", sizestr.c_str());
        return Action::ERR;
      }
      std::string centerstr = actionArgs.GetStringKey("center");
      ArgList centerargs = ArgList(centerstr, ",");
      double xcenter = centerargs.getNextDouble(0.0);
      double ycenter = centerargs.getNextDouble(0.0);
      double zcenter = centerargs.getNextDouble(0.0);
      // Allocate grid
      if ( grid_->Allocate_X_C_D(Vec3(xsize,ysize,zsize), 
                                 Vec3(xcenter,ycenter,zcenter), 
                                 Vec3(dx_,dy_,dz_)) ) return Action::ERR;
      Vec3 const& oxyz = grid_->GridOrigin();
      xmin_ = oxyz[0];
      ymin_ = oxyz[1];
      zmin_ = oxyz[2];
    } else {
      // Will generate grid around a mask. See if we have a center mask
      if (center.empty())
        centermask_.SetMaskString(reqmask);
      else
        centermask_.SetMaskString(center);
      buffer_ = actionArgs.getKeyDouble("buffer", 3.0);
      if (buffer_ < 0) {
        mprintf("Error: Volmap: The buffer must be non-negative.\n");
        return Action::ERR;
      }
      setupGridOnMask_ = true;
    }
  }

  // Setup output file
  DataFile* outfile = DFL->AddSetToFile(filename, (DataSet*)grid_);
  if (outfile == 0) {
    mprinterr("Error: volmap: Could not set up output file %s\n", filename.c_str());
    return Action::ERR;
  }

  // Info
  mprintf("    VOLMAP: Grid spacing will be %.2fx%.2fx%.2f Angstroms\n", dx_, dy_, dz_);
  if (sizestr.empty())
    mprintf("\tGrid centered around %s with %.2f Ang. clearance\n",
            centermask_.MaskExpression().c_str(), buffer_);
  else
    mprintf("\tGrid centered at origin.\n");
  mprintf("\tDensity will wrtten to %s\n", filename.c_str());
  if (!peakfilename_.empty())
    mprintf("\tDensity peaks above %.3f will be printed to %s in XYZ-format\n",
            peakcut_, peakfilename_.c_str());

  return Action::OK;
}

// Action_Volmap::Setup()
Action::RetType Action_Volmap::Setup(Topology* currentParm, Topology** parmAddress) {
  // Set up our density mask and make sure it's not empty
  if (currentParm->SetupIntegerMask( densitymask_ ))
    return Action::ERR;
  if (densitymask_.None()) {
    mprinterr("Error: Volmap: Density mask selection empty!\n");
    return Action::ERR;
  }
  mprintf("\tVolmap: Grid mask [%s] selects %d atoms.\n", densitymask_.MaskString(),
          densitymask_.Nselected());
  // If we did not specify a size, make sure we have a valid centermask_
  if (setupGridOnMask_) {
    if (currentParm->SetupIntegerMask( centermask_ ))
      return Action::ERR;
    // The masks must be populated
    if (centermask_.None()) {
      mprinterr("Error: Volmap: mask selection(s) empty!\n");
      return Action::ERR;
    }
    mprintf("\tVolmap: Centered mask [%s] selects %d atoms.\n", centermask_.MaskString(),
            centermask_.Nselected());
  }
  // Set up our radii_
  halfradii_.clear();
  halfradii_.reserve( currentParm->Natom() );
  for (int i = 0; i < currentParm->Natom(); i++)
    halfradii_.push_back( (float)(GetRadius_(*currentParm, i) * radscale_ / 2) );

  // DEBUG
//for (AtomMask::const_iterator it = densitymask_.begin(); it != densitymask_.end(); it++)
//  mprintf("Radius of atom %d is %f\n", *it, 2 * halfradii_[*it]);
  
  return Action::OK;
}

// Action_Volmap::GetRadius_()
/** Takes an input topology and gives back the VDW radius
  * \return vdW radius of the requested atom number from a Topology instance
  */
double Action_Volmap::GetRadius_(Topology const& top, int atom) {
  NonbondType const& LJ = top.GetLJparam(atom, atom);
  return 0.5 * pow(2 * LJ.A() / LJ.B(), one_over_6);
}

// Action_Volmap::DoAction()
Action::RetType Action_Volmap::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {
  // If this is our first frame, then we may have to set up our grid from our masks
  if (Nframes_ == 0) {
    if (setupGridOnMask_) {
      // Determine min/max coord values for atoms in centermask. Calculate
      // geometric center while doing this.
      double xmin, xmax, ymin, ymax, zmin, zmax;
      AtomMask::const_iterator it = centermask_.begin();
      Vec3 cxyz = Vec3(currentFrame->XYZ(*it));
      xmin = xmax = cxyz[0];
      ymin = ymax = cxyz[1];
      zmin = zmax = cxyz[2];
      ++it;
      for (; it != centermask_.end(); it++) {
        Vec3 pt = Vec3(currentFrame->XYZ(*it));
        cxyz += pt;
        xmin = std::min(xmin, pt[0]);
        xmax = std::max(xmax, pt[0]);
        ymin = std::min(ymin, pt[1]);
        ymax = std::max(ymax, pt[1]);
        zmin = std::min(zmin, pt[2]);
        zmax = std::max(zmax, pt[2]);
      }
      cxyz /= (double)centermask_.Nselected();
      // Extend min/max by buffer.
      xmin -= buffer_; 
      xmax += buffer_;
      ymin -= buffer_; 
      ymax += buffer_;
      zmin -= buffer_; 
      zmax += buffer_;
      // Allocate grid of given size centered on mask.
      if (grid_->Allocate_N_O_D( (xmax-xmin)/dx_, (ymax-ymin)/dy_, (zmax-zmin)/dz_,
                                 Vec3(xmin, ymin, zmin), Vec3(dx_, dy_, dz_) ))
        return Action::ERR;
      xmin_ = xmin;
      ymin_ = ymin;
      zmin_ = zmin;
      setupGridOnMask_ = false;
    }
  }
  // Now calculate the density for every point
  // TODO: Convert everything to size_t?
  int nX = (int)grid_->NX();
  int nY = (int)grid_->NY();
  int nZ = (int)grid_->NZ();
  for (AtomMask::const_iterator atom = densitymask_.begin();
                                atom != densitymask_.end(); ++atom) {
    Vec3 pt = Vec3(currentFrame->XYZ(*atom));
    int ix = (int) ( floor( (pt[0]-xmin_) / dx_ + 0.5 ) );
    int iy = (int) ( floor( (pt[1]-ymin_) / dy_ + 0.5 ) );
    int iz = (int) ( floor( (pt[2]-zmin_) / dz_ + 0.5 ) );
    /* See how many steps in each dimension we smear out our Gaussian. This
     * formula is taken to be consistent with VMD's volmap tool
     */
    int nxstep = (int) ceil(4.1 * halfradii_[*atom] / dx_);
    int nystep = (int) ceil(4.1 * halfradii_[*atom] / dy_);
    int nzstep = (int) ceil(4.1 * halfradii_[*atom] / dz_);
    // Calculate the gaussian normalization factor (in 3 dimensions with the
    // given half-radius)
    double norm = 1 / (sqrt_8_pi_cubed * 
                       halfradii_[*atom]*halfradii_[*atom]*halfradii_[*atom]);
    double exfac = -1.0 / (2.0 * halfradii_[*atom] * halfradii_[*atom]);
    if (ix < -nxstep || ix > nX + nxstep ||
        iy < -nystep || iy > nY + nystep ||
        iz < -nzstep || iz > nZ + nzstep)
      continue;

    int xend = std::min(ix+nxstep, nX);
    int yend = std::min(iy+nystep, nY);
    int zend = std::min(iz+nzstep, nZ);
    for (int xval = std::max(ix-nxstep, 0); xval < xend; xval++)
      for (int yval = std::max(iy-nystep, 0); yval < yend; yval++)
        for (int zval = std::max(iz-nzstep, 0); zval < zend; zval++) {
          Vec3 gridpt = Vec3(xmin_+xval*dx_, ymin_+yval*dy_, zmin_+zval*dz_) - pt;
          double dist2 = gridpt.Magnitude2();
          grid_->Increment(xval, yval, zval, norm * exp(exfac * dist2));
        }
  } // END loop over atoms in densitymask_
  
  // Increment frame counter
  Nframes_++;
  return Action::OK;
}

// Need this instead of MAX since size_t can never be negative
inline size_t setStart(size_t xIn) {
  if (xIn == 0)
    return 0UL;
  else
    return xIn - 1L;
}

// Action_Volmap::Print()
void Action_Volmap::Print() {

  // Divide our grid by the number of frames
  float nf = (float)Nframes_;
  for (DataSet_GridFlt::iterator gval = grid_->begin(); gval != grid_->end(); ++gval)
    *gval /= nf;
//    grid_.PrintXplor( filename_, "This line is ignored", 
//                      "rdparm generated grid density" );
  
  // See if we need to write the peaks out somewhere
  if (!peakfilename_.empty()) {
    // Extract peaks from the current grid, setup another Grid instance. This
    // works by taking every grid point and analyzing all grid points adjacent
    // to it (including diagonals). If any of those grid points have a higher 
    // value (meaning there is a direction towards "increased" density) then 
    // that value is _not_ a maximum. Any density peaks less than the minimum
    // filter are discarded. The result is a Grid instance with all non-peak 
    // grid points zeroed-out.
    Grid<float> peakgrid = grid_->InternalGrid();
    for (size_t i = 0; i < grid_->NX(); i++)
      for (size_t j = 0; j < grid_->NY(); j++)
        for (size_t k = 0; k < grid_->NZ(); k++) {
          float val = grid_->GridVal(i, j, k);
          if (val < peakcut_) {
            peakgrid.setGrid(i, j, k, 0.0f);
            continue;
          }
          size_t i_end = std::min(i+2, grid_->NX());
          size_t j_end = std::min(j+2, grid_->NY());
          size_t k_end = std::min(k+2, grid_->NZ()); 
          for (size_t ii = setStart(i); ii < i_end; ii++)
            for (size_t jj = setStart(j); jj < j_end; jj++)
              for (size_t kk = setStart(k); kk < k_end; kk++) {
                if (ii==i && jj==j && kk==k) continue;
                if (grid_->GridVal(ii, jj, kk) > val)
                  peakgrid.setGrid(i,j,k,0.0f); // TODO: break after this?
              }
        }
    int npeaks = 0;
    std::vector<double> peakdata;
    for (size_t i = 0; i < peakgrid.NX(); i++)
      for (size_t j = 0; j < peakgrid.NY(); j++)
        for (size_t k = 0; k < peakgrid.NZ(); k++) {
          double gval = peakgrid.element(i, j, k);
          if (gval > 0) {
            npeaks++;
            peakdata.push_back(xmin_+dx_*i);
            peakdata.push_back(ymin_+dy_*j);
            peakdata.push_back(zmin_+dz_*k);
            peakdata.push_back(gval);
          }
        }
    // If we have peaks, open up our peak data and print it
    if (npeaks > 0) {
      CpptrajFile outfile;
      if(outfile.OpenEnsembleWrite(peakfilename_, ensembleNum_)) {
        mprinterr("Error: Could not open %s for writing.\n", peakfilename_.c_str());
        return;
      }
      outfile.Printf("%d\n\n", npeaks);
      for (int i = 0; i < npeaks; i++)
        outfile.Printf("C %16.8f %16.8f %16.8f %16.8f\n", peakdata[4*i],
                       peakdata[4*i+1], peakdata[4*i+2], peakdata[4*i+3]);
      outfile.CloseFile();
      mprintf("Volmap: %d density peaks found with higher density than %.4lf\n",
              npeaks, peakcut_);
    }else{
      mprintf("No peaks found with a density greater than %.3lf\n", peakcut_);
    }
  }
}
