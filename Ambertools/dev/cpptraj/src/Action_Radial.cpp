// Action_Radial
#include <cmath> // sqrt
#include "Action_Radial.h"
#include "CpptrajStdio.h"
#include "Constants.h" // FOURTHIRDSPI
#ifdef _OPENMP
#  include "omp.h"
#endif

// CONSTRUCTOR
Action_Radial::Action_Radial() :
  RDF_(0),
  rdf_thread_(0),
  rmode_(NORMAL),
  currentParm_(0),
  intramol_distances_(0),
  useVolume_(false),
  volume_(0),
  maximum2_(0),
  spacing_(-1),
  one_over_spacing_(-1),
  numBins_(0),
  numthreads_(1),
  numFrames_(0),
  // Default particle density (molecules/Ang^3) for water based on 1.0 g/mL
  density_(0.033456),
  Dset_(0),
  intrdf_(0),
  rawrdf_(0),
  debug_(0)
{} 

void Action_Radial::Help() {
  mprintf("\t<outfilename> <spacing> <maximum> <solvent mask1> [<solute mask2>] [noimage]\n"
          "\t[density <density> | volume] [center1 | center2 | nointramol] [<name>]\n"
          "\t[intrdf <file>] [rawrdf <file>]\n"
          "  Calculate the radial distribution function (RDF) of atoms in <solvent mask1>.\n"
          "  If <solute mask2> is given calculate RDF of all atoms in <solvent mask1>\n"
          "  to each atom in <solute mask2>.\n");
}

// DESTRUCTOR
Action_Radial::~Action_Radial() {
  //fprintf(stderr,"Radial Destructor.\n");
  if (RDF_!=0) delete[] RDF_;
  if (rdf_thread_!=0) {
    for (int i=0; i < numthreads_; i++)
      delete[] rdf_thread_[i];
    delete[] rdf_thread_;
  }
}

inline Action::RetType RDF_ERR(const char* msg) {
  mprinterr("Error: %s\n", msg);
  return Action::ERR;
}

// Action_Radial::Init()
Action::RetType Action_Radial::Init(ArgList& actionArgs, TopologyList* PFL, DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  debug_ = debugIn;
  // Get Keywords
  InitImaging( !(actionArgs.hasKey("noimage")) );
  // Default particle density (mols/Ang^3) for water based on 1.0 g/mL
  density_ = actionArgs.getKeyDouble("density",0.033456);
  if (actionArgs.hasKey("center1"))
    rmode_ = CENTER1;
  else if (actionArgs.hasKey("center2"))
    rmode_ = CENTER2;
  else if (actionArgs.hasKey("nointramol"))
    rmode_ = NO_INTRAMOL;
  else
    rmode_ = NORMAL;
  useVolume_ = actionArgs.hasKey("volume");
  std::string intrdfname = actionArgs.GetStringKey("intrdf");
  std::string rawrdfname = actionArgs.GetStringKey("rawrdf");

  // Get required args
  std::string outfilename = actionArgs.GetStringNext();
  if (outfilename.empty()) {
    mprinterr("Error: Radial: No output filename given.\n");
    Help();
    return Action::ERR;
  }
  spacing_ = actionArgs.getNextDouble(-1.0);
  if (spacing_ < 0) {
    mprinterr("Error: Radial: No spacing argument or arg < 0.\n");
    Help();
    return Action::ERR;
  }
  double maximum = actionArgs.getNextDouble(-1.0);
  if (maximum < 0) {
    mprinterr("Error: Radial: No maximum argument or arg < 0.\n");
    Help();
    return Action::ERR;
  }
  // Store max^2, distances^2 greater than max^2 do not need to be
  // binned and therefore do not need a sqrt calc.
  maximum2_ = maximum * maximum;

  // Get First Mask
  std::string mask1 = actionArgs.GetMaskNext();
  if (mask1.empty()) {
    mprinterr("Error: Radial: No mask given.\n");
    return Action::ERR;
  }
  Mask1_.SetMaskString(mask1);

  // Check for second mask - if none specified use first mask
  std::string mask2 = actionArgs.GetMaskNext();
  if (!mask2.empty()) 
    Mask2_.SetMaskString(mask2);
  else
    Mask2_.SetMaskString(mask1);

  // Set up output dataset. 
  Dset_ = DSL->AddSet( DataSet::DOUBLE, actionArgs.GetStringNext(), "g(r)");
  if (Dset_ == 0) return RDF_ERR("Could not allocate RDF data set.");
  DataFile* outfile = DFL->AddSetToFile(outfilename, Dset_);
  if (outfile==0) {
    mprinterr("Error: Radial: Could not setup output file %s\n",outfilename.c_str());
    return Action::ERR;
  }
  // Make default precision a little higher than normal
  Dset_->SetPrecision(12,6);
  // Set DataSet legend from mask strings.
  Dset_->SetLegend(Mask1_.MaskExpression() + " => " + Mask2_.MaskExpression());
  // TODO: Set Yaxis label in DataFile
  // Calculate number of bins
  one_over_spacing_ = 1 / spacing_;
  double temp_numbins = maximum * one_over_spacing_;
  temp_numbins = ceil(temp_numbins);
  numBins_ = (int) temp_numbins;
  // Setup output datafile. Align on bin centers instead of left.
  // TODO: Use Rdim for binning?
  Dimension Rdim( spacing_ / 2.0, spacing_, numBins_, "Distance (Ang)" ); 
  Dset_->SetDim(Dimension::X, Rdim);
  // Set up output for integral of mask2 if specified.
  if (!intrdfname.empty()) {
    intrdf_ = DSL->AddSetAspect( DataSet::DOUBLE, Dset_->Name(), "int" );
    if (intrdf_ == 0) return RDF_ERR("Could not allocate RDF integral data set.");
    intrdf_->SetPrecision(12,6);
    intrdf_->SetLegend("Int[" + Mask2_.MaskExpression() + "]");
    intrdf_->SetDim(Dimension::X, Rdim);
    outfile = DFL->AddSetToFile( intrdfname, intrdf_ );
    if (outfile == 0) {
      mprinterr("Error: Could not add intrdf set to file %s\n", intrdfname.c_str());
      return Action::ERR;
    }
  } else
    intrdf_ = 0;
  // Set up output for raw rdf
  if (!rawrdfname.empty()) {
    rawrdf_ = DSL->AddSetAspect( DataSet::DOUBLE, Dset_->Name(), "raw" );
    if (rawrdf_ == 0) return RDF_ERR("Could not allocate raw RDF data set.");
    rawrdf_->SetPrecision(12,6);
    rawrdf_->SetLegend("Raw[" + Dset_->Legend() + "]");
    rawrdf_->SetDim(Dimension::X, Rdim);
    outfile = DFL->AddSetToFile( rawrdfname, rawrdf_ );
    if (outfile == 0) {
      mprinterr("Error: Could not add rawrdf set to file %s\n", rawrdfname.c_str());
      return Action::ERR;
    }
  } else
    rawrdf_ = 0;

  // Set up histogram
  RDF_ = new int[ numBins_ ];
  std::fill(RDF_, RDF_ + numBins_, 0);
# ifdef _OPENMP
  // Since RDF is shared by all threads and we cant guarantee that a given
  // bin in RDF wont be accessed at the same time by the same thread,
  // each thread needs its own bin space.
#pragma omp parallel
{
  if (omp_get_thread_num()==0)
    numthreads_ = omp_get_num_threads();
}
  rdf_thread_ = new int*[ numthreads_ ];
  for (int i=0; i < numthreads_; i++) {
    rdf_thread_[i] = new int[ numBins_ ];
    std::fill(rdf_thread_[i], rdf_thread_[i] + numBins_, 0);
  }
# endif
  
  mprintf("    RADIAL: Calculating RDF for atoms in mask [%s]",Mask1_.MaskString());
  if (!mask2.empty()) 
    mprintf(" to atoms in mask [%s]",Mask2_.MaskString());
  mprintf("\n            Output to %s.\n",outfilename.c_str());
  if (intrdf_ != 0)
    mprintf("            Integral of mask2 atoms will be output to %s\n", intrdfname.c_str());
  if (rawrdf_ != 0)
    mprintf("            Raw RDF bin values will be output to %s\n", rawrdfname.c_str());
  if (rmode_==CENTER1)
    mprintf("            Using center of atoms in mask1.\n");
  else if (rmode_==CENTER2)
    mprintf("            Using center of atoms in mask2.\n");
  mprintf("            Histogram max %f, spacing %f, bins %i.\n",maximum,
          spacing_,numBins_);
  if (useVolume_)
    mprintf("            Normalizing based on cell volume.\n");
  else
    mprintf("            Normalizing using particle density of %f molecules/Ang^3.\n",density_);
  if (!UseImage()) 
    mprintf("            Imaging disabled.\n");
  if (numthreads_ > 1)
    mprintf("            Parallelizing RDF calculation with %i threads.\n",numthreads_);

  return Action::OK;
}

// Action_Radial::Setup()
/** Determine what atoms each mask pertains to for the current parm file.
  * Also determine whether imaging should be performed.
  */
Action::RetType Action_Radial::Setup(Topology* currentParm, Topology** parmAddress) {

  if ( currentParm->SetupIntegerMask( Mask1_ ) ) return Action::ERR;
  if (Mask1_.None()) {
    mprintf("    Error: Radial: First mask has no atoms.\n");
    return Action::ERR;
  }
  if (currentParm->SetupIntegerMask( Mask2_ ) ) return Action::ERR;
  if (Mask2_.None()) {
    mprintf("    Error: Radial: Second mask has no atoms.\n");
    return Action::ERR;
  }
  SetupImaging( currentParm->BoxType() );

  // If not computing center for mask 1 or 2, make the outer loop for distance
  // calculation correspond to the mask with the most atoms.
  if (rmode_ == NORMAL || rmode_ == NO_INTRAMOL) {
    if (Mask1_.Nselected() > Mask2_.Nselected()) {
      OuterMask_ = Mask1_;
      InnerMask_ = Mask2_;
    } else {
      OuterMask_ = Mask2_;
      InnerMask_ = Mask1_;
    }
  } else if (rmode_ == CENTER1) {
    OuterMask_ = Mask1_;
    InnerMask_ = Mask2_;
  } else if (rmode_ == CENTER2) {
    OuterMask_ = Mask2_;
    InnerMask_ = Mask1_;
  }

  // If ignoring intra-molecular distances, need to count how many we
  // are ignoring.
  if (rmode_ == NO_INTRAMOL) {
    int ndist = 0;
    for (AtomMask::const_iterator atom1 = OuterMask_.begin(); 
                                  atom1 != OuterMask_.end(); ++atom1)
      for (AtomMask::const_iterator atom2 = InnerMask_.begin();
                                    atom2 != InnerMask_.end(); ++atom2)
        if ( (*currentParm)[*atom1].MolNum() == (*currentParm)[*atom2].MolNum() )
          ++ndist;
    if (currentParm_ != 0 && ndist != intramol_distances_)
      mprintf("Warning: # of intramolecular distances (%i) has changed from the last"
              " topology (%i).\nWarning: Normalization will not be correct.\n",
              ndist, intramol_distances_);
    intramol_distances_ = ndist;
    currentParm_ = currentParm;
    mprintf("\tIgnoring %i intra-molecular distances.\n", intramol_distances_);
  }

  // Check volume information
  if (useVolume_ && currentParm->BoxType()==Box::NOBOX) {
    mprintf("    Warning: Radial: 'volume' specified but no box information for %s, skipping.\n",
            currentParm->c_str());
    return Action::ERR;
  }

  // Print mask and imaging info for this parm
  mprintf("    RADIAL: %i atoms in Mask1, %i atoms in Mask2, ",
          Mask1_.Nselected(), Mask2_.Nselected());
  if (ImagingEnabled())
    mprintf("Imaging on.\n");
  else
    mprintf("Imaging off.\n");
  return Action::OK;  
}

// Action_Radial::DoAction()
/** Calculate distances from atoms in mask1 to atoms in mask 2 and
  * bin them.
  */
// NOTE: Because of maximum2 not essential to check idx>numBins?
Action::RetType Action_Radial::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {
  double D;
  Matrix_3x3 ucell, recip;
  int atom1, atom2;
  int nmask1, nmask2;
  int idx;
# ifdef _OPENMP
  int mythread;
# endif

  // Set imaging information and store volume if specified
  // NOTE: Ucell and recip only needed for non-orthogonal boxes.
  if (ImagingEnabled() || useVolume_) {
    D = currentFrame->BoxCrd().ToRecip(ucell,recip);
    if (useVolume_)  volume_ += D;
  }

  if ( rmode_ == NORMAL ) { 
    // Calculation of all atoms in Mask1 to all atoms in Mask2
    int outer_max = OuterMask_.Nselected();
    int inner_max = InnerMask_.Nselected();
#   ifdef _OPENMP
#   pragma omp parallel private(nmask1,nmask2,atom1,atom2,D,idx,mythread)
    {
    //mprintf("OPENMP: %i threads\n",omp_get_num_threads());
    mythread = omp_get_thread_num();
#   pragma omp for
#   endif
    for (nmask1 = 0; nmask1 < outer_max; nmask1++) {
      atom1 = OuterMask_[nmask1];
      for (nmask2 = 0; nmask2 < inner_max; nmask2++) {
        atom2 = InnerMask_[nmask2];
        if (atom1 != atom2) {
          D = DIST2( currentFrame->XYZ(atom1), currentFrame->XYZ(atom2),
                     ImageType(), currentFrame->BoxCrd(), ucell, recip);
          if (D <= maximum2_) {
            // NOTE: Can we modify the histogram to store D^2?
            D = sqrt(D);
            //mprintf("MASKLOOP: %10i %10i %10.4f\n",atom1,atom2,D);
            idx = (int) (D * one_over_spacing_);
            if (idx > -1 && idx < numBins_)
#             ifdef _OPENMP
              ++rdf_thread_[mythread][idx];
#             else
              ++RDF_[idx];
#             endif
          }
        }
      } // END loop over 2nd mask
    } // END loop over 1st mask
#   ifdef _OPENMP
    } // END pragma omp parallel
#   endif 
  } else if ( rmode_ == NO_INTRAMOL ) {
    // Calculation of all atoms in Mask1 to all atoms in Mask2, ignoring
    // intra-molecular distances.
    int outer_max = OuterMask_.Nselected();
    int inner_max = InnerMask_.Nselected();
#   ifdef _OPENMP
#   pragma omp parallel private(nmask1,nmask2,atom1,atom2,D,idx,mythread)
    {
    //mprintf("OPENMP: %i threads\n",omp_get_num_threads());
    mythread = omp_get_thread_num();
#   pragma omp for
#   endif
    for (nmask1 = 0; nmask1 < outer_max; nmask1++) {
      atom1 = OuterMask_[nmask1];
      for (nmask2 = 0; nmask2 < inner_max; nmask2++) {
        atom2 = InnerMask_[nmask2];
        if ( (*currentParm_)[atom1].MolNum() != (*currentParm_)[atom2].MolNum() ) {
          D = DIST2( currentFrame->XYZ(atom1), currentFrame->XYZ(atom2),
                     ImageType(), currentFrame->BoxCrd(), ucell, recip);
          if (D <= maximum2_) {
            // NOTE: Can we modify the histogram to store D^2?
            D = sqrt(D);
            //mprintf("MASKLOOP: %10i %10i %10.4f\n",atom1,atom2,D);
            idx = (int) (D * one_over_spacing_);
            if (idx > -1 && idx < numBins_)
#             ifdef _OPENMP
              ++rdf_thread_[mythread][idx];
#             else
              ++RDF_[idx];
#             endif
          }
        }
      } // END loop over 2nd mask
    } // END loop over 1st mask
#   ifdef _OPENMP
    } // END pragma omp parallel
#   endif
  } else { // CENTER1 || CENTER2
    // Calculation of center of one Mask to all atoms in other Mask
    Vec3 coord_center = currentFrame->VGeometricCenter(OuterMask_);
    int mask2_max = InnerMask_.Nselected();
#   ifdef _OPENMP
#   pragma omp parallel private(nmask2,atom2,D,idx,mythread)
    {
    mythread = omp_get_thread_num();
#   pragma omp for
#   endif
    for (nmask2 = 0; nmask2 < mask2_max; nmask2++) {
      atom2 = InnerMask_[nmask2];
      D = DIST2(coord_center.Dptr(), currentFrame->XYZ(atom2), ImageType(),
                currentFrame->BoxCrd(), ucell, recip);
      if (D <= maximum2_) {
        // NOTE: Can we modify the histogram to store D^2?
        D = sqrt(D);
        //mprintf("MASKLOOP: %10i %10i %10.4f\n",atom1,atom2,D);
        idx = (int) (D * one_over_spacing_);
        if (idx > -1 && idx < numBins_)
#         ifdef _OPENMP
          ++rdf_thread_[mythread][idx];
#         else
          ++RDF_[idx];
#         endif
      }
    } // END loop over 2nd mask
#   ifdef _OPENMP
    } // END pragma omp parallel
#   endif 
  }
  ++numFrames_;

  return Action::OK;
} 

// Action_Radial::Print()
/** Convert the histogram to a dataset, normalize, create datafile.
  */
// NOTE: Currently the normalization is based on number of atoms in each mask;
//       if multiple prmtops are loaded and number of atoms changes from 
//       prmtop to prmtop this will throw off normalization.
void Action_Radial::Print() {
  if (numFrames_==0) return;
# ifdef _OPENMP 
  // Combine results from each rdf_thread into rdf
  for (int thread=0; thread < numthreads_; thread++) 
    for (int bin = 0; bin < numBins_; bin++) 
      RDF_[bin] += rdf_thread_[thread][bin];
# endif

  mprintf("    RADIAL: %i frames,", numFrames_);
  double nmask1 = (double)Mask1_.Nselected();
  double nmask2 = (double)Mask2_.Nselected();
  int numSameAtoms = 0;
  if (rmode_ == NORMAL) {
    // If Mask1 and Mask2 have any atoms in common distances were not calcd
    // between them (because they are 0.0 of course); need to correct for this.
    numSameAtoms = Mask1_.NumAtomsInCommon( Mask2_ );
  } else if (rmode_ == NO_INTRAMOL) {
    // # of intra-molecular distances already counted. Any atoms in common
    // between mask1 and mask2 will have been included as an intramol dist.
    numSameAtoms = intramol_distances_ ;
  } else if (rmode_ == CENTER1) {
    // If the center1 option was specified only one distance was calcd
    // from mask 1. Assume COM of mask 1 != atom(s) in mask2.
    nmask1 = 1.0;
    numSameAtoms = 0;
  } else if (rmode_ == CENTER2) {
    // If the center2 option was specified only one distance was calcd
    // from mask 2. Assume COM of mask 2 != atom(s) in mask1.
    nmask2 = 1.0;
    numSameAtoms = 0;
  }
  mprintf(" # in mask1= %.0f, # in mask2 = %.0f, # in common = %i\n",
          nmask1, nmask2, numSameAtoms);
  
  // If useVolume, calculate the density from the average volume
  if (useVolume_) {
    double avgVol = volume_ / numFrames_;
    mprintf("            Average volume is %f Ang^3.\n",avgVol);
    density_ = (nmask1 * nmask2 - (double)numSameAtoms) / avgVol;
    mprintf("            Average density is %f distances / Ang^3.\n",density_);
  } else {
    density_ = density_ * 
               (nmask1 * nmask2 - (double)numSameAtoms) / nmask1;
    mprintf("            Density is %f distances / Ang^3.\n",density_);
  }
  // Need to normalize each bin, which holds the particle count at that
  // distance. Calculate the expected number of molecules for that 
  // volume slice. Expected # of molecules is particle density times volume 
  // of each slice:
  // Density * ( [(4/3)*PI*(R+dr)^3] - [(4/3)*PI*(dr)^3] )
  double sum = 0.0;
  for (int bin = 0; bin < numBins_; bin++) {
    //mprintf("DBG:\tNumBins= %i\n",rdf[bin]); 
    // Number of particles in this volume slice over all frames.
    double N = (double) RDF_[bin];
    if (rawrdf_ != 0)
      rawrdf_->Add(bin, &N);
    // r, r + dr
    double R = spacing_ * (double)bin;
    double Rdr = R + spacing_;
    // Volume of slice: 4/3_pi * [(r+dr)^3 - (dr)^3]
    double dv = Constants::FOURTHIRDSPI * ( (Rdr * Rdr * Rdr) - (R * R * R) );
    // Expected # distances in this volume slice
    double expectedD = dv * density_;
    if (debug_>0)
      mprintf("    \tBin %f->%f <Pop>=%f, V=%f, D=%f, norm %f distances.\n",
              R,Rdr,N/numFrames_,dv,density_,expectedD);
    // Divide by # frames
    double norm = expectedD * (double)numFrames_;
    N /= norm;
    Dset_->Add(bin, &N);
    // If specified, calc integral of # mask2 atoms as fn of distance
    if (intrdf_ != 0) {
      sum += N * expectedD / nmask2;
      intrdf_->Add(bin, &sum);
    }
  }
}
