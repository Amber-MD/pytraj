// Action_Rotdif
#include <cmath>
#include <cfloat> // DBL_MAX
#include <cstdio> //sscanf
#include <algorithm> // sort
#include "Action_Rotdif.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // NumberFilename
#include "Constants.h" // TWOPI
#include "DataSet_Mesh.h"
#include "ProgressBar.h"
#include "Corr.h"
#include "CurveFit.h"
#include "SimplexMin.h"

#ifndef NO_MATHLIB
// Definition of Fortran subroutines called from this class
extern "C" {
  // LAPACK
  void dgesvd_(char*, char*, int&, int&, double*,
               int&, double*, double*, int&, double*, int&,
               double*, int&, int& );
  void dsyev_(char*, char*, int&, double*, int&, double*,double*,int&,int&);
  // TODO: Use internal diagonalizer
  // DEBUG
  //void dgemm_(char*,char*,int&,int&,int&,double&,
  //            double*,int&,double*,int&,double&,double*,int&);
  // DEBUG
  //void dgemv_(char*,int&,int&,double&,double*,int&,double*,int&,double&,double*,int&);
}
#endif

// CONSTRUCTOR
Action_Rotdif::Action_Rotdif() :
  debug_(0),
  rseed_( 1 ),
  nvecs_( 0 ),
  tfac_( 0.0 ),
  ti_( 0.0 ),
  tf_( 0.0 ),
  NmeshPoints_( -1 ),
  itmax_( 0 ),
  delmin_( 0.0 ),
  d0_( 0.0 ),
  olegendre_( 2 ),
  ncorr_( 0 ),
  delqfrac_( 0 ),
  amoeba_ftol_(0.0000001 ),
  amoeba_itmax_(10000),
  amoeba_nsearch_(1),
  do_gridsearch_( false ),
  useMass_(false),
  usefft_( true )
{ } 
// TODO: MAKE ANALYSIS
void Action_Rotdif::Help() {
  mprintf("\t[outfile <outfilename>] [usefft]\n"
          "  Options for creating RMS best-fit rotation matrices:\n"
          "\t[<mask>] {%s}\n"
          "\t[rmout <rmOut>]\n"
          "  Options for generating random vectors:\n"
          "\t[nvecs <nvecs>] [rvecin <randvecIn>] [rseed <random seed>]\n"
          "\t[rvecout <randvecOut>]\n"
          "  Options for calculating vector time correlation functions:\n"
          "\t[order <olegendre>] [ncorr <ncorr>] [corrout <corrOut>]\n"
          "  *** The options below only apply if 'usefft' IS NOT specified. ***\n"
          "  Options for calculating local effective D, small anisotropy:\n"
          "\t[deffout <deffOut>] [itmax <itmax>] [tol <tolerance>] [d0 <d0>]\n"
          "\t[nmesh <NmeshPoints>] dt <tfac> [ti <ti>] tf <tf>\n"
          "  Options for calculating D with full anisotropy:\n"
          "\t[amoeba_tol <tolerance>] [amoeba_itmax <iterations>]\n"
          "\t[amoeba_nsearch <n>] [scalesimplex <scale>] [gridsearch]\n"
          "  *** The options below only apply if 'usefft' IS specified. ***\n"
          "  Options for curve-fitting:\n"
          "\t[fit_tol <tolerance>] [fit_itmax <max # iterations>]\n"
          "\n  Calculate rotational diffusion tensor. By default the procedure of\n"
          "  Wong & Case (2008) is used.\n"
          "  The 'usefft' option is currently an EXPERIMENTAL option, where time\n"
          "  correlation functions for each vector are calculated using spherical\n"
          "  harmonics. The average over all time correlation functions is then fit\n"
          "  to the following equation to obtain principal values of the diffusion\n"
          "  tensor:\n    C(t) = SUM[l=-2,...,+2]( cl * exp(-t / Tl)\n"
          "  (see equation 2.10b and related eqs. in Korzhnev et al., Prog. Nuc. Mag.\n"
          "   Res. Spec. 38 (2001) 197-266 for details).\n", DataSetList::RefArgs);
}

// Action_Rotdif::Init()
Action::RetType Action_Rotdif::Init(ArgList& actionArgs, TopologyList* PFL, DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  if (DSL->EnsembleNum() > -1) {
    mprinterr("Error: Rotational Diffusion calc. currently cannot be used in ensemble mode.\n");
    return Action::ERR;
  }
  debug_ = debugIn;
  // Get Keywords
  usefft_ = actionArgs.hasKey("usefft");
  nvecs_ = actionArgs.getKeyInt("nvecs",1000);
  rseed_ = actionArgs.getKeyInt("rseed",80531);
  ncorr_ = actionArgs.getKeyInt("ncorr",0);
  tfac_ = actionArgs.getKeyDouble("dt",0);
  if (tfac_<=0) {
    mprinterr("Error: 'dt <timestep>' must be specified and > 0.\n");
    return Action::ERR;
  }
  ti_ = actionArgs.getKeyDouble("ti",0);
  tf_ = actionArgs.getKeyDouble("tf",0);
  if (tf_ <= ti_) {
    mprinterr("Error: Initial time ti (%f) must be < final time tf (%f).\n",
              ti_, tf_);
    return Action::ERR;
  }
  NmeshPoints_ = actionArgs.getKeyInt("nmesh", -1);
  itmax_ = actionArgs.getKeyInt("itmax",500);
  delmin_ = actionArgs.getKeyDouble("tol",0.000001);
  d0_ = actionArgs.getKeyDouble("d0",0.03);
  olegendre_ = actionArgs.getKeyInt("order",2);
  if (olegendre_!=1 && olegendre_!=2) {
    mprinterr("Error: Order of legendre polynomial (%i) must be 1 or 2.\n",
              olegendre_);
    return Action::ERR;
  }
  delqfrac_ = actionArgs.getKeyDouble("delqfrac", 0.5);
  delqfrac_ = actionArgs.getKeyDouble("scalesimplex", delqfrac_);
  randvecOut_ = actionArgs.GetStringKey("rvecout");
  randvecIn_ = actionArgs.GetStringKey("rvecin");
  rmOut_ = actionArgs.GetStringKey("rmout");
  deffOut_ = actionArgs.GetStringKey("deffout");
  std::string outfilename = actionArgs.GetStringKey("outfile");
  if (outfilename.empty()) outfilename = actionArgs.GetStringKey("out");
  corrOut_ = actionArgs.GetStringKey("corrout");
  do_gridsearch_ = actionArgs.hasKey("gridsearch");
  amoeba_ftol_ = actionArgs.getKeyDouble("amoeba_tol", amoeba_ftol_);
  amoeba_itmax_ = actionArgs.getKeyInt("amoeba_itmax", amoeba_itmax_);
  amoeba_nsearch_ = actionArgs.getKeyInt("amoeba_nsearch", 1);
  if (usefft_) {
    amoeba_ftol_ = actionArgs.getKeyDouble("fit_tol", amoeba_ftol_);
    amoeba_itmax_ = actionArgs.getKeyInt("fit_itmax", amoeba_itmax_);
  }
  // Reference Keywords
  ReferenceFrame REF = DSL->GetReferenceFrame( actionArgs );
  if (REF.error()) return Action::ERR;
  if (REF.empty()) {
    mprinterr("Error: Must specify a reference structure.\n");
    return Action::ERR;
  }
  // Get Masks
  AtomMask RefMask( actionArgs.GetMaskNext() );
  TargetMask_.SetMaskString( RefMask.MaskExpression() );
  
  // Initialize random number generator
  RNgen_.rn_set( rseed_ );

  // Set up reference for RMSD
  // Setup reference mask
  if (REF.Parm().SetupIntegerMask( RefMask )) return Action::ERR;
  if (RefMask.None()) {
    mprinterr("Error: No atoms in reference mask.\n");
    return Action::ERR;
  }
  // Allocate frame for selected reference atoms
  SelectedRef_.SetupFrameFromMask(RefMask, REF.Parm().Atoms());
  // Set reference frame coordinates
  SelectedRef_.SetCoordinates(REF.Coord(), RefMask);
  // Always fitting; Pre-center reference frame
  SelectedRef_.CenterOnOrigin(useMass_); 

  // Open output file. Defaults to stdout if no name specified
  if (outfile_.OpenEnsembleWrite(outfilename, DSL->EnsembleNum())) {
    mprinterr("Error: Could not open Rotdif output file %s.\n", outfilename.c_str());
    return Action::ERR;
  }

  mprintf("    ROTDIF: Rotational diffusion tensor calculation.\n");
  mprintf("\tRotation matrices for rotating vectors will be generated by RMS fitting\n"
          "\t  to atoms in mask '%s', reference '%s'\n", TargetMask_.MaskString(),
          REF.FrameName().base());
  if (!rmOut_.empty())
    mprintf("\tRotation matrices will be written to file '%s'\n", rmOut_.c_str());
  if (randvecIn_.empty())
    mprintf("\tGenerating %i random vectors,", nvecs_);
  else
    mprintf("\tReading %i vectors from file '%s',", nvecs_, randvecIn_.c_str());
  mprintf(" random seed is %i.\n", rseed_); 
  if (!randvecOut_.empty())
    mprintf("\tWriting vectors to file '%s'\n", randvecOut_.c_str());
  mprintf("\tMax length to compute vector time correlation functions:");
  if (ncorr_ == 0)
    mprintf(" Total # of frames.\n");
  else
    mprintf(" %i frames.\n",ncorr_);
  mprintf("\tVector time correlation function order: %i\n", olegendre_);
  if (usefft_) {
    mprintf("Warning: 'usefft' is an EXPERIMENTAL option. Use at your own risk.\n");
    mprintf("\tVector time correlation functions will be calculated using spherical harmonics.\n");
    mprintf("\tVector time correlation time step is %.4g ps\n", tfac_);
    if (!corrOut_.empty())
      mprintf("\tAveraged vector time correlation function and fit curves"
              " will be written to '%s'\n", corrOut_.c_str());
    mprintf("\tCurve fit tolerance= %g, %i iterations.\n", amoeba_ftol_, amoeba_itmax_);
    if (!outfilename.empty())
      mprintf("\tDiffusion constants output to %s\n", outfilename.c_str());
    else
      mprintf("\tDiffusion constants output to STDOUT\n");
  } else {
    mprintf("\tVector time correlation functions will be calculated directly.\n");
    if (!corrOut_.empty())
      mprintf("\tVector time correlation functions will be written to '%s.X'\n",
              corrOut_.c_str());
    mprintf("\tVector time correlation functions assumed to fit single exponential\n"
            "\t  in the limit of small anisotropy.\n");
    mprintf("\tVector time correlation functions will be integrated from\n"
            "\t  %.4g to %.4g, time step %.4g\n", ti_, tf_, tfac_);
    mprintf("\tVector time correlation functions will be smoothed using cubic spline \n"
            "\t  interpolation. Data points will be increased");
    if (NmeshPoints_ != -1)
      mprintf(" by a factor of %i.\n", NmeshPoints_);
    else
      mprintf(" by a factor of 2.\n");
    mprintf("\tIntegral of single exponential iterative solver:\n"
            "\t  iterations= %i, tolerance= %g, initial guess= %g\n",
            itmax_, delmin_, d0_);
    mprintf("\tNelder Mead (downhill simplex) minimizer will be used to determine\n"
            "\t  Q with full anisotropy.\n");
    mprintf("\t  searches= %i, iterations= %i, tolerance= %g, simplex scaling= %g\n",
            amoeba_nsearch_, amoeba_itmax_, amoeba_ftol_, delqfrac_);
    if (do_gridsearch_)
      mprintf("\tGrid search will be performed for Q with full anisotropy (time consuming)\n");
#   ifdef NO_MATHLIB
    mprintf("------------------------------------------------------\n");
    mprintf("Warning: Cpptraj was compiled with -DNO_MATHLIB.\n");
    mprintf("         The final tensor fit cannot be performed.\n");
    mprintf("         Only Deffs will be calculated.\n");
    mprintf("------------------------------------------------------\n");
#   else
    if (!outfilename.empty())
      mprintf("\tDiffusion constants and tau will be written to %s\n",
              outfilename.c_str());
    else
      mprintf("\tDiffusion constants and tau will be written to STDOUT.\n");
#   endif
    if (!usefft_)
      mprintf("# Citation: Wong V.; Case, D. A.; \"Evaluating rotational diffusion from\n"
              "#           protein MD simulations.\"\n"
              "#           J. Phys. Chem. B (2008) V.112 pp.6013-6024.\n");
  }
  return Action::OK;
}

// Action_Rotdif::Setup()
/** Determine what atoms each mask pertains to for the current parm file.
  * Also determine whether imaging should be performed.
  */
Action::RetType Action_Rotdif::Setup(Topology* currentParm, Topology** parmAddress) {

  if ( currentParm->SetupIntegerMask( TargetMask_ ) ) return Action::ERR;
  if ( TargetMask_.None() ) {
    mprintf("Warning: No atoms in mask.\n");
    return Action::ERR;
  }
  // Allocate space for selected atoms in the frame. This will also put the
  // correct masses in based on the mask.
  SelectedTgt_.SetupFrameFromMask(TargetMask_, currentParm->Atoms());
  // Check that num atoms in frame mask from this parm match ref parm mask
  if ( SelectedRef_.Natom() != TargetMask_.Nselected() ) {
    mprinterr("Error: Number of atoms in RMS mask (%i) does not \n",TargetMask_.Nselected());
    mprinterr("Error:   equal number of atoms in Ref mask (%i).\n",SelectedRef_.Natom());
    return Action::ERR;
  }
  
  // Print info for this parm
  mprintf("    ROTDIF: %i atoms selected for RMS fit.\n",TargetMask_.Nselected());
        
  return Action::OK;  
}

// Action_Rotdif::DoAction()
/** Calculate and store the rotation matrix for frame to reference.
  */
Action::RetType Action_Rotdif::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {
  Matrix_3x3 U;
  Vec3 Trans; // Unused

  // Set selected frame atoms. Masses have already been set.
  SelectedTgt_.SetCoordinates(*currentFrame, TargetMask_);
  SelectedTgt_.RMSD_CenteredRef(SelectedRef_, U, Trans, useMass_);
  Rmatrices_.push_back( U );

  return Action::OK;
} 

// ========== ROTATIONAL DIFFUSION CALC ROUTINES ===============================
// Q_to_D()
/** Convert column vector Q with format:
  *   {Qxx, Qyy, Qzz, Qxy, Qyz, Qxz}
  * to diffusion tensor D via the formula:
  *   D = 3*Diso*I - 2Q
  * where
  *   Diso = trace(Q) / 3
  *
  * J. Biomol. NMR 9, 287 (1997) L. K. Lee, M. Rance, W. J. Chazin, A. G. Palmer
  * it has been assumed that the Q(l=1) has the same relationship to
  * D(l=1) as in the l=2 case; we have not yet proved this for the
  * non-symmetric case, but it is true for the axially symmetric case
  * also, Deff(l=1) = 1/(2*tau(l=1)) = e^T * Q(l=1) * e yields the 
  * correct values for tau(l=1) (known from Woessner model type 
  * correlation functions) along principal axe
  */
static void Q_to_D(SimplexMin::Darray const& Q_, Matrix_3x3& D) {
  double tq = Q_[0] + Q_[1] + Q_[2];
  D[0] = tq - (2 * Q_[0]); // tq-2Qxx
  D[1] = -2 * Q_[3];       // -2Qxy
  D[2] = -2 * Q_[5];       // -2Qxz
  D[3] = D[1];            // -2Qyx
  D[4] = tq - (2 * Q_[1]); // tq-2Qyy
  D[5] = -2 * Q_[4];       // -2Qyz
  D[6] = D[2];            // -2Qzx
  D[7] = D[5];            // -2Qzy
  D[8] = tq - (2 * Q_[2]); // tq-2Qzz
}

// D_to_Q()
/** Given diffusion tensor D[9], calculate Q[6] where Q is:
  *   Q = {Qxx, Qyy, Qzz, Qxy, Qyz, Qxz}
  * from
  *   Q = (3*Dav*I - D) / 2
  * where Dav = trace(D). See Q_to_D for more discussion.
  */
static void D_to_Q(Matrix_3x3 const& D, SimplexMin::Darray& Q_) {
  double td = D[0] + D[4] + D[8];
  Q_[0] = (td - D[0]) / 2; // Qxx
  Q_[1] = (td - D[4]) / 2; // Qyy
  Q_[2] = (td - D[8]) / 2; // Qzz
  Q_[3] = -D[1] / 2; // Qxy
  Q_[4] = -D[5] / 2; // Qyz
  Q_[5] = -D[2] / 2; // Qxz
}

// printMatrix()
static void printMatrix(const char *Title, const double *U, int mrows, int ncols) {
  mprintf("    %s",Title);
  int usize = mrows * ncols;
  for (int i = 0; i < usize; i++) {
    if ( (i%ncols)==0 ) mprintf("\n");
    mprintf(" %10.5g",U[i]);
  }
  mprintf("\n");
}

// chi_squared()
/** Calculate chi-squared of residual between Yvals and newY = fxn(Xvals).
  */
static double chi_squared(SimplexMin::SimplexFunctionType fxn,
                          SimplexMin::Darray const& Qin, DataSet* Xvals,
                                  std::vector<double> const& Yvals,
                                  std::vector<double>& newY)
{
#ifdef NO_MATHLIB
  return -1;
#else
  fxn(Xvals, Qin, newY);

  double chisq = 0.0;
  for (unsigned int i = 0; i != Yvals.size(); i++) {
    double diff = Yvals[i] - newY[i];
    chisq += (diff * diff); 
  }

  return chisq;  
#endif
}

// calculate_D_properties()
/** Given the principal components of D, calculate the isotropic diffusion
  * constant (Dav = (Dx + Dy + Dz) / 3), anisotropy (Dan = 2Dz / (Dx + Dy)),
  * and rhombicity (Drh = (3/2)(Dy - Dx) / (Dz - 0.5(Dx + Dy)).
  */
static Vec3 calculate_D_properties(Vec3 const& Dxyz) {
  double Dx = Dxyz[0];
  double Dy = Dxyz[1];
  double Dz = Dxyz[2];
  return Vec3( (Dx + Dy + Dz) / 3,
               (2 * Dz) / (Dx + Dy),
               (1.5 * (Dy - Dx)) / (Dz - (0.5 * (Dx + Dy))) );
}

/** Diagonalize 3x3 matrix Mat; eigenvectors are returned in columns due to
  * the fortran call. Workspace is fixed at 102, optimal value returned for
  * 3x3 matrix by LAPACK. Eigenvalues returned in Vec.
  */
static void Diagonalize(Matrix_3x3& Mat, Vec3& Vec) {
  int n_cols = 3, lwork = 102, info;
  double work[102];
# ifndef NO_MATHLIB
  dsyev_((char*)"Vectors", (char*)"Upper", n_cols, Mat.Dptr(), n_cols,
         Vec.Dptr(), work, lwork, info);
  if (info > 0)
    mprinterr("Error: The algorithm computing the eigenvalues/eigenvectors of D failed to converge.\n");
# endif
}

// =============================================================================
// Action_Rotdif::RandomVectors()
/** If no input file is specified by randvecIn, generate nvecs vectors of length 
  * 1.0 centered at the coordinate origin that are randomly oriented. The x,
  * y, and z components of each vector are generated from a polar coordinate 
  * system in which phi is randomly chosen in the range 0 to 2*PI and theta
  * is randomly chosen in the range 0 to PI/2 (single hemisphere).
  *   R = 1.0
  *   x = R * sin(theta) * cos(phi)
  *   y = R * sin(theta) * sin(phi)
  *   z = R * cos(theta)
  * \return Array containing random vectors; array will be empty on error.
  */
// NOTE: Theta could also be generated in the same way as phi. Currently done
//       to be consistent with the original implementation in randvec.F90
DataSet_Vector Action_Rotdif::RandomVectors() {
  DataSet_Vector XYZ;
  XYZ.ReserveVecs( nvecs_ );
  // ----- Read nvecs vectors from a file
  if (!randvecIn_.empty()) {
    CpptrajFile vecIn;
    if (vecIn.OpenRead(randvecIn_)) {
      mprinterr("Error: Could not open random vectors input file %s",randvecIn_.c_str());
      return XYZ;
    }
    Vec3 xIn;
    for (int i = 0; i < nvecs_; i++) {
      const char* buffer = vecIn.NextLine();
      if ( buffer == 0 ) {
        mprinterr("Error: Could not read vector %i from file %s\n", i+1,randvecIn_.c_str());
        XYZ.reset();
        return XYZ;
      }
      sscanf(buffer,"%*i %lf %lf %lf", xIn.Dptr(), xIn.Dptr()+1, xIn.Dptr()+2);
      xIn.Normalize();
      XYZ.AddVxyz( xIn );
    }
    vecIn.CloseFile();
  // ----- Generate nvecs normalized vectors
  } else {
    for (int i = 0; i < nvecs_; i++) {
      double phi = Constants::TWOPI * RNgen_.rn_gen();
      // XYZ[i+2] is cos(theta)
      double costheta = 1 - RNgen_.rn_gen();
      double theta = acos( costheta );
      double sintheta = sin( theta );
      XYZ.AddVxyz( Vec3(sintheta*cos(phi), sintheta*sin(phi), costheta) );
    }
  }
  // Print vectors
  if (!randvecOut_.empty()) {
    CpptrajFile rvout;
    if (rvout.OpenWrite(randvecOut_)) {
      mprinterr("Error: Could not set up %s for writing vectors.\n",randvecOut_.c_str());
    } else {
      int idx = 1;
      for (DataSet_Vector::const_iterator vec = XYZ.begin(); vec != XYZ.end(); ++vec)
        rvout.Printf("%6i  %15.8f  %15.8f  %15.8f\n",
                     idx++, (*vec)[0], (*vec)[1], (*vec)[2]);
      rvout.CloseFile();
    }
  }

  return XYZ;
}

// =============================================================================
// Action_Rotdif::PrintMatrix()
void Action_Rotdif::PrintMatrix(CpptrajFile& outfile, const char* Title, Matrix_3x3 const& U)
{
  outfile.Printf("    %s\n",Title);
  outfile.Printf(" %12.5e %12.5e %12.5e\n %12.5e %12.5e %12.5e\n %12.5e %12.5e %12.5e\n",
                 U[0], U[1], U[2], U[3], U[4], U[5], U[6], U[7], U[8]);
}

// Action_Rotdif::PrintVector()
void Action_Rotdif::PrintVector(CpptrajFile& outfile, const char* Title, Vec3 const& V)
{
  outfile.Printf("    %s\n",Title);
  outfile.Printf(" %12.5e %12.5e %12.5e\n", V[0], V[1], V[2]);
}

// Action_Rotdif::PrintVec6()
void Action_Rotdif::PrintVec6(CpptrajFile& outfile, const char* Title, SimplexMin::Darray const& V)
{
  outfile.Printf("    %s\n",Title);
  outfile.Printf(" %12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n", 
                 V[0], V[1], V[2], V[3], V[4], V[5]);
}

// Action_Rotdif::PrintTau()
void Action_Rotdif::PrintTau(std::vector<double> const& Tau)
{
  outfile_.Printf("     taueff(obs) taueff(calc)\n");
  for (int i = 0; i < nvecs_; i++)
    outfile_.Printf("%5i %12.5e %12.5e\n", i+1, D_eff_[i], Tau[i]);
}

// =============================================================================
/** Compute tau(l=2) with full anisotropy given vectors and Q. Q is converted
  * to the diffusion tensor D, and the eigenvectors/eigenvalues are computed.
  * n*ez = cos(theta), n*ey = sin(theta)*sin(phi), n*ex = sin(theta)*cos(phi)
  * It is expected that the principal components Dxyz and principal axes
  * D_matrix are in ascending order, i.e. Dx <= Dy <= Dz. It is also expected
  * that the eigenvectors are orthonormal. This is the default behavior when 
  * the LAPACK routines are used. Also, since LAPACK routines are being used
  * from fortran it is expected the matrix of eigenvectors will be in column
  * major order, i.e.:
  *   x0 y0 z0
  *   x1 y1 z1
  *   x2 y2 z2
  */
int AsymmetricFxn_L2(DataSet* Xvals, SimplexMin::Darray const& Qin, SimplexMin::Darray& Yvals) 
{
  Matrix_3x3 d_tensor;
  Vec3 d_xyz;

  // Convert input Q to Diffusion tensor D
  Q_to_D(Qin, d_tensor);
  // Diagonalize D; it is assumed workspace (work, lwork) set up prior 
  Diagonalize( d_tensor, d_xyz );

  double lambda[5];
  double dx = d_xyz[0];
  double dy = d_xyz[1];
  double dz = d_xyz[2];

  // DEBUG
//  printMatrix("Dxyz",d_xyz.Dptr(),1,3);
//  printMatrix("Dvec",d_tensor.Dptr(),3,3);

  // tau = sum(m){amp(l,m)/lambda(l,m)}
  // See Korzhnev DM, Billeter M, Arseniev AS, Orekhov VY; 
  // Prog. Nuc. Mag. Res. Spec., 38, 197 (2001) for details, Table 3 in 
  // particular. Only weights need to be computed for each vector, decay 
  // constants can be computed once. Decay constants correspond 
  // to l=2, m=-2,-1,0,+1,+2.

  // l=2, m=-2 term: lambda(2,-2) = Dx + Dy + 4*Dz
  lambda[0] = dx + dy + (4*dz);
  // l=2, m=-1 term: lambda(2,-1) = Dx + 4*Dy + Dz
  lambda[1] = dx + (4*dy) + dz;
  // l=2, m=0 term:  lambda(2, 0) = 6*[Dav - sqrt(Dav*Dav - Dpr*Dpr)]
  //                 Dav = (Dx + Dy + Dz)/3 
  //                 Dpr = sqrt[(Dx*Dy + Dy*Dz + Dx*Dz)/3 ]
  //                 Dpr2 = Dpr*Dpr = (Dx*Dy + Dy*Dz + Dx*Dz)/3
  double Dav = (dx + dy + dz) / 3;
  double Dpr2 = ((dx*dy) + (dy*dz) + (dx*dz)) / 3;
  if (Dpr2 < 0) {
    //mprinterr("Error: Rotdif::calc_Asymmetric: Cannot calculate Dpr (Dpr^2 < 0)\n");
    // NOTE: Original code just set Dpr to 0 and continued. 
    Dpr2 = 0;
  }
  double delta = (Dav*Dav) - Dpr2;
  if (delta < 0) {
    mprinterr("Error: calc_Asymmetric: Cannot calculate lambda l=2, m=0\n");
    return 1;
  }
  double sqrt_Dav_Dpr = sqrt( delta );
  lambda[2] = 6 * (Dav - sqrt_Dav_Dpr);
  // l=2, m=+1 term: lambda(2,+1) = 4*Dx + Dy + Dz
  lambda[3] = (4*dx) + dy + dz;
  // l=2, m=+2 term: lambda(2,+2) = 6*[Dav + sqrt(Dav*Dav - Dpr*Dpr)]
  lambda[4] = 6 * (Dav + sqrt_Dav_Dpr);

  // Check all normalization factors for 0.0, which can lead to inf values
  // propogating throughout the calc. Set to SMALL instead; will get 
  // very large numbers but should still not overflow.
//  mprintf("lambda=");
  for (int i = 0; i < 5; i++) {
    if (lambda[i] < Constants::SMALL) lambda[i] = Constants::SMALL;
//    mprintf(" %g", lambda[i]);
  }
//  mprintf("\n");

  // Loop over all random vectors
  int nvec = 0; // index into Yvals, sumc2
  DataSet_Vector const& random_vectors = static_cast<DataSet_Vector const&>( *Xvals );
  for (DataSet_Vector::const_iterator randvec = random_vectors.begin();
                                      randvec != random_vectors.end();
                                    ++randvec, ++nvec)
  {
    // Rotate vector i into D frame
    // This is an inverse rotation. Since matrix_D is in column major order
    // however, do a normal rotation instead.
    Vec3 rotated_vec = d_tensor * (*randvec);
    double dot1 = rotated_vec[0];
    double dot2 = rotated_vec[1];
    double dot3 = rotated_vec[2];
    //mprintf("DBG: dot1-3= %10.5g %10.5g %10.5g\n",dot1,dot2,dot3);
    
    // ----- Compute correlation time for l=2
    // pre-calculate squares
    double dot1_2 = dot1*dot1;
    double dot2_2 = dot2*dot2;
    double dot3_2 = dot3*dot3;
    //     m=-2 term:
    //     lambda(2,-2) = Dx + Dy + 4*Dz
    //     amp(2,-2) = 0.75*sin(theta)^4*sin(2*phi)^2 = 3*(l^2)*(m^2)
    double m_m2 = 3*dot1_2*dot2_2;
    //     m=-1 term:
    //     lambda(2,-1) = Dx + 4*Dy + Dz
    //     amp(2,-1) = 0.75*sin(2*theta)^2*cos(phi)^2 = 3*(l^2)*(n^2)
    double m_m1 = 3*dot1_2*dot3_2;
    //     m=0 term:
    //     lambda(2,0) = 6*[Dav - sqrt(Dav*Dav - Dpr*Dpr)]
    //         Dav = (Dx + Dy + Dz)/3 
    //         Dpr = sqrt[(Dx*Dy + Dy*Dz + Dx*Dz)/3 ]
    //     amp(2,0) = (w/N)^2*0.25*[3*cos(theta)^2-1]^2 +
    //                (u/N)^2*0.75*sin(theta)^4*cos(2*phi)^2 -
    //                (u/delta)*[sqrt(3)/8]*[3*cos(theta)^2-1]*sin(theta)^2*cos(2*phi)
    //         u = sqrt(3)*(Dx - Dy)
    //         delta = 3*sqrt(Dav*Dav - Dpr*Dpr)
    //         w = 2*Dz - Dx - Dy + 2*delta
    //         N = 2*sqrt(delta*w)
    delta = 3 * sqrt_Dav_Dpr;
    double dot1_4 = dot1_2 * dot1_2;
    double dot2_4 = dot2_2 * dot2_2;
    double dot3_4 = dot3_2 * dot3_2;
    double da = 0.25 * ( 3*(dot1_4 + dot2_4 + dot3_4) - 1);
    double ea = 0;
    if ( delta > Constants::SMALL) { 
      double epsx = 3*(dx-Dav)/delta;
      double epsy = 3*(dy-Dav)/delta;
      double epsz = 3*(dz-Dav)/delta;
      double d2d3 = dot2*dot3;
      double d1d3 = dot1*dot3;
      double d1d2 = dot1*dot2;
      ea = epsx*(3*dot1_4 + 6*(d2d3*d2d3) - 1) + 
           epsy*(3*dot2_4 + 6*(d1d3*d1d3) - 1) + 
           epsz*(3*dot3_4 + 6*(d1d2*d1d2) - 1);
      ea /= 12;
    }
    double m_0 = da + ea;
    //     m=+1 term:
    //     lambda(2,+1) = 4*Dx + Dy + Dz
    //     amp(2,+1) = 0.75*sin(2*theta)^2*sin(phi)^2 
    double m_p1 = 3*dot2_2*dot3_2;
    //     m=+2 term:
    //     lambda(2,+2) = 6*[Dav + sqrt(Dav*Dav = Dpr*Dpr)]
    //     amp(2,+2) = (u/n)^2*0.25*[3*cos(theta)^2-1]^2 +
    //                 (w/n)^2*0.75*sin(theta)^4*cos(2*phi)^2 +
    //                 (u/delta)*[sqrt(3)/8]*[3*cos(theta)^2-1]*sin(theta)^2*cos(2*phi)
    double m_p2 = da - ea;
    // Sum decay constants and weights for l=2 m=-2...2
    Yvals[nvec] = (m_m2 / lambda[0]) + (m_m1 / lambda[1]) + (m_0 / lambda[2]) +
                  (m_p1 / lambda[3]) + (m_p2 / lambda[4]);
    //sumc2_[nvec] = m_m2 + m_m1 + m_0 + m_p1 + m_p2;
//    mprintf("Yvals[%i]= (%g / %g) + (%g / %g) + (%g / %g) + (%g / %g) + (%g / %g) = %g\n", nvec,
//            m_m2, lambda[0], m_m1, lambda[1], m_0, lambda[2], m_p1, lambda[3], m_p2, lambda[4], Yvals[nvec]);
  }

  return 0;
}

// -----------------------------------------------
/** Compute tau(l=1) with full anisotropy given vectors and Q. See 
  * AsymmetricFxn_L2 for more details.
  */
int AsymmetricFxn_L1(DataSet* Xvals, SimplexMin::Darray const& Qin, SimplexMin::Darray& Yvals) 
{
  Matrix_3x3 d_tensor;
  Vec3 d_xyz;

  // Convert input Q to Diffusion tensor D
  Q_to_D(Qin, d_tensor);
  // Diagonalize D; it is assumed workspace (work, lwork) set up prior 
  Diagonalize( d_tensor, d_xyz );
  // to this call.
  // NOTE: Due to the fortran call, the eigenvectors are returned in COLUMN
  //       MAJOR order.

  double lambda[3];
  double dx = d_xyz[0];
  double dy = d_xyz[1];
  double dz = d_xyz[2];

  // DEBUG
  //printMatrix("Dxyz",Dxyz,1,3);
  //printMatrix("Dvec",matrix_D,3,3);

  // tau = sum(m){amp(l,m)/lambda(l,m)}
  // See Korzhnev DM, Billeter M, Arseniev AS, Orekhov VY; 
  // Prog. Nuc. Mag. Res. Spec., 38, 197 (2001) for details, Table 3 in 
  // particular. Only weights need to be computed for each vector, decay 
  // constants can be computed once. Three decay constants correspond 
  // to l=1, m=-1,0,+1

  // l=1, m=-1 term: lambda(1,-1) = Dy + Dz
  lambda[0] = dy + dz;
  // l=1, m=0 term:  lambda(1, 0) = Dx + Dy
  lambda[1] = dx + dy;
  // l=1, m=+1 term: lambda(1,+1) = Dx + Dz
  lambda[2] = dx + dz;

  // Check all normalization factors for 0.0, which can lead to inf values
  // propogating throughout the calc. Set to SMALL instead; will get 
  // very large numbers but should still not overflow.
  for (int i = 0; i < 3; i++) 
    if (lambda[i] < Constants::SMALL) lambda[i] = Constants::SMALL;

  // Loop over all random vectors
  int nvec = 0; // index into Yvals 
  DataSet_Vector const& random_vectors = static_cast<DataSet_Vector const&>( *Xvals );
  for (DataSet_Vector::const_iterator randvec = random_vectors.begin();
                                      randvec != random_vectors.end();
                                    ++randvec, ++nvec)
  {
    // Rotate vector i into D frame
    // This is an inverse rotation. Since matrix_D is in column major order
    // however, do a normal rotation instead.
    Vec3 rotated_vec = d_tensor * (*randvec);
    double dot1 = rotated_vec[0];
    double dot2 = rotated_vec[1];
    double dot3 = rotated_vec[2];
    //mprintf("DBG: dot1-3= %10.5g%10.5g%10.5g\n",dot1,dot2,dot3);
    
    // assuming e(3)*n = cos(theta), e(1)*n = sin(theta)*cos(phi),
    // e(2)*n = sin(theta)*sin(phi), theta >= 0;
    // sin(theta) = sqrt(1 - cos(theta)^2) and theta = tan^-1[sin(theta)/cos(theta)];
    // phi = tan^-1[sin(phi)/cos(phi)] = tan^-1[(e(2)*n)/(e(1)*n)]
    double theta = atan2( sqrt( 1 - (dot3*dot3) ), dot3 );
    double phi = atan2( dot2, dot1 );

    // Pre-calculate sines and cosines
    double sintheta = sin(theta);
    double costheta = cos(theta);
    double sinphi = sin(phi);
    double cosphi = cos(phi);

    // ----- Compute correlation time for l=1
    // tau(l=1) = [sin(theta)^2*cos(phi)^2]/(Dy+Dz) +
    //            [sin(theta)^2*sin(phi)^2]/(Dx+Dz) +
    //            [cos(theta)^2           ]/(Dx+Dy)
    double sintheta2 = sintheta * sintheta;
    //     cth2=cth*cth
    //     sphi2=sphi*sphi
    //     cphi2=cphi*cphi
    Yvals[nvec] = ((sintheta2 * (cosphi*cosphi)) / lambda[0]) +
                  ((sintheta2 * (sinphi*sinphi)) / lambda[2]) +
                  ((costheta*costheta)           / lambda[1]);
  }

  return 0;
}

// =============================================================================
// Action_Rotdif::Tensor_Fit()
/** Based on random_vectors and effective diffusion constants D_eff previously
  * calculated, first find the tensor Q (and therefore D) in the small
  * anisotropic limit by solving:
  *   D_eff(n) = At(n) * Q
  * where At(n) is composed of D_eff vector components:
  *   { x^2, y^2, z^2, 2xy, 2yz, 2xz }
  * \param vector_q Will be set with Q tensor for small anisotropic limit.
  */
int Action_Rotdif::Tensor_Fit(SimplexMin::Darray& vector_q) {
#ifdef NO_MATHLIB
  return 1;
#else
  //double cut_ratio = 0.000001; // threshold ratio for removing small singular values in SVD

  mprintf("\tDetermining diffusion tensor with small anisotropy.\n");
  // Generate matrix At
  // NOTE: The LAPACK fortran routines are COLUMN MAJOR, so m and n must be 
  //       flipped, i.e. matrix At must be tranposed before passing it in
  //       so we actually need to construct matrix A for SVD.
  // NOTE: LAPACK SVD routine destroys matrix A which is needed later, so 
  //       create a non-flipped version (i.e. matrix At) as well.
  int m_rows = nvecs_;
  int n_cols = 6;
  double *matrix_A = new double[ m_rows * n_cols ];
  double *matrix_At = new double[ m_rows * n_cols ];
  // Pre-compute array offsets of matrix_A for performing transpose
  int A1 = m_rows;
  int A2 = m_rows * 2;
  int A3 = m_rows * 3;
  int A4 = m_rows * 4;
  int A5 = m_rows * 5;
  double* At = matrix_At;
  int nvec = 0;
  for (DataSet_Vector::const_iterator randvec = random_vectors_.begin();
                                      randvec != random_vectors_.end(); ++randvec)
  {
    // Transpose of matrix A
    double *A = matrix_A + nvec;
    A[ 0] = (*randvec)[0] * (*randvec)[0];     // x^2
    A[A1] = (*randvec)[1] * (*randvec)[1];     // y^2
    A[A2] = (*randvec)[2] * (*randvec)[2];     // z^2
    A[A3] = 2*((*randvec)[0] * (*randvec)[1]); // 2xy
    A[A4] = 2*((*randvec)[1] * (*randvec)[2]); // 2yz
    A[A5] = 2*((*randvec)[0] * (*randvec)[2]); // 2xz
    // Matrix A
    At[0] = A[ 0];
    At[1] = A[A1];
    At[2] = A[A2];
    At[3] = A[A3];
    At[4] = A[A4];
    At[5] = A[A5];
    At += 6;
    ++nvec;
  }
  if (debug_>1) {
    printMatrix("matrix_A",matrix_A,n_cols,m_rows);
    printMatrix("matrix_At",matrix_At,m_rows,n_cols);
  }

  // Perform SVD on matrix At to generate U, Sigma, and Vt
  int lda = m_rows;
  int ldu = m_rows;
  int ldvt = n_cols;
  // NOTE: Final dimension of matrix_S is min(m,n)
  int s_dim = m_rows;
  if (n_cols < m_rows) s_dim = n_cols;
  double *matrix_S = new double[ s_dim ];
  double *matrix_U = new double[ m_rows * m_rows ];
  double *matrix_Vt = new double[ n_cols * n_cols ];
  int lwork = -1, info;
  // Allocate Workspace
  double wkopt = 0.0;
  dgesvd_((char*)"All",(char*)"All",m_rows, n_cols, matrix_A, lda, 
          matrix_S, matrix_U, ldu, matrix_Vt, ldvt, &wkopt, lwork, info );
  lwork = (int)wkopt;
  std::vector<double> workspace(lwork);
  // Compute SVD
  dgesvd_((char*)"All",(char*)"All", m_rows, n_cols, matrix_A, lda, 
          matrix_S, matrix_U, ldu, matrix_Vt, ldvt, &workspace[0], lwork, info );
  // matrix_A and work no longer needed
  delete[] matrix_A;
  workspace.clear();
  // DEBUG - Print Sigma
  if (debug_>1) {
    for (int i = 0; i < s_dim; i++) 
      mprintf("Sigma %6i %12.6g\n",i+1,matrix_S[i]);
  }
  // Check for convergence
  if ( info > 0 ) {
    mprinterr("Error: The algorithm computing SVD of At failed to converge.\n" );
    delete[] matrix_At;
    delete[] matrix_U;
    delete[] matrix_S;
    delete[] matrix_Vt;
    return 1;
  }

  // DEBUG: Print U and Vt
  // NOTE: U and Vt are in column-major order from fortran routine
  //       so are currently implicitly transposed.
  if (debug_>1) {
    printMatrix("matrix_Ut",matrix_U,m_rows,m_rows);
    printMatrix("matrix_V",matrix_Vt,n_cols,n_cols);
  }

  // Remove small singular values from SVD
/*  double wmax = 0;
  for (int i=0; i < n; i++)
    if (matrix_S[i] > wmax) wmax=matrix_S[i];
  double wmin=wmax*cut_ratio;
  for (int i=0; i < n; i++)
    if (matrix_S[i] < wmin) matrix_S[i]=0;*/

  // Calculate x = V * Sigma^-1 * Ut * b
  // Take advantage of the fact that Vt and U have already been implicitly 
  // transposed in dgesvd to do everything in one step.
  // First invert Sigma.
  for (int s = 0; s < s_dim; s++) 
    if (matrix_S[s] > 0) matrix_S[s] = 1 / matrix_S[s];
  int v_idx = 0;
  for (int n = 0; n < 6; n++) {
    vector_q[n] = 0.0;
    for (int m = 0; m < m_rows; m++) {
      double vsu_sum = 0.0;
      //mprintf("Row %i (",m);
      for (int vs = 0; vs < s_dim; vs++) {
        vsu_sum += (matrix_Vt[v_idx+vs] * matrix_S[vs] * matrix_U[(vs*m_rows)+m]);
        //mprintf(" VSU v=%i s=%i u=%i",v_idx+vs,vs,(vs*m_rows)+m);
      }
      //mprintf(") m=%i\n",m);
      vector_q[n] += (vsu_sum * D_eff_[m]);
    }
    v_idx += 6;
  }
  // matrix_U, matrix_Vt and matrix_S no longer needed
  delete[] matrix_S;
  delete[] matrix_Vt;
  delete[] matrix_U;

  outfile_.Printf("Results of small anisotropy (SVD) analysis:\n");
  // Print Q vector
  PrintVec6(outfile_,"Qxx Qyy Qzz Qxy Qyz Qxz",vector_q);
  // ---------------------------------------------

  // Convert vector Q to diffusion tensor D
  Q_to_D(vector_q, D_tensor_);
  PrintMatrix(outfile_,"D_tensor",D_tensor_);

  // Save D for later use in back calculating deff from Q
  Matrix_3x3 matrix_D_local = D_tensor_;

  // Diagonalize D tensor to find eigenvalues and eigenvectors
  Diagonalize( D_tensor_, D_XYZ_ );

  // eigenvectors are stored in columns due to implicit transpose from fortran,
  // i.e. Ex = {D_tensor[0], D_tensor[3], D_tensor[6]} etc.
  // Transpose back.
  //matrix_transpose_3x3( D_tensor );

  // Print eigenvalues/eigenvectors
  PrintVector(outfile_,"D eigenvalues",D_XYZ_);
  PrintMatrix(outfile_,"D eigenvectors (in columns)",D_tensor_);

  // Calculate Dav, Daniso, Drhomb
  Vec3 d_props = calculate_D_properties(D_XYZ_);
  PrintVector(outfile_,"Dav, Daniso, Drhomb",d_props);

  // Back-calculate the local diffusion constants via At*Q=Deff
  // First convert original D back to Q
  SimplexMin::Darray vector_q_local(6);
  D_to_Q(matrix_D_local, vector_q_local);
  if (debug_>1) {
    mprintf("    D_to_Q\n %10.5g %10.5g %10.5g %10.5g %10.5g %10.5g\n",
            vector_q_local[0], vector_q_local[1], vector_q_local[2],
            vector_q_local[3], vector_q_local[4], vector_q_local[4]);
  }
  std::vector<double> deff_local;
  deff_local.reserve( nvecs_ );
  At = matrix_At;
  // At*Q
  for (int i=0; i < nvecs_; i++) {
    deff_local.push_back( (At[0] * vector_q_local[0]) +
                          (At[1] * vector_q_local[1]) +
                          (At[2] * vector_q_local[2]) +
                          (At[3] * vector_q_local[3]) +
                          (At[4] * vector_q_local[4]) +
                          (At[5] * vector_q_local[5])   );
    At += 6;
  }
  // Convert deff to tau, Output
  double sgn = 0;
  for (int i = 0; i < nvecs_; i++) {
    // For the following chisq fits, convert deff to taueff
    D_eff_[i] = 1 / (6 * D_eff_[i]);
    deff_local[i] = 1 / (6 * deff_local[i]);
    // NOTE: in rotdif code, sig is 1.0 for all nvecs 
    double diff = deff_local[i] - D_eff_[i];
    sgn += (diff * diff);
  }
  PrintTau( deff_local );
  outfile_.Printf("  chisq for above is %15.5g\n",sgn);

  // Cleanup
  delete[] matrix_At;

  // DEBUG -------------------------------------------------
/* 
  // CHECK: Ensure U is orthogonal (U * Ut) = I
  double *Ut = matrix_transpose(matrix_U,m_rows,m_rows);
  printMatrix("Ut",Ut,m_rows,m_rows);
  double *UtU = new double[ m_rows * m_rows ];
  dgemm_((char*)"N",(char*)"N",m_rows,m_rows,m_rows,alpha,Ut,m_rows,
         matrix_U,m_rows,beta,UtU,m_rows);
  printMatrix("UtU",UtU,m_rows,m_rows);
  delete[] UtU;
  delete[] Ut;
*/
/*
  // CHECK: Ensure that U * S * Vt = A
  double *diag_S = new double[ m_rows * n_cols ];
  int sidx=0;
  int diag_sidx=0;
  for (int m=0; m < m_rows; m++) {
    for (int n=0; n < n_cols; n++) {
      if (m==n) 
        diag_S[diag_sidx]=matrix_S[sidx++];
      else
        diag_S[diag_sidx]=0.0;
      ++diag_sidx;
    }
  }
  printMatrix("diag_S",diag_S,m_rows,n_cols);
  // Transpose diag_S for fortran; U and Vt are already transposed (result from fortran SVD)
  double *diag_St = matrix_transpose( diag_S, m_rows, n_cols );
  printMatrix("diag_St",diag_St,n_cols,m_rows);
  // U * S
  double *US = new double[ m_rows * n_cols ];
  dgemm_((char*)"N",(char*)"N",m_rows,n_cols,m_rows,alpha,matrix_U,m_rows,
         diag_St,m_rows,beta,US,m_rows);
  // US * Vt
  double *USVt = new double[ m_rows * n_cols ];
  dgemm_((char*)"N",(char*)"N",m_rows,n_cols,n_cols,alpha,US,m_rows,
         matrix_Vt,n_cols,beta,USVt,m_rows);
  // The result from fortran will actually be (USVt)t
  //double *USVtt = matrix_transpose( USVt, m_rows, n_cols );
  printMatrix("matrix_A",matrix_A,m_rows,n_cols);
  printMatrix("USVt",USVt,n_cols,m_rows);
  // ----- Cleanup
  delete[] diag_S;
  delete[] diag_St;
  delete[] US;
  delete[] USVt;
*/
/*  
  // ---------------------------------------------
  // CHECK: Calculate x = V * Sigma^-1 * Ut * b with BLAS routines
  double alpha = 1.0;
  double beta = 0.0;
  // First create transposed (for fortran) n x m sigma^-1
  double *matrix_St = new double[ n_cols * m_rows];
  int ndiag = n_cols+1;
  int sidx = 0;
  for (int i = 0; i < n_cols*m_rows; i++) {
    if ( (i%ndiag)==0 && sidx<s_dim) 
      matrix_St[i] = (1 / matrix_S[sidx++]);
    else
      matrix_St[i] = 0;
  }
  printMatrix("matrix_St",matrix_St,m_rows,n_cols);
  // 1) V * Sigma^-1
  //    Call dgemm with first arg "T" to indicate we want matrix_Vt to
  //    be transposed. 
  double *matrix_VS = new double[ n_cols * m_rows ];
  dgemm_((char*)"T",(char*)"N",n_cols,m_rows,n_cols,alpha,matrix_Vt,n_cols,
         matrix_St,n_cols,beta,matrix_VS,n_cols);
  // matrix_S, matrix_Vt, and matrix_St no longer needed
  delete[] matrix_S;
  delete[] matrix_Vt;
  delete[] matrix_St;
  // 2) VSigma^-1 * Ut
  //    Call dgemm with second arg "T" to indicate we want matrix_U to
  //    be transposed.
  double *matrix_VSUt = new double[ n_cols * m_rows ];
  dgemm_((char*)"N",(char*)"T",n_cols,m_rows,m_rows,alpha,matrix_VS,n_cols,
         matrix_U,m_rows,beta,matrix_VSUt,n_cols);
  // matrix_VS and matrix_U no longer needed
  delete[] matrix_VS;
  delete[] matrix_U;
  // 3) VSigma^-1*Ut * b
  int incx = 1;
  dgemv_((char*)"N",n_cols,m_rows,alpha,matrix_VSUt,n_cols,deff,incx,beta,vector_q,incx);
  // matrix_VSUt no longer needed
  delete[] matrix_VSUt;
  // ---------------------------------------------
*/
  // END DEBUG ---------------------------------------------

  return 0;
#endif
}

// =============================================================================
// Action_Rotdif::fft_compute_corr()
int Action_Rotdif::fft_compute_corr(DataSet_Vector const& rotated_vectors, int nsteps, 
                                    std::vector<double>& pY)
{
  int n_of_vecs = rotated_vectors.Size();
  // Zero output array
  pY.assign(nsteps, 0.0);

  // Calculate correlation fn
  CorrF_FFT pubfft( n_of_vecs );
  ComplexArray data1 = pubfft.Array();
  // Loop over m = -olegendre, ..., +olegendre
  for (int midx = -olegendre_; midx <= olegendre_; ++midx) { 
    data1.Assign(rotated_vectors.SphericalHarmonics(midx));
    // Pad with zeros at the end
    // TODO: Does this always have to be done or can it just be done once
    data1.PadWithZero( n_of_vecs );
    pubfft.AutoCorr(data1);
    // Sum into pX
    for (int i = 0; i < nsteps; ++i)
      pY[i] += data1[i*2];
  }
  // Normalize correlation fn
  // 4/3*PI and 4/5*PI due to spherical harmonics addition theorem
  double norm = DataSet_Vector::SphericalHarmonicsNorm( olegendre_ );
  for (int i = 0; i < nsteps; ++i)
    pY[i] *= (norm / (n_of_vecs - i));

  return 0;
}

// Single exponential function
int ExpFxn(CurveFit::Darray const& Tvals, CurveFit::Darray const& Params,
              CurveFit::Darray& Ct)
{
  for (unsigned int n = 0; n != Tvals.size(); n++)
    Ct[n] = ( exp( -Tvals[n] * Params[0] ) );
  return 0;
}

// Integral of exponential with constant
//double IntExpFxn(double tau, double a0) {
//  return ( exp( a0 * tau ) / a0 );
//}

// Sum of 5 exponentials
int Sum5Exp(CurveFit::Darray const& Tvals, CurveFit::Darray const& Params,
            CurveFit::Darray& Ct)
{
  double penalty1 = Params[0] + Params[2] + Params[4] +
                   Params[6] + Params[8];
  penalty1 = 1000.0 * (1.0 - penalty1);

  double penalty2 = 0.0;
  for (int i = 1; i != 11; i += 2)
    if (Params[i] < 0.0) penalty2 += 200.0;

  for (unsigned int n = 0; n != Tvals.size(); n++) {
    double tau = -Tvals[n];
    Ct[n]= ( Params[0] * exp( tau * Params[1] ) +
             Params[2] * exp( tau * Params[3] ) +
             Params[4] * exp( tau * Params[5] ) +
             Params[6] * exp( tau * Params[7] ) +
             Params[8] * exp( tau * Params[9] ) ) + penalty1 + penalty2;
  }
  return 0;
}

/// Determine if penalty should be used during fit.
static bool USE_PENALTY = false;

// C(tau) = SUM(-l...l)[ c * exp(-tau / T) ]
// e-2 = (3*l^2*m^2) * exp( -tau * (x + y + 4z) )
// e-1 = (3*l^2*n^2) * exp( -tau * (x + 4y + z) )
// e0  = (d + e)     * exp( -tau * 6*(D-sqrt(D^2 - D'^2)) )
// e+1 = (3*m^2*n^2) * exp( -tau * (4x + y + z) )
// e+2 = (d - e)     * exp( -tau * 6*(D+sqrt(D^2 - D'^2)) )
// Params={l, m, n, x, y, z}
int Ctau_L2(CurveFit::Darray const& Tau, CurveFit::Darray const& Params,
               CurveFit::Darray& Ct)
{
  double lambda[5];
  double dot1 = Params[0]; // l
  double dot2 = Params[1]; // m
  double dot3 = Params[2]; // n
  double dx = Params[3];
  double dy = Params[4];
  double dz = Params[5];
  // pre-calculate squares
  double dot1_2 = dot1*dot1;
  double dot2_2 = dot2*dot2;
  double dot3_2 = dot3*dot3;

  // tau = sum(m){amp(l,m)/lambda(l,m)}
  // See Korzhnev DM, Billeter M, Arseniev AS, Orekhov VY; 
  // Prog. Nuc. Mag. Res. Spec., 38, 197 (2001) for details, Table 3 in 
  // particular. Only weights need to be computed for each vector, decay 
  // constants can be computed once. Decay constants correspond 
  // to l=2, m=-2,-1,0,+1,+2.

  // l=2, m=-2 term: lambda(2,-2) = Dx + Dy + 4*Dz
  lambda[0] = dx + dy + (4*dz);
  // l=2, m=-1 term: lambda(2,-1) = Dx + 4*Dy + Dz
  lambda[1] = dx + (4*dy) + dz;
  // l=2, m=0 term:  lambda(2, 0) = 6*[Dav - sqrt(Dav*Dav - Dpr*Dpr)]
  //                 Dav = (Dx + Dy + Dz)/3 
  //                 Dpr = sqrt[(Dx*Dy + Dy*Dz + Dx*Dz)/3 ]
  //                 Dpr2 = Dpr*Dpr = (Dx*Dy + Dy*Dz + Dx*Dz)/3
  double Dav = (dx + dy + dz) / 3;
  double Dpr2 = ((dx*dy) + (dy*dz) + (dx*dz)) / 3;
  if (Dpr2 < 0) {
    //mprinterr("Error: Rotdif::calc_Asymmetric: Cannot calculate Dpr (Dpr^2 < 0)\n");
    // NOTE: Original code just set Dpr to 0 and continued. 
    Dpr2 = 0;
  }
  double delta = (Dav*Dav) - Dpr2;
  // Do not allow square root of negative number. Return a large value to indicate
  // this is a poor fit.
  if (delta < 0) {
    //mprinterr("Error: Ctau_L2: Cannot calculate lambda l=2, m=0, tau=%g\n", tau);
    //mprinterr("\tl= %g  m= %g  n= %g  dx= %g  dy= %g  dz= %g\n",
    //          Params[0], Params[1], Params[2], Params[3], Params[4], Params[5]);
    Ct.assign(Ct.size(), 1000.0);
    return 1;
  }
  double sqrt_Dav_Dpr = sqrt( delta );
  lambda[2] = 6 * (Dav - sqrt_Dav_Dpr);
  // l=2, m=+1 term: lambda(2,+1) = 4*Dx + Dy + Dz
  lambda[3] = (4*dx) + dy + dz;
  // l=2, m=+2 term: lambda(2,+2) = 6*[Dav + sqrt(Dav*Dav - Dpr*Dpr)]
  lambda[4] = 6 * (Dav + sqrt_Dav_Dpr);

  // Check all normalization factors for 0.0, which can lead to inf values
  // propogating throughout the calc. Set to SMALL instead; will get 
  // very large numbers but should still not overflow.
//  mprintf("lambda=");
  for (int i = 0; i < 5; i++) {
    if (lambda[i] < Constants::SMALL) lambda[i] = Constants::SMALL;
//    mprintf(" %g", lambda[i]);
  }
//  mprintf("\n");

  //mprintf("DBG: dot1-3= %10.5g %10.5g %10.5g\n",dot1,dot2,dot3);
  
  // ----- Compute correlation time for l=2
  //     m=-2 term:
  //     lambda(2,-2) = Dx + Dy + 4*Dz
  //     amp(2,-2) = 0.75*sin(theta)^4*sin(2*phi)^2 = 3*(l^2)*(m^2)
  double m_m2 = 3*dot1_2*dot2_2;
  //     m=-1 term:
  //     lambda(2,-1) = Dx + 4*Dy + Dz
  //     amp(2,-1) = 0.75*sin(2*theta)^2*cos(phi)^2 = 3*(l^2)*(n^2)
  double m_m1 = 3*dot1_2*dot3_2;
  //     m=0 term:
  //     lambda(2,0) = 6*[Dav - sqrt(Dav*Dav - Dpr*Dpr)]
  //         Dav = (Dx + Dy + Dz)/3 
  //         Dpr = sqrt[(Dx*Dy + Dy*Dz + Dx*Dz)/3 ]
  //     amp(2,0) = (w/N)^2*0.25*[3*cos(theta)^2-1]^2 +
  //                (u/N)^2*0.75*sin(theta)^4*cos(2*phi)^2 -
  //                (u/delta)*[sqrt(3)/8]*[3*cos(theta)^2-1]*sin(theta)^2*cos(2*phi)
  //         u = sqrt(3)*(Dx - Dy)
  //         delta = 3*sqrt(Dav*Dav - Dpr*Dpr)
  //         w = 2*Dz - Dx - Dy + 2*delta
  //         N = 2*sqrt(delta*w)
  delta = 3 * sqrt_Dav_Dpr;
  double dot1_4 = dot1_2 * dot1_2;
  double dot2_4 = dot2_2 * dot2_2;
  double dot3_4 = dot3_2 * dot3_2;
  double da = 0.25 * ( 3*(dot1_4 + dot2_4 + dot3_4) - 1);
  double ea = 0;
  if ( delta > Constants::SMALL) { 
    double epsx = 3*(dx-Dav)/delta;
    double epsy = 3*(dy-Dav)/delta;
    double epsz = 3*(dz-Dav)/delta;
    double d2d3 = dot2*dot3;
    double d1d3 = dot1*dot3;
    double d1d2 = dot1*dot2;
    ea = epsx*(3*dot1_4 + 6*(d2d3*d2d3) - 1) + 
         epsy*(3*dot2_4 + 6*(d1d3*d1d3) - 1) + 
         epsz*(3*dot3_4 + 6*(d1d2*d1d2) - 1);
    ea /= 12;
  }
  double m_0 = da + ea;
  //     m=+1 term:
  //     lambda(2,+1) = 4*Dx + Dy + Dz
  //     amp(2,+1) = 0.75*sin(2*theta)^2*sin(phi)^2 
  double m_p1 = 3*dot2_2*dot3_2;
  //     m=+2 term:
  //     lambda(2,+2) = 6*[Dav + sqrt(Dav*Dav = Dpr*Dpr)]
  //     amp(2,+2) = (u/n)^2*0.25*[3*cos(theta)^2-1]^2 +
  //                 (w/n)^2*0.75*sin(theta)^4*cos(2*phi)^2 +
  //                 (u/delta)*[sqrt(3)/8]*[3*cos(theta)^2-1]*sin(theta)^2*cos(2*phi)
  double m_p2 = da - ea;

  double penalty = 0.0;
  if (USE_PENALTY) {
    // Calculate penalty functions. Require sqrt(l^2 + m^2 + n^2) = 1.0,
    // dx|dy|dz > 0.0, sum of all prefactors = 1.0.
    penalty = 1000.0 * (1.0 - sqrt(dot1_2 + dot2_2 + dot3_2));
    if (dx < 0.0) penalty += 400.0;
    if (dy < 0.0) penalty += 400.0;
    if (dz < 0.0) penalty += 400.0;
    penalty += (1000.0 * (1.0 - (m_m2+m_m1+m_0+m_p1+m_p2)));
  }

  // Calculate Y values
  for (unsigned int n = 0; n != Tau.size(); n++) {
    double tau = Tau[n]; 
    // Sum decay constants and weights for l=2 m=-2...2
    Ct[n]= ( (m_m2 * exp(-tau * lambda[0])) +
             (m_m1 * exp(-tau * lambda[1])) +
             (m_0  * exp(-tau * lambda[2])) +
             (m_p1 * exp(-tau * lambda[3])) +
             (m_p2 * exp(-tau * lambda[4])) ) + penalty;
  }
    //Yvals[nvec] = (m_m2 / lambda[0]) + (m_m1 / lambda[1]) + (m_0 / lambda[2]) +
    //              (m_p1 / lambda[3]) + (m_p2 / lambda[4]);
    //sumc2_[nvec] = m_m2 + m_m1 + m_0 + m_p1 + m_p2;
//    mprintf("Yvals[%i]= (%g / %g) + (%g / %g) + (%g / %g) + (%g / %g) + (%g / %g) = %g\n", nvec,
//            m_m2, lambda[0], m_m1, lambda[1], m_0, lambda[2], m_p1, lambda[3], m_p2, lambda[4], Yvals[nvec]);
  //mprintf("DEBUG: cl={ %g %g %g %g %g }, Tl={ %g %g %g %g %g }\n",
  //        m_m2, m_m1, m_0, m_p1, m_p2, 1.0/lambda[0], 1.0/lambda[1],
  //        1.0/lambda[2], 1.0/lambda[3], 1.0/lambda[4]);
  //mprintf("DEBUG: Sum of prefactors= %g\n", m_m2+m_m1+m_0+m_p1+m_p2);
  //mprintf("DEBUG: Penalty is %g\n", penalty);
  return 0;
}

// Action_Rotdif::DetermineDeffsAlt()
/** For each randomly generated vector:
  *   1) Rotate the vector according to saved rotation matrices.
  *   2) Calculate the autocorrelation function of the vector.
  * Then take the average of the autocorrelation functions and
  *   fit the autocorrelation to a single exponential and multi exp.
  */
int Action_Rotdif::DetermineDeffsAlt() {
  if (olegendre_ != 2) {
    mprintf("Warning: This calculation currently only works for order=2. Setting order to 2.\n");
    olegendre_ = 2;
  }
  // Determine max length of autocorrelation fxn.
   int vLength = (int)Rmatrices_.size() + 1;
   int ctMax = ncorr_;
   if (ncorr_ == 0)
     ctMax = vLength;
   else if (ncorr_ > vLength)
     ctMax = vLength;
  //int ctMax = (int)((tf_ - ti_) / tfac_) + 1;
  mprintf("DEBUG: Npoints for autocorrelation fxn= %i  vLength=%i  ncorr= %i\n",
          ctMax, vLength, ncorr_);

  // Reserve space for holding effective D values
  D_eff_.reserve( random_vectors_.Size() );
  // Hold vectors after rotation with Rmatrices
  DataSet_Vector rotated_vectors;
  rotated_vectors.ReserveVecs( vLength );
  // Hold single vector autocorrelation
  CurveFit::Darray Ct;
  Ct.reserve( ctMax );
  // Hold averaged vector autocorrelations
  CurveFit::Darray CtTotal( ctMax, 0.0 );
  // LOOP OVER RANDOM VECTORS
  for (DataSet_Vector::const_iterator rndvec = random_vectors_.begin();
                                      rndvec != random_vectors_.end(); ++rndvec)
  {
    // Reset rotated_vectors to the beginning and clear spherical harmonics 
    rotated_vectors.reset();
    // Normalize vector
    //rndvec->Normalize(); // FIXME: Should already be normalized
    // Assign normalized vector to rotated_vectors position 0
    rotated_vectors.AddVxyz( *rndvec );
    // Loop over rotation matrices
    for (std::vector<Matrix_3x3>::iterator rmatrix = Rmatrices_.begin();
                                           rmatrix != Rmatrices_.end();
                                           ++rmatrix)
      // Rotate normalized vector
      rotated_vectors.AddVxyz( *rmatrix * (*rndvec) );
    // Calculate spherical harmonics of given order for this vector.
    rotated_vectors.CalcSphericalHarmonics( olegendre_ );
    // Calculate autocorrelation for rotated vectors using SH.
    fft_compute_corr(rotated_vectors, ctMax, Ct);
    // Sum into CtTotal
    for (unsigned int i = 0; i != Ct.size(); i++)
      CtTotal[i] += Ct[i];
  }
  // Average curve over all vectors.
  double norm_nvecs = 1.0 / (double)random_vectors_.Size();
  for (CurveFit::Darray::iterator ctot = CtTotal.begin(); ctot != CtTotal.end(); ++ctot)
    *ctot *= norm_nvecs; 

  // Set up X and Y values for curve fitting.
  CurveFit::Darray Xvals;
  Xvals.reserve( ctMax );
  double xv = ti_;
  for (int n = 0; n != ctMax; n++) {
    Xvals.push_back( xv );
    xv += tfac_;
  }

  // Parameter array for curve fitting (1 for single exp)
  CurveFit::Darray Params(1, 1.0);
  // Parameter array for curve fitting (10 for multi exp)
  //CurveFit::Darray Mparams(10, 1.0);
  // Parameter array for curve fitting (6 for multi exp, l,m,n,x,y,z)
  CurveFit::Darray Cparams(6, 0.5);
  CurveFit fit;

  double stat_corr, stat_chisq, stat_theilu, stat_rpe;
  // Fit averaged autocorrelation to C(tau) = exp[-l(l+1)D * tau] 
  int info = fit.LevenbergMarquardt( ExpFxn, Xvals, CtTotal, Params, amoeba_ftol_, amoeba_itmax_ );
  mprintf("\tSingleExp: %s\n", fit.Message(info));
  if (info == 0) {
    mprinterr("Error: Single exp fit: %s\n", fit.ErrorMessage());
    return 1;
  }
  fit.Statistics(CtTotal, stat_corr, stat_chisq, stat_theilu, stat_rpe);
  CurveFit::Darray Ct_single = fit.FinalY();
  // We now in principal have A0 = -l(l+1)D, D = A0 / -l(l+1)
  double A0 = Params[0];
  double Deff = A0 / (double)(olegendre_ * (olegendre_ + 1));
  outfile_.Printf("# Results from single-exponential fit: <C(t)> = exp(-t * A0)\n");
  outfile_.Printf("# Corr= %g  ChiSq= %g  TheilU= %g  RMS_PE= %g\n",
                  stat_corr, stat_chisq, stat_theilu, stat_rpe);
  outfile_.Printf("%-12s %12s %12s\n", "#A0", "D" ,"T");
  outfile_.Printf("%12.5e %12.5e %12.5e\n", A0, Deff, 1.0 / A0);
//  mprintf("\t\tA0= %12.5e    Deff= %12.5e    Vxyz={%12.5e, %12.5e, %12.5e}\n",
//          A0, Deff, AvgVec[0], AvgVec[1], AvgVec[2]);
//  D_eff_.push_back( Deff );

  // Fit averaged autocorrelation to C(tau) = SUM(-l...l)[ c * exp(-tau / T) ], T = 1/E
//  for (unsigned int i = 0; i < 10; i += 2) {
//    Mparams[i  ] = 0.2;
//    Mparams[i+1] = A0;
//  }
////  CurveFit::Darray weights( CtTotal.size() );
////  for (unsigned int i = 0; i != CtTotal.size(); i++)
////    weights[i] = ( (double)i / (double)(CtTotal.size() - i) );
//  info = fit.LevenbergMarquardt( Sum5Exp, Xvals, CtTotal, Mparams,  
//                                 amoeba_ftol_, amoeba_itmax_ );
  Vec3 AvgVec = random_vectors_.Back();
  double rv_norm = 1.0 / sqrt( AvgVec[0]*AvgVec[0] +
                               AvgVec[1]*AvgVec[1] +
                               AvgVec[2]*AvgVec[2] );
  Cparams[0] = AvgVec[0] * rv_norm; // l
  Cparams[1] = AvgVec[1] * rv_norm; // m
  Cparams[2] = AvgVec[2] * rv_norm; // n
  Cparams[3] = A0;                  // dx (A0?)
  Cparams[4] = A0 + (A0 * 0.1);     // dy
  Cparams[5] = A0 - (A0 * 0.1);     // dz
  USE_PENALTY = true;
  info = fit.LevenbergMarquardt( Ctau_L2, Xvals, CtTotal, Cparams, amoeba_ftol_, amoeba_itmax_ ); 
  mprintf("\tMultiExp: %s\n", fit.Message(info));
  if (info == 0) {
    mprinterr("Error: Multi exp fit: %s\n", fit.ErrorMessage());
    return 1;
  }
  // Calculate final curve without applied penalty
  USE_PENALTY = false;
  fit.LevenbergMarquardt( Ctau_L2, Xvals, CtTotal, Cparams, amoeba_ftol_, 0 );
  fit.Statistics(CtTotal, stat_corr, stat_chisq, stat_theilu, stat_rpe);
//  for (unsigned int i = 0; i < 10; i += 2)
//    mprintf("\t\t# %2i  c= %12.5f  T= %12.5e\n", i/2, Mparams[i], Mparams[i+1]);
  // Sort Dx <= Dy <= Dz
  std::sort( Cparams.begin() + 3, Cparams.end() );
  outfile_.Printf("# Results from multi-exponential fit:"
                  " <C(t)> = SUM(l=-2,...,2)[ c(l) * exp(-t * E2(l)) ]\n");
  outfile_.Printf("# Corr= %g  ChiSq= %g  TheilU= %g  RMS_PE= %g\n",
                  stat_corr, stat_chisq, stat_theilu, stat_rpe);
  outfile_.Printf("%-12s %12s %12s %12s %12s %12s\n", "#l", "m", "n", "dx", "dy", "dz");
  outfile_.Printf("%12.5f %12.5f %12.5f %12.5e %12.5e %12.5e\n",
                  Cparams[0], Cparams[1], Cparams[2],
                  Cparams[3], Cparams[4], Cparams[5]);
//  mprintf("\t  l= %12.5f  m= %12.5f  n= %12.5f  L= %12.5f"
//          "  dx= %12.5e  dy= %12.5e  dz= %12.5e\n",
//          Cparams[0], Cparams[1], Cparams[2], 
//          sqrt(Cparams[0]*Cparams[0] + Cparams[1]*Cparams[1] + Cparams[2]*Cparams[2]),
//          Cparams[3], Cparams[4], Cparams[5]);
  double DX = Cparams[3];
  double DY = Cparams[4];
  double DZ = Cparams[5];
  Vec3 d_props = calculate_D_properties(Vec3(DX, DY, DZ));
  outfile_.Printf("%-12s %12s %12s\n%12.5e %12.5e %12.5e\n", "#Dav",
                  "anisotropy", "rhombicity", d_props[0], d_props[1], d_props[2]);
  double t_2_p1 = 1.0 / (4.0 * DX + DY + DZ);
  double t_2_m1 = 1.0 / (4.0 * DY + DX + DZ);
  double t_2_m2 = 1.0 / (4.0 * DZ + DX + DY);
  double Dav = d_props[0];
  double Dav2 = Dav * Dav;
  double Dp2 = (DX*DY + DY*DZ + DX*DZ) / 3.0;
  double Dm2 = sqrt(Dav2 - Dp2);
  mprintf("DEBUG: Dav2= %12.5e  Dp2= %12.5e  Dm2= %12.5e\n", Dav2, Dp2, Dm2);
  double t_2_p2 = 1.0 / (6.0 * (Dav + Dm2));
  double t_2_0  = 1.0 / (6.0 * (Dav - Dm2));
  double tR = 1.0 / (2.0 * (DX + DY + DZ));
  outfile_.Printf("%-12s %12s %12s %12s %12s %12s\n",
                  "#t2,-2", "t2,-1", "t2,0", "t2,+1", "t2,+2", "tR");
  outfile_.Printf("%12.5e %12.5e %12.5e %12.5e %12.5e %12.5e\n", 
                  t_2_m2, t_2_m1, t_2_0, t_2_p1, t_2_p2, tR);

  // Write out Ct and fit curves
  if (!corrOut_.empty() || debug_ > 3) {
    CpptrajFile outfile;
    std::string namebuffer;
    if (!corrOut_.empty())
      namebuffer = corrOut_;
    else
      namebuffer = "CtFit.dat";
    outfile.OpenWrite(namebuffer);
    outfile.Printf("%-12s %20s %20s %20s\n", "#Time", "<Ct>", "SingleExp", "MultiExp");
    CurveFit::Darray const& Ct_multi = fit.FinalY();
    for (int n = 0; n != ctMax; n++)
      outfile.Printf("%12.6g %20.8e %20.8e %20.8e\n",
                     Xvals[n], CtTotal[n], Ct_single[n], Ct_multi[n]);
    outfile.CloseFile();
  }

  return 0;
}

// =============================================================================
// Action_Rotdif::direct_compute_corr()
/** Given a normalized vector that has been randomly rotated itotframes 
  * times, compute the time correlation functions of the vector.
  * \param rotated_vectors array of vector coords for each frame 
  * \param maxdat Maximum length to compute time correlation functions (units of 'frames')
  * \param pY Will be set with values for correlation function, l=olegendre_
  */
// TODO: Make rotated_vectors const&
int Action_Rotdif::direct_compute_corr(DataSet_Vector const& rotated_vectors, int maxdat,
                                       std::vector<double>& pY)
{
  // Initialize output array 
  pY.assign(maxdat, 0.0);
  int itotframes = rotated_vectors.Size();
  // i loop:  each value of i is a value of delay (correlation function argument)
  for (int i = 0; i < maxdat; i++) {
    int jmax = itotframes - i;
    for (int j = 0; j < jmax; j++) {
      // Dot vector j with vector j+i 
      double dot = rotated_vectors[j] * rotated_vectors[j+i];
      //printf("DBG i=%6i j=%6i k=%i  %f\n",i,j,i+j,dot);
      if (olegendre_ == 2)
        pY[i] += ((1.5*dot*dot) - 0.5);
      else
        pY[i] += dot;
    }
    double one_jmax = (double) jmax;
    one_jmax = 1 / one_jmax;
    pY[i] = pY[i] * one_jmax;
  }
 
  return 0; 
}

// Action_Rotdif::calcEffectiveDiffusionConst()
/** computes effect diffusion constant for a vector using the integral over
  * its correlation function as input. Starting with definition:
  *
  *   6*D=integral[0,inf;C(t)] 
  *
  * C(t) has already been integrated from ti -> tf yielding F(ti,tf).
  * Iteratively solves the equation 
  *
  *   D(i+1)=[exp(6*D(i)*ti)-exp(6*D(i)*tf)]/[6*F(ti,tf)]
  *
  * (numerator obtained by integrating exp(6*D*t) from ti -> tf)
  *
  * Modified so that itsolv now solves
  *
  * F(ti,tf;C(t)]=integral[ti,tf;C(t)]
  *
  * D(i+1)={exp[l*(l+1)*D(i)*ti]-exp[l*(l+1)*D(i)*tf)]}/[l*(l+1)*F(ti,tf)]
  * /param f Integral of Cl(t) from ti to tf
  * /return Effective value of D
  */
double Action_Rotdif::calcEffectiveDiffusionConst(double f ) {
// Class variables used:
//   ti,tf: Integration limits.
//   itmax: Maximum number of iterations in subroutine itsolv.
//   delmin: convergence criterion used in subroutine itsolv;
//           maximum accepted fractional change in successive 
//           iterations
//   d0: initial guess for diffusion constant; accurate estimate not needed
//   olegendre: order of Legendre polynomial in the correlation function <P(l)>
//
// Solves the equation 6*D=[exp(-6*D*ti)-exp(-6*D*tf)]/F(ti,tf) iteratively, 
// by putting 6*D(i+1)=[exp(-6*D(i)*ti)-exp(-6*D(i)*tf)]/F(ti,tf)
// where F(ti,tf) is input (integral[dt*C(t)] from ti->tf).
  double l, d, del, fac, di; 
  int i;

  // Always use d0 as initial guess
  di = d0_;
  l = (double) olegendre_;
  fac = (l*(l+1));
  i=1;
  d = 0;
  del = DBL_MAX;
  while ( i<=itmax_ && del>delmin_) {
     d = ( exp(-fac*di*ti_) - exp(-fac*di*tf_) );
     d = d / (fac*f);
     del = (d-di)/di;
     if (del < 0) del = -del; // Abs value
     if (debug_>2)
       mprintf("ITSOLV: %6i  %15.8g  %15.8g  %15.8g\n", i,di,d,del);
     di = d;
     ++i;
  }
  if ( i>itmax_ && del>delmin_) {
     mprintf("\tWarning, itsolv did not converge: # iterations=%i, fractional change=%lf\n",
             i, del);
  } else {
    if (debug_>1) mprintf("\tITSOLV Converged: # iterations=%i\n",i);
  }
  //mprintf("DEBUG: Final D= %g\n", d);
  return d; 
}

// Action_Rotdif::DetermineDeffs()
/** Calculate effective diffusion constant for each random vector. 
  * Vectors will be normalized during this phase. First rotate the vector 
  * by all rotation matrices, storing the resulting vectors. The first entry 
  * of rotated_vectors is the original random vector. Then calculate the time 
  * correlation function of the vector. Finally compute the area under the 
  * time correlation function curve and estimate the diffusion constant.
  * Sets D_Eff, normalizes random_vectors.
  */
// TODO: OpenMP Parallelize
int Action_Rotdif::DetermineDeffs() {
  int itotframes;                 // Total number of frames (rotation matrices) 
  DataSet_Vector rotated_vectors; // Hold vectors after rotation with Rmatrices
  int maxdat;                     // Length of C(t) 
  std::vector<double> pX;         // Hold X values of C(t)
  std::vector<double> pY;         // Hold Y values of C(t) for p(olegendre_)
  int meshSize;                   // Total mesh size, maxdat * NmeshPoints

  mprintf("\tDetermining local diffusion constants for each vector.\n");
  ProgressBar progress( nvecs_ );

  itotframes = (int) Rmatrices_.size();
  if (ncorr_ == 0) ncorr_ = itotframes;
  maxdat = ncorr_ + 1;
  // Allocate memory to hold calcd effective D values
  D_eff_.reserve( nvecs_ );
  // Allocate memory to hold rotated vectors. Need +1 since the original
  // vector is stored at position 0. 
  rotated_vectors.ReserveVecs( itotframes + 1 );
  // Allocate memory for C(t)
  pY.reserve( maxdat );
  pX.reserve( maxdat );
  // Set X values of C(t) based on tfac
  for (int i = 0; i < maxdat; i++)
    pX.push_back( (double)i * tfac_ );
  // Allocate mesh to hold interpolated C(t)
  if (NmeshPoints_ < 1)
    meshSize = maxdat * 2;
  else
    meshSize = maxdat * NmeshPoints_;
  // Cubic splines will be used to interpolate C(t) for smoother integration.
  DataSet_Mesh spline( meshSize, ti_, tf_ );
  // LOOP OVER RANDOM VECTORS
  int nvec = 0;
  for (DataSet_Vector::const_iterator rndvec = random_vectors_.begin();
                                      rndvec != random_vectors_.end();
                                    ++rndvec, ++nvec)
  {
    progress.Update( nvec );
    // Reset rotated_vectors to the beginning 
    rotated_vectors.reset();
    // Normalize vector
    //rndvec->Normalize(); // FIXME: Should already be normalized
    // Assign normalized vector to rotated_vectors position 0
    rotated_vectors.AddVxyz( *rndvec );
    // Loop over rotation matrices
    for (std::vector<Matrix_3x3>::iterator rmatrix = Rmatrices_.begin();
                                           rmatrix != Rmatrices_.end();
                                           ++rmatrix)
    {
      // Rotate normalized vector
      rotated_vectors.AddVxyz( *rmatrix * (*rndvec) );
      // DEBUG
      //Vec3 current = rotated_vectors.CurrentVec();
      //mprintf("DBG:Rotated %6u: %15.8f%15.8f%15.8f\n", rmatrix - Rmatrices_.begin(),
      //        current[0], current[1], current[2]); 
    }
    // Calculate time correlation function for this vector
//    if (usefft_) {
//      // Calculate spherical harmonics for each vector.
//      rotated_vectors.CalcSphericalHarmonics( olegendre_ );
//      fft_compute_corr(rotated_vectors, maxdat, pY);
//    } else
      direct_compute_corr(rotated_vectors, maxdat, pY);
    // Calculate mesh Y values
    spline.SetSplinedMeshY(pX, pY);
    // Integrate
    double integral = spline.Integrate_Trapezoid();
    //mprintf("DEBUG: Vec %i integral= %g\n", nvec, integral);
    // Solve for deff
    D_eff_.push_back( calcEffectiveDiffusionConst(integral) );

    // DEBUG: Write out p1 and p2 ------------------------------------
    if (!corrOut_.empty() || debug_ > 3) {
      CpptrajFile outfile;
      std::string namebuffer;
      if (!corrOut_.empty())
        namebuffer = NumberFilename( corrOut_, nvec );
      else
        namebuffer = NumberFilename( "p1p2.dat", nvec );
      outfile.OpenWrite(namebuffer);
      for (int i = 0; i < maxdat; i++) 
        //outfile.Printf("%lf %lf %lf\n",pX[i], p2[i], p1[i]);
        outfile.Printf("%12.6g %20.8e\n",pX[i], pY[i]);
      outfile.CloseFile();
      //    Write Mesh
      if (debug_>3) {
        namebuffer = NumberFilename( "mesh.dat", nvec );
        outfile.OpenWrite(namebuffer);
        for (int i=0; i < (int)spline.Size(); i++)
          outfile.Printf("%12.6g %20.8e\n", spline.X(i), spline.Y(i));
        outfile.CloseFile();
      }
    }
    if (debug_ > 1) {
      mprintf("DBG: Vec %i Spline integral= %12.4g\n",nvec,integral);
      mprintf("DBG: deff is %g\n",D_eff_[nvec]);
    }
    // END DEBUG -----------------------------------------------------
  }

  return 0;
}

// =============================================================================
void Action_Rotdif::PrintDeffs(std::string const& nameIn) const {
  // Print deffs
  if (!nameIn.empty()) {
    CpptrajFile dout;
    if (dout.SetupWrite(nameIn, debug_)) {
      mprinterr("Error: Could not set up Deff file %s\n",nameIn.c_str());
    } else {
      dout.OpenFile();
      for (int vec = 0; vec < nvecs_; vec++)
        dout.Printf("%6i %15.8e\n", vec+1, D_eff_[vec]);
      dout.CloseFile();
    }
  }
}

// Action_Rotdif::Print()
/** Main tensorfit calculation.
  * - Read/generate random vectors.
  * - For each random vector:
  *   a) Normalize it.
  *   b) Rotate it by each of the rotation matrices.
  *   c) Calculate the time correlation function of the rotated vector.
  *   d) Integrate the time correlation function to obtain F.
  *   e) Given f, iteratively solve eq 18  in Wong&Case 2008 for
  *      effective value of diffusion constant (deff) for that vector
  * - Given deff for each vector, solve for Q assuming small anisotropy
  *   via SVD using eq 13 from Wong&Case 2008
  * - Based on Q from small anisotropic limit, use downhill simplex
  *   minimizer to optimize Q in full anisotropic limit
  */
void Action_Rotdif::Print() {
  mprintf("    ROTDIF:\n");
  // Read/Generate nvecs random vectors
  random_vectors_ = RandomVectors();
  if (random_vectors_.Size() < 1) return;
  // ---------------------------------------------
  // If no rotation matrices generated, exit
  if (Rmatrices_.empty()) return;
  // HACK: To match results from rmscorr.f (where rotation matrices are
  //       implicitly transposed), transpose each rotation matrix.
  // NOTE: Is this actually correct? Want inverse rotation?
  for (std::vector<Matrix_3x3>::iterator rmatrix = Rmatrices_.begin();
                                         rmatrix != Rmatrices_.end(); ++rmatrix)
    rmatrix->Transpose();
  // Print rotation matrices
  if (!rmOut_.empty()) {
    CpptrajFile rmout;
    if (rmout.SetupWrite(rmOut_,debug_)) {
      mprinterr("Error: Could not set up %s for writing rotation matrices.\n",rmOut_.c_str());
    } else {
      rmout.OpenFile();
      int rmframe=1;
      for (std::vector<Matrix_3x3>::const_iterator rmatrix = Rmatrices_.begin();
                                                   rmatrix != Rmatrices_.end();
                                                 ++rmatrix, ++rmframe) 
      {
        rmout.Printf("%13i %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f %12.9f\n",
            rmframe,
            (*rmatrix)[0], (*rmatrix)[1], (*rmatrix)[2],
            (*rmatrix)[3], (*rmatrix)[4], (*rmatrix)[5],
            (*rmatrix)[6], (*rmatrix)[7], (*rmatrix)[8]);
      }
      rmout.CloseFile();
    }
  }
  mprintf("\t%i vectors, %u rotation matrices.\n",nvecs_,Rmatrices_.size());
  if (usefft_) {
    // ---------------------------------------------
    // Test calculation; determine constants directly with SH and curve fitting.
    DetermineDeffsAlt();
    //PrintDeffs( deffOut_ );
  } else {
    // ---------------------------------------------
    // Original Wong & Case method using Legendre polynomial Ct and SVD.
    // Determine effective D for each vector.
    DetermineDeffs( );
    PrintDeffs( deffOut_ );
    // All remaining functions require LAPACK
#   ifndef NO_MATHLIB
    SimplexMin::Darray Q_isotropic(6);
    if (Tensor_Fit( Q_isotropic )) return;
    // Using Q (small anisotropy) as a guess, calculate Q with full anisotropy
    mprintf("\tDetermining diffusion tensor with full anisotropy.\n");
    SimplexMin::Darray Q_anisotropic = Q_isotropic;
    // Set up simplex minimizer and function to use.
    SimplexMin minimizer;
    SimplexMin::SimplexFunctionType fxn;
    if (olegendre_ == 1)
      fxn = AsymmetricFxn_L1;
    else
      fxn = AsymmetricFxn_L2;
    std::vector<double> Tau(nvecs_, 0.0);
    // First, back-calculate with the SVD tensor, but with the full anisotropy
    // chi_squared performs diagonalization. The workspace for dsyev should
    // already have been set up in Tensor_Fit.
    double initial_chisq = chi_squared(fxn, Q_anisotropic, &random_vectors_, D_eff_, Tau);
    outfile_.Printf("\nSame diffusion tensor, but full anisotropy:\n");
    outfile_.Printf("  chi_squared for SVD tensor is %15.5g\n", initial_chisq);
    PrintTau( Tau );
    // Use the SVD solution as one of the initial points for the minimizer.
    minimizer.Minimize(fxn, Q_anisotropic, &random_vectors_, D_eff_, delqfrac_,
                       amoeba_itmax_, amoeba_ftol_, amoeba_nsearch_, RNgen_);
    outfile_.Printf("\nOutput from amoeba:\n");
    // Print Q vector
    PrintVec6(outfile_,"Qxx Qyy Qzz Qxy Qyz Qxz", Q_anisotropic);
    Q_to_D( Q_anisotropic, D_tensor_ );
    Diagonalize( D_tensor_, D_XYZ_ );
    // Copy final Tau values from minimizer in case grid search needed.
    Tau = minimizer.FinalY();
    //mprintf("Output from amoeba - average at cycle %i\n",i+1);
    Vec3 d_props = calculate_D_properties(D_XYZ_);
    outfile_.Printf("    Final chisq = %15.5g\n", minimizer.FinalChiSquared());
    PrintVector(outfile_,"Dav, aniostropy, rhombicity:",d_props);
    PrintVector(outfile_,"D tensor eigenvalues:",D_XYZ_);
    PrintMatrix(outfile_,"D tensor eigenvectors (in columns):",D_tensor_);
    PrintTau( Tau );
    // Brute force grid search
    if (do_gridsearch_) {
      // Store initial solution
      SimplexMin::Darray best = Q_anisotropic;
      SimplexMin::Darray xsearch(6);
      bool success = false;
      const int gridsize = 5;
      int gridmax = gridsize + 1;
      int gridmin = -gridsize;
      // Calc initial chi squared
      double sgn0 = chi_squared(fxn, Q_anisotropic, &random_vectors_, D_eff_, Tau);
      mprintf("Grid search: Starting chisq is %15.5g\n",sgn0);
      //mprintf("_Grid q0 %10.5g\n",Q_anisotropic[0]);
      // Loop over all grid points
      ProgressBar progress( gridmax );
      for (int i = gridmin; i < gridmax; i++) {
        progress.Update( i );
        xsearch[0] = Q_anisotropic[0] + (i*delqfrac_/100.0);
        for (int j = gridmin; j < gridmax; j++) {
          xsearch[1] = Q_anisotropic[1] + (j*delqfrac_/100.0);
          for (int k = gridmin; k < gridmax; k++) {
            xsearch[2] = Q_anisotropic[2] + (k*delqfrac_/100.0);
            for (int l = gridmin; l < gridmax; l++) {
              xsearch[3] = Q_anisotropic[3] + (l*delqfrac_/100.0);
              for (int m = gridmin; m < gridmax; m++) {
                xsearch[4] = Q_anisotropic[4] + (m*delqfrac_/100.0);
                for (int n = gridmin; n < gridmax; n++) {
                  xsearch[5] = Q_anisotropic[5] + (n*delqfrac_/100.0);
  
                  double sgn = chi_squared(fxn, xsearch, &random_vectors_, D_eff_, Tau);
                  //mprintf("_Grid x0 %10.5g\n",xsearch[0]);
                  //mprintf("_Grid %3i%3i%3i%3i%3i%3i%10.5g%10.5g\n",
                  //        i,j,k,l,m,n,sgn,sgn0);
                  if (sgn < sgn0) {
                    best = xsearch;
                    sgn0 = sgn;
                    success = true;
                  }
                } // n
              } // m
            } // l
          } // k
        } // j
      } // i
      // Set best solution
      if (success) {
        mprintf("  Grid search succeeded.\n");
        Q_anisotropic = best;
        PrintVec6(outfile_,"Qxx Qyy Qzz Qxy Qyz Qxz", Q_anisotropic);
      } else
        mprintf("  Grid search could not find a better solution.\n");
    } // END grid search
#   endif
  }
  outfile_.CloseFile();
}
