#include <cmath>
#include "ClusterDist.h"
#include "Constants.h" // RADDEG, DEGRAD
#include "ProgressBar.h"
#ifdef _OPENMP
#  include "omp.h"
#endif
// TODO: All DataSet stuff const&
/// Calculate smallest difference between two angles (in degrees).
static double DistCalc_Dih(double d1, double d2) {
  double diff = fabs(d1 - d2);
  if (diff > 180.0)
    return (360.0 - diff);
  else
    return diff;
}

/// Calculate basic difference.
static double DistCalc_Std(double d1, double d2) {
  return fabs(d1 - d2);
}

/* Update centroid value for adding/removing a frame.
 * \param fval value of frame being added/removed.
 * \param cval current centroid value.
 * \param isTorsion data is periodic.
 * \param oldSize Previous size of the centroid.
 * \param OP Operation being performed.
 */
static double DistCalc_FrameCentroid(double fval, double cval, bool isTorsion,
                                     double oldSize, ClusterDist::CentOpType OP,
                                     double& sumx, double& sumy)
{
  double newcval;
  if (isTorsion) {
    double ftheta = fval * Constants::DEGRAD;
    if (OP == ClusterDist::ADDFRAME) {
      sumy += sin( ftheta );
      sumx += cos( ftheta );
    } else { // SUBTRACTFRAME
      sumy -= sin( ftheta );
      sumx -= cos( ftheta );
    }
    newcval = atan2(sumy, sumx) * Constants::RADDEG;
  } else {
    newcval = cval * oldSize;
    if (OP == ClusterDist::ADDFRAME) {
      newcval += fval;
      newcval /= ( oldSize + 1 );
    } else { // SUBTRACTFRAME
      newcval -= fval;
      newcval /= ( oldSize - 1 );
    }
  }
  return newcval;
}

// -----------------------------------------------------------------------------
/** Calculate unambiguous average dihedral angle (in degrees) by converting to 
  * cartesian coords using x = cos(theta), y = sin(theta), and:
  *   tan(avgtheta) = avgy / avgx = SUM[sin(theta)] / SUM[cos(theta)]
  * See Eq. 2 from Altis et al., J. Chem. Phys., 126 p. 244111 (2007).
  */
static double AvgCalc_Dih( DataSet_1D const& dsIn, ClusterDist::Cframes const& cframesIn,
                           double& sumx, double& sumy ) {
  sumy = 0.0;
  sumx = 0.0;
  // TODO: Convert angles to radians prior to this call?
  for (ClusterDist::Cframes_it frm = cframesIn.begin(); frm != cframesIn.end(); ++frm) {
    double theta = dsIn.Dval( *frm ) * Constants::DEGRAD;
    sumy += sin( theta );
    sumx += cos( theta );
  }
  return atan2(sumy, sumx) * Constants::RADDEG; 
}

static double AvgCalc_Std( DataSet_1D const& dsIn, ClusterDist::Cframes const& cframesIn ) {
  double val = 0.0;
  for (ClusterDist::Cframes_it frm = cframesIn.begin(); frm != cframesIn.end(); ++frm)
    val += dsIn.Dval( *frm );
  return (val / (double)cframesIn.size());
}

// ---------- Distance calc routines for single DataSet ------------------------
ClusterDist_Num::ClusterDist_Num( DataSet* dsIn ) :
  data_((DataSet_1D*)dsIn)
{
  if ( data_->IsTorsionArray() )
    dcalc_ = DistCalc_Dih;
  else
    dcalc_ = DistCalc_Std;
}

void ClusterDist_Num::PairwiseDist(ClusterMatrix& frameDistances, 
                                   ClusterSieve::SievedFrames const& frames)
{
  int f1, f2;
  int f2end = (int)frames.size();
  int f1end = f2end - 1;
#ifdef _OPENMP
#pragma omp parallel private(f1, f2)
{
#pragma omp for schedule(dynamic)
#endif
  for (f1 = 0; f1 < f1end; f1++) {
    for (f2 = f1 + 1; f2 < f2end; f2++)
      frameDistances.SetElement( f1, f2, dcalc_(data_->Dval(frames[f1]), 
                                                data_->Dval(frames[f2])) );
  }
#ifdef _OPENMP
}
#endif
}

double ClusterDist_Num::FrameDist(int f1, int f2) {
  return dcalc_(data_->Dval(f1), data_->Dval(f2));
}

double ClusterDist_Num::CentroidDist(Centroid* c1, Centroid* c2) {
  return dcalc_(((Centroid_Num*)c1)->cval_, ((Centroid_Num*)c2)->cval_);
}

double ClusterDist_Num::FrameCentroidDist(int f1, Centroid* c1) {
  return dcalc_(data_->Dval(f1), ((Centroid_Num*)c1)->cval_);
}

/** Calculate avg value of given frames. */
void ClusterDist_Num::CalculateCentroid(Centroid* centIn, Cframes const& cframesIn) {
  Centroid_Num* cent = (Centroid_Num*)centIn;
  if (data_->IsTorsionArray())
    cent->cval_ = AvgCalc_Dih(*data_, cframesIn, cent->sumx_, cent->sumy_);
  else
    cent->cval_ = AvgCalc_Std(*data_, cframesIn);
}

/** \return A new centroid of the given frames. */
Centroid* ClusterDist_Num::NewCentroid( Cframes const& cframesIn ) {
  Centroid_Num* cent = new Centroid_Num();
  CalculateCentroid( cent, cframesIn );
  return cent;
}

void ClusterDist_Num::FrameOpCentroid(int frame, Centroid* centIn, double oldSize,
                                      CentOpType OP)
{
  Centroid_Num* cent = (Centroid_Num*)centIn;
  cent->cval_ = DistCalc_FrameCentroid(data_->Dval(frame), cent->cval_,
                                       data_->IsTorsionArray(), oldSize, OP,
                                       cent->sumx_, cent->sumy_);
}
 
// ---------- Distance calc routines for multiple DataSets (Euclid) ------------
ClusterDist_Euclid::ClusterDist_Euclid(DsArray const& dsIn)
{
  for (DsArray::const_iterator ds = dsIn.begin(); ds != dsIn.end(); ++ds) {
    dsets_.push_back( (DataSet_1D*)*ds );
    if ( dsets_.back()->IsTorsionArray() )
      dcalcs_.push_back( DistCalc_Dih );
    else
      dcalcs_.push_back( DistCalc_Std );
  }
}

void ClusterDist_Euclid::PairwiseDist(ClusterMatrix& frameDistances,
                                      ClusterSieve::SievedFrames const& frames)
{
  int f1, f2;
  double dist, diff;
  DcArray::iterator dcalc;
  D1Array::iterator ds;
  int f2end = (int)frames.size();
  int f1end = f2end - 1;
#ifdef _OPENMP
#pragma omp parallel private(f1, f2, dist, diff, dcalc, ds)
{
#pragma omp for schedule(dynamic)
#endif
  for (f1 = 0; f1 < f1end; f1++) {
    for (f2 = f1 + 1; f2 < f2end; f2++) {
      dist = 0.0;
      dcalc = dcalcs_.begin();
      for (ds = dsets_.begin(); ds != dsets_.end(); ++ds, ++dcalc) {
        diff = (*dcalc)((*ds)->Dval(frames[f1]), (*ds)->Dval(frames[f2]));
        dist += (diff * diff);
      }
      frameDistances.SetElement( f1, f2, sqrt(dist) );
    }
  }
#ifdef _OPENMP
}
#endif
}

double ClusterDist_Euclid::FrameDist(int f1, int f2) {
  double dist = 0.0;
  DcArray::iterator dcalc = dcalcs_.begin();
  for (D1Array::iterator ds = dsets_.begin(); ds != dsets_.end(); ++ds, ++dcalc) {
    double diff = (*dcalc)((*ds)->Dval(f1), (*ds)->Dval(f2));
    dist += (diff * diff);
  }
  return sqrt(dist);
}

double ClusterDist_Euclid::CentroidDist(Centroid* c1, Centroid* c2) {
  double dist = 0.0;
  std::vector<double>::iterator c2val = ((Centroid_Multi*)c2)->cvals_.begin();
  DcArray::iterator dcalc = dcalcs_.begin();
  for (std::vector<double>::iterator c1val = ((Centroid_Multi*)c1)->cvals_.begin();
                                     c1val != ((Centroid_Multi*)c1)->cvals_.end();
                                     ++c1val, ++dcalc)
  {
    double diff = (*dcalc)(*c1val, *(c2val++));
    dist += (diff * diff);
  }
  return sqrt(dist);
}

double ClusterDist_Euclid::FrameCentroidDist(int f1, Centroid* c1) {
  double dist = 0.0;
  std::vector<double>::iterator c1val = ((Centroid_Multi*)c1)->cvals_.begin();
  DcArray::iterator dcalc = dcalcs_.begin();
  for (D1Array::iterator ds = dsets_.begin(); ds != dsets_.end(); ++ds) {
    double diff = (*dcalc)((*ds)->Dval(f1), *(c1val++));
    dist += (diff * diff);
  }
  return sqrt(dist);
}

void ClusterDist_Euclid::CalculateCentroid(Centroid* centIn, Cframes const& cframesIn) {
  Centroid_Multi* cent = (Centroid_Multi*)centIn;
  cent->cvals_.resize( dsets_.size(), 0.0 );
  cent->Sumx_.resize( dsets_.size(), 0.0 );
  cent->Sumy_.resize( dsets_.size(), 0.0 );
  for (unsigned int idx = 0; idx != dsets_.size(); ++idx) {
    if (dsets_[idx]->IsTorsionArray())
      cent->cvals_[idx] = AvgCalc_Dih(*dsets_[idx], cframesIn, cent->Sumx_[idx], cent->Sumy_[idx]);
    else
      cent->cvals_[idx] = AvgCalc_Std(*dsets_[idx], cframesIn);
  }
//  mprintf("DEBUG: Centroids:");
//  for (unsigned int i = 0; i != cent->cvals_.size(); i++)
//    mprintf("   %f (sumy=%f sumx=%f)", cent->cvals_[i], cent->Sumy_[i], cent->Sumx_[i]);
//  mprintf("\n");
}

Centroid* ClusterDist_Euclid::NewCentroid(Cframes const& cframesIn) {
  Centroid_Multi* cent = new Centroid_Multi();
  CalculateCentroid(cent, cframesIn);
  return cent;
}

//static const char* OPSTRING[] = {"ADD", "SUBTRACT"}; // DEBUG

void ClusterDist_Euclid::FrameOpCentroid(int frame, Centroid* centIn, double oldSize,
                                         CentOpType OP)
{
  Centroid_Multi* cent = (Centroid_Multi*)centIn;
//  mprintf("DEBUG: Old Centroids:");
//  for (unsigned int i = 0; i != cent->cvals_.size(); i++)
//    mprintf("   sumy=%f sumx=%f", cent->Sumy_[i], cent->Sumx_[i]);
//    //mprintf(" %f", cent->cvals_[i]);
//  mprintf("\n");
  for (unsigned int i = 0; i != dsets_.size(); ++i)
    cent->cvals_[i] = DistCalc_FrameCentroid(dsets_[i]->Dval(frame), 
                          cent->cvals_[i], dsets_[i]->IsTorsionArray(), oldSize, OP,
                          cent->Sumx_[i], cent->Sumy_[i]);
//  mprintf("DEBUG: New Centroids after %s frame %i:", OPSTRING[OP], frame);
//  for (unsigned int i = 0; i != cent->cvals_.size(); i++)
//    mprintf(" %f", cent->cvals_[i]);
//  mprintf("\n");
}

// ---------- Distance calc routines for COORDS DataSet using DME --------------
ClusterDist_DME::ClusterDist_DME(DataSet* dIn, AtomMask const& maskIn) :
  coords_((DataSet_Coords*)dIn),
  mask_(maskIn)
{
  frm1_.SetupFrameFromMask(mask_, coords_->Top().Atoms());
  frm2_ = frm1_;
}

void ClusterDist_DME::PairwiseDist(ClusterMatrix& frameDistances,
                                   ClusterSieve::SievedFrames const& frames)
{
  int f1, f2;
  Frame frm2 = frm1_;
  int f2end = (int)frames.size();
  int f1end = f2end - 1;
#ifdef _OPENMP
  Frame frm1 = frm1_;
# define frm1_ frm1
#pragma omp parallel private(f1, f2) firstprivate(frm1, frm2)
{
#pragma omp for schedule(dynamic)
#endif
  for (f1 = 0; f1 < f1end; f1++) { 
    coords_->GetFrame( frames[f1], frm1_, mask_ );
    for (f2 = f1 + 1; f2 < f2end; f2++) {
      coords_->GetFrame( frames[f2], frm2,  mask_ );
      frameDistances.SetElement( f1, f2, frm1_.DISTRMSD( frm2 ) );
    }
  }
#ifdef _OPENMP
# undef frm1_
} // END pragma omp parallel
#endif
}

double ClusterDist_DME::FrameDist(int f1, int f2) {
  coords_->GetFrame( f1, frm1_, mask_ );
  coords_->GetFrame( f2, frm2_, mask_ );
  return frm1_.DISTRMSD( frm2_ );
}

double ClusterDist_DME::CentroidDist(Centroid* c1, Centroid* c2) {
  return ((Centroid_Coord*)c1)->cframe_.DISTRMSD( ((Centroid_Coord*)c2)->cframe_ );
}

double ClusterDist_DME::FrameCentroidDist(int f1, Centroid* c1) {
  coords_->GetFrame( f1, frm1_, mask_ );
  return frm1_.DISTRMSD( ((Centroid_Coord*)c1)->cframe_ );
}

/** Compute the centroid (avg) coords for each atom from all frames in this
  * cluster. NOTE: For DME the centroid should probably be calculated via
  * internal coordinates; use RMS best-fit as a cheat.
  */
void ClusterDist_DME::CalculateCentroid(Centroid* centIn, Cframes const& cframesIn) {
  Matrix_3x3 Rot;
  Vec3 Trans;
  Centroid_Coord* cent = (Centroid_Coord*)centIn;
  // Reset atom count for centroid.
  cent->cframe_.ClearAtoms();
  for (Cframes_it frm = cframesIn.begin(); frm != cframesIn.end(); ++frm)
  {
    coords_->GetFrame( *frm, frm1_, mask_ );
    if (cent->cframe_.empty()) {
      cent->cframe_ = frm1_;
      cent->cframe_.CenterOnOrigin(false);
    } else {
      frm1_.RMSD_CenteredRef( cent->cframe_, Rot, Trans, false );
      frm1_.Rotate( Rot );
      cent->cframe_ += frm1_;
    }
  }
  cent->cframe_.Divide( (double)cframesIn.size() );
  //mprintf("\t\tFirst 3 centroid coords (of %i): %f %f %f\n", cent->cframe_.Natom(), 
  //        cent->cent->cframe_[0], cent->cframe_[1],cent->cframe_[2]);
}

Centroid* ClusterDist_DME::NewCentroid( Cframes const& cframesIn ) {
  // TODO: Incorporate mass?
  Centroid_Coord* cent = new Centroid_Coord( mask_.Nselected() );
  CalculateCentroid( cent, cframesIn );
  return cent;
}

void ClusterDist_DME::FrameOpCentroid(int frame, Centroid* centIn, double oldSize,
                                      CentOpType OP)
{
  Matrix_3x3 Rot;
  Vec3 Trans;
  Centroid_Coord* cent = (Centroid_Coord*)centIn;
  coords_->GetFrame( frame, frm1_, mask_ );
  frm1_.RMSD_CenteredRef( cent->cframe_, Rot, Trans, false );
  frm1_.Rotate( Rot );
  cent->cframe_.Multiply( oldSize );
  if (OP == ADDFRAME) {
    cent->cframe_ += frm1_;
    cent->cframe_.Divide( oldSize + 1 );
  } else { // SUBTRACTFRAME
    cent->cframe_ -= frm1_;
    cent->cframe_.Divide( oldSize - 1 );
  }
}

// ---------- Distance calc routines for COORDS DataSets using RMSD ------------
ClusterDist_RMS::ClusterDist_RMS(DataSet* dIn, AtomMask const& maskIn, 
                                 bool nofit, bool useMass) :
  coords_((DataSet_Coords*)dIn),
  mask_(maskIn),
  nofit_(nofit),
  useMass_(useMass)
{
  frm1_.SetupFrameFromMask(mask_, coords_->Top().Atoms());
  frm2_ = frm1_;
}

void ClusterDist_RMS::PairwiseDist(ClusterMatrix& frameDistances,
                                   ClusterSieve::SievedFrames const& frames)
{
  double rmsd;
  int f1, f2;
  Frame frm2 = frm1_;
  int f2end = (int)frames.size();
  int f1end = f2end - 1;
  ParallelProgress progress(f1end);
#ifdef _OPENMP
  Frame frm1 = frm1_;
# define frm1_ frm1
#pragma omp parallel private(f1, f2, rmsd) firstprivate(frm1, frm2, progress)
{
  progress.SetThread(omp_get_thread_num());
#pragma omp for schedule(dynamic)
#endif
  for (f1 = 0; f1 < f1end; f1++) {
    progress.Update(f1);
    coords_->GetFrame( frames[f1], frm1_, mask_ );
    for (f2 = f1 + 1; f2 < f2end; f2++) {
      coords_->GetFrame( frames[f2], frm2,  mask_ );
      if (nofit_) 
        rmsd = frm1_.RMSD_NoFit( frm2, useMass_ );
      else
        rmsd = frm1_.RMSD( frm2, useMass_ );
      frameDistances.SetElement( f1, f2, rmsd );
    }
  }
#ifdef _OPENMP
# undef frm1_
} // END pragma omp parallel
#endif
  progress.Finish();
}

double ClusterDist_RMS::FrameDist(int f1, int f2) {
  coords_->GetFrame( f1, frm1_, mask_ );
  coords_->GetFrame( f2, frm2_, mask_ );
  if (nofit_)
    return frm1_.RMSD_NoFit( frm2_, useMass_ );
  else
    return frm1_.RMSD( frm2_, useMass_ );
}

double ClusterDist_RMS::CentroidDist(Centroid* c1, Centroid* c2) {
  if (nofit_)
    return ((Centroid_Coord*)c1)->cframe_.RMSD_NoFit( ((Centroid_Coord*)c2)->cframe_, useMass_ );
  else // Centroid is already at origin.
    return ((Centroid_Coord*)c1)->cframe_.RMSD_CenteredRef( ((Centroid_Coord*)c2)->cframe_, 
                                                            useMass_ );
}

double ClusterDist_RMS::FrameCentroidDist(int f1, Centroid* c1) {
  coords_->GetFrame( f1, frm1_, mask_ );
  if (nofit_)
    return frm1_.RMSD_NoFit( ((Centroid_Coord*)c1)->cframe_, useMass_ );
  else // Centroid is already at origin.
    return frm1_.RMSD_CenteredRef( ((Centroid_Coord*)c1)->cframe_, useMass_ );
}

/** Compute the centroid (avg) coords for each atom from all frames in this
  * cluster. If fitting,  RMS fit to centroid as it is being built.
  */
void ClusterDist_RMS::CalculateCentroid(Centroid* centIn,  Cframes const& cframesIn) {
  Matrix_3x3 Rot;
  Vec3 Trans;
  Centroid_Coord* cent = (Centroid_Coord*)centIn;
  // Reset atom count for centroid.
  cent->cframe_.ClearAtoms();
  for (Cframes_it frm = cframesIn.begin(); frm != cframesIn.end(); ++frm)
  {
    coords_->GetFrame( *frm, frm1_, mask_ );
    if (cent->cframe_.empty()) {
      cent->cframe_ = frm1_;
      if (!nofit_)
        cent->cframe_.CenterOnOrigin(useMass_);
    } else {
      if (!nofit_) {
        frm1_.RMSD_CenteredRef( cent->cframe_, Rot, Trans, useMass_ );
        frm1_.Rotate( Rot );
      }
      cent->cframe_ += frm1_;
    }
  }
  cent->cframe_.Divide( (double)cframesIn.size() );
  //mprintf("\t\tFirst 3 centroid coords (of %i): %f %f %f\n", cent->cframe_.Natom(), 
  //        cent->cent->cframe_[0], cent->cframe_[1],cent->cframe_[2]);
}

Centroid* ClusterDist_RMS::NewCentroid( Cframes const& cframesIn ) {
  // TODO: Incorporate mass?
  Centroid_Coord* cent = new Centroid_Coord( mask_.Nselected() );
  CalculateCentroid( cent, cframesIn );
  return cent;
}

// Subtract Notes
// FIXME: Handle single frame
// FIXME: Check if frame is in cluster?
void ClusterDist_RMS::FrameOpCentroid(int frame, Centroid* centIn, double oldSize,
                                      CentOpType OP)
{
  Matrix_3x3 Rot;
  Vec3 Trans;
  Centroid_Coord* cent = (Centroid_Coord*)centIn;
  coords_->GetFrame( frame, frm1_, mask_ );
  if (!nofit_) {
    frm1_.RMSD_CenteredRef( cent->cframe_, Rot, Trans, useMass_ );
    frm1_.Rotate( Rot );
  }
  cent->cframe_.Multiply( oldSize );
  if (OP == ADDFRAME) {
    cent->cframe_ += frm1_;
    cent->cframe_.Divide( oldSize + 1 );
  } else { // SUBTRACTFRAME
    cent->cframe_ -= frm1_;
    cent->cframe_.Divide( oldSize - 1 );
  }
}

// ---------- Distance calc routines for COORDS DataSets using SRMSD -----------
ClusterDist_SRMSD::ClusterDist_SRMSD(DataSet* dIn, AtomMask const& maskIn, 
                                     bool nofit, bool useMass, int debugIn) :
  coords_((DataSet_Coords*)dIn),
  mask_(maskIn),
  SRMSD_(mask_, !nofit, useMass, coords_->Top(), debugIn)
{
  frm1_.SetupFrameFromMask(mask_, coords_->Top().Atoms());
  frm2_ = frm1_;
}

void ClusterDist_SRMSD::PairwiseDist(ClusterMatrix& frameDistances,
                                   ClusterSieve::SievedFrames const& frames)
{
  double rmsd;
  int f1, f2;
  Frame frm2 = frm1_;
  int f2end = (int)frames.size();
  int f1end = f2end - 1;
  ParallelProgress progress(f1end);
#ifdef _OPENMP
  Frame frm1 = frm1_;
# define frm1_ frm1
  SymmetricRmsdCalc SRMSD_OMP = SRMSD_;
# define SRMSD_ SRMSD_OMP
#pragma omp parallel private(f1, f2, rmsd) firstprivate(SRMSD_OMP, frm1, frm2, progress)
{
  progress.SetThread(omp_get_thread_num());
#pragma omp for schedule(dynamic)
#endif
  for (f1 = 0; f1 < f1end; f1++) {
    progress.Update(f1);
    coords_->GetFrame( frames[f1], frm1_, mask_ );
    for (f2 = f1 + 1; f2 < f2end; f2++) {
      coords_->GetFrame( frames[f2], frm2,  mask_ );
      rmsd = SRMSD_.SymmRMSD(frm1_, frm2);
      frameDistances.SetElement( f1, f2, rmsd );
    }
  }
#ifdef _OPENMP
# undef frm1_
# undef SRMSD_
} // END pragma omp parallel
#endif
  progress.Finish();
}

double ClusterDist_SRMSD::FrameDist(int f1, int f2) {
  coords_->GetFrame( f1, frm1_, mask_ );
  coords_->GetFrame( f2, frm2_, mask_ );
  return SRMSD_.SymmRMSD(frm1_, frm2_);
}

double ClusterDist_SRMSD::CentroidDist(Centroid* c1, Centroid* c2) {
  // Centroid is already at origin.
  return SRMSD_.SymmRMSD_CenteredRef( ((Centroid_Coord*)c1)->cframe_,
                                      ((Centroid_Coord*)c2)->cframe_ );
}

double ClusterDist_SRMSD::FrameCentroidDist(int f1, Centroid* c1) {
  coords_->GetFrame( f1, frm1_, mask_ );
  // Centroid is already at origin.
  return SRMSD_.SymmRMSD_CenteredRef( frm1_, ((Centroid_Coord*)c1)->cframe_ );
}

/** Compute the centroid (avg) coords for each atom from all frames in this
  * cluster. If fitting,  RMS fit to centroid as it is being built.
  */
void ClusterDist_SRMSD::CalculateCentroid(Centroid* centIn,  Cframes const& cframesIn) {
  Matrix_3x3 Rot;
  Vec3 Trans;
  Centroid_Coord* cent = (Centroid_Coord*)centIn;
  // Reset atom count for centroid.
  cent->cframe_.ClearAtoms();
  for (Cframes_it frm = cframesIn.begin(); frm != cframesIn.end(); ++frm)
  {
    coords_->GetFrame( *frm, frm1_, mask_ );
    if (cent->cframe_.empty()) {
      cent->cframe_ = frm1_;
      if (SRMSD_.Fit())
        cent->cframe_.CenterOnOrigin(SRMSD_.UseMass());
    } else {
      SRMSD_.SymmRMSD_CenteredRef( frm1_, cent->cframe_ );
      // Remap atoms
      frm2_.SetCoordinatesByMap( frm1_, SRMSD_.AMap() );
      if (SRMSD_.Fit()) {
        frm2_.Translate( SRMSD_.TgtTrans() );
        frm2_.Rotate( SRMSD_.RotMatrix() );
      }
      cent->cframe_ += frm2_;
    }
  }
  cent->cframe_.Divide( (double)cframesIn.size() );
  //mprintf("\t\tFirst 3 centroid coords (of %i): %f %f %f\n", cent->cframe_.Natom(), 
  //        cent->cent->cframe_[0], cent->cframe_[1],cent->cframe_[2]);
}

Centroid* ClusterDist_SRMSD::NewCentroid( Cframes const& cframesIn ) {
  // TODO: Incorporate mass?
  Centroid_Coord* cent = new Centroid_Coord( mask_.Nselected() );
  CalculateCentroid( cent, cframesIn );
  return cent;
}

void ClusterDist_SRMSD::FrameOpCentroid(int frame, Centroid* centIn, double oldSize,
                                        CentOpType OP)
{
  Matrix_3x3 Rot;
  Vec3 Trans;
  Centroid_Coord* cent = (Centroid_Coord*)centIn;
  coords_->GetFrame( frame, frm1_, mask_ );
  SRMSD_.SymmRMSD_CenteredRef( frm1_, cent->cframe_ );
  // Remap atoms
  frm2_.SetCoordinatesByMap( frm1_, SRMSD_.AMap() );
  if (SRMSD_.Fit()) {
    frm2_.Translate( SRMSD_.TgtTrans() );
    frm2_.Rotate( SRMSD_.RotMatrix() );
  }
  cent->cframe_.Multiply( oldSize );
  if (OP == ADDFRAME) {
    cent->cframe_ += frm2_;
    cent->cframe_.Divide( oldSize + 1 );
  } else { // SUBTRACTFRAME
    cent->cframe_ -= frm2_;
    cent->cframe_.Divide( oldSize - 1 );
  }
}
