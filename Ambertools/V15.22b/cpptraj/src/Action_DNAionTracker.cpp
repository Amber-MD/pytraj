#include <cfloat> // DBL_MAX
#include <cmath> // sqrt
#include "Action_DNAionTracker.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_DNAionTracker::Action_DNAionTracker() :
  distance_(0),
  bintype_(COUNT),
  poffset_(0),
  useMass_(true)
{ }

void Action_DNAionTracker::Help() {
  mprintf("\tname mask_p1 mask_p2 mask_base mask_ions\n"
          "\t[poffset <value>] [out <filename>] [time <interval>] [noimage]\n"
          "\t[shortest | counttopcone| countbottomcone | count]\n");
}

Action::RetType Action_DNAionTracker::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  // Get keywords
  std::string filename_ = actionArgs.GetStringKey("out");
  poffset_ = actionArgs.getKeyDouble("poffset", 5.0);
  InitImaging( !actionArgs.hasKey("noimage") );
  if (actionArgs.hasKey("shortest"))
    bintype_ = SHORTEST;
  else if (actionArgs.hasKey("counttopcone"))
    bintype_ = TOPCONE;
  else if (actionArgs.hasKey("countbottomcone"))
    bintype_ = BOTTOMCONE;
  else if (actionArgs.hasKey("count"))
    bintype_ = COUNT;

  // Get masks - 4 must be specified
  std::string m1 = actionArgs.GetMaskNext();
  std::string m2 = actionArgs.GetMaskNext();
  std::string m3 = actionArgs.GetMaskNext();
  std::string m4 = actionArgs.GetMaskNext();
  if (m1.empty() || m2.empty() || m3.empty() || m4.empty()) {
    mprinterr("Error: dnaiontracker requires 4 masks.\n");
    return Action::ERR;
  }
  p1_.SetMaskString(m1);
  p2_.SetMaskString(m2);
  base_.SetMaskString(m3);
  ions_.SetMaskString(m4);

  // Add dataset to dataset list (and datafile list if filename specified)
  distance_ = DSL->AddSet(DataSet::DOUBLE, actionArgs.GetStringNext(), "DNAion");
  // NOTE: Set to mode distance in PTRAJ
  distance_->SetScalar( DataSet::M_DISTANCE );
  if (distance_==0) return Action::ERR;
  if (!filename_.empty())
    DFL->AddSetToFile( filename_, distance_ );

  // INFO
  mprintf("    DNAIONTRACKER: Data representing the ");
  switch (bintype_) {
    case COUNT : 
      mprintf("count within the cone will be\n"); break;
    case SHORTEST: 
      mprintf("shortest distance to a phosphate or base centroid will be\n"); break;
    case TOPCONE: 
      mprintf("count in the top half of the cone (and sort-of bound) will be\n"); break;
    case BOTTOMCONE: 
      mprintf("count in the bottom half of the cone will be\n"); break;
  }
  mprintf("      saved to array named %s\n", distance_->Legend().c_str());
  mprintf("      Perpendicular offset for cone is %5.2f angstroms\n", poffset_);
  if (!UseImage())
    mprintf("      Imaging has been disabled\n");
  mprintf("\tPhosphate1 Mask [%s]\n", p1_.MaskString());
  mprintf("\tPhosphate2 Mask [%s]\n", p2_.MaskString());
  mprintf("\tBase Mask       [%s]\n", base_.MaskString());
  mprintf("\tIons Mask       [%s]\n", ions_.MaskString());
  if (!filename_.empty())
    mprintf("      Data will be printed to a file named %s\n", filename_.c_str());

  return Action::OK;
}

Action::RetType Action_DNAionTracker::Setup(Topology* currentParm, Topology** parmAddress) {
  // Setup masks
  if (currentParm->SetupIntegerMask( p1_ )) return Action::ERR; 
  if ( p1_.None() ) {
    mprinterr("Error: dnaiontracker: No atoms in mask1\n");
    return Action::ERR;
  }
  if (currentParm->SetupIntegerMask( p2_ )) return Action::ERR;  
  if ( p2_.None() ) {
    mprinterr("Error: dnaiontracker: No atoms in mask2\n");
    return Action::ERR;
  }
  if (currentParm->SetupIntegerMask( base_ )) return Action::ERR;  
  if ( base_.None() ) {
    mprinterr("Error: dnaiontracker: No atoms in mask3\n");
    return Action::ERR;
  }
  if (currentParm->SetupIntegerMask( ions_ )) return Action::ERR;  
  if ( ions_.None() ) {
    mprinterr("Error: dnaiontracker: No atoms in mask4\n");
    return Action::ERR;
  }
  SetupImaging( currentParm->BoxType() );
  mprintf("\tPhosphate1 Mask [%s] %i atoms.\n", p1_.MaskString(), p1_.Nselected());
  mprintf("\tPhosphate2 Mask [%s] %i atoms.\n", p2_.MaskString(), p2_.Nselected());
  mprintf("\t      Base Mask [%s] %i atoms.\n", base_.MaskString(), base_.Nselected());
  mprintf("\t      Ions Mask [%s] %i atoms.\n", ions_.MaskString(), ions_.Nselected());

  return Action::OK;
}

Action::RetType Action_DNAionTracker::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {
  Matrix_3x3 ucell, recip;
  double d_tmp, dval;
  Vec3 P1, P2, BASE;
  // Setup imaging info if necessary
  if (ImageType()==NONORTHO) 
    currentFrame->BoxCrd().ToRecip(ucell,recip);

  // Get center for P1, P2, and Base
  if (useMass_) {
    P1 = currentFrame->VCenterOfMass( p1_ );
    P2 = currentFrame->VCenterOfMass( p2_ );
    BASE = currentFrame->VCenterOfMass( base_ );
  } else {
    P1 = currentFrame->VGeometricCenter( p1_ );
    P2 = currentFrame->VGeometricCenter( p2_ );
    BASE = currentFrame->VGeometricCenter( base_ );
  }
 
  // Calculate P -- P distance and centroid
  double d_pp = DIST2(P1.Dptr(), P2.Dptr(), ImageType(), currentFrame->BoxCrd(), ucell, recip);
  Vec3 pp_centroid = (P1 + P2) / 2.0;

  // Cutoff^2
  double d_cut = d_pp*0.25 + (poffset_*poffset_); // TODO: precalc offset^2

  // Calculate P -- base centroid to median point
  double d_pbase = DIST2(pp_centroid.Dptr(), BASE.Dptr(), ImageType(), currentFrame->BoxCrd(), 
                         ucell, recip);

  //double d_min = DBL_MAX;
  if (bintype_ == SHORTEST)
    dval = DBL_MAX; //d_min;
  else
    dval = 0;
  // Loop over ion positions
  for (AtomMask::const_iterator ion = ions_.begin(); ion != ions_.end(); ++ion)
  {
    const double* ionxyz = currentFrame->XYZ(*ion);
    double d_p1ion =   DIST2(P1.Dptr(),   ionxyz, ImageType(), currentFrame->BoxCrd(), 
                             ucell, recip);
    double d_p2ion =   DIST2(P2.Dptr(),   ionxyz, ImageType(), currentFrame->BoxCrd(), 
                             ucell, recip);
    double d_baseion = DIST2(BASE.Dptr(), ionxyz, ImageType(), currentFrame->BoxCrd(), 
                             ucell, recip);
    //mprintf("DEBUG: ion atom %i to P1 is %f\n", *ion+1, sqrt(d_p1ion));
    //mprintf("DEBUG: ion atom %i to P2 is %f\n", *ion+1, sqrt(d_p2ion));
    //mprintf("DEBUG: ion atom %i to base is %f\n", *ion+1, sqrt(d_baseion));
    //mprintf("DEBUG: d_pp is %f, poffset is %f, d_cut is %f\n", sqrt(d_pp), poffset_, sqrt(d_cut));

    int bound = 0;
    int boundLower = 0;
    int boundUpper = 0;

    if (d_p1ion < d_cut && d_p2ion < d_cut)
      bound = 1;
    if (d_baseion < d_pbase)
      boundLower = 1;
    if (bound && boundLower == 0)
      boundUpper = 1;
    //if (d_tmp > d_min)
    //  d_min = d_tmp;

    switch (bintype_) {
      case COUNT: 
        dval += (double)bound; break;
      case SHORTEST:
        if (d_p1ion < d_p2ion)
          d_tmp = d_p1ion;
        else
          d_tmp = d_p2ion;
        if (d_tmp > d_baseion)
          d_tmp = d_baseion;
        if (d_tmp < dval) 
          dval = d_tmp;
        break;
      case TOPCONE: 
        dval += (double)boundUpper; break;
      case BOTTOMCONE: 
        dval += (double)boundLower; break;
    }
  }
  if (bintype_ == SHORTEST)
    dval = sqrt(dval);
  distance_->Add(frameNum, &dval);
  
  return Action::OK;
}
