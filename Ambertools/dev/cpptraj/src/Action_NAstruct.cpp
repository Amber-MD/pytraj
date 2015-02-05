#include <cmath>
#include "Action_NAstruct.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // integerToString
#include "DistRoutines.h"
#include "TorsionRoutines.h"
#include "Constants.h" // RADDEG
#ifdef NASTRUCTDEBUG
#include "PDBfile.h"
#endif

// CONSTRUCTOR
Action_NAstruct::Action_NAstruct() :
  puckerMethod_(ALTONA),
  HBdistCut2_(12.25),  // Hydrogen Bond distance cutoff^2: 3.5^2
  // NOTE: Is this too big?
  originCut2_(6.25),   // Origin cutoff^2 for base-pairing: 2.5^2
  maxResSize_(0),
  debug_(0),
  ensembleNum_(-1),
  printheader_(true),
  useReference_(false),
  masterDSL_(0)
# ifdef NASTRUCTDEBUG
  ,calcparam_(true)
# endif
{}

void Action_NAstruct::Help() {
  mprintf("\t[<dataset name>] [resrange <range>] [naout <suffix>]\n"
          "\t[noheader] [resmap <ResName>:{A,C,G,T,U} ...]\n"
          "\t[hbcut <hbcut>] [origincut <origincut>] [altona | cremer]\n"
          "\t[ %s ]\n", DataSetList::RefArgs);
  mprintf("  Perform nucleic acid structure analysis. Base pairing is determined\n"
          "  from specified reference or first frame.\n"
          "  Base pair parameters are written to BP.<suffix>, base pair step parameters\n"
          "  are written to BPstep.<suffix>, and helix parameters are written to\n"
          "  Helix.<suffix>\n");
}

// DESTRUCTOR
Action_NAstruct::~Action_NAstruct() { 
  ClearLists();
  // NOTE: Since BasePairAxes are set up to correspond with SHEAR etc they 
  //       are only freed at the very end.
  BasePairAxes_.clear();
}

// Output Format Strings
static const char* BP_OUTPUT_FMT = "%8i %8i %8i %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %2i %10.4f %10.4f\n";
static const char* NA_OUTPUT_FMT = "%8i %4i-%-4i %4i-%-4i %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f\n";

// ------------------------- PRIVATE FUNCTIONS ---------------------------------
// Action_NAstruct::ClearLists()
/** Clear all parm-dependent lists */
void Action_NAstruct::ClearLists() {
  Bases_.clear();
  BaseAxes_.clear();
}

// -----------------------------------------------------------------------------
#ifdef NASTRUCTDEBUG
/// Write given NA_Axis to a PDB file.
static void WriteAxes(PDBfile& outfile, int resnum, const char* resname, NA_Axis const& axis)
{
  // Origin
  Vec3 oxyz = axis.Oxyz();
  outfile.WriteATOM("Orig", resnum, oxyz[0], oxyz[1], oxyz[2], resname, 0.0);
  // X vector
  Vec3 vec = axis.Rx() + oxyz;
  outfile.WriteATOM("X", resnum, vec[0], vec[1], vec[2], resname, 0.0);
  // Y vector
  vec = axis.Ry() + oxyz;
  outfile.WriteATOM("Y", resnum, vec[0], vec[1], vec[2], resname, 0.0);
  // Z vector
  vec = axis.Rz() + oxyz;
  outfile.WriteATOM("Z", resnum, vec[0], vec[1], vec[2], resname, 0.0);
}
#endif

// Action_NAstruct::setupBaseAxes()
/** For each residue in Bases (set up in Setup()), get the corresponding input
  * coords and fit the reference coords on top of input coords. This sets up 
  * the reference axes for each base.
  */
int Action_NAstruct::setupBaseAxes(Frame const& InputFrame) {
  Frame refFrame(maxResSize_); // Hold copy of base reference coords for RMS fit
  Frame inpFrame(maxResSize_); // Hold copy of input base coords for RMS fit
# ifdef NASTRUCTDEBUG
  PDBfile baseaxesfile;
  baseaxesfile.OpenWrite("baseaxes.pdb");
  PDBfile basesfile;
  basesfile.OpenWrite("bases.pdb");
  mprintf("\n=================== Setup Base Axes ===================\n");
# endif
  std::vector<NA_Axis>::iterator baseAxis = BaseAxes_.begin();
  for (std::vector<NA_Base>::iterator base = Bases_.begin(); 
                                      base != Bases_.end(); ++base, ++baseAxis)
  {
    // Set input coords for entire NA residue. 
    (*base).SetInputFrame( InputFrame );
    // Set input coords for RMS fit.
    inpFrame.SetCoordinates( (*base).Input(), (*base).InputFitMask());
    // Set ref coords for RMS fit. 
    refFrame.SetCoordinates( (*base).Ref(), (*base).RefFitMask());
#   ifdef NASTRUCTDEBUG
    mprintf("Base %i:%4s\n", (*base).ResNum()+1, (*base).ResName()); 
    (*base).InputFitMask().PrintMaskAtoms("InpMask");
    (*base).RefFitMask().PrintMaskAtoms("RefMask");
    mprintf("%-2s %4s %8s %8s %8s %2s %8s %8s %8s\n","#","Atom","Ex","Ey","Ez","#","Rx","Ry","Rz");
    AtomMask::const_iterator refatom = (*base).RefFitMask().begin();
    for (AtomMask::const_iterator inpatom = (*base).InputFitMask().begin();
                                  inpatom != (*base).InputFitMask().end(); ++inpatom)
    {
      const double* XYZ = (*base).Input().XYZ(*inpatom);
      mprintf("%-2i %4s %8.3f %8.3f %8.3f", *inpatom+1, (*base).AtomName(*inpatom),
              XYZ[0], XYZ[1], XYZ[2]);
      XYZ = (*base).Ref().XYZ(*refatom);
      mprintf(" %2i %8.3f %8.3f %8.3f\n", *refatom+1, XYZ[0], XYZ[1], XYZ[2]);
      ++refatom;
    }
#   endif 
    /* Now that we have a set of reference coords and the corresponding input
     * coords, RMS fit the reference coords to the input coords to obtain the
     * appropriate rotation and translations that will put the reference coords 
     * on top of input (experimental) coords. Per 3DNA procedure, not all 
     * reference atoms are used in the RMS fit; only ring atoms are used. 
     */ 
    Matrix_3x3 RotMatrix;
    Vec3 TransVec, refTrans;
    double rmsd = refFrame.RMSD( inpFrame, RotMatrix, TransVec, refTrans, false);
    /* RotMatrix and TransVec now contain rotation and translation
     * that will orient refcoord to expframe. The first translation is that of
     * the reference frame to the absolute origin, the second translation is
     * that of the reference frame to the exp. coords after rotation.
     * The rotation matrix contains the coordinates of the X, Y, and Z unit 
     * vectors of the base axes.
     */
    // Store the Rotation matrix and the rotated and translated origin.
    (*baseAxis).SetupBaseAxis( RotMatrix, (RotMatrix*TransVec)+refTrans, (*base).ResNum() );
    if (debug_>0) { 
      mprintf("Base %u: RMS of RefCoords from ExpCoords is %f\n",base-Bases_.begin(),rmsd);
      (*baseAxis).PrintAxisInfo("BaseAxes");
    }
#   ifdef NASTRUCTDEBUG
    // DEBUG - Write base axis to file
    WriteAxes(baseaxesfile, (*base).ResNum()+1, (*base).ResName(), *baseAxis);
     // Overlap ref coords onto input coords.
    Frame reftemp = (*base).Ref(); 
    reftemp.Trans_Rot_Trans(TransVec, RotMatrix, refTrans);
    // DEBUG - Write reference frame to file
    for (int i = 0; i < reftemp.Natom(); i++) {
      const double* XYZ = reftemp.XYZ(i);
      basesfile.WriteATOM( (*base).RefName(i), (*base).ResNum()+1, XYZ[0], XYZ[1], XYZ[2], 
                           (*base).ResName(), 0.0 );
    }
#   endif
  } // END loop over bases
  return 0;
}

// -----------------------------------------------------------------------------
// Action_NAstruct::GCpair()
/** Look for 3 HB based on heavy atom distances:
  * 1. G:O6 -- C:N4  6 -- 6
  * 2. G:N1 -- C:N3  7 -- 4
  * 3. G:N2 -- C:O2  9 -- 3
  */
int Action_NAstruct::GCpair(NA_Base const& DG, NA_Base const& DC) {
  int Nhbonds = 0;
  for (int hb = 0; hb < 3; hb++) {
    double dist2 = DIST2_NoImage(DG.HBxyz(hb), DC.HBxyz(hb));
    if ( dist2 < HBdistCut2_ ) {
      ++Nhbonds;
#     ifdef NASTRUCTDEBUG
      int dg_hbatom = DG.HBidx(hb);
      int dc_hbatom = DC.HBidx(hb);
      mprintf("\t\t%s:%s -- %s:%s = %f\n",
              DG.ResName(), DG.AtomName(dg_hbatom),
              DC.ResName(), DC.AtomName(dc_hbatom), sqrt(dist2));
#     endif
    }
  }
  return Nhbonds;
}

// Action_NAstruct::ATpair()
/** Look for 2 HB based on heavy atom distances
  * 1. A:N6 -- T:O4  6 -- 6
  * 2. A:N1 -- T:N3  7 -- 4
  */
int Action_NAstruct::ATpair(NA_Base const& DA, NA_Base const& DT) {
  int Nhbonds = 0;
  for (int hb = 0; hb < 2; hb++) {
    double dist2 = DIST2_NoImage(DA.HBxyz(hb), DT.HBxyz(hb)); 
    if ( dist2 < HBdistCut2_ ) {
      ++Nhbonds;
#     ifdef NASTRUCTDEBUG
      int da_hbatom = DA.HBidx(hb);
      int dt_hbatom = DT.HBidx(hb);
      mprintf("\t\t%s:%s -- %s:%s = %f\n",
              DA.ResName(), DA.AtomName(da_hbatom),
              DT.ResName(), DT.AtomName(dt_hbatom), sqrt(dist2));
#     endif
    }
  }
  return Nhbonds;
}

/** Given two NA_Bases for which IDs have been given and input coords set,
  * calculate the number of hydrogen bonds between them.
  * NOTE: Currently only set up for WC detection
  */
int Action_NAstruct::CalcNumHB(NA_Base const& base1, NA_Base const& base2) {
  // G C
  if      ( base1.Type()==NA_Base::GUA && base2.Type()==NA_Base::CYT ) return GCpair(base1,base2);
  else if ( base1.Type()==NA_Base::CYT && base2.Type()==NA_Base::GUA ) return GCpair(base2,base1);
  // A T
  else if ( base1.Type()==NA_Base::ADE && base2.Type()==NA_Base::THY ) return ATpair(base1,base2);
  else if ( base1.Type()==NA_Base::THY && base2.Type()==NA_Base::ADE ) return ATpair(base2,base1);
  // A U
  else if ( base1.Type()==NA_Base::ADE && base2.Type()==NA_Base::URA ) return ATpair(base1,base2);
  else if ( base1.Type()==NA_Base::URA && base2.Type()==NA_Base::ADE ) return ATpair(base2,base1);
//  else {
//    mprintf("Warning: NAstruct: Unrecognized pair: %s - %s\n",NAbaseName[base1->ID],
//             NAbaseName[base2->ID]);
//  }
  return 0;
}

// Action_NAstruct::determineBasePairing()
/** Determine which bases are paired from the individual base axes. Also 
  * sets up BP and BP step parameter DataSets. This routine should only
  * be called once.
  */
int Action_NAstruct::determineBasePairing() {
# ifdef NASTRUCTDEBUG  
  mprintf("\n=================== Setup Base Pairing ===================\n");
# endif
  if (!BasePairAxes_.empty()) {
    mprinterr("Internal Error: Base pairing has already been determined.\n");
    return 1;
  }
  // For each unpaired base, find the closest potential pairing base 
  // determined by the distance between their axis origins.
  NumberOfHbonds_.clear();
  int Nbases = (int)Bases_.size();
  int Nbases1 = Nbases - 1;
  std::vector<bool> isPaired( Nbases, false);
  for (int base1 = 0; base1 < Nbases1; base1++) {
#   ifdef NASTRUCTDEBUG
    mprintf("Base %i:%s\n",Bases_[base1].ResNum()+1, Bases_[base1].ResName());
#   endif
    if (isPaired[base1]) continue;
    int minBaseNum = -1;
    double minDistance = 0;
    for (int base2 = base1+1; base2 < Nbases; base2++) {
      if (isPaired[base2]) continue;
      double dist2 = DIST2_NoImage(BaseAxes_[base1].Oxyz(), BaseAxes_[base2].Oxyz());
      if (dist2 < originCut2_) {
#       ifdef NASTRUCTDEBUG
        mprintf("  Axes distance for %i:%s -- %i:%s is %f\n",
                Bases_[base1].ResNum()+1, Bases_[base1].ResName(), 
                Bases_[base2].ResNum()+1, Bases_[base2].ResName(), sqrt(dist2));
#       endif
        if (minBaseNum == -1) {
          minDistance = dist2;
          minBaseNum = base2;
        } else if (dist2 < minDistance) {
            minDistance = dist2;
            minBaseNum = base2;
        }
      }
    }
    if (minBaseNum != -1) {
      // Closest base to base1 is minBaseNum
#     ifdef NASTRUCTDEBUG
      mprintf("    Closest base is %i, %f Ang.\n", 
              Bases_[minBaseNum].ResNum()+1, sqrt(minDistance));
#         endif
      // Figure out if z vectors point in same (<90 deg) or opposite (>90 deg) direction
      bool AntiParallel;
      double dist2 = BaseAxes_[base1].Rz().Angle( BaseAxes_[minBaseNum].Rz() );
      //mprintf("    Dot product of Z vectors: %f\n", dist2);
      if (dist2 > Constants::PIOVER2) { // If theta(Z) > 90 deg.
#       ifdef NASTRUCTDEBUG
        mprintf("      %s is anti-parallel to %s\n", Bases_[base1].ResName(),
                Bases_[minBaseNum].ResName());
#       endif
        AntiParallel = true;
      } else {
#       ifdef NASTRUCTDEBUG
        mprintf("      %s is parallel to %s\n", Bases_[base1].ResName(),
                Bases_[minBaseNum].ResName());
#       endif
        AntiParallel = false;
      }
      int NHB = CalcNumHB(Bases_[base1], Bases_[minBaseNum]);
      if (NHB > 0) {
        BasePairAxes_.push_back( NA_Axis(base1, minBaseNum, AntiParallel) );
        NumberOfHbonds_.push_back( NHB );
        isPaired[base1] = true;
        isPaired[minBaseNum] = true;
      }
    } // END if minBaseNum!=-1
  } // END Loop over base1

  if (debug_>0) mprintf("    NAstruct: Detected %u base pairs.\n", BasePairAxes_.size());
  // Print Base Pair info
  if (debug_>1) {
    int nBP = 1; //  Base pair #
    for (std::vector<NA_Axis>::iterator BP = BasePairAxes_.begin();
                                        BP != BasePairAxes_.end(); ++BP)
    {
      int bp_1 = (*BP).Res1();
      int bp_2 = (*BP).Res2();
      mprintf("        BP %i: Res %i:%c to %i:%c", nBP++,
              bp_1+1, Bases_[bp_1].BaseChar(), 
              bp_2+1, Bases_[bp_2].BaseChar());
      if ( (*BP).IsAnti() )
        mprintf(" AntiParallel.\n");
      else
        mprintf(" Parallel.\n");
    }
  }
  // For each BP, set up a dataset for each structural parameter
  if (dataname_.empty())
    dataname_ = masterDSL_->GenerateDefaultName("NA");
  int dsidx = 1; // Base pair # and DataSet idx
  for (std::vector<NA_Axis>::iterator BP = BasePairAxes_.begin();
                                      BP != BasePairAxes_.end(); ++BP)
  {
    // Create legend
    int bp_1 = (*BP).Res1();
    int bp_2 = (*BP).Res2();
    std::string b1name = integerToString( Bases_[bp_1].ResNum()+1 ) +
                         Bases_[bp_1].BaseChar();
    std::string b2name = integerToString( Bases_[bp_2].ResNum()+1 ) + 
                         Bases_[bp_2].BaseChar();
    std::string bpname = b1name + b2name;
    // Create pucker Data Sets
    PUCKER_.push_back( 
      (DataSet_1D*)masterDSL_->AddSetIdxAspect(DataSet::FLOAT, dataname_,
                                               Bases_[bp_1].ResNum()+1, "pucker",
                                               b1name) );
    PUCKER_.back()->SetScalar( DataSet::M_PUCKER, DataSet::PUCKER );
    PUCKER_.push_back( 
      (DataSet_1D*)masterDSL_->AddSetIdxAspect(DataSet::FLOAT, dataname_,
                                               Bases_[bp_2].ResNum()+1, "pucker",
                                               b2name) );
    PUCKER_.back()->SetScalar( DataSet::M_PUCKER, DataSet::PUCKER );
    // Create sets
    SHEAR_.push_back( (DataSet_1D*)masterDSL_->AddSetIdxAspect(DataSet::FLOAT,dataname_,dsidx,"shear",bpname) );
    STRETCH_.push_back( (DataSet_1D*)masterDSL_->AddSetIdxAspect(DataSet::FLOAT,dataname_,dsidx,"stretch",bpname));
    STAGGER_.push_back( (DataSet_1D*)masterDSL_->AddSetIdxAspect(DataSet::FLOAT,dataname_,dsidx,"stagger",bpname));
    BUCKLE_.push_back( (DataSet_1D*)masterDSL_->AddSetIdxAspect(DataSet::FLOAT,dataname_,dsidx,"buckle",bpname) );
    PROPELLER_.push_back( (DataSet_1D*)masterDSL_->AddSetIdxAspect(DataSet::FLOAT,dataname_,dsidx,"prop",bpname) );
    OPENING_.push_back( (DataSet_1D*)masterDSL_->AddSetIdxAspect(DataSet::FLOAT,dataname_,dsidx,"open",bpname) );
    BPHBONDS_.push_back( (DataSet_1D*)masterDSL_->AddSetIdxAspect(DataSet::INTEGER,dataname_,dsidx,"hb",bpname) );
    MAJOR_.push_back( (DataSet_1D*)masterDSL_->AddSetIdxAspect(DataSet::FLOAT,dataname_,dsidx,"major",bpname) );
    MINOR_.push_back( (DataSet_1D*)masterDSL_->AddSetIdxAspect(DataSet::FLOAT,dataname_,dsidx,"minor",bpname) );
    ++dsidx;
  }
  // For each BP step, set up a dataset for each structural parameter. 
  // One less than total # BP.
  if (BasePairAxes_.size() > 1) {
    dsidx = 1; // Base pair step # and DataSet idx
    std::vector<NA_Axis>::iterator NBPstep = BasePairAxes_.end() - 1;
    for (std::vector<NA_Axis>::iterator BP1 = BasePairAxes_.begin();
                                        BP1 != NBPstep; ++BP1)
    {
      std::vector<NA_Axis>::iterator BP2 = BP1 + 1;
      // Create legend
      int bp_1 = (*BP1).Res1();
      int bp_2 = (*BP1).Res2();
      int bp_3 = (*BP2).Res1();
      int bp_4 = (*BP2).Res2();
      std::string sname = integerToString( Bases_[bp_1].ResNum()+1 ) +
                          Bases_[bp_1].BaseChar() +
                          integerToString( Bases_[bp_2].ResNum()+1 ) +
                          Bases_[bp_2].BaseChar() + "-" +
                          integerToString( Bases_[bp_3].ResNum()+1 ) +
                          Bases_[bp_3].BaseChar() +
                          integerToString( Bases_[bp_4].ResNum()+1 ) +
                          Bases_[bp_4].BaseChar();
      // Create Sets
      SHIFT_.push_back( (DataSet_1D*)masterDSL_->AddSetIdxAspect(DataSet::FLOAT,dataname_,dsidx,"shift",sname) );
      SLIDE_.push_back( (DataSet_1D*)masterDSL_->AddSetIdxAspect(DataSet::FLOAT,dataname_,dsidx,"slide",sname) );
      RISE_.push_back( (DataSet_1D*)masterDSL_->AddSetIdxAspect(DataSet::FLOAT,dataname_,dsidx,"rise",sname) );
      TILT_.push_back( (DataSet_1D*)masterDSL_->AddSetIdxAspect(DataSet::FLOAT,dataname_,dsidx,"tilt",sname) );
      ROLL_.push_back( (DataSet_1D*)masterDSL_->AddSetIdxAspect(DataSet::FLOAT,dataname_,dsidx,"roll",sname) );
      TWIST_.push_back( (DataSet_1D*)masterDSL_->AddSetIdxAspect(DataSet::FLOAT,dataname_,dsidx,"twist",sname) );
      XDISP_.push_back( (DataSet_1D*)masterDSL_->AddSetIdxAspect(DataSet::FLOAT,dataname_,dsidx,"xdisp",sname) );
      YDISP_.push_back( (DataSet_1D*)masterDSL_->AddSetIdxAspect(DataSet::FLOAT,dataname_,dsidx,"ydisp",sname) );
      HRISE_.push_back( (DataSet_1D*)masterDSL_->AddSetIdxAspect(DataSet::FLOAT,dataname_,dsidx,"hrise",sname) );
      INCL_.push_back( (DataSet_1D*)masterDSL_->AddSetIdxAspect(DataSet::FLOAT,dataname_,dsidx,"incl",sname) );
      TIP_.push_back( (DataSet_1D*)masterDSL_->AddSetIdxAspect(DataSet::FLOAT,dataname_,dsidx,"tip",sname) );
      HTWIST_.push_back( (DataSet_1D*)masterDSL_->AddSetIdxAspect(DataSet::FLOAT,dataname_,dsidx,"htwist",sname) );
      ++dsidx;
    }
  }
  return 0;
}

// -----------------------------------------------------------------------------
// AverageMatrices()
static Matrix_3x3 AverageMatrices(Matrix_3x3 const& RotatedR1, Matrix_3x3 const& RotatedR2) {
  Matrix_3x3 R;
  // Average R1 and R2 to get the middle frame
  for (int i = 0; i < 9; i++)
    R[i] = (RotatedR1[i] + RotatedR2[i]) / 2;
  // Normalize X, Y and Z vectors
  double r2 = sqrt( R[0]*R[0] + R[3]*R[3] + R[6]*R[6] );
  R[0] /= r2;
  R[3] /= r2;
  R[6] /= r2;
  r2 = sqrt( R[1]*R[1] + R[4]*R[4] + R[7]*R[7] );
  R[1] /= r2;
  R[4] /= r2;
  R[7] /= r2;
  r2 = sqrt( R[2]*R[2] + R[5]*R[5] + R[8]*R[8] );
  R[2] /= r2;
  R[5] /= r2;
  R[8] /= r2;
  return R;
}

// Action_NAstruct::calculateParameters()
/** Given two base axes, calculate translational and rotational parameters
  * between them. Store base pair axis info in BPaxis if given.
  */
int Action_NAstruct::calculateParameters(NA_Axis const& Axis1, NA_Axis const& Axis2, 
                                         NA_Axis* BPaxis, double *Param) 
{
# ifdef NASTRUCTDEBUG
  NA_Axis tempAxis;
  PDBfile paramfile;
  if (calcparam_)
    paramfile.OpenWrite("Param.pdb");
  Axis1.Oxyz().Print("O1");
  Axis1.Rot().Print("R1");
  Axis2.Oxyz().Print("O2");
  Axis2.Rot().Print("R2");
# endif
  // Hinge axis is cross product between Z1 and Z2
  Vec3 hingeAxis = Axis1.Rz().Cross( Axis2.Rz() );
# ifdef NASTRUCTDEBUG
  hingeAxis.Print("hinge");
# endif
  // Normalize hinge axis
  hingeAxis.Normalize();
# ifdef NASTRUCTDEBUG
  hingeAxis.Print("norm(hinge)");
# endif
  // Roll/Tilt is Angle between Z1 and Z2
  double rolltilt = Axis1.Rz().Angle( Axis2.Rz() );
# ifdef NASTRUCTDEBUG
  mprintf("\tAngle between Z1 and Z2= %f\n",rolltilt*Constants::RADDEG);
# endif
  // Calculate forward and backwards half rolltilt rotation around
  // hinge axis.
  Matrix_3x3 R;
  R.CalcRotationMatrix(hingeAxis, -0.5*rolltilt);
# ifdef NASTRUCTDEBUG
  R.Print("Rhalf");
# endif
  // Rotate R2 by -0.5 * rolltilt degrees around the hinge
  Matrix_3x3 RotatedR2 = R * Axis2.Rot();
  // Rotate R1 by 0.5 * rolltilt degrees around the hinge (inverse rotation)
  R.Transpose();
  Matrix_3x3 RotatedR1 = R * Axis1.Rot();
# ifdef NASTRUCTDEBUG
  // Print rotated R1 and R2
  RotatedR1.Print("Rotated R1");
  RotatedR1.Print("Rotated R2");
  if (calcparam_) {
    tempAxis.StoreRotMatrix(RotatedR1, Axis1.Oxyz()); 
    WriteAxes(paramfile, 1, "R1'", tempAxis);
    tempAxis.StoreRotMatrix(RotatedR2, Axis2.Oxyz());
    WriteAxes(paramfile, 2, "R2'", tempAxis);
  }
# endif
  // Average rotated R1 and R2 to get the middle frame
  Matrix_3x3 Rm = AverageMatrices(RotatedR1, RotatedR2);
  // Take average of origins
  Vec3 OM = (Axis1.Oxyz() + Axis2.Oxyz()) / 2.0;
# ifdef NASTRUCTDEBUG
  OM.Print("Origin Mean");
  // Print Rm and hinge axis
  Rm.Print("Rm");
  if (calcparam_) {
    // Use R to store hinge axis in Z
    R.Zero();
    R[2] = hingeAxis[0]; 
    R[5] = hingeAxis[1]; 
    R[8] = hingeAxis[2];
    tempAxis.StoreRotMatrix(R, OM);
    WriteAxes(paramfile, 3, "Hng", tempAxis);
    // Store middle frame
    tempAxis.StoreRotMatrix(Rm, OM);
    WriteAxes(paramfile, 4, "Rm", tempAxis);
  }
# endif

  // If BPaxis is not null, store Rm and OM as BP axis.
  if (BPaxis != 0)  
    BPaxis->StoreRotMatrix(Rm, OM);

  // Shift Slide Rise / Shear Stretch Stagger
  OM = Axis2.Oxyz() - Axis1.Oxyz();
  // Since this is really vector times matrix, use matrix transpose times vec
  Vec3 Vec = Rm.TransposeMult( OM );
# ifdef NASTRUCTDEBUG
  OM.Print("O21");
  Vec.Print("Vec");
# endif
  Param[0] = Vec[0];
  Param[1] = Vec[1];
  Param[2] = Vec[2];

  // Set Z1 to Z from middle frame
  Vec3 Z1 = Rm.Col3();
  // Twist / Opening
  // Angle between rotated Y1 and rotated Y2
  // Sign of twistopen related to (Y1'xY2') dot Z of middle frame
  Vec3 Y1 = RotatedR1.Col2();
  Vec3 Y2 = RotatedR2.Col2();
  double twistopen = Y1.SignedAngle(Y2, Z1);
# ifdef NASTRUCTDEBUG
  mprintf("\tFinal Twist/Opening is %10.4f\n",twistopen*Constants::RADDEG);
# endif
  Param[3] = twistopen;

  // Phase angle
  // Angle between hinge axis and middle frame Y axis
  // Sign of phi related to (hingeAxis x Ym) dot Z of middle frame
  Y1 = Rm.Col2();
  double phi = hingeAxis.SignedAngle(Y1, Z1);
  double sinphi = sin( phi );
  double cosphi = cos( phi );
# ifdef NASTRUCTDEBUG
  mprintf("\tPhase angle is %f, sinphi is %f, cosphi is %f\n",phi*Constants::RADDEG,sinphi,cosphi);
# endif

  // Roll / Propeller
  double rollprop = rolltilt * cosphi;
  Param[4] = rollprop;

  // Tilt / Buckle
  double tiltbuck = rolltilt * sinphi;
  Param[5] = tiltbuck;

# ifdef NASTRUCTDEBUG
  mprintf("\tRoll/Propeller %10.4f\n",rollprop*Constants::RADDEG);
  mprintf("\tTilt/Buckle %10.4f\n",tiltbuck*Constants::RADDEG);
  if (calcparam_) calcparam_=false;
# endif
  return 0;
}

// Action_NAstruct::helicalParameters()
int Action_NAstruct::helicalParameters(NA_Axis const& Axis1, NA_Axis const& Axis2, double *Param) 
{
  // O1 = X2 - X1
  Vec3 O1 = Axis2.Rx() - Axis1.Rx();
  // O2 = Y2 - Y1
  Vec3 O2 = Axis2.Ry() - Axis1.Ry();
  // Local helical axis: (X2-X1) x (Y2-Y1)
  Vec3 helicalAxis = O1.Cross( O2 );
# ifdef NASTRUCTDEBUG
  O1.Print("X2 - X1");
  O2.Print("Y2 - Y1");
  helicalAxis.Print("(X2-X1) x (Y2-Y1)");
# endif
  helicalAxis.Normalize( );
# ifdef NASTRUCTDEBUG
  helicalAxis.Print("NORM[(X2-X1)x(Y2-Y1)]");
# endif

  // Tip/inclination is angle between helical axis and z1
  double tipinc = helicalAxis.Angle( Axis1.Rz() );
  // Hinge axis is normalized cross product of helical axis to z1
  Vec3 hingeAxis = helicalAxis.Cross( Axis1.Rz() );
  hingeAxis.Normalize();
  // Rotate R1 around hinge axis by -tipinc
  Matrix_3x3 R;
  R.CalcRotationMatrix(hingeAxis, -tipinc);
  Matrix_3x3 RotatedR1 = R * Axis1.Rot();
# ifdef NASTRUCTDEBUG
  mprintf("\tTip/Inclination: %f\n",tipinc*Constants::RADDEG);
  hingeAxis.Print("Hinge axis 1");
  RotatedR1.Print("Rotated R1");
# endif

  // Tip/inclination should be same for z2
  //mprintf("\tTipCheck= %f\n",dot_product_angle(helicalAxis, Z2)*Constants::RADDEG);
  // Hinge axis (Vec) is normalized cross product from h to z2
  Vec3 Vec = helicalAxis.Cross( Axis2.Rz() );
  Vec.Normalize();
  // Rotate R2 around hinge axis by -tipinc
  R.CalcRotationMatrix(Vec, -tipinc); 
  Matrix_3x3 RotatedR2 = R * Axis2.Rot();
# ifdef NASTRUCTDEBUG
  Vec.Print("Hinge axis 2");
  RotatedR2.Print("Rotated R2");
# endif

  // Average Rotated R1 and R2 to get middle helical frame
  R = AverageMatrices(RotatedR1, RotatedR2);

  // Helical twist is angle from Rotated Y1 to Rotated Y2
  // Sign is given by (Y1'xY2' dot helicalAxis)
  Vec3 Y1 = RotatedR1.Col2();
  Vec3 Y2 = RotatedR2.Col2();
  double Twist = Y1.SignedAngle(Y2, helicalAxis);
  Param[5] = Twist;

  // Calc Vec = O2 - O1
  Vec = Axis2.Oxyz() - Axis1.Oxyz();
  // Project (O2-O1) onto helical axis
  double Rise = Vec * helicalAxis;
  Param[2] = Rise;
# ifdef NASTRUCTDEBUG
  R.Print("Hm");
  mprintf("\tTwist is %f\n",Twist*Constants::RADDEG);
  mprintf("\tRise is %f\n",Rise);
# endif

  // Phase angle is angle from hinge Axis 1 to RotatedR1 Y
  // Sign is given by (hingeAxis x Y1') dot helicalAxis
  double phase = hingeAxis.SignedAngle(Y1, helicalAxis);

  // Tip is tipinc * cos( phase )
  double Tip = tipinc * cos( phase );
  Param[4] = Tip;
  // Inclination is tipinc * sin( phase )
  double Inc = tipinc * sin( phase );
  Param[3] = Inc;
# ifdef NASTRUCTDEBUG
  mprintf("\tPhase angle is %f\n",phase*Constants::RADDEG);
  mprintf("\tTip is %f\n",Tip*Constants::RADDEG);
  mprintf("\tInclination is %f\n",Inc*Constants::RADDEG);
# endif

  Vec3 Z1 = helicalAxis * Rise;
  // Calc vector AB (store in O1) 
  // Vec contains O2-O1
  O1 = Vec - Z1; 

  // Calc vector AD
  double AD_angle = Constants::PIOVER2 - (0.5 * Twist);
  // rotation of AD_angle around helicalAxis
  // NOTE: Assuming we dont need RotatedR2 anymore
  RotatedR2.CalcRotationMatrix(helicalAxis, AD_angle);
  // rotate AB, = AD (store in O2)
  O2 = RotatedR2 * O1;
  O2.Normalize();
# ifdef NASTRUCTDEBUG
  O1.Print("AB");
  mprintf("\tAD_angle is %f\n",AD_angle*Constants::RADDEG);
  O2.Print("AD");
# endif

  // Calc magnitude of AD; 0.5 * |AB| / sin( 0.5 * Twist )
  double AD_mag = (0.5 * sqrt(O1.Magnitude2())) / sin( 0.5 * Twist );

  // Calc origin of local helical frame for BP 1
  // O1 = Origin1 + (AD_mag * AD)
  O1 = Axis1.Oxyz() + (O2 * AD_mag); 

  // Calc origin of local helical frame for BP 2
  // O2 = O1 + (Rise * helicalAxis)
  // Z1 contains helicalAxis * Rise
  O2 = O1 + Z1;

  // Calculate origin of middle helical frame, store in Vec 
  Vec = (O2 + O1) / 2.0;
# ifdef NASTRUCTDEBUG
  mprintf("\t|AD| = %f\n",AD_mag);
  O1.Print("o1_h");
  O2.Print("o2_h");
  Vec.Print("Om_h");
# endif

  // Calc vector from O1 to Origin1
  Vec = Axis1.Oxyz() - O1;

  // X-disp is projection of vector from O1 to Origin1 onto 
  // X axis of RotatedR1.
  // TODO: X1 could just be O1
  Vec3 X1 = RotatedR1.Col1();
  double X_disp = Vec * X1;
  Param[0] = X_disp;

  // Y-disp is projection of vector from O1 to Origin1 onto 
  // Y axis of RotatedR1.
  X1 = RotatedR1.Col2();
  double Y_disp = Vec * X1;
  Param[1] = Y_disp;
# ifdef NASTRUCTDEBUG
  mprintf("\tX-displacement= %f\n",X_disp);
  mprintf("\tY-displacement= %f\n",Y_disp);
# endif

  return 0;
}

/// Calculate nucleic acid sugar pucker.
void Action_NAstruct::CalcPucker( NA_Base const& base, int framenum, int nbase ) {
  double pval=0.0, aval, tval;
  switch (puckerMethod_) {
    case ALTONA:
      pval = Pucker_AS( base.C1xyz(), base.C2xyz(), base.C3xyz(),
                        base.C4xyz(), base.O4xyz(), aval );
      break;
    case CREMER:
      pval = Pucker_CP( base.C1xyz(), base.C2xyz(), base.C3xyz(),
                        base.C4xyz(), base.O4xyz(), 0,
                        5, aval, tval );
      break;
  }
  float fval = (float)(pval * Constants::RADDEG);
  PUCKER_[nbase]->Add(framenum, &fval);
}

// Action_NAstruct::determineBaseParameters()
/** For each base in a base pair, get the values of buckle, propeller twist,
  * opening, shear, stretch, and stagger. Also determine the origin and 
  * rotation matrix for each base pair reference frame.
  */
int Action_NAstruct::determineBaseParameters(int frameNum) {
  double Param[6];
# ifdef NASTRUCTDEBUG
  PDBfile basepairaxesfile;
  basepairaxesfile.OpenWrite("basepairaxes.pdb");
  mprintf("\n=================== Determine BP Parameters ===================\n");
# endif

  int nbasepair = 0;
  for (std::vector<NA_Axis>::iterator BP = BasePairAxes_.begin();
                                      BP != BasePairAxes_.end(); ++BP)
  {
    int b1 = (*BP).Res1();
    int b2 = (*BP).Res2();
    const NA_Base& base1 = Bases_[b1];
    const NA_Base& base2 = Bases_[b2]; 
#   ifdef NASTRUCTDEBUG
    mprintf("BasePair %i:%s to %i:%s", b1+1, base1.ResName(), b2+1, base2.ResName());
    if ((*BP).IsAnti())
      mprintf(" Anti-parallel.\n");
    else
      mprintf(" Parallel.\n");
#   endif
    // Check Antiparallel / Parallel
    // Flip YZ (rotate around X) for antiparallel
    // Flip XY (rotate around Z) for parallel
    if ((*BP).IsAnti())
      BaseAxes_[b2].FlipYZ();
    else
      BaseAxes_[b2].FlipXY();
    // TEST - calc P--P distance
    float dPtoP = 0.0;
    //mprintf("\tDEBUG: %i %i:", BaseAxes[base1].Pidx(), BaseAxes[base2].Pidx() );
    if ( base1.HasPatom() && base2.HasPatom() ) {
      double DP = DIST2_NoImage( base1.Pxyz(), base2.Pxyz() );
      //mprintf(" %i to %i P--P D= %f", BaseAxes[base1].ResNum()+1, BaseAxes[base2].ResNum()+1,
      //        sqrt(dPtoP) );
      DP = sqrt(DP);
      dPtoP = (float)DP;
    }
    //mprintf("\n");
    float dOtoO = 0.0;
    //mprintf("\tDEBUG: %i %i:", BaseAxes[base1].O4idx(), BaseAxes[base2].O4idx() );
    if ( base1.HasO4atom() && base2.HasO4atom() ) {
      double DO4 = DIST2_NoImage( base1.O4xyz(), base2.O4xyz() );
      //mprintf(" %i to %i O4'--O4' D= %f", BaseAxes[base1].ResNum()+1, BaseAxes[base2].ResNum()+1,
      //        sqrt(dOtoO) );
      DO4 = sqrt(DO4);
      dOtoO = (float)DO4;
    }
    if (base1.HasSugarAtoms())
      CalcPucker( base1, frameNum, nbasepair*2 );
    if (base2.HasSugarAtoms())
      CalcPucker( base2, frameNum, (nbasepair*2) + 1 );
    //mprintf("\n");
    // Calc BP parameters, set up basepair axes
    //calculateParameters(BaseAxes[base1],BaseAxes[base2],&BasePairAxes[nbasepair],Param);
    calculateParameters(BaseAxes_[b2], BaseAxes_[b1], &(*BP), Param);
    // Store data
    Param[3] *= Constants::RADDEG;
    Param[4] *= Constants::RADDEG;
    Param[5] *= Constants::RADDEG;
    //mprintf("DBG: BP %i # hbonds = %i\n", nbasepair+1, NumberOfHbonds_[nbasepair]);
    // Convert everything to float to save space
    float shear = (float)Param[0];
    float stretch = (float)Param[1];
    float stagger = (float)Param[2];
    float opening = (float)Param[3];
    float prop = (float)Param[4];
    float buckle = (float)Param[5];
    int n_of_hb = NumberOfHbonds_[nbasepair];
    // Add to DataSets
    SHEAR_[nbasepair]->Add(frameNum, &shear);
    STRETCH_[nbasepair]->Add(frameNum, &stretch);
    STAGGER_[nbasepair]->Add(frameNum, &stagger);
    OPENING_[nbasepair]->Add(frameNum, &opening);
    PROPELLER_[nbasepair]->Add(frameNum, &prop);
    BUCKLE_[nbasepair]->Add(frameNum, &buckle);
    BPHBONDS_[nbasepair]->Add(frameNum, &n_of_hb);
    MAJOR_[nbasepair]->Add(frameNum, &dPtoP);
    MINOR_[nbasepair]->Add(frameNum, &dOtoO);
#   ifdef NASTRUCTDEBUG
    // DEBUG - write base pair axes
    WriteAxes(basepairaxesfile, b1+1, base1.ResName(), *BP);
#   endif
    ++nbasepair; 
  }

  return 0;
}

// Action_NAstruct::determineBasepairParameters() 
/** For each base pair step, determine values of Tilt, Roll, Twist, Shift,
  * Slide, and Rise.
  */
int Action_NAstruct::determineBasepairParameters(int frameNum) {
  double Param[6];
# ifdef NASTRUCTDEBUG
  mprintf("\n=================== Determine BPstep Parameters ===================\n");
# endif
  if (BasePairAxes_.size() < 2) return 0;
  int bpi = 0;
  std::vector<NA_Axis>::iterator NBPstep = BasePairAxes_.end() - 1;
  for (std::vector<NA_Axis>::iterator BP1 = BasePairAxes_.begin();
                                      BP1 != NBPstep; ++BP1)
  {
    std::vector<NA_Axis>::iterator BP2 = BP1 + 1;
#   ifdef NASTRUCTDEBUG
    mprintf("BasePair step %i to %i\n", bpi+1, bpi+2);
#   endif
    // Calc step parameters
    calculateParameters(*BP1, *BP2, 0, Param);
    // Store data
    Param[3] *= Constants::RADDEG;
    Param[4] *= Constants::RADDEG;
    Param[5] *= Constants::RADDEG;
    // Convert everything to float to save space
    float shift = (float)Param[0];
    float slide = (float)Param[1];
    float rise = (float)Param[2];
    float twist = (float)Param[3];
    float roll = (float)Param[4];
    float tilt = (float)Param[5];
    SHIFT_[bpi]->Add(frameNum, &shift);
    SLIDE_[bpi]->Add(frameNum, &slide);
    RISE_[bpi]->Add(frameNum, &rise);
    TWIST_[bpi]->Add(frameNum, &twist);
    ROLL_[bpi]->Add(frameNum, &roll);
    TILT_[bpi]->Add(frameNum, &tilt);
    // Calc helical parameters
    helicalParameters(*BP1, *BP2, Param);
    Param[3] *= Constants::RADDEG;
    Param[4] *= Constants::RADDEG;
    Param[5] *= Constants::RADDEG;
    // Convert to float
    float xdisp = (float)Param[0];
    float ydisp = (float)Param[1];
    float hrise = (float)Param[2];
    float incl = (float)Param[3];
    float tip = (float)Param[4];
    float htwist = (float)Param[5];
    XDISP_[bpi]->Add(frameNum, &xdisp);
    YDISP_[bpi]->Add(frameNum, &ydisp);
    HRISE_[bpi]->Add(frameNum, &hrise);
    INCL_[bpi]->Add(frameNum, &incl);
    TIP_[bpi]->Add(frameNum, &tip);
    HTWIST_[bpi]->Add(frameNum, &htwist);
    ++bpi; 
  }

  return 0;
}
// ----------------------------------------------------------------------------

// Action_NAstruct::Init()
Action::RetType Action_NAstruct::Init(ArgList& actionArgs, TopologyList* PFL, DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  debug_ = debugIn;
  masterDSL_ = DSL;
  ensembleNum_ = DSL->EnsembleNum();
  // Get keywords
  outputsuffix_ = actionArgs.GetStringKey("naout");
  double hbcut = actionArgs.getKeyDouble("hbcut", -1);
  if (hbcut > 0) 
    HBdistCut2_ = hbcut * hbcut;
  double origincut = actionArgs.getKeyDouble("origincut", -1);
  if (origincut > 0)
    originCut2_ = origincut * origincut;
  if      (actionArgs.hasKey("altona")) puckerMethod_=ALTONA;
  else if (actionArgs.hasKey("cremer")) puckerMethod_=CREMER;
  // Get residue range
  resRange_.SetRange(actionArgs.GetStringKey("resrange"));
  if (!resRange_.Empty())
    resRange_.ShiftBy(-1); // User res args start from 1
  printheader_ = !actionArgs.hasKey("noheader");
  // Reference for setting up basepairs
  ReferenceFrame REF = DSL->GetReferenceFrame( actionArgs );
  if (REF.error()) return Action::ERR;
  if (!REF.empty()) 
    useReference_ = true;

  // Get custom residue maps
  ArgList maplist;
  NA_Base::NAType mapbase;
  while ( actionArgs.Contains("resmap") ) {
    // Split maparg at ':'
    maplist.SetList( actionArgs.GetStringKey("resmap"), ":" );
    // Expect only 2 args
    if (maplist.Nargs()!=2) {
      mprinterr("Error: resmap format should be '<ResName>:{A,C,G,T,U}' (%s)\n",
                maplist.ArgLine());
      return Action::ERR;
    }
    // Check that second arg is A,C,G,T,or U
    if      (maplist[1] == "A") mapbase = NA_Base::ADE;
    else if (maplist[1] == "C") mapbase = NA_Base::CYT;
    else if (maplist[1] == "G") mapbase = NA_Base::GUA;
    else if (maplist[1] == "T") mapbase = NA_Base::THY;
    else if (maplist[1] == "U") mapbase = NA_Base::URA;
    else {
      mprinterr("Error: resmap format should be '<ResName>:{A,C,G,T,U}' (%s)\n",
                maplist.ArgLine());
      return Action::ERR;
    }
    // Check that residue name is <= 4 chars
    std::string resname = maplist[0]; 
    if (resname.size() > 4) {
      mprinterr("Error: resmap resname > 4 chars (%s)\n",maplist.ArgLine());
      return Action::ERR;
    }
    // Format residue name
    // TODO: Use NameType in map
    NameType mapresname = resname;
    resname.assign( *mapresname );
    mprintf("\tCustom Map: [%s]\n",resname.c_str());
    //maplist.PrintList();
    // Add to CustomMap
    ResMapType::iterator customRes = CustomMap_.find(resname);
    if (customRes != CustomMap_.end()) {
      mprintf("Warning: resmap: %s already mapped.\n",resname.c_str());
    } else {
      CustomMap_.insert( std::pair<std::string,NA_Base::NAType>(resname,mapbase) );
    }
  }
  // Get Masks
  // Dataset
  dataname_ = actionArgs.GetStringNext();
  // DataSets are added to data file list in print()

  mprintf("    NAstruct: ");
  if (resRange_.Empty())
    mprintf("Scanning all NA residues");
  else
    mprintf("Scanning residues %s",resRange_.RangeArg());
  if (!outputsuffix_.empty()) {
      mprintf(", formatted output using file suffix %s",outputsuffix_.c_str());
    if (!printheader_) mprintf(", no header");
  }
  mprintf(".\n");
  mprintf("\tHydrogen bond cutoff for determining base pairs is %.2f Angstroms.\n",
          sqrt( HBdistCut2_ ) );
  mprintf("\tBase reference axes origin cutoff for determining base pairs is %.2f Angstroms.\n",
          sqrt( originCut2_ ) );

  // Use reference to determine base pairing
  if (useReference_) {
    mprintf("\tUsing reference %s to determine base-pairing.\n", REF.FrameName().base());
    if (Setup((Topology*)(&REF.Parm()), 0)) return Action::ERR;
    // Set up base axes
    if ( setupBaseAxes(REF.Coord()) ) return Action::ERR;
    // Determine Base Pairing
    if ( determineBasePairing() ) return Action::ERR;
    mprintf("\tSet up %zu base pairs.\n", BasePairAxes_.size() ); 
  } else
    mprintf("\tUsing first frame to determine base pairing.\n");
  if (puckerMethod_==ALTONA)
    mprintf("\tCalculating sugar pucker using Altona & Sundaralingam method.\n");
  else if (puckerMethod_==CREMER)
    mprintf("\tCalculating sugar pucker using Cremer & Pople method.\n");
  DSL->SetDataSetsPending(true);
  return Action::OK;
}

// Action_NAstruct::Setup()
/** Determine the number of NA bases that will be analyzed, along with 
  * the masks that correspond to the reference frame atoms.
  */
Action::RetType Action_NAstruct::Setup(Topology* currentParm, Topology** parmAddress) {
  // Clear Bases and BaseAxes
  ClearLists();
  // If range is empty (i.e. no resrange arg given) look through all 
  // solute residues.
  Range actualRange;
  if (resRange_.Empty()) 
    actualRange = currentParm->SoluteResidues();
  else 
    actualRange = resRange_;
  // Exit if no residues specified
  if (actualRange.Empty()) {
    mprinterr("Error: NAstruct::setup: No residues specified for %s\n",currentParm->c_str());
    return Action::ERR;
  }

  // DEBUG - print all residues
  //if (debug>0)
  //  actualRange.PrintRange("    NAstruct: NA res:",1);

  // Set up NA_base for each selected NA residue 
  for (Range::const_iterator resnum = actualRange.begin();
                             resnum != actualRange.end(); ++resnum)
  {
    NA_Base::NAType baseType = NA_Base::UNKNOWN_BASE;
#   ifdef NASTRUCTDEBUG
    mprintf(" ----- Setting up %i:%s -----\n", *resnum+1, currentParm->Res(*resnum).c_str());
#   endif
    // Check if the residue at resnum matches any of the custom maps
    if (!CustomMap_.empty()) {
      std::string resname( currentParm->Res(*resnum).c_str() );
      ResMapType::iterator customRes = CustomMap_.find( resname );
      if (customRes != CustomMap_.end()) {
        mprintf("\tCustom map found: %i [%s]\n",*resnum+1,(*customRes).first.c_str());
        baseType = (*customRes).second;
      }
    }
    // If not in custom map, attempt to identify base from name
    if (baseType == NA_Base::UNKNOWN_BASE)
      baseType = NA_Base::ID_BaseFromName( currentParm->Res(*resnum).Name() );
    // If still unknown skip to the next base
    if (baseType == NA_Base::UNKNOWN_BASE) {
      // Print a warning if the user specified this range.
      if (!resRange_.Empty()) {
        mprintf("Warning: Residue %i:%s not recognized as NA residue.\n",
                *resnum+1, currentParm->Res(*resnum).c_str());
      }
      continue;
    }
    // Set up ref coords for this base type. If there is a problem
    // the type will be reset to UNKNOWN.
    Bases_.push_back( NA_Base( *currentParm, *resnum, baseType ) );
    if (Bases_.back().Type() == NA_Base::UNKNOWN_BASE) {
      mprinterr("Error: NAstruct::setup: Could not get ref coords for %i:%s\n",
                *resnum+1, currentParm->Res(*resnum).c_str());
      return Action::ERR;
    }
    // Determine the largest residue for setting up frames for RMS fit later.
    if (Bases_.back().InputFitMask().Nselected() > maxResSize_)
      maxResSize_ = Bases_.back().InputFitMask().Nselected();
    if (debug_>1) {
      mprintf("\tNAstruct: Res %i:%s ", *resnum+1, currentParm->Res(*resnum).c_str());
      Bases_.back().PrintAtomNames();
      Bases_.back().InputFitMask().PrintMaskAtoms("InpMask");
      Bases_.back().RefFitMask().PrintMaskAtoms("RefMask");
    }
  } // End Loop over NA residues
  // Allocate space for base axes
  BaseAxes_.resize( Bases_.size() );
  mprintf("\tSet up %zu bases.\n", Bases_.size());

  return Action::OK;  
}

// Action_NAstruct::DoAction()
Action::RetType Action_NAstruct::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {
  // Set up base axes
  if ( setupBaseAxes(*currentFrame) ) return Action::ERR;

  if (!useReference_) {
    // Determine Base Pairing based on first frame
    if ( determineBasePairing() ) return Action::ERR;
    useReference_ = true;
  } else {
    // Base pairing determined from ref. Just calc # hbonds for each pair.
    int bp = 0;
    for (std::vector<NA_Axis>::iterator BP = BasePairAxes_.begin();
                                        BP != BasePairAxes_.end(); ++BP)
      NumberOfHbonds_[bp++] = CalcNumHB(Bases_[(*BP).Res1()], Bases_[(*BP).Res2()]);
  }

  // Determine base parameters
  determineBaseParameters(frameNum);

  // Determine base pair parameters
  determineBasepairParameters(frameNum);

  return Action::OK;
} 

// Action_NAstruct::print()
void Action_NAstruct::Print() {
  CpptrajFile outfile;
  int nframes;
  if (outputsuffix_.empty()) return;

  // ---------- Base pair parameters ----------
  std::string outfilename = "BP." + outputsuffix_;
  // Check that there is actually data
  if ( SHEAR_.empty() || SHEAR_[0]->Empty() )
    mprinterr("Error: Could not write BP file %s: No BP data.\n",outfilename.c_str()); 
  else {
    if (outfile.OpenEnsembleWrite( outfilename, ensembleNum_ ) == 0) {
      // Determine number of frames from SHEAR[0] DataSet
      nframes = SHEAR_[0]->Size();
      mprintf("\tBase pair output file %s; %i frames, %zu base pairs.\n", 
              outfilename.c_str(), nframes, BasePairAxes_.size());
      //  File header
      if (printheader_)
        outfile.Printf("%-8s %8s %8s %10s %10s %10s %10s %10s %10s %2s %10s %10s\n",
                       "#Frame","Base1","Base2", "Shear","Stretch","Stagger",
                       "Buckle","Propeller","Opening", "HB", "Major", "Minor");
      // Loop over all frames
      for (int frame = 0; frame < nframes; ++frame) {
        int nbp = 0;
        for (std::vector<NA_Axis>::iterator BP = BasePairAxes_.begin();
                                        BP != BasePairAxes_.end(); ++BP) {
          int bp_1 = (*BP).Res1();
          int bp_2 = (*BP).Res2();
          // FIXME: Hack for integer
          int n_of_hb = (int)BPHBONDS_[nbp]->Dval(frame);
          outfile.Printf(BP_OUTPUT_FMT, frame+1, 
                         Bases_[bp_1].ResNum()+1, Bases_[bp_2].ResNum()+1,
                         SHEAR_[nbp]->Dval(frame), STRETCH_[nbp]->Dval(frame),
                         STAGGER_[nbp]->Dval(frame), BUCKLE_[nbp]->Dval(frame),
                         PROPELLER_[nbp]->Dval(frame), OPENING_[nbp]->Dval( frame),
                         n_of_hb, MAJOR_[nbp]->Dval(frame), MINOR_[nbp]->Dval(frame));
          ++nbp;
        }
        outfile.Printf("\n");
      }
      outfile.CloseFile();
    } else {
      mprinterr("Error: Could not open %s for writing.\n", outfilename.c_str());
    }
  }

  // ---------- Base pair step parameters ----------
  CpptrajFile outfile2;
  outfilename = "BPstep." + outputsuffix_;
  std::string outfilename2 = "Helix." + outputsuffix_;
  // Check that there is actually data
  // TODO: Check helix data as well
  if ( SHIFT_.empty() || SHIFT_[0]->Empty() )
    mprinterr("Error: Could not write BPstep / helix files: No data.\n"); 
  else {
    int err = 0;
    err += outfile.OpenEnsembleWrite( outfilename, ensembleNum_ );
    err += outfile2.OpenEnsembleWrite( outfilename2, ensembleNum_ );
    if (err == 0) {
      // Determine number of frames from SHIFT[0] DataSet. Should be same as SHEAR.
      nframes = SHIFT_[0]->Size();
      mprintf("\tBase pair step output file %s; ",outfilename.c_str());
      mprintf("Helix output file %s; %i frames, %zu base pair steps.\n", outfilename2.c_str(),
              nframes, BasePairAxes_.size() - 1);
      //  File headers
      if (printheader_) {
        outfile.Printf("%-8s %-9s %-9s %10s %10s %10s %10s %10s %10s\n","#Frame","BP1","BP2",
                       "Shift","Slide","Rise","Tilt","Roll","Twist");
        outfile2.Printf("%-8s %-9s %-9s %10s %10s %10s %10s %10s %10s\n","#Frame","BP1","BP2",
                        "X-disp","Y-disp","Rise","Incl.","Tip","Twist");
      }
      // Loop over all frames
      for (int frame = 0; frame < nframes; ++frame) {
        int nstep = 0;
        std::vector<NA_Axis>::iterator NBPstep = BasePairAxes_.end() - 1;
        for (std::vector<NA_Axis>::iterator BP1 = BasePairAxes_.begin();
                                            BP1 != NBPstep; ++BP1)
        {
          std::vector<NA_Axis>::iterator BP2 = BP1 + 1;
          int bp1_1 = Bases_[(*BP1).Res1()].ResNum() + 1;
          int bp1_2 = Bases_[(*BP1).Res2()].ResNum() + 1;
          int bp2_1 = Bases_[(*BP2).Res1()].ResNum() + 1;
          int bp2_2 = Bases_[(*BP2).Res2()].ResNum() + 1;
          // BPstep write
          outfile.Printf(NA_OUTPUT_FMT, frame+1, 
                         bp1_1, bp1_2, bp2_1, bp2_2,
                         SHIFT_[nstep]->Dval(frame), SLIDE_[nstep]->Dval(frame),
                         RISE_[nstep]->Dval(frame), TILT_[nstep]->Dval(frame),
                         ROLL_[nstep]->Dval(frame), TWIST_[nstep]->Dval(frame));
          // Helix write
          outfile2.Printf(NA_OUTPUT_FMT, frame+1,
                          bp1_1, bp1_2, bp2_1, bp2_2,
                          XDISP_[nstep]->Dval(frame), YDISP_[nstep]->Dval(frame),
                          HRISE_[nstep]->Dval(frame), INCL_[nstep]->Dval(frame),
                          TIP_[nstep]->Dval(frame), HTWIST_[nstep]->Dval(frame));
          ++nstep;
        }
        outfile.Printf("\n");
        outfile2.Printf("\n");
      }
      outfile.CloseFile();
      outfile2.CloseFile();
    } else {
      mprinterr("Error: Could not open BPstep files for writing.\n");
    }
  }
}

