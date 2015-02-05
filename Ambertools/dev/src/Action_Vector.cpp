#include <cmath>
#include <algorithm> // sort
#include "Constants.h" // PI
#include "Action_Vector.h"
#include "CpptrajStdio.h"
#include "DistRoutines.h" // MinImagedVec, includes Matrix_3x3 for principal

// CONSTRUCTOR
Action_Vector::Action_Vector() :
  ensembleNum_(-1),
  Vec_(0),
  Magnitude_(0),
  vcorr_(0),
  ptrajoutput_(false),
  needBoxInfo_(false),
  CurrentParm_(0)
{}

void Action_Vector::Help() {
  mprintf("\t[<name>] <Type> [out <filename> [ptrajoutput]] [<mask1>] [<mask2>]\n"
          "\t[magnitude] [ired]\n"
          "\t<Type> = { mask | principal [x|y|z] | dipole | box | center | corrplane |\n"
          "             ucellx | ucelly | ucellz }\n"
          "  Calculate the specified coordinate vector.\n"
          "    mask: (Default) Vector from <mask1> to <mask2>.\n"
          "    principal [x|y|z]: X, Y, or Z principal axis vector for atoms in <mask1>.\n"
          "    dipole: Dipole and center of mass of the atoms specified in <mask1>\n"
          "    box: (No mask needed) Store the box lengths of the trajectory.\n"
          "    center: Store the center of mass of atoms in <mask1>.\n"
          "    corrplane: Vector perpendicular to plane through the atoms in <mask1>.\n"
          "    ucell{x|y|z}: (No mask needed) Store specified unit cell vector.\n");
}

// DESTRUCTOR
Action_Vector::~Action_Vector() {
  if (vcorr_!=0) delete[] vcorr_;
}

const char* Action_Vector::ModeString[] = {
  "NO_OP", "Principal X", "Principal Y", "Principal Z",
  "Dipole", "Box", "Mask", "Ired",
  "CorrPlane", "Center", "Unit cell X", "Unit cell Y", "Unit cell Z",
  "Box Center", "MinImage"
};

static Action::RetType WarnDeprecated() {
  mprinterr("Error: Vector: 'corrired' and 'corr' are deprecated.\n"
            "Error: 'corrired' functionality is now part of the\n"
            "Error: IRED analysis. 'corr' can now be done with a normal 2-mask\n"
            "Error: vector and TIMECORR analysis.\n");
  return Action::ERR;
}

static inline Action::RetType DeprecatedErr(const char* key) {
  mprinterr("Error: '%s' is deprecated. To write vector pseudo-traj use:\n"
            "\twritedata <filename> <vector set> [trajfmt <format>] [parmout <file>]\n", key);
  return Action::ERR;
}

// Action_Vector::Init()
Action::RetType Action_Vector::Init(ArgList& actionArgs, TopologyList* PFL, DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  ensembleNum_ = DSL->EnsembleNum();
  // filename is saved in case ptraj-compatible output is desired
  filename_ = actionArgs.GetStringKey("out");
  if (actionArgs.hasKey("trajout")) return DeprecatedErr( "trajout" );
  if (actionArgs.hasKey("trajfmt")) return DeprecatedErr( "trajfmt" );
  if (actionArgs.hasKey("parmout")) return DeprecatedErr( "parmout" );
  ptrajoutput_ = actionArgs.hasKey("ptrajoutput");
  if (ptrajoutput_ && filename_.empty()) {
    mprinterr("Error: 'ptrajoutput' specified but no 'out <filename>' arg given.\n");
    return Action::ERR;
  }
  bool calc_magnitude = actionArgs.hasKey("magnitude");
  if (calc_magnitude && ptrajoutput_) {
    mprinterr("Error: 'ptrajoutput' and 'magnitude' are incompatible.\n");
    return Action::ERR;
  }
  needBoxInfo_ = false;
  // Deprecated: corrired, corr, ired
  if ( actionArgs.hasKey("principal") ) {
    mode_ = PRINCIPAL_X;
    if ( actionArgs.hasKey("x") ) mode_ = PRINCIPAL_X;
    if ( actionArgs.hasKey("y") ) mode_ = PRINCIPAL_Y;
    if ( actionArgs.hasKey("z") ) mode_ = PRINCIPAL_Z;
  } else if (actionArgs.hasKey("center"))
    mode_ = CENTER;
  else if (actionArgs.hasKey("dipole"))
    mode_ = DIPOLE;
  else if (actionArgs.hasKey("box"))
    mode_ = BOX;
  else if (actionArgs.hasKey("corrplane"))
    mode_ = CORRPLANE;
  else if (actionArgs.hasKey("corrired"))
    return WarnDeprecated();
  else if (actionArgs.hasKey("corr"))
    return WarnDeprecated();
  else if (actionArgs.hasKey("mask"))
    mode_ = MASK;
  else if (actionArgs.hasKey("ucellx"))
    mode_ = BOX_X;
  else if (actionArgs.hasKey("ucelly"))
    mode_ = BOX_Y;
  else if (actionArgs.hasKey("ucellz"))
    mode_ = BOX_Z;
  else if (actionArgs.hasKey("boxcenter"))
    mode_ = BOX_CTR;
  else if (actionArgs.hasKey("minimage"))
    mode_ = MINIMAGE; 
  else
    mode_ = MASK;
  if (mode_ == BOX || mode_ == BOX_X || mode_ == BOX_Y || mode_ == BOX_Z ||
      mode_ == BOX_CTR || mode_ == MINIMAGE)
    needBoxInfo_ = true;
  // Check if IRED vector
  bool isIred = actionArgs.hasKey("ired"); 
  // Vector Mask
  if (mode_ != BOX && mode_ != BOX_X && mode_ != BOX_Y && mode_ != BOX_Z)
    mask_.SetMaskString( actionArgs.GetMaskNext() );
  // Get second mask if necessary
  if (mode_ == MASK || mode_ == MINIMAGE) {
    std::string maskexpr = actionArgs.GetMaskNext();
    if (maskexpr.empty()) {
      mprinterr("Error: Specified vector mode (%s) requires a second mask.\n",
                ModeString[ mode_ ]);
      return Action::ERR;
    }
    mask2_.SetMaskString( maskexpr );
  }
  // Set up vector dataset and IRED status
  Vec_ = (DataSet_Vector*)DSL->AddSet(DataSet::VECTOR, actionArgs.GetStringNext(), "Vec");
  if (Vec_ == 0) return Action::ERR;
  if (isIred)
    Vec_->SetIred( );
  // Add set to output file if not doing ptraj-compatible output
  if (!ptrajoutput_)
    DFL->AddSetToFile(filename_, Vec_);
  // Set up magnitude data set.
  if (calc_magnitude) {
    Magnitude_ = DSL->AddSetAspect(DataSet::FLOAT, Vec_->Name(), "Mag");
    if (Magnitude_ == 0) return Action::ERR;
    DFL->AddSetToFile(filename_, Magnitude_);
  }
  
  mprintf("    VECTOR: Type %s", ModeString[ mode_ ]);
  if (calc_magnitude)
    mprintf(" (with magnitude)");
  if (isIred)
    mprintf(", IRED");
  if (mask_.MaskStringSet())
     mprintf(", mask [%s]", mask_.MaskString());
  if (mask2_.MaskStringSet())
    mprintf(", second mask [%s]", mask2_.MaskString());
  if (!filename_.empty()) {
    if (ptrajoutput_)
      mprintf(", ptraj-compatible output to");
    else
      mprintf(", output to");
    mprintf(" %s", filename_.c_str());
  }
  mprintf("\n");

  return Action::OK;
}

// Action_Vector::Setup()
Action::RetType Action_Vector::Setup(Topology* currentParm, Topology** parmAddress) {
  if (needBoxInfo_) {
    // Check for box info
    if (currentParm->BoxType() == Box::NOBOX) {
      mprinterr("Error: vector box: Parm %s does not have box information.\n",
                currentParm->c_str());
      return Action::ERR;
    }
  }
  if (mask_.MaskStringSet()) {
    // Setup mask 1
    if (currentParm->SetupIntegerMask(mask_)) return Action::ERR;
    mprintf("\tVector mask [%s] corresponds to %i atoms.\n",
            mask_.MaskString(), mask_.Nselected());
  }

  // Allocate space for CORRPLANE. 
  if (mode_ == CORRPLANE) {
    if (vcorr_!=0) delete[] vcorr_;
    vcorr_ = new double[ 3 * mask_.Nselected() ];
  }

  // Setup mask 2
  if (mask2_.MaskStringSet()) {
    if (currentParm->SetupIntegerMask(mask2_)) return Action::ERR;
    mprintf("\tVector mask [%s] corresponds to %i atoms.\n",
            mask2_.MaskString(), mask2_.Nselected());
  }
  CurrentParm_ = currentParm;
  return Action::OK;
}

// -----------------------------------------------------------------------------
// Action_Vector::solve_cubic_eq()
/** Solves a cubic equation: ax^3 + bx^2 + cx + d = 0 using "Cardan's formula" 
  * (see: Bronstein, S.131f)
  */
double Action_Vector::solve_cubic_eq(double a, double b, double c, double d) {
  const double one3 = 1.0 / 3.0;
  const double one27 = 1.0 / 27.0;
  double droot = 0;

  double r, s, t;
  double p, q, rho, phi;
  double D, u, v;
  std::vector<double> dtmp(3);

  /* Coeff. for normal form x^3 + rx^2 + sx + t = 0 */
  r = b / a;
  s = c / a;
  t = d / a;

  /* Coeff. for red. eq. y^3 + py + q = 0 with y = x + r/3 bzw. (x = y - r/3) */
  p = s - r * r * one3;
  q = 2.0 * r * r * r * one27 - r * s * one3 + t;

  /* Dummy variables */
  rho = sqrt(-p * p * p * one27);
  phi = acos(-q / (2.0 * rho));

  /* Discriminante(?) */
  D = pow((p * one3),3) + q * q * 0.25;

  if(D > 0){ /* x real -> one real solution */
    u = pow(-q * 0.5 + sqrt(D), one3);
    v = -p / u * one3;
    droot = (u + v) - r * one3;
  } else if(D <= 0){
  /* three real solutions (d < 0) | one real solution + one real double solution or 
                                                     one real triple solution (d = 0) */
    dtmp[0] = 2.0 * pow(rho, one3) * cos(phi * one3) - r * one3;
    dtmp[1] = 2.0 * pow(rho, one3) * cos((phi + Constants::TWOPI ) * one3) - r * one3;
    dtmp[2] = 2.0 * pow(rho, one3) * cos((phi + Constants::FOURPI) * one3) - r * one3;

    sort(dtmp.begin(), dtmp.end());

    //qsort((void *) dtmp, (size_t) 3, sizeof(double), cmpdouble);
    droot = dtmp[0];
  }
  return droot;
}

// Action_Vector::leastSquaresPlane()
/** Calcs (least-squares best) plane through a series of points
  * relative to their center of geom. (the latter has to be done outside this routine), 
  * returns (normalized) coeff. for plane eq. ax + by + cz = 0
  * following: Crystal Structure Analysis for Chem. and Biol.,
  * Glusker, Lewis, Rossi, S. 460ff
  */
Vec3 Action_Vector::leastSquaresPlane(int n, const double* vcorr) {
  double Xout, Yout, Zout;
  if (n == 9) {  // Special case, only 3 coords
    double x1 = vcorr[3] - vcorr[0];
    double y1 = vcorr[4] - vcorr[1];
    double z1 = vcorr[5] - vcorr[2];
    double x2 = vcorr[6] - vcorr[3];
    double y2 = vcorr[7] - vcorr[4];
    double z2 = vcorr[8] - vcorr[5];

    Xout = y1 * z2 - z1 * y2;
    Yout = z1 * x2 - x1 * z2;
    Zout = x1 * y2 - y1 * x2;
  } else { // General case
    // Calc Var.
    double dSumXX = 0.0;
    double dSumYY = 0.0;
    double dSumZZ = 0.0;
    double dSumXY = 0.0;
    double dSumXZ = 0.0;
    double dSumYZ = 0.0;

    for (int i = 0; i < n; i+=3) {
      dSumXX += vcorr[i  ] * vcorr[i  ];
      dSumYY += vcorr[i+1] * vcorr[i+1];
      dSumZZ += vcorr[i+2] * vcorr[i+2];

      dSumXY += vcorr[i  ] * vcorr[i+1];
      dSumXZ += vcorr[i  ] * vcorr[i+2];
      dSumYZ += vcorr[i+1] * vcorr[i+2];
    }

    // Calc coeff. for -l^3 + o * l^2 + p * l + q = 0
    double o = dSumXX + dSumYY + dSumZZ;
    double p = pow(dSumXY,2) + pow(dSumXZ,2) + pow(dSumYZ,2) -
               (dSumXX * dSumYY + dSumXX * dSumZZ + dSumYY * dSumZZ);
    double q = dSumXX * dSumYY * dSumZZ + 2.0 * dSumXY * dSumXZ * dSumYZ -
               (dSumXX * dSumYZ * dSumYZ + dSumYY * dSumXZ * dSumXZ + dSumZZ * dSumXY * dSumXY);

    // Solve cubic eq.
    double root = Action_Vector::solve_cubic_eq(-1.0, o, p, q);

    // Calc determinents
    Xout = (dSumYY - root) * dSumXZ - dSumXY * dSumYZ;
    Yout = (dSumXX - root) * dSumYZ - dSumXY * dSumXZ;
    Zout =  dSumXY         * dSumXY - (dSumYY - root) * (dSumXX - root);
  }
  // Normalize
  double dnorm = 1.0 / sqrt((Xout * Xout) + (Yout * Yout) + (Zout * Zout));
  return Vec3(Xout * dnorm, Yout * dnorm, Zout * dnorm);
}

// -----------------------------------------------------------------------------
void Action_Vector::Mask(Frame const& currentFrame) {
  Vec3 CXYZ = currentFrame.VCenterOfMass(mask_);
  Vec3 VXYZ = currentFrame.VCenterOfMass(mask2_);
  VXYZ -= CXYZ;
  Vec_->AddVxyz(VXYZ, CXYZ);
}

void Action_Vector::Dipole(Frame const& currentFrame) {
  Vec3 VXYZ(0.0, 0.0, 0.0);
  Vec3 CXYZ(0.0, 0.0, 0.0);
  double total_mass = 0;
  for (AtomMask::const_iterator atom = mask_.begin();
                                atom != mask_.end(); ++atom)
  {
    double mass = (*CurrentParm_)[*atom].Mass();
    total_mass += mass;
    Vec3 XYZ = currentFrame.XYZ( *atom );
    CXYZ += ( XYZ * mass );
    double charge = (*CurrentParm_)[*atom].Charge();
    XYZ *= charge;
    VXYZ += ( XYZ );
  }
  CXYZ /= total_mass;
  Vec_->AddVxyz( VXYZ, CXYZ );
}

void Action_Vector::Principal(Frame const& currentFrame) {
  Matrix_3x3 Inertia;
  Vec3 Eval;

  // Origin is center of atoms in mask_ 
  Vec3 OXYZ = currentFrame.CalculateInertia( mask_, Inertia );
  // NOTE: Diagonalize_Sort_Chirality places sorted eigenvectors in rows.
  Inertia.Diagonalize_Sort_Chirality( Eval, 0 );
  // Eval.Print("PRINCIPAL EIGENVALUES");
  // Inertia.Print("PRINCIPAL EIGENVECTORS (Rows)");
  if ( mode_ == PRINCIPAL_X ) 
    Vec_->AddVxyz( Inertia.Row1(), OXYZ ); // First row = first eigenvector
  else if ( mode_ == PRINCIPAL_Y )
    Vec_->AddVxyz( Inertia.Row2(), OXYZ ); // Second row = second eigenvector
  else // PRINCIPAL_Z
    Vec_->AddVxyz( Inertia.Row3(), OXYZ ); // Third row = third eigenvector
}

void Action_Vector::CorrPlane(Frame const& currentFrame) {
  Vec3 CXYZ = currentFrame.VCenterOfMass(mask_);
  int idx = 0;
  for (AtomMask::const_iterator atom = mask_.begin();
                              atom != mask_.end(); ++atom)
  {
    Vec3 XYZ = currentFrame.XYZ( *atom );
    XYZ -= CXYZ;
    vcorr_[idx++] = XYZ[0];
    vcorr_[idx++] = XYZ[1];
    vcorr_[idx++] = XYZ[2];
  }
  Vec3 VXYZ = leastSquaresPlane(idx, vcorr_);
  Vec_->AddVxyz(VXYZ, CXYZ);
}

void Action_Vector::UnitCell(Box const& box) {
  Matrix_3x3 ucell, recip;
  box.ToRecip( ucell, recip );
  switch ( mode_ ) {
    case BOX_X: Vec_->AddVxyz( ucell.Row1() ); break;
    case BOX_Y: Vec_->AddVxyz( ucell.Row2() ); break;
    case BOX_Z: Vec_->AddVxyz( ucell.Row3() ); break;
    case BOX_CTR: Vec_->AddVxyz( ucell.TransposeMult(Vec3(0.5)) ); break;
    default: return;
  }
}

void Action_Vector::MinImage(Frame const& frm) {
  Matrix_3x3 ucell, recip;
  frm.BoxCrd().ToRecip( ucell, recip );
  Vec3 com1 = frm.VCenterOfMass(mask_);
  Vec_->AddVxyz( MinImagedVec(com1, frm.VCenterOfMass(mask2_), ucell, recip), com1 );
}
                 
// Action_Vector::DoAction()
Action::RetType Action_Vector::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {
  switch ( mode_ ) {
    case MASK        : Mask(*currentFrame); break;
    case CENTER      : Vec_->AddVxyz( currentFrame->VCenterOfMass(mask_) ); break; 
    case DIPOLE      : Dipole(*currentFrame); break;
    case PRINCIPAL_X :
    case PRINCIPAL_Y :
    case PRINCIPAL_Z : Principal(*currentFrame); break;
    case CORRPLANE   : CorrPlane(*currentFrame); break;
    case BOX         : Vec_->AddVxyz( currentFrame->BoxCrd().Lengths() ); break;
    case BOX_X       : 
    case BOX_Y       : 
    case BOX_Z       : 
    case BOX_CTR     : UnitCell( currentFrame->BoxCrd() ); break;
    case MINIMAGE    : MinImage( *currentFrame ); break; 
    default          : return Action::ERR; // NO_OP
  } // END switch over vectorMode
  if (Magnitude_ != 0) {
    float mag = (float)(sqrt(Vec_->Back().Magnitude2()));
    Magnitude_->Add(frameNum, &mag);
  }
  return Action::OK;
}

// Action_Vector::Print()
void Action_Vector::Print() {
  if (ptrajoutput_) {
    CpptrajFile outfile;
    if (outfile.OpenEnsembleWrite(filename_, ensembleNum_)) return;
    mprintf("    VECTOR: writing ptraj-style vector information for %s\n", Vec_->Legend().c_str());
    outfile.Printf("# FORMAT: frame vx vy vz cx cy cz cx+vx cy+vy cz+vz\n"
                   "# FORMAT where v? is vector, c? is center of mass...\n");
    int totalFrames = Vec_->Size();
    for (int i=0; i < totalFrames; ++i) {
      Vec3 const& vxyz = (*Vec_)[i];
      Vec3 const& cxyz = Vec_->OXYZ(i);
      Vec3 txyz  = cxyz + vxyz;
      outfile.Printf("%i %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n",
              i+1, vxyz[0], vxyz[1], vxyz[2], cxyz[0], cxyz[1], cxyz[2],
              txyz[0], txyz[1], txyz[2]);
    }
    outfile.CloseFile();
  }
}
