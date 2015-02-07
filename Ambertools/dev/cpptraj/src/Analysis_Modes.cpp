#include <cmath> // tanh, sqrt
#include "Analysis_Modes.h"
#include "CpptrajStdio.h"
#include "Constants.h" // TWOPI

// CONSTRUCTOR
Analysis_Modes::Analysis_Modes() :
  debug_(0),
  type_(FLUCT),
  beg_(0),
  end_(0),
  bose_(false),
  factor_(0),
  modinfo_(0),
  modinfo2_(0),
  results_(0),
  tOutParm_(0),
  tMode_(0),
  pcmin_(0.0),
  pcmax_(0.0)
{}

void Analysis_Modes::Help() {
  mprintf("\t{fluct|displ|corr|eigenval|trajout|rmsip} name <modesname> [name2 <modesname>]\n"
          "\t[beg <beg>] [end <end>] [bose] [factor <factor>]\n"
          "\t[out <outfile>] [maskp <mask1> <mask2> [...]]\n"
          "    Options for 'trajout': (Generate pseudo-trajectory)\n"
          "\t[trajout <name> [<parm arg>] [trajoutfmt <format>] [trajoutmask <mask>]\n"
          "\t  [pcmin <pcmin>] [pcmax <pcmax>] [tmode <mode>]]\n"
          "  Perform one of the following analysis on calculated Eigenmodes.\n"
          "    fluct: rms fluctations from normal modes\n"
          "    displ: displacement of cartesian coordinates along normal mode directions\n"
          "    eigenval: Calculate eigenvalue fractions.\n"
          "    rmsip: Root mean square inner product.\n"
          "  Results vector usage:\n"
          "    fluct:\n"
          "\t[rmsx(at1), rmsy(at1), rmsz(at1), rms(at1), ..., rmsx(atN), ..., rms(atN)]\n"
          "    displ:\n"
          "\t[displx(at1), disply(at1), displz(at1), ..., displx(atN), ..., displz(atN)]\n"
          "    corr:\n"
          "\t[corr(pair1, vec1), ..., corr(pair1, vecN), ..., corr(pairM, vec1), ..., corr(pairM, vecN)\n");
}

/// hc/2kT in cm, with T=300K; use for quantum Bose statistics)
const double Analysis_Modes::CONSQ = 2.39805E-3;
/// kT/c^2 in cgs units (grams), with T=300K
const double Analysis_Modes::TKBC2 = 0.46105E-34;
/// Avogadros number
const double Analysis_Modes::AVO   = 6.023E23;
const double Analysis_Modes::CNST  = TKBC2 * AVO;
// cm to angstroms
const double Analysis_Modes::CMTOA = 1.000E8;
const double Analysis_Modes::CONT  = CMTOA / Constants::TWOPI;
// NOTE: Original TWOPI was 6.2832, results in small roundoff diffs from ptraj

/// Strings describing modes calculation types.
const char* Analysis_Modes::analysisTypeString[] = {
  "rms fluctuations",
  "displacements",
  "correlation functions",
  "coordinate projection",
  "eigenvalue fraction",
  "root mean square inner product"
};

// DESTRUCTOR
Analysis_Modes::~Analysis_Modes() {
  if (results_!=0)
    delete[] results_;
  if (tOutParm_ != 0)
    delete tOutParm_;
}

// Analysis_Modes::CheckDeprecated()
void Analysis_Modes::CheckDeprecated(ArgList& analyzeArgs, std::string& modesname, 
                                     const char* key) {
  std::string arg = analyzeArgs.GetStringKey( key );
  if (!arg.empty()) {
    mprintf("Warning: Argument '%s' is deprecated, use 'name <modes>' instead.\n",
            key);
    if (modesname.empty()) modesname = arg;
  }
}

// Analysis_Modes::Setup()
Analysis::RetType Analysis_Modes::Setup(ArgList& analyzeArgs, DataSetList* DSLin,
                            TopologyList* PFLin, DataFileList* DFLin, int debugIn)
{
  debug_ = debugIn;
  // Analysis type
  if (analyzeArgs.hasKey("fluct"))
    type_ = FLUCT;
  else if (analyzeArgs.hasKey("displ"))
    type_ = DISPLACE;
  else if (analyzeArgs.hasKey("corr"))
    type_ = CORR;
  else if (analyzeArgs.Contains("trajout"))
    type_ = TRAJ;
  else if (analyzeArgs.hasKey("eigenval"))
    type_ = EIGENVAL;
  else if (analyzeArgs.hasKey("rmsip"))
    type_ = RMSIP;
  else {
    mprinterr("Error: No analysis type specified.\n");
    return Analysis::ERR;
  }

  // Get modes name
  std::string modesfile = analyzeArgs.GetStringKey("name");
  if (modesfile.empty()) {
    // Check for deprecated args
    CheckDeprecated(analyzeArgs, modesfile, "file");
    CheckDeprecated(analyzeArgs, modesfile, "stack");
    if (modesfile.empty()) {
      mprinterr("Error: No 'name <modes data set name>' argument given.\n");
      return Analysis::ERR;
    }
  }
  // Get second modes name for RMSIP
  std::string modesfile2 = analyzeArgs.GetStringKey("name2");
  if (type_ == RMSIP) {
    if (modesfile2.empty()) {
      mprinterr("Error: 'rmsip' requires second modes data 'name2 <modes>'\n");
      return Analysis::ERR;
    }
  } else
    modesfile2.clear(); 

  // Get trajectory format args for projected traj
  if (type_ == TRAJ ) {
    beg_ = analyzeArgs.getKeyInt("beg",1) - 1; // Args start at 1
    std::string tOutName = analyzeArgs.GetStringKey("trajout");
    if (tOutName.empty()) {
      mprinterr("Error: Require output trajectory filename, 'trajout <name>'\n");
      return Analysis::ERR;
    }
    TrajectoryFile::TrajFormatType tOutFmt = 
      TrajectoryFile::GetFormatFromString( analyzeArgs.GetStringKey("trajoutfmt") );
    Topology* parm = PFLin->GetParm( analyzeArgs );
    if (parm == 0) {
      mprinterr("Error: Could not get topology for output trajectory.\n");
      return Analysis::ERR;
    }
    AtomMask tOutMask( analyzeArgs.GetStringKey("trajoutmask") );
    if ( parm->SetupIntegerMask( tOutMask ) || tOutMask.None() ) {
      mprinterr("Error: Could not setup output trajectory mask.\n");
      return Analysis::ERR;
    }
    tOutMask.MaskInfo();
    // Strip topology to match mask.
    if (tOutParm_ != 0) delete tOutParm_;
    tOutParm_ = parm->modifyStateByMask( tOutMask );
    if (tOutParm_ == 0) {
      mprinterr("Error: Could not create topology to match mask.\n");
      return Analysis::ERR;
    }
    // Setup output traj
    if (trajout_.InitTrajWrite( tOutName, tOutParm_, tOutFmt ) != 0) {
      mprinterr("Error: Could not setup output trajectory.\n");
      return Analysis::ERR;
    }
    // Get min and max for PC
    pcmin_ = analyzeArgs.getKeyDouble("pcmin", -10.0);
    pcmax_ = analyzeArgs.getKeyDouble("pcmax",  10.0);
    if (pcmax_ < pcmin_ || pcmax_ - pcmin_ < Constants::SMALL) {
      mprinterr("Error: pcmin must be less than pcmax\n");
      return Analysis::ERR;
    }
    tMode_ = analyzeArgs.getKeyInt("tmode", 1);
  } else {
    beg_ = analyzeArgs.getKeyInt("beg",7) - 1; // Args start at 1
    // Get factor, bose
    bose_ = analyzeArgs.hasKey("bose");
    factor_ = analyzeArgs.getKeyDouble("factor",1.0);
  }
  end_ = analyzeArgs.getKeyInt("end", 50);

  // Check if modes name exists on the stack
  modinfo_ = (DataSet_Modes*)DSLin->FindSetOfType( modesfile, DataSet::MODES );
  if (modinfo_ == 0) {
    mprinterr("Error: '%s' not found: %s\n", modesfile.c_str(), DataSet_Modes::DeprecateFileMsg);
    return Analysis::ERR;
  }
  if (!modesfile2.empty()) {
    modinfo2_ = (DataSet_Modes*)DSLin->FindSetOfType( modesfile2, DataSet::MODES );
    if (modinfo2_ == 0) {
      mprinterr("Error: Set %s not found.\n", modesfile2.c_str());
      return Analysis::ERR;
    }
  }

  // Check modes type for specified analysis
  if (type_ == FLUCT || type_ == DISPLACE || type_ == CORR || type_ == TRAJ) {
    if (modinfo_->Type() != DataSet_2D::COVAR && 
        modinfo_->Type() != DataSet_2D::MWCOVAR)
    {
      mprinterr("Error: Modes must be of type COVAR or MWCOVAR for %s.\n",
                analysisTypeString[type_]);
      return Analysis::ERR;
    }
  }

  // Get output filename
  filename_ = analyzeArgs.GetStringKey("out");

  // Get mask pair info for ANALYZEMODES_CORR option and build the atom pair stack
  if ( type_ == CORR ) {
    Topology* analyzeParm = PFLin->GetParm( analyzeArgs );
    if (analyzeParm == 0) {
      mprinterr("Error: 'corr' requires topology (parm <file>, parmindex <#>).\n");
      return Analysis::ERR;
    }
    while (analyzeArgs.hasKey("maskp")) {
      // Next two arguments should be one-atom masks
      std::string a1mask = analyzeArgs.GetMaskNext();
      std::string a2mask = analyzeArgs.GetMaskNext();
      if (a1mask.empty() || a2mask.empty()) {
        mprinterr("Error: For 'corr' two 1-atom masks are expected.\n");
        return Analysis::ERR;
      }
      // Check that each mask is just 1 atom
      AtomMask m1( a1mask );
      AtomMask m2( a2mask );
      analyzeParm->SetupIntegerMask( m1 ); 
      analyzeParm->SetupIntegerMask( m2 );
      if ( m1.Nselected()==1 && m2.Nselected()==1 )
        // Store atom pair
        atompairStack_.push_back( std::pair<int,int>( m1[0], m2[0] ) );
      else {
        mprinterr("Error: For 'corr', masks should specify only one atom.\n"
                  "\tM1[%s]=%i atoms, M2[%s]=%i atoms.\n", m1.MaskString(), m1.Nselected(),
                  m2.MaskString(), m2.Nselected());
        return Analysis::ERR;
      }
    }
    if ( atompairStack_.empty() ) {
      mprinterr("Error: No atom pairs found (use 'maskp <mask1> <mask2>')\n");
      return Analysis::ERR;
    }
  }

  // Status
  mprintf("    ANALYZE MODES: Calculating %s using modes from %s", 
          analysisTypeString[type_], modinfo_->Legend().c_str());
  if ( type_ != TRAJ ) {
    if (type_ != EIGENVAL)
      mprintf(", modes %i to %i", beg_+1, end_);
    mprintf("\n\tResults are written to");
    if (filename_.empty())
      mprintf(" STDOUT\n");
    else
      mprintf(" %s\n", filename_.c_str());
    if (type_ != EIGENVAL && type_ != RMSIP) {
      if (bose_)
        mprintf("\tBose statistics used.\n");
      else
        mprintf("\tBoltzmann statistics used.\n");
    }
    if (type_ == DISPLACE)
      mprintf("\tFactor for displacement: %lf\n", factor_);
    if (type_ == CORR) {
      mprintf("\tUsing the following atom pairs:");
      for (modestack_it apair = atompairStack_.begin();
                        apair != atompairStack_.end(); ++apair)
        mprintf(" (%i,%i)", (*apair).first+1, (*apair).second+1 );
      mprintf("\n");
    }
    if (type_ == RMSIP)
      mprintf("\tRMSIP calculated to modes in %s\n", modinfo2_->Legend().c_str());
  } else {
    mprintf("\n\tCreating trajectory for mode %i\n"
              "\tWriting to trajectory %s\n"
              "\tPC range: %f to %f\n", tMode_, 
            trajout_.TrajFilename().full(), pcmin_, pcmax_);
  }

  return Analysis::OK;
}

// Analysis_Modes::Analyze()
Analysis::RetType Analysis_Modes::Analyze() {
  CpptrajFile outfile;
  // Check # of modes
  if (type_ != TRAJ && type_ != EIGENVAL) {
    if (beg_ < 0 || beg_ >= modinfo_->Nmodes()) {
      mprinterr("Error: 'beg %i' is out of bounds.\n", beg_+1);
      return Analysis::ERR;
    }
    if (end_ > modinfo_->Nmodes()) {
      mprintf("Warning: 'end %i' is > # of modes, setting to %i\n",
              end_, modinfo_->Nmodes());
      end_ = modinfo_->Nmodes();
    }
    if (end_ <= beg_) {
      mprinterr("Warning: beg must be <= end, (%i -- %i)\n", beg_+1, end_);
      return Analysis::ERR;
    }
  }

  // ----- FLUCT PRINT -----
  if (type_ == FLUCT) {
    // Calc rms atomic fluctuations
    int natoms = modinfo_->NavgCrd() / 3; // Possible because COVAR/MWCOVAR
    results_ = new double[ natoms * 4 ];
    std::fill(results_, results_ + natoms * 4, 0);
    double* Ri = results_;
    for (int vi = 0; vi < natoms; ++vi) {
      double sumx = 0.0;
      double sumy = 0.0;
      double sumz = 0.0;
      // Loop over all eigenvector elements for this atom
      const double* Vec = modinfo_->Eigenvector(beg_) + vi*3;
      for (int mode = beg_; mode < end_; ++mode) {
        double frq = modinfo_->Eigenvalue(mode);
        if (frq >= 0.5) {
          // Don't use eigenvectors associated with zero or negative eigenvalues
          double distx = Vec[0] * Vec[0];
          double disty = Vec[1] * Vec[1];
          double distz = Vec[2] * Vec[2];
          double fi = 1.0 / (frq*frq);
          if (bose_) {
            double argq = CONSQ * frq;
            fi *= (argq / tanh(argq));
          }
          sumx += distx * fi;
          sumy += disty * fi;
          sumz += distz * fi;
        }
        Vec += modinfo_->VectorSize();
      }
      sumx *= CNST;
      sumy *= CNST;
      sumz *= CNST;
      Ri[0] = sqrt(sumx) * CONT;
      Ri[1] = sqrt(sumy) * CONT;
      Ri[2] = sqrt(sumz) * CONT;
      Ri[3] = sqrt(sumx + sumy + sumz) * CONT;
      Ri += 4;
    }
    // Output
    if (outfile.OpenWrite( filename_ )) return Analysis::ERR;
    outfile.Printf("#Analysis of modes: RMS FLUCTUATIONS\n");
    outfile.Printf("%-10s %10s %10s %10s %10s\n", "#Atom_no.", "rmsX", "rmsY", "rmsZ", "rms");
    int anum = 1;
    for (int i4 = 0; i4 < modinfo_->NavgCrd()*4/3; i4+=4) 
      outfile.Printf("%10i %10.3f %10.3f %10.3f %10.3f\n", anum++, results_[i4], 
                     results_[i4+1], results_[i4+2], results_[i4+3]); 
  // ----- DISPLACE PRINT -----
  } else if (type_ == DISPLACE) {
    // Calc displacement of coordinates along normal mode directions
    results_ = new double[ modinfo_->NavgCrd() ];
    std::fill(results_, results_ + modinfo_->NavgCrd(), 0);
    double sqrtcnst = sqrt(CNST) * CONT * factor_;
    // Loop over all modes
    for (int mode = beg_; mode < end_; ++mode) {
      double frq = modinfo_->Eigenvalue(mode);
      if (frq >= 0.5) {
        // Don't use eigenvectors associated with zero or negative eigenvalues
        double fi = 1.0 / frq;
        if (bose_) {
          double argq = CONSQ * frq;
          fi *= (fi * argq / tanh(argq));
          fi = sqrt(fi);
        }
        fi *= sqrtcnst; // * CONT * factor
        // Loop over all vector elements
        const double* Vec = modinfo_->Eigenvector(mode);
        for (int vi = 0; vi < modinfo_->NavgCrd(); vi += 3) {
          results_[vi  ] += Vec[vi  ] * fi;
          results_[vi+1] += Vec[vi+1] * fi;
          results_[vi+2] += Vec[vi+2] * fi;
        }
      }
    }
    // Output
    if (outfile.OpenWrite( filename_ )) return Analysis::ERR;
    outfile.Printf("#Analysis of modes: DISPLACEMENT\n");
    outfile.Printf("%-10s %10s %10s %10s\n", "#Atom_no.", "displX", "displY", "displZ");
    int anum = 1;
    for (int i3 = 0; i3 < modinfo_->NavgCrd(); i3 += 3)
      outfile.Printf("%10i %10.3f %10.3f %10.3f\n", anum++, results_[i3], 
                     results_[i3+1], results_[i3+2]);
  // ----- CORR PRINT -----
  } else if (type_ == CORR) {
    // Calc dipole-dipole correlation functions
    CalcDipoleCorr();
    if (results_==0) return Analysis::ERR;
    if (outfile.OpenWrite( filename_ )) return Analysis::ERR;
    outfile.Printf("#Analysis of modes: CORRELATION FUNCTIONS\n");
    outfile.Printf("%-10s %10s %10s %10s %10s %10s\n", "#Atom1", "Atom2", "Mode", 
                   "Freq", "1-S^2", "P2(cum)");
    int ncnt = 0;
    for (modestack_it apair = atompairStack_.begin();
                      apair != atompairStack_.end(); ++apair)
    {
      outfile.Printf("%10i %10i\n", (*apair).first+1, (*apair).second+1);
      double val = 1.0;
      for (int mode = beg_; mode < end_; ++mode) {
        double frq = modinfo_->Eigenvalue(mode);
        if (frq >= 0.5) {
          val += results_[ncnt];
          outfile.Printf("%10s %10s %10i %10.5f %10.5f %10.5f\n",
                         "", "", mode, frq, results_[ncnt], val);
          ++ncnt;
        }
      }
    } // END loop over CORR atom pairs
  // ----- PSEUDO-TRAJECTORY -----
  } else if (type_ == TRAJ) {
    ProjectCoords();
  // ----- EIGENVALUE FRACTION -----
  } else if (type_ == EIGENVAL) {
    double sum = 0.0;
    for (unsigned int mode = 0; mode != modinfo_->Size(); mode++)
      sum += modinfo_->Eigenvalue( mode );
    mprintf("\t%zu eigenvalues, sum is %f\n", modinfo_->Size(), sum);
    double cumulative = 0.0;
    if (outfile.OpenWrite( filename_ )) return Analysis::ERR;
    outfile.Printf("%6s %12s %12s %12s\n", "#Mode", "Frac.", "Cumulative", "Eigenval");
    for (unsigned int mode = 0; mode != modinfo_->Size(); mode++) {
      double frac = modinfo_->Eigenvalue( mode ) / sum;
      cumulative += frac;
      outfile.Printf("%6u %12.6f %12.6f %12.6f\n", mode+1, frac, cumulative, 
                     modinfo_->Eigenvalue( mode ));
    }
  } else if (type_ == RMSIP) {
    if (outfile.OpenWrite( filename_ )) return Analysis::ERR;
    if (CalcRMSIP(outfile)) return Analysis::ERR; 
  } else // SANITY CHECK
    return Analysis::ERR;
  outfile.CloseFile();

  return Analysis::OK;
}

// Analysis_Modes::CalcDipoleCorr()
void Analysis_Modes::CalcDipoleCorr() {
  double qcorr = 1.0; // For when bose is false
  int rsize = atompairStack_.size() * (end_ - beg_ + 1);
  results_ = new double[ rsize ];
  std::fill(results_, results_ + rsize, 0);
  // Loop over atom pairs 
  double* Res = results_;
  DataSet_Modes::Darray const& Avg = modinfo_->AvgCrd();
  for (modestack_it apair = atompairStack_.begin(); apair != atompairStack_.end(); ++apair)
  {
    int idx1 = (*apair).first  * 3;
    int idx2 = (*apair).second * 3;
    // Check if coordinates are out of bounds
    if (idx1 >= modinfo_->NavgCrd() || idx2 >= modinfo_->NavgCrd()) {
      mprintf("Warning: Atom pair %i -- %i is out of bounds (# avg atoms in modes = %i)\n",
              (*apair).first + 1, (*apair).second + 1, modinfo_->NavgCrd() / 3);
      continue;
    }
    // Calc unit vector along at2->at1 bond
    Vec3 vec = Vec3(Avg[idx1], Avg[idx1+1], Avg[idx1+2]) -
               Vec3(Avg[idx2], Avg[idx2+1], Avg[idx2+2]);
    double dnorm = sqrt( vec.Magnitude2() );
    vec /= dnorm;
    // Precalc certain values
    dnorm = 3.0 / (dnorm * dnorm);
    double vx2 = vec[0] * vec[0];
    double vxy = vec[0] * vec[1];
    double vxz = vec[0] * vec[2];
    double vy2 = vec[1] * vec[1];
    double vyz = vec[1] * vec[2];
    double vz2 = vec[2] * vec[2]; 
    // Loop over desired modes
    for (int mode = beg_; mode < end_; ++mode) {
      double eval = modinfo_->Eigenvalue(mode);
      if (eval >= 0.5) {
        // Don't use eigenvectors associated with zero or negative eigenvalues
        // NOTE: Where is 11791.79 from? Should it be a const?
        double frq = eval * eval / 11791.79;
        if (bose_) {
          double argq = CONSQ * eval;
          qcorr = argq / tanh(argq);
        }
        /* Calc the correlation matrix for delta
         *    as in eq. 7.16 of lamm and szabo, J Chem Phys 1986, 85, 7334.
         *  Note that the rhs of this eq. should be multiplied by kT
         */
        double qcorrf = (qcorr / frq) * 0.6;
        const double* Evec = modinfo_->Eigenvector(mode);
        double dx = (Evec[idx1  ] - Evec[idx2  ]);
        double dy = (Evec[idx1+1] - Evec[idx2+1]);
        double dz = (Evec[idx1+2] - Evec[idx2+2]);
        double delx2 = qcorrf * dx * dx;
        double delxy = qcorrf * dx * dy;
        double delxz = qcorrf * dx * dz;
        double dely2 = qcorrf * dy * dy;
        double delyz = qcorrf * dy * dz;
        double delz2 = qcorrf * dz * dz;
        // Correlation in length, eq. 10.2 of lamm and szabo
        // NOTE: Commented out in PTRAJ
        /*****
        rtr0 = 0.0;
        for(j = 0; j < 3; j++)
          for(k = 0; k < 3; k++)
            rtr0 += e[j] * e[k] * del[j][k];
        *****/
        // Librational correlation function, using eq. 7.12 of lamm and szabo
        // (w/o beta on the lhs).
        *(Res++) = (-delx2 + (vx2 * delx2) + (vxy * delxy) + (vxz * delxz)
                    -dely2 + (vxy * delxy) + (vy2 * dely2) + (vyz * delyz)
                    -delz2 + (vxz * delxz) + (vyz * delyz) + (vz2 * delz2))
                   * dnorm; // Here dnorm is actually (3 / dnorm^2)
      } // END if positive definite eigenvalue
    } // END loop over modes
  } // END loop over atom pairs
}

// Calculate projection of coords along given mode.
void Analysis_Modes::CalculateProjection(int set, Frame const& Crd, int mode) {
  double proj = 0.0;
  DataSet_Modes::Darray const& Avg = modinfo_->AvgCrd();
  const double* Vec = modinfo_->Eigenvector(mode);
  for (int idx = 0; idx < Crd.size(); ++idx)
    proj += (Crd[idx] - Avg[idx]) * 1.0 * Vec[idx];
  mprintf("\tFrame %i mode %i projection = %f\n", set, mode, proj);
}

/** Project average coords along eigenvectors */
int Analysis_Modes::ProjectCoords() {
  double scale = 1.0; // TODO: read in or calculate - use factor?
  int max_it = (int)(pcmax_ - pcmin_);
  // Check that size of eigenvectors match # coords
  int ncoord = tOutParm_->Natom() * 3;
  if (ncoord != modinfo_->NavgCrd()) {
    mprinterr("Error: # selected coords (%i) != eigenvector size (%i)\n",
               ncoord, modinfo_->NavgCrd());
    return 1;
  }
  // Check that mode is valid.
  if (tMode_ < 1 || tMode_ > modinfo_->Nmodes() ) {
    mprinterr("Error: mode %i is out of bounds.\n", tMode_);
    return Analysis::ERR;
  }
  // Setup frame to hold output coords, initalized to avg coords.
  Frame outframe;
  outframe.SetupFrameXM( modinfo_->AvgCrd(), modinfo_->Mass() );
  // Point to correct eigenvector
  const double* Vec = modinfo_->Eigenvector(tMode_-1);
  // Initialize coords to pcmin
  for (int idx = 0; idx < ncoord; idx++)
    outframe[idx] += pcmin_ * Vec[idx];
  if (debug_>0) CalculateProjection(0, outframe, tMode_-1);
  // Write first frame with coords at pcmin.
  int set = 0;
  trajout_.WriteFrame(set++, tOutParm_, outframe);
  // Main loop
  for (int it = 0; it < max_it; it++) {
    double* crd = outframe.xAddress();
    // Move coordinates along eigenvector.
    for (int idx = 0; idx < ncoord; ++idx)
      crd[idx] += scale * Vec[idx];
    if (debug_>0) CalculateProjection(set, outframe, tMode_-1);
    // DEBUG: calc PC projection for first 3 modes
    //for (int m = 0; m < 3; m++)
    //  CalculateProjection(set, outframe, m);
    // Write frame
    trajout_.WriteFrame(set++, tOutParm_, outframe);
  }
  trajout_.EndTraj();
  return 0;
}

// Analysis_Modes::CalcRMSIP()
int Analysis_Modes::CalcRMSIP(CpptrajFile& outfile) {
  if (modinfo_->VectorSize() != modinfo2_->VectorSize()) {
    mprinterr("Error: '%s' vector size (%i) != '%s' vector size (%i)\n",
              modinfo_->Legend().c_str(), modinfo_->VectorSize(),
              modinfo2_->Legend().c_str(), modinfo2_->VectorSize());
    return 1;
  }
  if ( beg_ >= modinfo2_->Nmodes() || end_ > modinfo2_->Nmodes() ) {
    mprinterr("Error: beg/end out of range for %s (%i modes)\n",
              modinfo2_->Legend().c_str(), modinfo2_->Nmodes());
    return 1;
  }
  double sumsq = 0.0;
  for (int m1 = beg_; m1 < end_; m1++) {
    const double* ev1 = modinfo_->Eigenvector(m1);
    for (int m2 = beg_; m2 < end_; m2++) {
      const double* ev2 = modinfo2_->Eigenvector(m2);
      double dot = 0.0;
      for (int iv = 0; iv < modinfo_->VectorSize(); iv++)
        dot += ev1[iv] * ev2[iv];
      sumsq += (dot * dot);
    }
  }
  sumsq /= (double)(end_ - beg_);
  double rmsip = sqrt( sumsq );
  mprintf("\tRMSIP= %g\n", rmsip);
  if (!filename_.empty())
    outfile.Printf("%g\n", rmsip);
  return 0;
}
