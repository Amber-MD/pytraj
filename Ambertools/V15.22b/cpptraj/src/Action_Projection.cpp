#include <cmath> // sqrt
#include "Action_Projection.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // integerToString
#include "Constants.h" // DEGRAD

// CONSTRUCTOR
Action_Projection::Action_Projection() :
  modinfo_(0),
  beg_(0),
  end_(0)
{}

void Action_Projection::Help() {
  mprintf("\tevecs <evecs dataset> [out <outfile>] [beg <beg>] [end <end>] [<mask>]\n"
          "\t[dihedrals <dataset arg>]\n\t%s\n"
          "  Calculate projection along given eigenvectors.\n", ActionFrameCounter::HelpText);
}

// Action_Projection::Init()
Action::RetType Action_Projection::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  // Get ibeg, iend, start, stop, offset
  // NOTE: Must get 'end' before InitFrameCounter since the latter checks for 'end'
  beg_ = actionArgs.getKeyInt("beg", 1) - 1;
  end_ = actionArgs.getKeyInt("end", 2);
  if (InitFrameCounter(actionArgs)) return Action::ERR;

  std::string modesname = actionArgs.GetStringKey("modes"); // For backwards compat.
  if (modesname.empty()) modesname = actionArgs.GetStringKey("evecs");
  if (modesname.empty()) {
    mprinterr("Error: No eigenvectors data set specified ('evecs <name>'). To load\n"
              "Error:   eigenvectors from a file use 'readdata <file>' prior to this command.\n");
    return Action::ERR;
  }
  // Check if DataSet exists
  modinfo_ = (DataSet_Modes*)DSL->FindSetOfType( modesname, DataSet::MODES );
  if (modinfo_ == 0) {
    // To preserve backwards compat., if no modes data set specified try to
    // load the data file.
    DataFile dataIn;
    dataIn.SetDebug( debugIn );
    if (dataIn.ReadDataOfType( modesname, DataFile::EVECS, *DSL ))
      return Action::ERR;
    modinfo_ = (DataSet_Modes*)DSL->FindSetOfType( modesname, DataSet::MODES );
    if (modinfo_ == 0) return Action::ERR;
  }
  // Check if beg and end are in bounds.
  if (end_ > modinfo_->Nmodes()) {
    mprintf("Warning: 'end' %i is greater than # evecs (%i); setting end to %i\n",
            end_, modinfo_->Nmodes(), modinfo_->Nmodes());
    end_ = modinfo_->Nmodes();
  }
  if (beg_ < 0 || beg_ >= end_) {
    mprinterr("Error: 'beg' %i out of bounds.\n", beg_+1);
    return Action::ERR;
  }

  // Check modes type
  if (modinfo_->Type() != DataSet_2D::COVAR &&
      modinfo_->Type() != DataSet_2D::MWCOVAR &&
      modinfo_->Type() != DataSet_2D::DIHCOVAR &&
      modinfo_->Type() != DataSet_2D::IDEA)
  {
    mprinterr("Error: evecs type is not COVAR, MWCOVAR, DIHCOVAR, or IDEA.\n");
    return Action::ERR;
  }

  // Output Filename
  DataFile* DF = DFL->AddDataFile( actionArgs.GetStringKey("out"), actionArgs );

  // Get dihedral data sets 
  if (modinfo_->Type() == DataSet_2D::DIHCOVAR) {
    DihedralSets_.clear();
    DihedralSets_.AddTorsionSets( DSL->GetMultipleSets(actionArgs.GetStringKey("dihedrals")) );
    if ( DihedralSets_.empty() ) {
      mprinterr("Error: No valid data sets found.\n");
      return Action::ERR;
    } else if ((int)DihedralSets_.size() * 2 != modinfo_->VectorSize()) {
      mprinterr("Error: Number of dihedral data sets %u does not correspond to"
                " number of eigenvectors %i\n", DihedralSets_.size()*2, modinfo_->VectorSize());
      return Action::ERR;
    } else if ((int)DihedralSets_.size() * 2 != modinfo_->NavgCrd()) {
      mprinterr("Error: Number of dihedral data sets %u does not correspond to"
                " number of average elements %i\n", DihedralSets_.size()*2, modinfo_->NavgCrd());
      return Action::ERR;
    }
  } else
    // Get mask
    mask_.SetMaskString( actionArgs.GetMaskNext() );

  // Set up data sets
  std::string setname = actionArgs.GetStringNext();
  if (setname.empty())
    setname = DSL->GenerateDefaultName("Proj");
  for (int mode = beg_; mode < end_; ++mode) {
    int imode = mode + 1;
    if (modinfo_->Type() != DataSet_2D::IDEA) { // COVAR, MWCOVAR
      DataSet* dout = DSL->AddSetIdx( DataSet::FLOAT, setname, imode );
      if (dout == 0) {
        mprinterr("Error: Could not create output dataset for mode %i\n", imode);
        return Action::ERR;
      }
      dout->SetLegend("Mode"+integerToString(imode));
      project_.push_back( dout );
      if (DF != 0) DF->AddSet( dout );
    } else { // IDEA TODO: Error check
      project_.push_back( DSL->AddSetIdxAspect( DataSet::FLOAT, setname, imode, "X") );
      if (DF != 0) DF->AddSet( project_.back() );
      project_.push_back( DSL->AddSetIdxAspect( DataSet::FLOAT, setname, imode, "Y") );
      if (DF != 0) DF->AddSet( project_.back() );
      project_.push_back( DSL->AddSetIdxAspect( DataSet::FLOAT, setname, imode, "Z") );
      if (DF != 0) DF->AddSet( project_.back() );
      project_.push_back( DSL->AddSetIdxAspect( DataSet::FLOAT, setname, imode, "R") );
      if (DF != 0) DF->AddSet( project_.back() );
    }
  }
  // Set datafile args
  mprintf("    PROJECTION: Calculating projection using eigenvectors %i to %i of %s\n",
          beg_+1, end_, modinfo_->Legend().c_str());
  if (DF != 0)
    mprintf("\tResults are written to %s\n", DF->DataFilename().full());
  FrameCounterInfo();
  if (modinfo_->Type() == DataSet_2D::DIHCOVAR)
    mprintf("\t%zu dihedral data sets.\n", DihedralSets_.size());
  else
    mprintf("\tAtom Mask: [%s]\n", mask_.MaskString());

  return Action::OK;
}

// Action_Projection::Setup()
Action::RetType Action_Projection::Setup(Topology* currentParm, Topology** parmAddress) {
  if (modinfo_->Type() != DataSet_2D::DIHCOVAR) {
    // Setup mask
    if (currentParm->SetupIntegerMask( mask_ )) return Action::ERR;
    if (mask_.None()) {
      mprinterr("Error: No atoms selected.\n");
      return Action::ERR;
    }
    mask_.MaskInfo();
    // Check # of selected atoms against modes info
    if ( modinfo_->Type() == DataSet_2D::COVAR || 
         modinfo_->Type() == DataSet_2D::MWCOVAR)
    {
      // Check if (3 * number of atoms in mask) and nvectelem agree
      int natom3 = mask_.Nselected() * 3;
      if ( natom3 != modinfo_->NavgCrd() ) {
        mprinterr("Error: number selected coords (%i) != number avg coords (%i) in %s\n",
                  natom3, modinfo_->NavgCrd(), modinfo_->Legend().c_str());
        return Action::ERR;
      }
      if ( natom3 != modinfo_->VectorSize() ) {
        mprinterr("Error: number selected coords (%i) != eigenvector size (%i)\n",
                  natom3, modinfo_->VectorSize() );
        return Action::ERR;
      }
    } else if ( modinfo_->Type() == DataSet_2D::IDEA ) {
      // Check if (number of atoms in mask) and nvectelem agree
      if (//mask_.Nselected() != modinfo_.Navgelem() ||
          mask_.Nselected() != modinfo_->VectorSize()) 
      {
        mprinterr("Error: number selected atoms (%i) != eigenvector size (%i)\n",
                  mask_.Nselected(), modinfo_->VectorSize() );
        return Action::ERR;
      }
    }

    // Precalc sqrt of mass for each coordinate
    sqrtmasses_.clear();
    if ( modinfo_->Type() == DataSet_2D::MWCOVAR ) {
      sqrtmasses_.reserve( mask_.Nselected() );
      for (AtomMask::const_iterator atom = mask_.begin(); atom != mask_.end(); ++atom)
        sqrtmasses_.push_back( sqrt( (*currentParm)[*atom].Mass() ) );
    } else {
      // If not MWCOVAR no mass-weighting necessary
      sqrtmasses_.resize( mask_.Nselected(), 1.0 );
    }
  }
  return Action::OK;
}

// Action_Projection::DoAction()
Action::RetType Action_Projection::DoAction(int frameNum, Frame* currentFrame, 
                                            Frame** frameAddress)
{
  if ( CheckFrameCounter( frameNum ) ) return Action::OK;
  // Always start at first eigenvector element of first mode.
  const double* Vec = modinfo_->Eigenvector(beg_);
  // Project snapshots on modes
  if ( modinfo_->Type() == DataSet_2D::COVAR || 
       modinfo_->Type() == DataSet_2D::MWCOVAR ) 
  {
    for (int mode = beg_; mode < end_; ++mode) {
      DataSet_Modes::AvgIt Avg = modinfo_->AvgBegin();
      double proj = 0;
      std::vector<double>::const_iterator sqrtmass = sqrtmasses_.begin();
      for (AtomMask::const_iterator atom = mask_.begin(); atom != mask_.end(); ++atom)
      {
        const double* XYZ = currentFrame->XYZ( *atom );
        double mass = *(sqrtmass++);
        proj += (XYZ[0] - *(Avg++)) * mass * Vec[0]; 
        proj += (XYZ[1] - *(Avg++)) * mass * Vec[1]; 
        proj += (XYZ[2] - *(Avg++)) * mass * Vec[2]; 
        Vec += 3;
      }
      float fproj = (float)proj;
      project_[mode]->Add( frameNum, &fproj );
    }
  } else if (modinfo_->Type() == DataSet_2D::DIHCOVAR ) {
    for (int mode = beg_; mode < end_; ++mode) {
      DataSet_Modes::AvgIt Avg = modinfo_->AvgBegin();
      double proj = 0.0;
      for (Array1D::const_iterator dih = DihedralSets_.begin();
                                   dih != DihedralSets_.end(); ++dih)
      {
        double theta = (*dih)->Dval( frameNum ) * Constants::DEGRAD;
        proj += (cos(theta) - *(Avg++)) * Vec[0];
        proj += (sin(theta) - *(Avg++)) * Vec[1];
        Vec += 2;
      }
      // TODO: Convert to degrees?
      float fproj = (float)proj;
      project_[mode]->Add( frameNum, &fproj );
    }
  } else { // if modinfo_.Type() == IDEA
    int ip = 0;
    for (int mode = beg_; mode < end_; ++mode) {
      double proj1 = 0;
      double proj2 = 0;
      double proj3 = 0;
      for (AtomMask::const_iterator atom = mask_.begin(); atom != mask_.end(); ++atom)
      {
        const double* XYZ = currentFrame->XYZ(*atom);
        proj1 += XYZ[0] * *(Vec  );
        proj2 += XYZ[1] * *(Vec  );
        proj3 += XYZ[2] * *(Vec++);
      }
      float fproj1 = (float)proj1;
      float fproj2 = (float)proj2;
      float fproj3 = (float)proj3;
      double proj4 = sqrt(proj1*proj1 + proj2*proj2 + proj3*proj3);
      float fproj4 = (float)proj4;
      project_[ip++]->Add( frameNum, &fproj1 );
      project_[ip++]->Add( frameNum, &fproj2 );
      project_[ip++]->Add( frameNum, &fproj3 );
      project_[ip++]->Add( frameNum, &fproj4 );
    }
  }
  return Action::OK;
}
