#include <cmath> // sqrt
#include <cfloat> // DBL_MAX
#include <cstdlib> // abs, intel 11 compilers choke on std::abs
#include "Action_NativeContacts.h"
#include "CpptrajStdio.h"
#include "DistRoutines.h"

// CONSTRUCTOR
Action_NativeContacts::Action_NativeContacts() :
  distance_(7.0),
  debug_(0),
  ensembleNum_(-1),
  matrix_min_(0),
  resoffset_(1),
  nframes_(0),
  first_(false),
  byResidue_(false),
  includeSolvent_(false),
  numnative_(0),
  nonnative_(0),
  mindist_(0),
  maxdist_(0),
  map_(0),
  CurrentParm_(0)
{}
// TODO: mapout, avg contacts over traj, 1=native, -1=nonnative
void Action_NativeContacts::Help() {
  mprintf("\t[<mask1> [<mask2>]] [writecontacts <outfile>]\n"
          "\t[noimage] [distance <cut>] [out <filename>] [includesolvent]\n"
          "\t[ first | %s ] [resoffset <n>]\n"
          "\t[name <dsname>] [mindist] [maxdist] [byresidue]\n"
          "\t[map [mapout <mapfile>]]\n"
          "  Calculate number of contacts in <mask1>, or between <mask1> and <mask2>\n"
          "  if both are specified. Native contacts are determined based on the given\n"
          "  reference structure (or first frame if not specified) and the specified\n"
          "  distance cut-off (7.0 Ang. default). If [byresidue] is specified contacts\n"
          "  between two residues spaced <resoffset> residues apart are ignored, and\n"
          "  the map (if specified) is written per-residue.\n", FrameList::RefArgs);
}

/** Set up atom/residue indices corresponding to atoms selected in mask.
  * This is done to make creating an atom/residue contact map easier.
  */
Action_NativeContacts::Iarray Action_NativeContacts::SetupContactIndices(
                                AtomMask const& mask, Topology const& parmIn)
{
  Iarray contactIdx;
  for (AtomMask::const_iterator atom = mask.begin(); atom != mask.end(); ++atom)
    if (byResidue_)
      contactIdx.push_back( parmIn[*atom].ResNum() );
    else
      contactIdx.push_back( *atom );
  return contactIdx;
}

// DEBUG
static void DebugContactList(AtomMask const& mask, Topology const& parmIn)
{
  for (AtomMask::const_iterator atom = mask.begin(); atom != mask.end(); ++atom)
    mprintf("\tPotential Contact %u: %s\n", atom - mask.begin(),
            parmIn.AtomMaskName(*atom).c_str());
}

/** Remove any selected solvent atoms from mask. */
static void removeSelectedSolvent( Topology const& parmIn, AtomMask& mask ) {
  AtomMask newMask = mask;
  newMask.ClearSelected();
  for (AtomMask::const_iterator atom = mask.begin(); atom != mask.end(); ++atom) {
    int molnum = parmIn[*atom].MolNum();
    if (!parmIn.Mol(molnum).IsSolvent())
      newMask.AddSelectedAtom( *atom );
  }
  mask = newMask;
}

/** Based on selected atoms in Mask1 (and optionally Mask2), set up
  * potential contact lists. Also set up atom/residue indices corresponding
  * to each potential contact.
  */
int Action_NativeContacts::SetupContactLists(Topology const& parmIn, Frame const& frameIn)
{
  // Setup first contact list
  if ( parmIn.SetupIntegerMask( Mask1_, frameIn ) ) return 1;
  if (!includeSolvent_) removeSelectedSolvent( parmIn, Mask1_ );
  Mask1_.MaskInfo();
  if (Mask1_.None())
  {
    mprinterr("Warning: Nothing selected for '%s'\n", Mask1_.MaskString());
    return 1;
  }
  if (debug_ > 0) DebugContactList( Mask1_, parmIn );
  contactIdx1_ = SetupContactIndices( Mask1_, parmIn );
  // Setup second contact list if necessary
  if ( Mask2_.MaskStringSet() ) {
    if (parmIn.SetupIntegerMask( Mask2_, frameIn ) ) return 1;
    if (!includeSolvent_) removeSelectedSolvent( parmIn, Mask2_ );
    Mask2_.MaskInfo();
    if (Mask2_.None())
    {
      mprinterr("Warning: Nothing selected for '%s'\n", Mask2_.MaskString());
      return 1;
    }
    // Warn if masks overlap
    int nOverlap = Mask1_.NumAtomsInCommon( Mask2_ );
    if (nOverlap > 0) {
      mprintf("Warning: Masks '%s' and '%s' overlap by %i atoms.\n"
              "Warning: Some contacts may be double-counted.\n", 
              Mask1_.MaskString(), Mask2_.MaskString(), nOverlap);
      if (mindist_ != 0)
        mprintf("Warning: Minimum distance will always be 0.0\n");
    }
    if (debug_ > 0) DebugContactList( Mask2_, parmIn );
    contactIdx2_ = SetupContactIndices( Mask2_, parmIn );
  }
  return 0;
}

// Action_NativeContacts::DetermineNativeContacts()
/** Determine potential contacts for given Topology and Frame, then determine 
  * which pairs of contacts satisfy the cutoff and set those as native contacts.
  */
int Action_NativeContacts::DetermineNativeContacts(Topology const& parmIn, Frame const& frameIn)
{
  if ( SetupContactLists(parmIn, frameIn) ) return 1;
  // If specified, set up contacts map; base size on atom masks.
  if (map_ != 0) {
    int matrix_max;
    if (Mask2_.MaskStringSet()) {
      matrix_min_ = std::min( Mask1_[0], Mask2_[0] );
      matrix_max = std::max( Mask1_.back(), Mask2_.back() );
    } else {
      matrix_min_ = Mask1_[0];
      matrix_max = Mask1_.back();
    }
    std::string label("Atom");
    if (byResidue_) {
      matrix_min_ = parmIn[matrix_min_].ResNum();
      matrix_max = parmIn[matrix_max].ResNum();
      label.assign("Residue");
    }
    int matrix_cols = matrix_max - matrix_min_ + 1;
    map_->AllocateHalf( matrix_cols );
    Dimension matrix_dim( matrix_min_+1, 1, matrix_cols, label );
    map_->SetDim(Dimension::X, matrix_dim);
    map_->SetDim(Dimension::Y, matrix_dim);
  }
  double maxDist2 = 0.0;
  double minDist2 = DBL_MAX;
  nativeContacts_.clear();
  if ( Mask2_.MaskStringSet() ) {
    for (AtomMask::const_iterator c1 = Mask1_.begin(); c1 != Mask1_.end(); ++c1)
      for (AtomMask::const_iterator c2 = Mask2_.begin(); c2 != Mask2_.end(); ++c2)
      {
        if (ValidContact(*c1, *c2, parmIn)) {
          double Dist2 = DIST2(frameIn.XYZ(*c1), frameIn.XYZ(*c2), image_.ImageType(),
                               frameIn.BoxCrd(), ucell_, recip_);
          minDist2 = std::min( Dist2, minDist2 );
          maxDist2 = std::max( Dist2, maxDist2 );
          if (Dist2 < distance_)
            nativeContacts_.insert( contactType(*c1, *c2) );
        }
      }
  } else {
    for (AtomMask::const_iterator c1 = Mask1_.begin(); c1 != Mask1_.end(); ++c1)
      for (AtomMask::const_iterator c2 = c1 + 1; c2 != Mask1_.end(); ++c2)
      {
        if (ValidContact(*c1, *c2, parmIn)) {
          double Dist2 = DIST2(frameIn.XYZ(*c1), frameIn.XYZ(*c2), image_.ImageType(),
                               frameIn.BoxCrd(), ucell_, recip_);
          minDist2 = std::min( Dist2, minDist2 );
          maxDist2 = std::max( Dist2, maxDist2 );
          if (Dist2 < distance_)
            nativeContacts_.insert( contactType(*c1, *c2) );
        }
      }
  }
  //mprintf("\tMinimum observed distance= %f, maximum observed distance= %f\n",
  //        sqrt(minDist2), sqrt(maxDist2));
  // Print contacts
  CpptrajFile outfile;
  if (outfile.OpenEnsembleWrite(cfile_, ensembleNum_)) return 1;
  if (!cfile_.empty()) {
    outfile.Printf("#Native contacts determine from mask '%s'", Mask1_.MaskString());
    if (Mask2_.MaskStringSet())
      outfile.Printf(" and mask '%s'", Mask2_.MaskString());
    outfile.Printf("\n");
    mprintf("\tSetup %zu native contacts, written to '%s'\n", 
            nativeContacts_.size(), outfile.Filename().full());
  }
  outfile.Printf("\tSetup %zu native contacts:\n", nativeContacts_.size());
  for (contactListType::const_iterator contact = nativeContacts_.begin();
                                       contact != nativeContacts_.end(); ++contact)
  {
    int a1 = contact->first;
    int a2 = contact->second;
    outfile.Printf("\t\tAtom '%s' to '%s'\n", 
                   parmIn.AtomMaskName(a1).c_str(),
                   parmIn.AtomMaskName(a2).c_str());
  }
  outfile.CloseFile();
  return 0;  
}
// -----------------------------------------------------------------------------
// Action_NativeContacts::Init()
Action::RetType Action_NativeContacts::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  ensembleNum_ = DSL->EnsembleNum();
  debug_ = debugIn;
  // Get Keywords
  image_.InitImaging( !(actionArgs.hasKey("noimage")) );
  double dist = actionArgs.getKeyDouble("distance", 7.0);
  byResidue_ = actionArgs.hasKey("byresidue");
  resoffset_ = actionArgs.getKeyInt("resoffset", 0) + 1;
  if (resoffset_ < 1) {
    mprinterr("Error: Residue offset must be >= 0\n");
    return Action::ERR;
  }
  includeSolvent_ = actionArgs.hasKey("includesolvent");
  distance_ = dist * dist; // Square the cutoff
  first_ = actionArgs.hasKey("first");
  DataFile* outfile = DFL->AddDataFile( actionArgs.GetStringKey("out"), actionArgs );
  cfile_ = actionArgs.GetStringKey("writecontacts");
  // Get reference
  ReferenceFrame REF = FL->GetFrameFromArgs( actionArgs );
  if (!first_) {
    if (REF.error()) return Action::ERR;
    if (REF.empty()) {
      mprintf("Warning: No reference structure specified. Defaulting to first.\n");
      first_ = true;
    }
  }
  // Create data sets
  std::string name = actionArgs.GetStringKey("name");
  if (name.empty())
    name = DSL->GenerateDefaultName("Contacts");
  numnative_ = DSL->AddSetAspect(DataSet::INTEGER, name, "native");
  nonnative_ = DSL->AddSetAspect(DataSet::INTEGER, name, "nonnative");
  if (outfile != 0) {
    outfile->AddSet(numnative_);
    outfile->AddSet(nonnative_);
  }
  if (numnative_ == 0 || nonnative_ == 0) return Action::ERR;
  if (actionArgs.hasKey("mindist")) {
    mindist_ = DSL->AddSetAspect(DataSet::DOUBLE, name, "mindist");
    if (mindist_ == 0) return Action::ERR;
    if (outfile != 0) outfile->AddSet(mindist_);
  }
  if (actionArgs.hasKey("maxdist")) {
    maxdist_ = DSL->AddSetAspect(DataSet::DOUBLE, name, "maxdist");
    if (maxdist_ == 0) return Action::ERR;
    if (outfile != 0) outfile->AddSet(maxdist_);
  }
  if (actionArgs.hasKey("map")) {
    map_ = (DataSet_MatrixDbl*)DSL->AddSetAspect(DataSet::MATRIX_DBL, name, "map");
    if (map_ == 0) return Action::ERR;
    DFL->AddSetToFile( actionArgs.GetStringKey("mapout"), map_ );
  }
  // Get Masks
  Mask1_.SetMaskString( actionArgs.GetMaskNext() );
  std::string mask2 = actionArgs.GetMaskNext();
  if (!mask2.empty())
    Mask2_.SetMaskString( mask2 );
  mprintf("    NATIVECONTACTS: Mask1='%s'", Mask1_.MaskString());
  if (Mask2_.MaskStringSet())
    mprintf(" Mask2='%s'", Mask2_.MaskString());
  mprintf(", contacts set up based on");
  if (first_)
    mprintf(" first frame.\n");
  else
    mprintf("'%s'.\n", REF.FrameName().base());
  if (byResidue_) {
    mprintf("\tContacts will be ignored for residues spaced < %i apart.\n", resoffset_);
    if (map_ != 0)
      mprintf("\tMap will be printed by residue.\n");
  }
  mprintf("\tDistance cutoff is %g Angstroms,", sqrt(distance_));
  if (!image_.UseImage())
    mprintf(" imaging is off.\n");
  else
    mprintf(" imaging is on.\n");
  if (includeSolvent_)
    mprintf("\tMask selection will including solvent.\n");
  else
    mprintf("\tMask selection will not include solvent.\n");
  mprintf("\tData set name: %s\n", name.c_str());
  if (maxdist_ != 0)
    mprintf("\tSaving maximum observed distances in set '%s'\n", maxdist_->Legend().c_str());
  if (mindist_ != 0)
    mprintf("\tSaving minimum observed distances in set '%s'\n", mindist_->Legend().c_str());
  if (outfile != 0)
    mprintf("\tOutput to '%s'\n", outfile->DataFilename().full());
  if (map_ != 0)
    mprintf("\tContacts map will be saved as set '%s'\n", map_->Legend().c_str());
  // Set up reference if necessary.
  if (!first_) {
    // Set up imaging info for ref parm
    image_.SetupImaging( REF.Parm().BoxType() );
    if (image_.ImageType() == NONORTHO)
      REF.Coord().BoxCrd().ToRecip(ucell_, recip_);
    if (DetermineNativeContacts( REF.Parm(), REF.Coord() )) return Action::ERR;
  }

  return Action::OK;
}

// Action_NativeContacts::Setup()
Action::RetType Action_NativeContacts::Setup(Topology* currentParm, Topology** parmAddress) {
  // Setup potential contact lists for this topology
  if (SetupContactLists( *currentParm, Frame()))
    return Action::ERR;
  mprintf("\t%zu potential contact sites for '%s'\n", Mask1_.Nselected(), Mask1_.MaskString());
  if (Mask2_.MaskStringSet())
    mprintf("\t%zu potential contact sites for '%s'\n", Mask2_.Nselected(), Mask2_.MaskString());
  // Set up imaging info for this parm
  image_.SetupImaging( currentParm->BoxType() );
  if (image_.ImagingEnabled())
    mprintf("\tImaging enabled.\n");
  else
    mprintf("\tImaging disabled.\n");
  CurrentParm_ = currentParm;
  return Action::OK;
}

/// \return true if valid contact; when by residue atoms cannot be in same residue
bool Action_NativeContacts::ValidContact(int a1, int a2, Topology const& parmIn) const {
  if (byResidue_) {
    if ( abs(parmIn[a1].ResNum() - parmIn[a2].ResNum()) < resoffset_ )
      return false;
  }
  return true;
}

// Action_NativeContacts::DoAction()
Action::RetType Action_NativeContacts::DoAction(int frameNum, Frame* currentFrame,
                                                Frame** frameAddress)
{
  if (image_.ImageType() == NONORTHO) currentFrame->BoxCrd().ToRecip(ucell_, recip_);
  if (first_) {
    mprintf("\tUsing first frame to determine native contacts.\n");
    if (DetermineNativeContacts( *CurrentParm_, *currentFrame )) return Action::ERR;
    first_ = false;
  }
  nframes_++;
  // Loop over all potential contacts
  int Nnative = 0;
  int NnonNative = 0;
  double maxDist2 = 0.0;
  double minDist2 = DBL_MAX;
  if ( Mask2_.MaskStringSet() ) {
    for (int c1 = 0; c1 != Mask1_.Nselected(); c1++)
      for (int c2 = 0; c2 != Mask2_.Nselected(); c2++)
      {
        if (ValidContact(Mask1_[c1], Mask2_[c2], *CurrentParm_)) {
          double Dist2 = DIST2(currentFrame->XYZ(Mask1_[c1]), currentFrame->XYZ(Mask2_[c2]), 
                               image_.ImageType(), currentFrame->BoxCrd(), ucell_, recip_);
          minDist2 = std::min( Dist2, minDist2 );
          maxDist2 = std::max( Dist2, maxDist2 );
          if (Dist2 < distance_) { // Potential contact
            if (nativeContacts_.find( contactType(Mask1_[c1], Mask2_[c2]) ) != 
                nativeContacts_.end())
            {
              ++Nnative;    // Native contact
              if (map_ != 0) map_->Element(contactIdx1_[c1] - matrix_min_, 
                                           contactIdx2_[c2] - matrix_min_) += 1;
            } else {
              ++NnonNative; // Non-native contact
              if (map_ != 0) map_->Element(contactIdx1_[c1] - matrix_min_,
                                           contactIdx2_[c2] - matrix_min_) -= 1;
            }
          }
        }
      }
  } else {
    for (int c1 = 0; c1 != Mask1_.Nselected(); c1++)
      for (int c2 = c1 + 1; c2 != Mask1_.Nselected(); c2++)
      {
        if (ValidContact(Mask1_[c1], Mask1_[c2], *CurrentParm_)) {
          double Dist2 = DIST2(currentFrame->XYZ(Mask1_[c1]), currentFrame->XYZ(Mask1_[c2]),
                               image_.ImageType(), currentFrame->BoxCrd(), ucell_, recip_);
          minDist2 = std::min( Dist2, minDist2 );
          maxDist2 = std::max( Dist2, maxDist2 );
          if (Dist2 < distance_) { // Potential contact
            if (nativeContacts_.find( contactType(Mask1_[c1], Mask1_[c2]) ) != 
                nativeContacts_.end())
            {
              ++Nnative;    // Native contact
              if (map_ != 0) map_->Element(contactIdx1_[c1] - matrix_min_,
                                           contactIdx1_[c2] - matrix_min_) += 1;
            } else {
              ++NnonNative; // Non-native contact
              if (map_ != 0) map_->Element(contactIdx1_[c1] - matrix_min_, 
                                           contactIdx1_[c2] - matrix_min_) -= 1;
            }
          }
        }
      }
  }
  numnative_->Add(frameNum, &Nnative);
  nonnative_->Add(frameNum, &NnonNative);
  if (mindist_ != 0) {
    minDist2 = sqrt(minDist2);
    mindist_->Add(frameNum, &minDist2);
  }
  if (maxdist_ != 0) {
    maxDist2 = sqrt(maxDist2);
    maxdist_->Add(frameNum, &maxDist2);
  }
  return Action::OK;
}

// Action_NativeContacts::Print()
void Action_NativeContacts::Print() {
  if (map_ != 0) {
    // Normalize map by number of frames.
    double norm = 1.0 / (double)nframes_;
    for (DataSet_MatrixDbl::iterator m = map_->begin(); m != map_->end(); ++m)
      *m *= norm;
  }
}
