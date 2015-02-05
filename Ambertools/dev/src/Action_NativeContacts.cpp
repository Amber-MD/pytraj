#include <cmath> // sqrt
#include <cfloat> // DBL_MAX
#include <cstdlib> // abs, intel 11 compilers choke on std::abs
#include <algorithm> // std::max, std::sort
#include "Action_NativeContacts.h"
#include "CpptrajStdio.h"
#include "DistRoutines.h"
#include "PDBfile.h"

// CONSTRUCTOR
Action_NativeContacts::Action_NativeContacts() :
  distance_(7.0),
  pdbcut_(0.0),
  debug_(0),
  ensembleNum_(-1),
  matrix_min_(0),
  resoffset_(1),
  nframes_(0),
  first_(false),
  byResidue_(false),
  includeSolvent_(false),
  series_(false),
  usepdbcut_(false),
  numnative_(0),
  nonnative_(0),
  mindist_(0),
  maxdist_(0),
  nativeMap_(0),
  nonnatMap_(0),
  CurrentParm_(0),
  refParm_(0),
  masterDSL_(0)
{}
// TODO: mapout, avg contacts over traj, 1=native, -1=nonnative
void Action_NativeContacts::Help() {
  mprintf("\t[<mask1> [<mask2>]] [writecontacts <outfile>] [resout <resfile>]\n"
          "\t[noimage] [distance <cut>] [out <filename>] [includesolvent]\n"
          "\t[ first | %s ]\n"
          "\t[resoffset <n>] [contactpdb <file>] [pdbcut <cut>] [mindist] [maxdist]\n"
          "\t[name <dsname>] [byresidue] [map [mapout <mapfile>]] [series]\n"
          "  Calculate number of contacts in <mask1>, or between <mask1> and <mask2>\n"
          "  if both are specified. Native contacts are determined based on the given\n"
          "  reference structure (or first frame if not specified) and the specified\n"
          "  distance cut-off (7.0 Ang. default). If [byresidue] is specified contacts\n"
          "  between two residues spaced <resoffset> residues apart are ignored, and\n"
          "  the map (if specified) is written per-residue.\n", DataSetList::RefArgs);
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

/** This macro is used by DetermineNativeContacts to set up a new contact
  * if it is valid.
  */
#define SetNativeContact() { \
        if (ValidContact(*c1, *c2, parmIn)) { \
          double Dist2 = DIST2(frameIn.XYZ(*c1), frameIn.XYZ(*c2), image_.ImageType(), \
                               frameIn.BoxCrd(), ucell_, recip_); \
          minDist2 = std::min( Dist2, minDist2 ); \
          maxDist2 = std::max( Dist2, maxDist2 ); \
          if (Dist2 < distance_) { \
            std::string legend(parmIn.AtomMaskName(*c1) + "_" + parmIn.AtomMaskName(*c2)); \
            int r1 = parmIn[*c1].ResNum(); \
            int r2 = parmIn[*c2].ResNum(); \
            ret = nativeContacts_.insert( Mpair(Cpair(*c1,*c2), contactType(legend,r1,r2)) ); \
            if (ret.second && series_) \
              ret.first->second.SetData(masterDSL_->AddSetIdxAspect(DataSet::INTEGER, \
                                                numnative_->Name(), nativeContacts_.size(), \
                                                "NC", legend)); \
          } \
        } \
}

// Action_NativeContacts::DetermineNativeContacts()
/** Determine potential contacts for given Topology and Frame, then determine 
  * which pairs of contacts satisfy the cutoff and set those as native contacts.
  * Should only be called once.
  */
int Action_NativeContacts::DetermineNativeContacts(Topology const& parmIn, Frame const& frameIn)
{
  if (!pfile_.empty()) {
    refFrame_ = frameIn; // Save frame for later PDB output.
    refParm_ = &parmIn;  // Save parm for later PDB output.
  }
  if ( SetupContactLists(parmIn, frameIn) ) return 1;
  // If specified, set up contacts maps; base size on atom masks.
  if (nativeMap_ != 0) {
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
    if (nativeMap_->AllocateHalf( matrix_cols )) return 1;
    if (nonnatMap_->AllocateHalf( matrix_cols )) return 1;
    Dimension matrix_dim( matrix_min_+1, 1, matrix_cols, label );
    nativeMap_->SetDim(Dimension::X, matrix_dim);
    nativeMap_->SetDim(Dimension::Y, matrix_dim);
    nonnatMap_->SetDim(Dimension::X, matrix_dim);
    nonnatMap_->SetDim(Dimension::Y, matrix_dim);
  }
  double maxDist2 = 0.0;
  double minDist2 = DBL_MAX;
  nativeContacts_.clear();
  std::pair<contactListType::iterator, bool> ret; 
  if ( Mask2_.MaskStringSet() ) {
    for (AtomMask::const_iterator c1 = Mask1_.begin(); c1 != Mask1_.end(); ++c1)
      for (AtomMask::const_iterator c2 = Mask2_.begin(); c2 != Mask2_.end(); ++c2)
      {
        SetNativeContact();
      }
  } else {
    for (AtomMask::const_iterator c1 = Mask1_.begin(); c1 != Mask1_.end(); ++c1)
      for (AtomMask::const_iterator c2 = c1 + 1; c2 != Mask1_.end(); ++c2)
      {
        SetNativeContact();
      }
  }
  //mprintf("\tMinimum observed distance= %f, maximum observed distance= %f\n",
  //        sqrt(minDist2), sqrt(maxDist2));
  // Print contacts
  mprintf("\tSetup %zu native contacts:\n", nativeContacts_.size());
  for (contactListType::const_iterator contact = nativeContacts_.begin();
                                       contact != nativeContacts_.end(); ++contact)
  {
    int a1 = contact->first.first;
    int a2 = contact->first.second;
    mprintf("\t\tAtom '%s' to '%s'\n", parmIn.AtomMaskName(a1).c_str(),
            parmIn.AtomMaskName(a2).c_str());
  }
  return 0;  
}
// -----------------------------------------------------------------------------
// Action_NativeContacts::Init()
Action::RetType Action_NativeContacts::Init(ArgList& actionArgs, TopologyList* PFL, DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  masterDSL_ = DSL;
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
  series_ = actionArgs.hasKey("series");
  distance_ = dist * dist; // Square the cutoff
  first_ = actionArgs.hasKey("first");
  DataFile* outfile = DFL->AddDataFile( actionArgs.GetStringKey("out"), actionArgs );
  cfile_ = actionArgs.GetStringKey("writecontacts");
  pfile_ = actionArgs.GetStringKey("contactpdb");
  rfile_ = actionArgs.GetStringKey("resout");
  pdbcut_ = (float)actionArgs.getKeyDouble("pdbcut", -1.0);
  usepdbcut_ = (pdbcut_ > -1.0);
  // Get reference
  ReferenceFrame REF = DSL->GetReferenceFrame( actionArgs );
  if (!first_) {
    if (REF.error()) return Action::ERR;
    if (REF.empty()) {
      mprintf("Warning: No reference structure specified. Defaulting to first.\n");
      first_ = true;
    }
  } else {
    if (!REF.empty()) {
      mprinterr("Error: Must only specify 'first' or a reference structure, not both.\n");
      return Action::ERR;
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
  DataFile *natmapfile = 0, *nonmapfile = 0;
  if (actionArgs.hasKey("map")) {
    nativeMap_ = (DataSet_MatrixDbl*)DSL->AddSetAspect(DataSet::MATRIX_DBL, name, "nativemap");
    if (nativeMap_ == 0) return Action::ERR;
    nonnatMap_ = (DataSet_MatrixDbl*)DSL->AddSetAspect(DataSet::MATRIX_DBL, name, "nonnatmap");
    if (nonnatMap_ == 0) return Action::ERR;
    FileName mapFilename;
    mapFilename.SetFileName( actionArgs.GetStringKey("mapout") );
    if (!mapFilename.empty()) {
      natmapfile = DFL->AddDataFile(mapFilename.DirPrefix() + "native." + mapFilename.Base());
      if (natmapfile != 0) natmapfile->AddSet(nativeMap_);
      nonmapfile = DFL->AddDataFile(mapFilename.DirPrefix() + "nonnative." + mapFilename.Base());
      if (nonmapfile != 0) nonmapfile->AddSet(nonnatMap_);
    }
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
    if (nativeMap_ != 0)
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
  if (!cfile_.empty()) mprintf("\tContact stats will be written to '%s'\n", cfile_.c_str());
  if (!rfile_.empty()) mprintf("\tContact res pairs will be written to '%s'\n", rfile_.c_str());
  if (!pfile_.empty()) {
    mprintf("\tContact PDB will be written to '%s'\n", pfile_.c_str());
    if (usepdbcut_) mprintf("\tOnly atoms with values > %g will be written to PDB.\n", pdbcut_);
  }
  if (nativeMap_ != 0) {
    mprintf("\tNative contacts map will be saved as set '%s'\n"
            "\tNon-native contacts map will be saved as set '%s'\n",
            nativeMap_->Legend().c_str(), nonnatMap_->Legend().c_str());
    if (natmapfile!=0) mprintf("\tNative map output to '%s'\n",natmapfile->DataFilename().full());
    if (nonmapfile!=0) mprintf("\tNative map output to '%s'\n",nonmapfile->DataFilename().full());
  }
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

/** This macro is used by DoAction to check if a contact is valid, formed,
  * if it is native, and if so update it.
  */
#define UpdateNativeContact(M1_, M2_, CI1_, CI2_) { \
        if (ValidContact(M1_[c1], M2_[c2], *CurrentParm_)) { \
          double Dist2 = DIST2(currentFrame->XYZ(M1_[c1]), currentFrame->XYZ(M2_[c2]), \
                               image_.ImageType(), currentFrame->BoxCrd(), ucell_, recip_); \
          minDist2 = std::min( Dist2, minDist2 ); \
          maxDist2 = std::max( Dist2, maxDist2 ); \
          if (Dist2 < distance_) { \
            contactListType::iterator it = nativeContacts_.find( Cpair(M1_[c1], M2_[c2]) ); \
            if (it != nativeContacts_.end()) \
            { \
              ++Nnative; \
              it->second.Increment(frameNum, sqrt(Dist2), Dist2); \
              if (nativeMap_ != 0) nativeMap_->Element(CI1_[c1] - matrix_min_, \
                                                       CI2_[c2] - matrix_min_) += 1; \
            } else { \
              ++NnonNative; \
              if (nonnatMap_ != 0) nonnatMap_->Element(CI1_[c1] - matrix_min_, \
                                                       CI2_[c2] - matrix_min_) += 1; \
            } \
          } \
        } \
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
        UpdateNativeContact(Mask1_, Mask2_, contactIdx1_, contactIdx2_);
      }
  } else {
    for (int c1 = 0; c1 != Mask1_.Nselected(); c1++)
      for (int c2 = c1 + 1; c2 != Mask1_.Nselected(); c2++)
      {
        UpdateNativeContact(Mask1_, Mask1_, contactIdx1_, contactIdx1_);
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
  if (nativeMap_ != 0) {
    // Normalize maps by number of frames.
    double norm = 1.0 / (double)nframes_;
    for (DataSet_MatrixDbl::iterator m = nativeMap_->begin(); m != nativeMap_->end(); ++m)
      *m *= norm;
    for (DataSet_MatrixDbl::iterator m = nonnatMap_->begin(); m != nonnatMap_->end(); ++m)
      *m *= norm;
  }
  if (series_) {
    // Ensure all series have been updated for all frames.
    for (contactListType::iterator it = nativeContacts_.begin();
                                   it != nativeContacts_.end(); ++it)
      if (it->second.Data().Size() < nframes_)
        it->second.Data().AddVal( nframes_ - 1, 0 );
  }
  CpptrajFile outfile;
  if (outfile.OpenEnsembleWrite(cfile_, ensembleNum_)) {
    mprinterr("Error: Could not open file '%s' for writing.\n", cfile_.c_str());
    return;
  }
  if (!cfile_.empty()) {
    mprintf("    CONTACTS: %s: Writing native contacts to file '%s'\n",
            numnative_->Name().c_str(), cfile_.c_str());
    outfile.Printf("# Contacts: %s\n", numnative_->Name().c_str());
    outfile.Printf("# Native contacts determined from mask '%s'", Mask1_.MaskString());
    if (Mask2_.MaskStringSet())
      outfile.Printf(" and mask '%s'", Mask2_.MaskString());
    outfile.Printf("\n");
  } else
    mprintf("    CONTACTS: %s\n", numnative_->Name().c_str());
  // Map of residue pairs to total contact values.
  typedef std::map<Cpair, resContact> resContactMap;
  resContactMap ResContacts;
  std::pair<resContactMap::iterator, bool> ret;
  // Normalize native contacts. Place them into an array where they will
  // be sorted. Sum up total contact over residue pairs.
  std::vector<contactType> sortedList;
  for (contactListType::iterator it = nativeContacts_.begin();
                                 it != nativeContacts_.end(); ++it)
  {
    it->second.Finalize();
    sortedList.push_back( it->second );
    ret = ResContacts.insert( Rpair(Cpair(it->second.Res1(),it->second.Res2()),
                                    resContact(it->second.Nframes())) );
    if (!ret.second) // residue pair exists, update it.
      ret.first->second.Increment( it->second.Nframes() );
  }
  std::sort( sortedList.begin(), sortedList.end() );
  // Place residue pairs into an array to be sorted.
  std::vector<Rpair> ResList;
  for (resContactMap::const_iterator it = ResContacts.begin(); it != ResContacts.end(); ++it)
    ResList.push_back( *it );
  std::sort( ResList.begin(), ResList.end(), res_cmp() );
  // Print out total fraction frames for residue pairs.
  CpptrajFile resout;
  if (resout.OpenWrite(rfile_)==0) {
    resout.Printf("%-8s %8s %10s %10s\n", "#Res1", "#Res2", "TotalFrac", "Contacts");
    //for (resContactMap::const_iterator it = ResContacts.begin(); it != ResContacts.end(); ++it)
    for (std::vector<Rpair>::const_iterator it = ResList.begin();
                                            it != ResList.end(); ++it)
      resout.Printf("%-8i %8i %10g %10i\n", it->first.first+1, it->first.second+1,
                    (double)it->second.Nframes()/(double)nframes_,
                    it->second.Ncontacts());
  }
  resout.CloseFile();
  // Print out sorted atom contacts.
  outfile.Printf("%-8s %20s %8s %8s %8s %8s\n", "#", "Contact", "Nframes", "Frac.", "Avg", "Stdev");
  unsigned int num = 1;
  for (std::vector<contactType>::const_iterator NC = sortedList.begin();
                                                NC != sortedList.end(); ++NC, ++num)
  { 
    double fracPresent = (double)NC->Nframes() / (double)nframes_;
    outfile.Printf("%8u %20s %8i %8.3g %8.3g %8.3g\n", num, NC->id(),
                   NC->Nframes(), fracPresent, NC->Avg(), NC->Stdev());
  }
  outfile.CloseFile();
  // Break down contacts by atom, write to PDB.
  if (!pfile_.empty()) {
    std::vector<double> atomContactFrac(refParm_->Natom(), 0.0);
    double norm = 1.0 / ((double)nframes_ * 2.0);
    for (contactListType::const_iterator it = nativeContacts_.begin();
                                         it != nativeContacts_.end(); ++it)
    {
      int a1 = it->first.first;
      int a2 = it->first.second;
      contactType const& NC = it->second;
      double fracShared = (double)NC.Nframes() * norm;
      atomContactFrac[a1] += fracShared;
      atomContactFrac[a2] += fracShared;
    }
    // Normalize so the strongest contact value is 100.00
    norm = (double)*std::max_element(atomContactFrac.begin(), atomContactFrac.end());
    norm = 100.00 / norm;
    PDBfile contactPDB;
    if (contactPDB.OpenWrite(pfile_))
      mprinterr("Error: Could not open contact PDB for write.\n");
    else {
      mprintf("\tWriting contacts PDB to '%s'\n", pfile_.c_str());
      contactPDB.WriteTITLE( numnative_->Name() + " " + Mask1_.MaskExpression() + " " +
                             Mask2_.MaskExpression() );
      int cidx = 0;
      for (int aidx = 0; aidx != refParm_->Natom(); aidx++, cidx += 3) {
        float bfac = (float)(atomContactFrac[aidx] * norm);
        if (!usepdbcut_ || (bfac > pdbcut_)) {
          int resnum = (*refParm_)[aidx].ResNum();
          contactPDB.WriteCoord(PDBfile::ATOM, aidx+1, (*refParm_)[aidx].Name(),
                                refParm_->Res(resnum).Name(), ' ', resnum+1,
                                refFrame_[cidx], refFrame_[cidx+1], refFrame_[cidx+2],
                                1.0, bfac, (*refParm_)[aidx].ElementName(), 0, false);
        }
      }
      contactPDB.CloseFile();
    }
  }
}
// -----------------------------------------------------------------------------
void Action_NativeContacts::contactType::Finalize() {
  if (nframes_ > 0) {
    dist_ /= (double)nframes_;
    dist2_ /= (double)nframes_;
    dist2_ -= (dist_ * dist_);
    if (dist2_ > 0)
      dist2_ = sqrt(dist2_);
    else
      dist2_ = 0.0;
  }
}
