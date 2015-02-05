// AxisType
#include <map>
#include "AxisType.h"
#include "CpptrajStdio.h"

// ---------- NA REFERENCE BASE ATOM NAMES AND COORDS --------------------------
struct NA_RefAtom {
  double x;
  double y;
  double z;
  int hb_index;
  int rms_fit;
  const char* aname;
};

typedef const NA_RefAtom* RefPtr;

// x y z hb_idx rms_fit name
static const NA_RefAtom R_ADE[] = {
  { -2.479000,  5.346000,  0.000000, -1, 0, "C1' "},
  { -1.291000,  4.498000,  0.000000, -1, 1, "N9  "},
  {  0.024000,  4.897000,  0.000000, -1, 1, "C8  "},
  {  0.877000,  3.902000,  0.000000, -1, 1, "N7  "},
  {  0.071000,  2.771000,  0.000000, -1, 1, "C5  "},
  {  0.369000,  1.398000,  0.000000, -1, 1, "C6  "},
  {  1.611000,  0.909000,  0.000000,  0, 0, "N6  "},
  { -0.668000,  0.532000,  0.000000,  1, 1, "N1  "},
  { -1.912000,  1.023000,  0.000000, -1, 1, "C2  "},
  { -2.320000,  2.290000,  0.000000, -1, 1, "N3  "},
  { -1.267000,  3.124000,  0.000000, -1, 1, "C4  "},
  {  0.000000,  0.000000,  0.000000,  0, 0, 0     }
};

static const NA_RefAtom R_CYT[] = {
  { -2.477000,  5.402000,  0.000000, -1, 0, "C1' "},
  { -1.285000,  4.542000,  0.000000, -1, 1, "N1  "},
  { -1.472000,  3.158000,  0.000000, -1, 1, "C2  "},
  { -2.628000,  2.709000,  0.000000,  2, 0, "O2  "},
  { -0.391000,  2.344000,  0.000000,  1, 1, "N3  "},
  {  0.837000,  2.868000,  0.000000, -1, 1, "C4  "},
  {  1.875000,  2.027000,  0.000000,  0, 0, "N4  "},
  {  1.056000,  4.275000,  0.000000, -1, 1, "C5  "},
  { -0.023000,  5.068000,  0.000000, -1, 1, "C6  "},
  {  0.000000,  0.000000,  0.000000,  0, 0, 0     }
};

static const NA_RefAtom R_GUA[] = {
  { -2.477000,  5.399000,  0.000000, -1, 0, "C1' "},
  { -1.289000,  4.551000,  0.000000, -1, 1, "N9  "},
  {  0.023000,  4.962000,  0.000000, -1, 1, "C8  "},
  {  0.870000,  3.969000,  0.000000, -1, 1, "N7  "},
  {  0.071000,  2.833000,  0.000000, -1, 1, "C5  "},
  {  0.424000,  1.460000,  0.000000, -1, 1, "C6  "},
  {  1.554000,  0.955000,  0.000000,  0, 0, "O6  "},
  { -0.700000,  0.641000,  0.000000,  1, 1, "N1  "},
  { -1.999000,  1.087000,  0.000000, -1, 1, "C2  "},
  { -2.949000,  0.139000, -0.001000,  2, 0, "N2  "},
  { -2.342000,  2.364000,  0.001000, -1, 1, "N3  "},
  { -1.265000,  3.177000,  0.000000, -1, 1, "C4  "},
  {  0.000000,  0.000000,  0.000000,  0, 0, 0     }
};

static const NA_RefAtom R_THY[] = {
  { -2.481000,  5.354000,  0.000000, -1, 0, "C1' "},
  { -1.284000,  4.500000,  0.000000, -1, 1, "N1  "},
  { -1.462000,  3.135000,  0.000000, -1, 1, "C2  "},
  { -2.562000,  2.608000,  0.000000, -1, 0, "O2  "},
  { -0.298000,  2.407000,  0.000000,  1, 1, "N3  "},
  {  0.994000,  2.897000,  0.000000, -1, 1, "C4  "},
  {  1.944000,  2.119000,  0.000000,  0, 0, "O4  "},
  {  1.106000,  4.338000,  0.000000, -1, 1, "C5  "},
  {  2.466000,  4.961000,  0.001000, -1, 0, "C7  "},
  { -0.024000,  5.057000,  0.000000, -1, 1, "C6  "},
  {  0.000000,  0.000000,  0.000000,  0, 0, 0     }
};

static const NA_RefAtom R_URA[] = {
  { -2.481000,  5.354000,  0.000000, -1, 0, "C1' "},
  { -1.284000,  4.500000,  0.000000, -1, 1, "N1  "},
  { -1.462000,  3.131000,  0.000000, -1, 1, "C2  "},
  { -2.563000,  2.608000,  0.000000, -1, 0, "O2  "},
  { -0.302000,  2.397000,  0.000000,  1, 1, "N3  "},
  {  0.989000,  2.884000,  0.000000, -1, 1, "C4  "},
  {  1.935000,  2.094000, -0.001000,  0, 0, "O4  "},
  {  1.089000,  4.311000,  0.000000, -1, 1, "C5  "},
  { -0.024000,  5.053000,  0.000000, -1, 1, "C6  "},
  {  0.000000,  0.000000,  0.000000,  0, 0, 0     }
};

// UNKNOWN_BASE, ADE, CYT, GUA, THY, URA
/// 1 character base names corresponding to NAbaseType
static const char NAbaseChar[] = { '?', 'A', 'C', 'G', 'T', 'U' }; 
#ifdef NASTRUCTDEBUG
/// Base names corresponding to NAbaseType
static const char* NAbaseName[] = { "UNK", "ADE", "CYT", "GUA", "THY", "URA" };
#endif

// ---------- NA_Base ----------------------------------------------------------
NA_Base::NA_Base() :
  rnum_(0),
  bchar_('?'),
  type_(UNKNOWN_BASE)
{
  hbidx_[0] = -1;
  hbidx_[1] = -1;
  hbidx_[2] = -1;
  std::fill( atomIdx_, atomIdx_+6, -1 );
}

// COPY CONSTRUCTOR
NA_Base::NA_Base(const NA_Base& rhs) :
  rnum_(rhs.rnum_),
  bchar_(rhs.bchar_),
  type_(rhs.type_),
  Ref_(rhs.Ref_),
  anames_(rhs.anames_),
# ifdef NASTRUCTDEBUG
  rname_(rhs.rname_),
  refnames_(rhs.refnames_),
# endif
  Inp_(rhs.Inp_),
  parmMask_(rhs.parmMask_),
  inpFitMask_(rhs.inpFitMask_),
  refFitMask_(rhs.refFitMask_)
{
  hbidx_[0] = rhs.hbidx_[0];
  hbidx_[1] = rhs.hbidx_[1];
  hbidx_[2] = rhs.hbidx_[2];
  std::copy( rhs.atomIdx_, rhs.atomIdx_+6, atomIdx_ );
}

// ASSIGNMENT
NA_Base& NA_Base::operator=(const NA_Base& rhs) {
  if (this == &rhs) return *this;
  rnum_ = rhs.rnum_;
  bchar_ = rhs.bchar_;
  type_ = rhs.type_;
  Ref_ = rhs.Ref_;
  anames_ = rhs.anames_;
# ifdef NASTRUCTDEBUG
  rname_ = rhs.rname_;
  refnames_ = rhs.refnames_;
# endif
  Inp_ = rhs.Inp_;
  hbidx_[0] = rhs.hbidx_[0];
  hbidx_[1] = rhs.hbidx_[1];
  hbidx_[2] = rhs.hbidx_[2];
  std::copy( rhs.atomIdx_, rhs.atomIdx_+6, atomIdx_ );
  parmMask_ = rhs.parmMask_;
  inpFitMask_ = rhs.inpFitMask_;
  refFitMask_ = rhs.refFitMask_;
  return *this;
}

// NA_Base::ID_BaseFromName()
/** Identify NA base typ from residue name. */
NA_Base::NAType NA_Base::ID_BaseFromName(NameType const& resname) {
  if (resname[0]=='D') {
    // If residue name begins with D, assume AMBER DNA residue
    switch (resname[1]) {
      case 'A': return ADE;
      case 'C': return CYT;
      case 'G': return GUA;
      case 'T': return THY;
    }
  } else if (resname[0]=='R') {
    // If residue name beings with R, assume AMBER RNA residue
    switch (resname[1]) {
      case 'A': return ADE;
      case 'C': return CYT;
      case 'G': return GUA;
      case 'U': return URA;
    }
  } else if (resname[2] == ' ' && (resname[1] == '3' || resname[1] == '5')) {
    // Look for 1 letter terminal NA residue names
    if (resname[0] == 'A') return ADE;
    if (resname[0] == 'C') return CYT;
    if (resname[0] == 'G') return GUA;
    if (resname[0] == 'T') return THY;
    if (resname[0] == 'U') return URA;
  } else {
    // Look for standard 3 letter/1 letter NA residue names
    if ( resname == "ADE " ) return ADE;
    if ( resname == "CYT " ) return CYT;
    if ( resname == "GUA " ) return GUA;
    if ( resname == "THY " ) return THY;
    if ( resname == "URA " ) return URA;
    if ( resname == "A   " ) return ADE;
    if ( resname == "C   " ) return CYT;
    if ( resname == "G   " ) return GUA;
    if ( resname == "T   " ) return THY;
    if ( resname == "U   " ) return URA;
  } 
  return UNKNOWN_BASE;
}

int NA_Base::FindAtom(NameType const& atname) const {
  int atom = 0;
  for (std::vector<NameType>::const_iterator Name = anames_.begin();
                                             Name != anames_.end(); ++Name)
  {
    if (*Name == atname) return atom;
    ++atom;
  }
  return -1;
}

/** Set NA residue reference coordinates for given NA base. Ensure that
  * the atom ordering in the reference matches that in the given parm.
  * If an error occurs the type will be set to UNKNOWN_BASE. 
  */
NA_Base::NA_Base(Topology const& currentParm, int resnum, NA_Base::NAType baseType) 
{
  type_ = UNKNOWN_BASE;
  // Set up reference info for this base type 
  RefPtr REF = 0;
  switch (baseType) {
    case ADE : REF = R_ADE; break;
    case CYT : REF = R_CYT; break;
    case GUA : REF = R_GUA; break;
    case THY : REF = R_THY; break;
    case URA : REF = R_URA; break;
    case UNKNOWN_BASE: // Sanity check; should never be called with UNKNOWN 
      REF = 0;
      mprinterr("Internal Error: Residue %i is not a recognized NA residue.\n", resnum+1);
  }
  if ( REF != 0 ) {
    int resstart = currentParm.Res(resnum).FirstAtom();
    int resstop = currentParm.Res(resnum).LastAtom();
    // Create mask for all input coords for this residue
    parmMask_.AddAtomRange(resstart, resstop);
    // Allocate space to hold input coords
    Inp_.SetupFrame( parmMask_.Nselected() );
    // Save atom names for input coords. Look for specific atom names for
    // calculating things like groove width and pucker.
    int inpatom = 0;
    std::fill( atomIdx_, atomIdx_+6, -1 );
    for (int atom = resstart; atom < resstop; ++atom) {
      anames_.push_back( currentParm[atom].Name() );
      // Is this atom P?
      if (anames_.back() == "P   ")
        atomIdx_[PHOS] = inpatom;
      // Is this atom O4'/O4*?
      else if (anames_.back() == "O4' " || anames_.back() == "O4* ")
        atomIdx_[O4p] = inpatom;
      else if (anames_.back() == "C1' " || anames_.back() == "C1* ")
        atomIdx_[C1p] = inpatom;
      else if (anames_.back() == "C2' " || anames_.back() == "C2* ")
        atomIdx_[C2p] = inpatom;
      else if (anames_.back() == "C3' " || anames_.back() == "C3* ")
        atomIdx_[C3p] = inpatom;
      else if (anames_.back() == "C4' " || anames_.back() == "C4* ")
        atomIdx_[C4p] = inpatom;
      inpatom++;
    }
    // For each atom defined as a reference atom for this base, find the
    // corresponding atom in the parm.
    std::map<int,int> BaseMap;
    int refatom = 0;
    for (RefPtr ref = REF; ref->aname != 0; ++ref) {
      NameType atomName(ref->aname);
      inpatom = FindAtom(atomName);
      // Sometimes C1' is listed as C1*; if search for C1' fails look for C1*.
      if (inpatom < 0 && atomName == "C1' ")
        inpatom = FindAtom("C1* ");
      if (inpatom < 0) {
        mprinterr("Error: Ref Atom [%s] not found in NA base [%s].\n",
                  ref->aname, currentParm.Res(resnum).c_str());
        BaseMap.clear();
        break;
      } else {
        BaseMap.insert( std::pair<int,int>(inpatom, refatom) );
#       ifdef NASTRUCTDEBUG
        mprintf("Ref atom %i:%s found in parm (%i:%s)\n",refatom+1,ref->aname,
                inpatom+1,*anames_[inpatom]);
#       endif
      }
      ++refatom;
    }
    if (!BaseMap.empty()) {
      // Now create reference frame with same order as parm. Create RMS fit
      // masks for Ref and Inp frames. If atom is indicated as H-bonding store
      // its index in Inp.
      hbidx_[0] = -1;
      hbidx_[1] = -1;
      hbidx_[2] = -1;
      int refidx = 0; // Will correspond to Ref frame
      for (std::map<int,int>::iterator atom = BaseMap.begin(); 
                                       atom != BaseMap.end(); atom++, refidx++) {
        inpatom = (*atom).first;
        refatom = (*atom).second;
        // Check if this is an H bonding atom. If so, store the memory address
        // of the appropriate place in the coordinate array for future Hbond
        // calculations.
        int hb_index = REF[refatom].hb_index;
        if ( hb_index != -1 ) 
          hbidx_[hb_index] = inpatom;
        // Store coords
        Ref_.AddVec3( Vec3(REF[refatom].x, REF[refatom].y, REF[refatom].z) );
#       ifdef NASTRUCTDEBUG
        // Store reference atom names
        refnames_.push_back( REF[refatom].aname );
#       endif
        // Will this atom be used for RMS fitting?
        if (REF[refatom].rms_fit == 1) {
          inpFitMask_.AddAtom( inpatom );
          refFitMask_.AddAtom( refidx );
        }
      }
      // Make sure all masks have atoms
      if (parmMask_.None() || inpFitMask_.None() || refFitMask_.None()) 
        mprinterr("Error: One or more masks for NA residue %i has no atoms.\n", resnum+1);
      else {
        rnum_ = resnum;
        type_ = baseType;
        bchar_ = NAbaseChar[type_];
#       ifdef NASTRUCTDEBUG
        rname_ = currentParm.Res(resnum).Name();
        mprintf("\tSet up residue %i:%s as %s (%c)\n", rnum_+1, *rname_, NAbaseName[type_], bchar_);
        mprintf("\tReference Atoms:\n");
        for (int atom = 0; atom < Ref_.Natom(); ++atom) {
          mprintf("\t\t%s: ", *(refnames_[atom]));
          Ref_.printAtomCoord(atom);
        }
        mprintf("\tResidue is %i atoms:\n", Inp_.Natom());
        for (int atom = 0; atom < (int)anames_.size(); ++atom)
          mprintf("\t\t%s: %i\n", *(anames_[atom]), atom+1);
        mprintf("\tHBidxs={%i, %i, %i}  P=%i  O4'=%i\n",
                hbidx_[0]+1, hbidx_[1]+1, hbidx_[2]+1, atomIdx_[PHOS]+1, atomIdx_[O4p]+1);
        parmMask_.PrintMaskAtoms("ParmMask");
        inpFitMask_.PrintMaskAtoms("InputFitMask");
        refFitMask_.PrintMaskAtoms("RefFitMask");
#       endif
      }
    } else
      mprinterr("Error: Could not set up reference for residue %i\n", resnum+1);
  }
}

void NA_Base::SetInputFrame(Frame const& inputFrame) {
  Inp_.SetCoordinates( inputFrame, parmMask_ );
}

void NA_Base::PrintAtomNames() const {
  mprintf("\tInp Atoms:");
  for (std::vector<NameType>::const_iterator aname = anames_.begin();
                                             aname != anames_.end(); ++aname)
    mprintf(" %s", *(*aname));
  mprintf("\n");
}

bool NA_Base::HasSugarAtoms() const {
  for (int i = 1; i < 6; i++)
    if (atomIdx_[i] < 0) return false;
  return true;
}
// ---------- NA_Axis ----------------------------------------------------------
// CONSTRUCTOR
NA_Axis::NA_Axis() : residue_number_(0), second_resnum_(-1), isAnti_(false) {}

void NA_Axis::SetupBaseAxis(Matrix_3x3 const& Rin, Vec3 const& oIn, int rnum) {
  StoreRotMatrix(Rin, oIn);
  residue_number_ = rnum;
}

NA_Axis::NA_Axis(int rnum1, int rnum2, bool isAntiIn) :
  residue_number_(rnum1),
  second_resnum_(rnum2),
  isAnti_(isAntiIn)
{}

void NA_Axis::StoreRotMatrix(Matrix_3x3 const& Rin, Vec3 const& vIn) {
  R_ = Rin;
  RX_ = R_.Col1();
  RY_ = R_.Col2();
  RZ_ = R_.Col3();
  origin_ = vIn;
}
 
void NA_Axis::PrintAxisInfo(const char *title) const {
  mprintf("         %s origin: %8.4f %8.4f %8.4f\n",title,origin_[0],origin_[1],origin_[2]);
  mprintf("         %s R_x vec: %8.4f %8.4f %8.4f\n",title,R_[0],R_[3],R_[6]);
  mprintf("         %s R_y vec: %8.4f %8.4f %8.4f\n",title,R_[1],R_[4],R_[7]);
  mprintf("         %s R_z vec: %8.4f %8.4f %8.4f\n",title,R_[2],R_[5],R_[8]);
}

/** Flip the Z and Y axes. Equivalent to rotation around the X axis.
  * Done for antiparallel stranded DNA.
  */
void NA_Axis::FlipYZ() {
  R_[1] = -R_[1]; // -Yx
  R_[4] = -R_[4]; // -Yy
  R_[7] = -R_[7]; // -Yz
  R_[2] = -R_[2]; // -Zx
  R_[5] = -R_[5]; // -Zy
  R_[8] = -R_[8]; // -Zz
  RY_.Neg();
  RZ_.Neg();
}

/** Flip the X and Y axes. Equivalent to rotation around the Z axis.
  * Done for parallel stranded DNA.
  */
void NA_Axis::FlipXY() {
  R_[0] = -R_[0]; // -Xx
  R_[3] = -R_[3]; // -Xy
  R_[6] = -R_[6]; // -Xz
  R_[1] = -R_[1]; // -Yx
  R_[4] = -R_[4]; // -Yy
  R_[7] = -R_[7]; // -Yz
  RX_.Neg();
  RY_.Neg();
}
