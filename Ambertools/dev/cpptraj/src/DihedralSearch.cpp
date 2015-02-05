#include "DihedralSearch.h"
#include "CpptrajStdio.h"
#include "AxisType.h" // FIXME: For pyrimidine chi

/// Function to search for atoms within a residue by name.
static int ByName(Topology const& topIn, int res, NameType const& tgtname)
{
  if (res < 0 || res >= topIn.Nres()) return -1;
  for (int i = topIn.Res(res).FirstAtom(); i < topIn.Res(res).LastAtom(); i++)
    if (topIn[i].Name() == tgtname) return i;
  return -1;
}

/// Function to search for atoms within a residue by type.
static int ByType(Topology const& topIn, int res, NameType const& tgttype)
{
  if (res < 0 || res >= topIn.Nres()) return -1;
  for (int i = topIn.Res(res).FirstAtom(); i < topIn.Res(res).LastAtom(); i++)
    if (topIn[i].Type() == tgttype) return i;
  return -1;
}

/// Token to store pre-defined dihedral types.
struct DihedralSearch::DihedralToken::DIH_TYPE {
  int off;
  const char* an0;
  const char* an1;
  const char* an2;
  const char* an3;
  const char* name;
  AtomSearchFxn search0;
  AtomSearchFxn search1;
  AtomSearchFxn search2;
  AtomSearchFxn search3;
};

/** Recognized dihedral types go here. Must correspond to 
  * DihedralSearch::DihedralType, except for NDIHTYPE which does not get 
  * an entry.
  */
const DihedralSearch::DihedralToken::DIH_TYPE DihedralSearch::DihedralToken::DIH[] = {
  {-1, "C"  , "N"  , "CA" , "C"  , "phi"    ,ByName,ByName,ByName,ByName}, // PHI: C0-N1-CA1-C1
  { 1, "N"  , "CA" , "C"  , "N"  , "psi"    ,ByName,ByName,ByName,ByName}, // PSI: N0-CA0-C0-N1
  { 0, "N"  , "CA" , "CB" , "CG",  "chip"   ,ByName,ByName,ByName,ByName}, // Protein CHI:
  {-2, "CA" , "C"  , "N"  , "CA" , "omega"  ,ByName,ByName,ByName,ByName}, // OMEGA: CA0-C0-N1-CA1
  {-1, "O3'", "P"  , "O5'", "C5'", "alpha"  ,ByName,ByName,ByName,ByName}, // ALPHA: 
  { 0, "P"  , "O5'", "C5'", "C4'", "beta"   ,ByName,ByName,ByName,ByName}, // BETA:
  { 0, "O5'", "C5'", "C4'", "C3'", "gamma"  ,ByName,ByName,ByName,ByName}, // GAMMA:
  { 0, "C5'", "C4'", "C3'", "O3'", "delta"  ,ByName,ByName,ByName,ByName}, // DELTA:
  { 1, "C4'", "C3'", "O3'", "P"  , "epsilon",ByName,ByName,ByName,ByName}, // EPSILON:
  { 2, "C3'", "O3'", "P"  , "O5'", "zeta"   ,ByName,ByName,ByName,ByName}, // ZETA:
  { 0, "O4'", "C1'", "C2'", "C3'", "nu1"    ,ByName,ByName,ByName,ByName}, // NU1:
  { 0, "C1'", "C2'", "C3'", "C4'", "nu2"    ,ByName,ByName,ByName,ByName}, // NU2:
  { 0, "O4'", "C1'", "N*",  "C4" , "chin"   ,ByName,ByName,ByType,ByName}  // Nucleic CHI:
  // NOTE: pyrimidine chi last atom name manually changed to C2 below.
};

// DihedralSearch::ListKnownTypes()
void DihedralSearch::ListKnownTypes() {
  for (int i = 0; i < (int)NDIHTYPE; ++i)
    mprintf(" %s", DihedralSearch::DihedralToken::DIH[i].name);
  mprintf("\n");
}

// TODO: const char*
void DihedralSearch::OffsetHelp() {
  mprintf("\t\tOffset -2=<a0><a1> in previous res, -1=<a0> in previous res,\n"
          "\t\t        0=All <aX> in single res,\n"
          "\t\t        1=<a3> in next res, 2=<a2><a3> in next res.\n");
}

// DihedralSearch::GetType()
DihedralSearch::DihedralType DihedralSearch::GetType(std::string const& typeIn) {
  for (int i = 0; i < (int)NDIHTYPE; ++i)
    if (typeIn.compare(DihedralSearch::DihedralToken::DIH[i].name)==0)
      return (DihedralType)i;
  return NDIHTYPE;
}

// -----------------------------------------------------------------------------
// CONSTRUCTOR - DihedralMask
DihedralSearch::DihedralMask::DihedralMask() : 
  a0_(-1), a1_(-1), a2_(-1), a3_(-1), res_(-1), type_(NDIHTYPE) {}

// CONSTRUCTOR - DihedralMask
DihedralSearch::DihedralMask::DihedralMask(int a0, int a1, int a2, int a3, 
                                           int res, std::string const& n,
                                           DihedralType t) :
  a0_(a0), a1_(a1), a2_(a2), a3_(a3), res_(res), name_(n), type_(t) {}

// -----------------------------------------------------------------------------
// CONSTRUCTOR - Custom type 
DihedralSearch::DihedralToken::DihedralToken(int off, 
                                             NameType const& an0, NameType const& an1,
                                             NameType const& an2, NameType const& an3,
                                             std::string const& name) :
  offset_(off),
  name_(name),
  type_(NDIHTYPE)
{
  aname_[0] = an0;
  aname_[1] = an1;
  aname_[2] = an2;
  aname_[3] = an3;
  // Default to ByName search
  search_[0] = ByName;
  search_[1] = ByName;
  search_[2] = ByName;
  search_[3] = ByName;
}

// CONSTRUCTOR - Recognized type 
DihedralSearch::DihedralToken::DihedralToken(DIH_TYPE const& dih, DihedralType dt) :
  offset_(dih.off),
  name_(dih.name),
  type_(dt)
{
  aname_[0] = dih.an0;
  aname_[1] = dih.an1;
  aname_[2] = dih.an2;
  aname_[3] = dih.an3;
  search_[0] = dih.search0;
  search_[1] = dih.search1;
  search_[2] = dih.search2;
  search_[3] = dih.search3;
}

// DihedralSearch::DihedralToken::FindDihedralAtoms()
DihedralSearch::DihedralMask 
  DihedralSearch::DihedralToken::FindDihedralAtoms(Topology const& topIn, int resIn)
{
  int resnum[4] = {resIn, resIn, resIn, resIn};
  switch (offset_) {
    case -2: --resnum[1];        // -1 a2 and a1
    case -1: --resnum[0]; break; // -1 a1 only
    case  2: ++resnum[2];        // +1 a3 and a4
    case  1: ++resnum[3]; break; // +1 a4 only
  }
  int atnum[4];
  for (int i = 0; i < 4; i++) {
    atnum[i] = search_[i](topIn, resnum[i], aname_[i]);
    if (atnum[i] == -1) return DihedralMask();
    // Ensure this atom is bonded to previous atom
    if (i > 0) {
      if ( !topIn[ atnum[i] ].IsBondedTo( atnum[i-1] ) ) {
        mprintf("Warning: Atom %s is not bonded to atom %s\n",
                topIn.AtomMaskName(atnum[i]).c_str(),
                topIn.AtomMaskName(atnum[i-1]).c_str());
        return DihedralMask();
      }
    }
  }
  // All atoms found at this point.
  return DihedralMask(atnum[0], atnum[1], atnum[2], atnum[3], resIn, name_, type_);
}

// -----------------------------------------------------------------------------
// CONSTRUCTOR
DihedralSearch::DihedralSearch() {}

// DihedralSearch::SearchFor()
int DihedralSearch::SearchFor(DihedralType typeIn) {
  dihedralTokens_.push_back( DihedralToken(DihedralSearch::DihedralToken::DIH[typeIn],
                                           typeIn) );
  return 0;
}

// DihedralSearch::SearchForArgs()
/** See if ArgList has any recognized dihedral type keywords. */
void DihedralSearch::SearchForArgs(ArgList& argIn) {
  for (int i = 0; i < (int)NDIHTYPE; ++i) {
    if (argIn.hasKey( DihedralSearch::DihedralToken::DIH[i].name ))
      SearchFor( (DihedralType)i );
  }
}

// DihedralSearch::SearchForNewType()
/** Add new type to search for. */
int DihedralSearch::SearchForNewType(int off, std::string const& an0, std::string const& an1,
                                     std::string const& an2, std::string const& an3,
                                     std::string const& name)
{
  for (std::vector<DihedralToken>::iterator tkn = dihedralTokens_.begin();
                                            tkn != dihedralTokens_.end(); ++tkn)
    if ( tkn->Name() == name ) {
      mprintf("Warning: Dihedral type %s already defined.\n", name.c_str());
      return 1;
    }
  dihedralTokens_.push_back( DihedralToken(off, an0, an1, an2, an3, name) );
  return 0;
}

// DihedralSearch::SearchForAll()
/** If no dihedrals selected yet, select all. */
int DihedralSearch::SearchForAll() {
  if (!dihedralTokens_.empty()) return 0;
  for (int dih=0; dih < (int)NDIHTYPE; dih++)
    SearchFor((DihedralType)dih);
  return 0;
}

// DihedralSearch::FindDihedrals()
int DihedralSearch::FindDihedrals(Topology const& currentParm, Range const& rangeIn)
{
  dihedrals_.clear();
  for (Range::const_iterator res = rangeIn.begin(); res != rangeIn.end(); ++res)
  { // FIXME: Kludge for pyrimidine chi
    NA_Base::NAType baseType = NA_Base::ID_BaseFromName( currentParm.Res(*res).Name() );
    for (std::vector<DihedralToken>::iterator tkn = dihedralTokens_.begin();
                                              tkn != dihedralTokens_.end(); ++tkn)
    {
      // FIXME: Kludge for pyrimidine chi
      if (tkn->Type() == CHIN && (baseType == NA_Base::URA ||
                                  baseType == NA_Base::THY ||
                                  baseType == NA_Base::CYT)  )
      {
        DihedralToken pyrimidineChi = *tkn;
        pyrimidineChi.SetAtomName(3, "C2");
        dihedrals_.push_back( pyrimidineChi.FindDihedralAtoms(currentParm, *res) );
      } else
        dihedrals_.push_back( tkn->FindDihedralAtoms(currentParm, *res) );
      if (dihedrals_.back().None()) {
        mprintf("Warning: Dihedral %s not found for residue %i\n", 
                tkn->Name().c_str(), *res + 1);
        dihedrals_.pop_back();
      } 
    }
  }
  if (dihedrals_.empty()) {
    mprintf("Warning: No dihedrals selected for topology %s\n", currentParm.c_str());
    return 1;
  }
  //mprintf("\tFound %u dihedrals.\n", dihedrals_.size()); // DEBUG
  return 0;
}

// DihedralSearch::Clear()
void DihedralSearch::Clear() {
  dihedralTokens_.clear();
  dihedrals_.clear();
}

// DihedralSearch::ClearFound()
void DihedralSearch::ClearFound() {
  dihedrals_.clear();
}

// DihedralSearch::PrintTypes()
void DihedralSearch::PrintTypes() {
  for (std::vector<DihedralToken>::iterator tkn = dihedralTokens_.begin();
                                            tkn != dihedralTokens_.end(); ++tkn)
  {
    mprintf(" %s", tkn->Name().c_str());
  }
}

// VisitAtom()
static void VisitAtom( Topology const& topIn, int atm, std::vector<bool>& Visited )
{
  // If this atom has already been visited return
  if (Visited[atm]) return;
  // Mark this atom as visited
  Visited[atm] = true;
  // Visit each atom bonded to this atom
  for (Atom::bond_iterator bondedatom = topIn[atm].bondbegin();
                           bondedatom != topIn[atm].bondend(); ++bondedatom)
    VisitAtom(topIn, *bondedatom, Visited);
}

// DihedralSearch::MovingAtoms()
AtomMask DihedralSearch::MovingAtoms(Topology const& topIn, int atom0, int atom1) {
  std::vector<bool> Visited( topIn.Natom(), false );
  // Mark atom0 as already visited
  Visited[atom0] = true;
  for (Atom::bond_iterator bndatm = topIn[atom1].bondbegin();
                           bndatm != topIn[atom1].bondend(); ++bndatm)
  {
    if ( *bndatm != atom0 )
      VisitAtom( topIn, *bndatm, Visited );
  }
  // Everything marked T will move.
  AtomMask Rmask;
  // Needed for conversion to atom mask
  Rmask.SetNatom( topIn.Natom() );
  for (int maskatom = 0; maskatom < (int)Visited.size(); maskatom++) {
    if (Visited[maskatom])
      Rmask.AddAtom(maskatom);
  }
  return Rmask;
}
