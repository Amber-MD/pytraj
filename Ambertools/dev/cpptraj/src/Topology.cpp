#include <stack> // For ParseMask
#include <algorithm> // sort
#ifdef _OPENMP
#  include "omp.h"
#endif
#include "Topology.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // integerToString 
#include "DistRoutines.h"
#include "Constants.h"

const NonbondType Topology::LJ_EMPTY = NonbondType();

// CONSTRUCTOR
Topology::Topology() :
  offset_(0.20),
  debug_(0),
  ipol_(0),
  NsolventMolecules_(0),
  pindex_(0),
  nframes_(0),
  n_extra_pts_(0),
  n_atom_types_(0),
  hasVelInfo_(false),
  nRepDim_(0)
{ }

// Topology::SetParmName()
void Topology::SetParmName(std::string const& title, FileName const& filename) {
  parmName_ = title;
  fileName_ = filename;
}

// Topology::SetReferenceCoords()
void Topology::SetReferenceCoords( Frame const& frameIn ) {
  if (!frameIn.empty()) {
    if (frameIn.Natom() == Natom())
      refCoords_ = frameIn;
    else if (frameIn.Natom() > Natom()) {
      mprintf("Warning: Active reference has %i atoms, parm '%s' has only %i.\n"
              "Warning: Truncating reference coords for this parm (distance-based masks only).\n",
              frameIn.Natom(), c_str(), Natom());
      refCoords_.SetupFrame(Natom());
      std::copy(frameIn.xAddress(), frameIn.xAddress() + refCoords_.size(),
                refCoords_.xAddress());
    } else {
      mprintf("Warning: Active reference has only %i atoms, parm '%s' has %i.\n"
              "Warning: Parm will only have reference coordinates for the first %i atoms"
              " (distance-based masks only).\n",
              frameIn.Natom(), c_str(), Natom(), frameIn.Natom());
      refCoords_.SetupFrame(Natom());
      std::copy(frameIn.xAddress(), frameIn.xAddress() + frameIn.size(), refCoords_.xAddress());
      std::fill(refCoords_.xAddress() + frameIn.size(),
                refCoords_.xAddress() + refCoords_.size(), 0.0);
    }
  }
}

// -----------------------------------------------------------------------------
/** \return Range containing only solute residues. */
Range Topology::SoluteResidues() const {
  Range solute_res;
  atom_iterator atom = atoms_.begin();
  while (atom != atoms_.end()) {
    // If atom is in a solvent molecule skip molecule. Otherwise add res num
    // and skip to next residue.
    if (molecules_[atom->MolNum()].IsSolvent())
      atom += molecules_[atom->MolNum()].NumAtoms();
    else if (molecules_[atom->MolNum()].NumAtoms() == 1) // Assume ion.
      ++atom;
    else {
      solute_res.AddToRange( atom->ResNum() );
      if (debug_ > 0)
        mprintf("DEBUG:\t\tAdding solute residue %i\n", atom->ResNum()+1);
      atom += residues_[atom->ResNum()].NumAtoms();
    }
  }
  return solute_res;
}

// Topology::c_str()
/** Return a printf-compatible char* of the parm filename, or the parm
  * name (title) if the parm filename is empty.
  */
const char *Topology::c_str() const {
  if (!parmTag_.empty())
    return (parmTag_.c_str());
  else if (!fileName_.empty()) 
    return fileName_.base();
  return parmName_.c_str();
}

// -----------------------------------------------------------------------------

// Topology::TruncResAtomName()
/** Given an atom number, return a string containing the corresponding 
  * residue name and number (starting from 1) along with the atom name 
  * with format: 
  * "<resname><resnum>@<atomname>", e.g. "ARG_11@CA".
  * Truncate the residue and atom names so there are no blanks.
  */
std::string Topology::TruncResAtomName(int atom) const {
  std::string res_name;
  if (atom < 0 || atom >= (int)atoms_.size()) return res_name;
  // Atom name with no trailing spaces.
  std::string atom_name = atoms_[atom].Name().Truncated();
  int res = atoms_[atom].ResNum();
  // Residue name with no trailing spaces.
  // NOTE: ensure a residue size of 4?
  res_name = residues_[res].Name().Truncated();
  ++res; // want output as res+1
  res_name += "_";
  res_name += integerToString(res);
  res_name += "@";
  res_name += atom_name;
  return res_name;
}

// Topology::AtomMaskName()
/** \return A string of format :r@a where r is atoms residue number and
  *         a is atoms name.
  */
std::string Topology::AtomMaskName(int atom) const {
  if (atom < 0 || atom >= (int)atoms_.size()) return std::string(""); 
  std::string maskName = ":";
  maskName += integerToString( atoms_[atom].ResNum() + 1 );
  maskName += "@";
  maskName += atoms_[atom].Name().Truncated();
  return maskName;
}

std::string Topology::TruncAtomNameNum(int atom) const {
  if (atom < 0 || atom >= (int)atoms_.size()) return std::string("");
  std::string atom_name = atoms_[atom].Name().Truncated();
  atom_name += "_";
  atom_name += integerToString(atom + 1);
  return atom_name;
}

// Topology::TruncResNameNum()
/** Given a residue number (starting from 0), return a string containing 
  * residue name and number (starting from 1) with format: 
  * "<resname>:<resnum>", e.g. "ARG:11".
  * Truncate residue name so there are no blanks.
  */
// FIXME: Add residue bounds check.
std::string Topology::TruncResNameNum(int res) const {
  // Residue name with no trailing spaces.
  return residues_[res].Name().Truncated() + ":" + integerToString( res+1 );
}

// Topology::FindAtomInResidue()
/** Find the atom # of the specified atom name in the given residue.
  * \param res Residue number to search.
  * \param atname Atom name to find.
  * \return the atom number of the specified atom if found in the given residue.
  * \return -1 if atom not found in given residue.
  */
int Topology::FindAtomInResidue(int res, NameType const& atname) const {
  if (res < 0 || res >= (int)residues_.size()) return -1;
  for (int at = residues_[res].FirstAtom(); at < residues_[res].LastAtom(); ++at)
    if ( atoms_[at].Name() == atname )
      return at;
  return -1;
}

// Topology::FindResidueMaxNatom()
/** Return the # atoms in the largest residue. */
int Topology::FindResidueMaxNatom() const {
  if (residues_.size() <= 1)
    return (int)atoms_.size();
  int largest_natom = 0;
  for (std::vector<Residue>::const_iterator res = residues_.begin();
                                            res != residues_.end(); res++)
  {
    int diff = (*res).NumAtoms();
    if (diff > largest_natom) largest_natom = diff;
  }
  return largest_natom;
}

// -----------------------------------------------------------------------------
// Topology::Summary()
void Topology::Summary() const {
  mprintf("\tTopology %s contains %zu atoms.\n", c_str(), atoms_.size());
  if (!parmName_.empty())
    mprintf("\t\tTitle: %s\n", parmName_.c_str());
  mprintf("\t\t%zu residues.\n", residues_.size());
  mprintf("\t\t%zu molecules.\n", molecules_.size());
  size_t s1 = bondsh_.size();
  size_t s2 = bonds_.size();
  if (s1 + s2 > 0)
    mprintf("\t\t%zu bonds (%zu to H, %zu other).\n", s1+s2, s1, s2);
  s1 = anglesh_.size();
  s2 = angles_.size();
  if (s1 + s2 > 0)
    mprintf("\t\t%zu angles (%zu with H, %zu other).\n", s1+s2, s1 ,s2);
  s1 = dihedralsh_.size();
  s2 = dihedrals_.size();
  if (s1 + s2 > 0)
    mprintf("\t\t%zu dihedrals (%zu with H, %zu other).\n", s1+s2, s1, s2);
  mprintf("\t\tBox: %s\n",box_.TypeName());
  if (NsolventMolecules_>0) {
    mprintf("\t\t%i solvent molecules.\n", NsolventMolecules_);
  }
  if (!radius_set_.empty())
    mprintf("\t\tGB radii set: %s\n", radius_set_.c_str());
  if (chamber_.HasChamber()) {
    mprintf("\t\tCHAMBER: %zu Urey-Bradley terms, %zu Impropers\n",
            chamber_.UB().size(), chamber_.Impropers().size());
    if (chamber_.HasCmap())
      mprintf("\t\t         %zu CMAP grids, %zu CMAP terms.\n", 
              chamber_.CmapGrid().size(), chamber_.Cmap().size());
  }
  if (lesparm_.HasLES())
    mprintf("\t\tLES info: %i types, %i copies\n", lesparm_.Ntypes(), lesparm_.Ncopies());
  if (cap_.HasWaterCap())
    mprintf("\t\tCAP info: Last atom before cap = %s, Cut= %g, X= %g, Y= %g, Z= %g\n",
            AtomMaskName(cap_.NatCap()).c_str(), cap_.CutCap(), 
            cap_.xCap(), cap_.yCap(), cap_.zCap());
}

// Topology::Brief()
void Topology::Brief(const char* heading) const {
  if (heading != 0)
    mprintf("\t%s", heading);
  if (!parmTag_.empty())
    mprintf(" %s", parmTag_.c_str());
  if (!fileName_.empty())
    mprintf(" '%s',", fileName_.base());
  else if (!parmName_.empty())
    mprintf(" %s,", parmName_.c_str());
  mprintf(" %zu atoms, %zu res, box: %s, %zu mol", atoms_.size(), 
          residues_.size(), box_.TypeName(), molecules_.size());
  if (NsolventMolecules_>0)
    mprintf(", %i solvent", NsolventMolecules_);
  if (heading != 0)
    mprintf("\n");
}

// Topology::PrintAtomInfo()
/** Since this function may be called from command line with worldsilent
  * set to true, use loudPrintf and mprinterr.
  */
void Topology::PrintAtomInfo(std::string const& maskString) const {
  AtomMask mask( maskString );
  ParseMask(refCoords_, mask, true); // integer mask
  if ( mask.None() )
    mprinterr("\tSelection is empty.\n");
  else {
    int width = DigitWidth(atoms_.size());
    if (width < 5) width = 5;
    loudPrintf("%-*s %4s %*s %4s %*s %4s %8s %8s %8s %2s\n", 
               width, "#Atom", "Name", 
               width, "#Res",  "Name",
               width, "#Mol",  "Type", "Charge", "Mass", "GBradius", "El");
    for (AtomMask::const_iterator atnum = mask.begin(); atnum != mask.end(); atnum++) {
      const Atom& atom = atoms_[*atnum];
      int resnum = atom.ResNum();
      loudPrintf("%*i %4s %*i %4s %*i %4s %8.4f %8.4f %8.4f %2s\n", 
                 width, *atnum+1, atom.c_str(), 
                 width, resnum+1, residues_[resnum].c_str(),
                 width, atom.MolNum()+1, *(atom.Type()), atom.Charge(), 
                 atom.Mass(), atom.GBRadius(), atom.ElementName());
    }
  }
}

// Topology::PrintBonds()
/** \param maskIn AtomMask which should have already been set up as a char mask
  */
void Topology::PrintBonds(BondArray const& barray, AtomMask const& maskIn, int& nb) const
{
  int rwidth = DigitWidth(residues_.size()) + 7;
  for (BondArray::const_iterator batom = barray.begin();
                                 batom != barray.end(); ++batom)
  {
    int atom1 = (*batom).A1();
    int atom2 = (*batom).A2();
    if (maskIn.AtomInCharMask( atom1 ) || maskIn.AtomInCharMask( atom2 )) {
      mprintf("%8i:", nb);
      int bidx = (*batom).Idx();
      if ( bidx > -1 )
        mprintf(" %6.2f %6.3f", bondparm_[bidx].Rk(), bondparm_[bidx].Req());
      mprintf(" %-*s %-*s (%i,%i)",
              rwidth, AtomMaskName(atom1).c_str(), rwidth, AtomMaskName(atom2).c_str(),
              atom1+1, atom2+1);
      // Atom types
      const char* atype1 = *atoms_[atom1].Type();
      const char* atype2 = *atoms_[atom2].Type();
      mprintf(" %c%c-%c%c\n",atype1[0],atype1[1],atype2[0],atype2[1]);
    }
    nb++;
  }
  mprintf("\n");
}

// Topology::PrintBondInfo()
void Topology::PrintBondInfo(std::string const& maskString) const {
  AtomMask mask( maskString );
  if (ParseMask(refCoords_, mask, false)) return; // Char mask
  mprintf("#");
  mask.MaskInfo();
  if (mask.None()) return;
  mprintf("#   Bond     Kb     Req       atom names   (numbers)\n");
  int nb = 1;
  if (!bondsh_.empty())
    PrintBonds( bondsh_, mask, nb );
  if (!bonds_.empty())
    PrintBonds( bonds_, mask, nb );
}

// Topology::PrintAngles()
void Topology::PrintAngles(AngleArray const& aarray, AtomMask const& maskIn, int& na) const
{
  int rwidth = DigitWidth(residues_.size()) + 7;
  for (AngleArray::const_iterator aatom = aarray.begin();
                                  aatom != aarray.end(); ++aatom)
  {
    int atom1 = (*aatom).A1();
    int atom2 = (*aatom).A2();
    int atom3 = (*aatom).A3();
    if (maskIn.AtomInCharMask( atom1 ) || maskIn.AtomInCharMask( atom2 ) ||
        maskIn.AtomInCharMask( atom3 ))
    {
      mprintf("%8i:", na);
      int aidx = (*aatom).Idx();
      if ( aidx > -1 )
        mprintf(" %6.3f %6.2f", angleparm_[aidx].Tk(), angleparm_[aidx].Teq() * Constants::RADDEG);
      mprintf(" %-*s %-*s %-*s (%i,%i,%i)", rwidth, AtomMaskName(atom1).c_str(), 
              rwidth, AtomMaskName(atom2).c_str(), rwidth, AtomMaskName(atom3).c_str(),
              atom1+1, atom2+1, atom3+1); 
      // Atom types
      const char* atype1 = *atoms_[atom1].Type();
      const char* atype2 = *atoms_[atom2].Type();
      const char* atype3 = *atoms_[atom3].Type();
      mprintf(" %c%c-%c%c-%c%c\n",atype1[0],atype1[1],atype2[0],atype2[1],
              atype3[0],atype3[1]);
    }
    na++;
  }
  mprintf("\n");
}

// Topology::PrintAngleInfo()
void Topology::PrintAngleInfo(std::string const& maskString) const {
  AtomMask mask( maskString );
  if (ParseMask(refCoords_, mask, false)) return; // Char mask
  mprintf("#");
  mask.MaskInfo();
  if (mask.None()) return;
  mprintf("# Angle   Kthet  degrees        atom names        (numbers)\n");
  int na = 1;
  if (!anglesh_.empty())
    PrintAngles( anglesh_, mask, na );
  if (!angles_.empty())
    PrintAngles( angles_, mask, na );
}

// Topology::PrintDihedrals()
void Topology::PrintDihedrals(DihedralArray const& darray, AtomMask const& maskIn, 
                              int& nd) const
{
  int rwidth = DigitWidth(residues_.size()) + 7;
  for (DihedralArray::const_iterator datom = darray.begin();
                                     datom != darray.end(); ++datom)
  {
    int atom1 = (*datom).A1();
    int atom2 = (*datom).A2();
    int atom3 = (*datom).A3();
    int atom4 = (*datom).A4();
    if (maskIn.AtomInCharMask( atom1 ) || maskIn.AtomInCharMask( atom2 ) ||
        maskIn.AtomInCharMask( atom3 ) || maskIn.AtomInCharMask( atom4 )   )
    {
      // Determine dihedral type: 'E'nd, 'I'mproper, or 'B'oth
      char type = ' ';
      if ((*datom).Type() == DihedralType::END) type = 'E';
      else if ((*datom).Type() == DihedralType::IMPROPER) type = 'I';
      else if ((*datom).Type() == DihedralType::BOTH) type = 'B';
      mprintf("%c %8i:", type, nd);
      int didx = (*datom).Idx();
      if ( didx > -1 )
        mprintf(" %6.3f %4.2f %4.1f", dihedralparm_[didx].Pk(), dihedralparm_[didx].Phase(),
                 dihedralparm_[didx].Pn());
      mprintf(" %-*s %-*s %-*s %-*s (%i,%i,%i,%i)",
              rwidth, AtomMaskName(atom1).c_str(), rwidth, AtomMaskName(atom2).c_str(), 
              rwidth, AtomMaskName(atom3).c_str(), rwidth, AtomMaskName(atom4).c_str(),
              atom1+1, atom2+1, atom3+1, atom4+1);
      // Atom types
      const char* atype1 = *atoms_[atom1].Type();
      const char* atype2 = *atoms_[atom2].Type();
      const char* atype3 = *atoms_[atom3].Type();
      const char* atype4 = *atoms_[atom4].Type();
      mprintf(" %c%c-%c%c-%c%c-%c%c\n",atype1[0],atype1[1],atype2[0],atype2[1],
              atype3[0],atype3[1],atype4[0],atype4[1]);
    }
    nd++;
  }
  mprintf("\n");
}

// Topology::PrintDihedralInfo()
void Topology::PrintDihedralInfo(std::string const& maskString) const {
  AtomMask mask( maskString );
  if (ParseMask(refCoords_, mask, false)) return; // Char mask
  mprintf("#");
  mask.MaskInfo();
  if (mask.None()) return;
  mprintf("#Dihedral    pk     phase pn                atoms\n");
  int nd = 1;
  if (!dihedralsh_.empty())
    PrintDihedrals( dihedralsh_, mask, nd );
  if (!dihedrals_.empty())
    PrintDihedrals( dihedrals_, mask, nd );
}


// Topology::PrintMoleculeInfo()
void Topology::PrintMoleculeInfo(std::string const& maskString) const {
  if (molecules_.empty())
    mprintf("\t'%s' No molecule info.\n",c_str());
  else {
    AtomMask mask( maskString );
    ParseMask(refCoords_, mask, false); // Char mask
    if ( mask.None() )
      mprintf("\tSelection is empty.\n");
    else {
      int mwidth = DigitWidth(molecules_.size());
      if (mwidth < 5) mwidth = 5;
      int awidth = DigitWidth(atoms_.size());
      if (awidth < 5) awidth = 5;
      int rwidth = DigitWidth(residues_.size());
      if (rwidth < 5) rwidth = 5;
      mprintf("%-*s %*s %*s %4s\n", mwidth, "#Mol", awidth, "Natom", 
              rwidth, "#Res", "Name");
      unsigned int mnum = 1;
      for (std::vector<Molecule>::const_iterator mol = molecules_.begin(); 
                                                 mol != molecules_.end(); mol++)
      {
        if ( mask.AtomsInCharMask( mol->BeginAtom(), mol->EndAtom() ) ) {
          int firstres = atoms_[ mol->BeginAtom() ].ResNum();
          mprintf("%*u %*i %*i %4s %c", mwidth, mnum, awidth, mol->NumAtoms(),
                  rwidth, firstres+1, residues_[firstres].c_str(),
                  atoms_[mol->BeginAtom()].ChainID());
          if ( mol->IsSolvent() ) mprintf(" SOLVENT");
          mprintf("\n");
        }
        ++mnum;
      }
    }
  }
}

// Topology::PrintResidueInfo()
/** Since this function may be called from command line with worldsilent
  * set to true, use loudPrintf and mprinterr.
  */
void Topology::PrintResidueInfo(std::string const& maskString) const {
  AtomMask mask( maskString );
  ParseMask(refCoords_, mask, true); // Integer mask
  if ( mask.None() )
    mprinterr("\tSelection is empty.\n");
  else {
    int awidth = DigitWidth(atoms_.size());
    if (awidth < 5) awidth = 5;
    int rwidth = DigitWidth(residues_.size());
    if (rwidth < 5) rwidth = 5;
    int mwidth = DigitWidth(molecules_.size());
    if (mwidth < 5) mwidth = 5;
    loudPrintf("%-*s %4s %*s %*s %*s %*s %*s\n", rwidth, "#Res", "Name",
               awidth, "First", awidth, "Last", 
               awidth, "Natom", rwidth, "#Orig", mwidth, "#Mol");
    int rn = -1;
    for (AtomMask::const_iterator atom = mask.begin();
                                  atom != mask.end(); ++atom)
    {
      if (atoms_[*atom].ResNum() > rn) {
        rn = atoms_[*atom].ResNum();
        Residue const& res = residues_[rn];
        loudPrintf("%*i %4s %*i %*i %*i %*i %*i %c\n", rwidth, rn+1, res.c_str(),
                   awidth, res.FirstAtom()+1, awidth, res.LastAtom(),
                   awidth, res.NumAtoms(), rwidth, res.OriginalResNum(),
                   mwidth, atoms_[*atom].MolNum()+1, atoms_[*atom].ChainID());
      }
    }
  }
}

/** Print residue info using single char names. */
void Topology::PrintShortResInfo(std::string const& maskString, int maxChar) const {
  AtomMask mask( maskString );
  ParseMask(refCoords_, mask, true); // Integer mask
  if ( mask.None() )
    mprinterr("\tSelection is empty.\n");
  else {
    // Determine last selected residue.
    int max_res = atoms_[mask.back()].ResNum();
    int total = 0;
    int rn = -1, startRes = -1;
    std::string resLine;
    for (AtomMask::const_iterator atom = mask.begin();
                                  atom != mask.end(); ++atom)
    {
      int current_res = atoms_[*atom].ResNum();
      if (current_res > rn) {
        int n_res_skipped = 1;
        if (startRes == -1)
          startRes = current_res;
        else
          n_res_skipped = current_res - rn;
        // If we skipped any residues print last consective segment and start a new one.
        if (n_res_skipped > 1) {
          mprintf("%-8i %s\n", startRes+1, resLine.c_str());
          startRes = current_res;
          resLine = residues_[current_res].SingleCharName();
          total = 1;
        } else {
          // Convert residue name.
          resLine += residues_[current_res].SingleCharName();
          total++;
        }
        // Print if max line length reached or final res.
        if ((total%maxChar)==0 || current_res == max_res)
        {
          mprintf("%-8i %s\n", startRes+1, resLine.c_str());
          if (current_res == max_res) break;
          startRes = -1;
          resLine.clear();
        } else if ((total % 10) == 0)
          resLine += ' ';
        rn = current_res;
      }
    }
  }
}

// Topology::PrintChargeMassInfo()
int Topology::PrintChargeMassInfo(std::string const& maskString, int type) const {
  AtomMask mask( maskString );
  if (ParseMask(refCoords_, mask, true)) return 1; // Int mask
  if (type == 0 || type == 2) {
    mprintf("\tSum of charges in mask");
    mask.BriefMaskInfo();
    double sumq = 0.0;
    for (AtomMask::const_iterator aidx = mask.begin(); aidx != mask.end(); ++aidx)
      sumq += atoms_[*aidx].Charge();
    mprintf(" is %g\n", sumq);
  }
  if (type == 1 || type == 2) {
    mprintf("\tSum of masses in mask");
    mask.BriefMaskInfo();
    double summ = 0.0;
    for (AtomMask::const_iterator aidx = mask.begin(); aidx != mask.end(); ++aidx)
      summ += atoms_[*aidx].Mass();
    mprintf(" is %g\n", summ);
  }
  return 0; 
}

// -----------------------------------------------------------------------------
// Topology::AddTopAtom()
int Topology::AddTopAtom(Atom const& atomIn, int o_resnum, 
                         NameType const& resname, const double* XYZin)
{
  // If no residues or res num has changed, this is a new residue.
  if ( residues_.empty() || residues_.back().OriginalResNum() != o_resnum )
  {
    // Last atom of old residue is == current # atoms.
    if (!residues_.empty())
      residues_.back().SetLastAtom( atoms_.size() );
    // First atom of new residue is == current # atoms.
    residues_.push_back( Residue(o_resnum, resname, atoms_.size()) );
  }
  atoms_.push_back(atomIn);
  // Set this atoms internal residue number 
  atoms_.back().SetResNum( residues_.size()-1 );
  // Add coordinate if given
  refCoords_.AddXYZ( XYZin );
  return 0;
}

// Topology::StartNewMol()
void Topology::StartNewMol() {
  // If this is the first time this routine has been called, consider all
  // atoms to this point as belonging to first molecule. 
  if (molecules_.empty()) {
    //mprintf("DEBUG:\tFirst molecule, atoms 0 to %zu\n",atoms_.size());
    molecules_.push_back( Molecule(0, atoms_.size()) );
  } else {
    // The first atom of this molecule will be end atom of last molecule.
    int molFirstAtom = molecules_.back().EndAtom();
    // Only add a new molecule if #atoms > first atom of the molecule.
    if ((int)atoms_.size() > molFirstAtom) 
      molecules_.push_back( Molecule( molFirstAtom, atoms_.size()) );
    // First atom
    //mprintf("DEBUG:\tMolecule %zu, atoms %i to %zu\n",
    //       molecules_.size(), lastAtom, atoms_.size());
  }
}

// Topology::CommonSetup()
int Topology::CommonSetup(bool bondsearch) {
  // Set residue last atom (PDB/Mol2/PSF) 
  residues_.back().SetLastAtom( atoms_.size() );
  // Set up bond information if specified and necessary
  if (bondsearch) {
    if (bonds_.empty() && bondsh_.empty() && !refCoords_.empty()) {
      GetBondsFromAtomCoords();
      molecules_.clear();
    }
  }
  // Assign default lengths if necessary (for e.g. CheckStructure)
  if (bondparm_.empty())
    AssignBondParameters();
  // Determine molecule info
  if (molecules_.empty())  
    if (DetermineMolecules()) 
      mprinterr("Error: Could not determine molecule information for %s.\n", c_str());
  // Set up solvent information
  if (SetSolventInfo())
    mprinterr("Error: Could not determine solvent information for %s.\n", c_str());
  // Determine excluded atoms
  DetermineExcludedAtoms();
  // Determine # of extra points.
  DetermineNumExtraPoints();

  return 0;
}

/** For topology formats that do not contain residue info, base residues
  * on molecules.
  */
int Topology::Setup_NoResInfo() {
  mprintf("\tAttempting to determine residue info from molecules.\n");
  if (DetermineMolecules()) {
    mprintf("Warning: Could not determine molecule info. Not setting up residues.\n");
    return 0;
  }
  // Save residue name if its there at all.
  NameType default_res_name, res_name;
  if (!residues_.empty())
    default_res_name = residues_[0].Name();
  else
    default_res_name = "RES";
  // Set residue info to match molecule info.
  residues_.clear();
  int resnum = 0;
  for (std::vector<Molecule>::const_iterator mol = molecules_.begin();
                                             mol != molecules_.end();
                                           ++mol, ++resnum)
  {
    // Try to detect at least water as solvent. Assume CommonSetup will be
    // run after this to set up molecule solvent info.
    if (mol->NumAtoms() == 3) {
      int nH = 0;
      int nO = 0;
      for (int atnum = mol->BeginAtom(); atnum != mol->EndAtom(); atnum++)
      {
        if (atoms_[atnum].Element() == Atom::HYDROGEN) nH++;
        if (atoms_[atnum].Element() == Atom::OXYGEN)   nO++;
      }
      if (nO == 1 && nH == 2) res_name = "HOH";
    } else
      res_name = default_res_name;
    residues_.push_back( Residue(resnum+1, res_name, mol->BeginAtom()) );
    residues_.back().SetLastAtom( mol->EndAtom() );
    // Update atom residue numbers
    for (int atnum = residues_.back().FirstAtom(); 
             atnum != residues_.back().LastAtom(); ++atnum)
      atoms_[atnum].SetResNum( resnum );
  }
  return 0;
}

static inline int NoAtomsErr(const char* msg) {
  mprinterr("Error: Cannot set up %s, no atoms present.\n");
  return 1;
}

// Topology::SetBondInfo()
int Topology::SetBondInfo(BondArray const& bondsIn, BondArray const& bondshIn,
                          BondParmArray const& bondparmIn)
{
  if (atoms_.empty()) return NoAtomsErr("bonds");
  // Clear away bond info from atoms array.
  for (std::vector<Atom>::iterator atom = atoms_.begin(); atom != atoms_.end(); atom++)
    atom->ClearBonds();
  bonds_ = bondsIn;
  bondsh_ = bondshIn;
  // Create bonds in atom array.
  SetAtomBondInfo( bonds_ );
  SetAtomBondInfo( bondsh_ );
  bondparm_ = bondparmIn;

  return 0;
}

// Topology::SetAngleInfo()
int Topology::SetAngleInfo(AngleArray const& anglesIn, AngleArray const& angleshIn,
                           AngleParmArray const& angleparmIn)
{
  if (atoms_.empty()) return NoAtomsErr("angles"); 
  angles_ = anglesIn;
  anglesh_ = angleshIn;
  angleparm_ = angleparmIn;
  return 0;
}

// Topology::SetDihedralInfo()
int Topology::SetDihedralInfo(DihedralArray const& dihedralsIn, DihedralArray const& dihedralshIn,
                              DihedralParmArray const& dihedralparmIn)
{
  if (atoms_.empty()) return NoAtomsErr("dihedrals"); 
  dihedrals_ = dihedralsIn;
  dihedralsh_ = dihedralshIn;
  dihedralparm_ = dihedralparmIn;
  return 0;
}

/** This is for any extra information that is not necessarily pertinent to
  * all topologies, like Ambers ITREE or PDB chain ID etc.
  */
int Topology::SetExtraAtomInfo(int natyp, std::vector<AtomExtra> const& extraIn) 
{
  n_atom_types_ = natyp;
  if (!extraIn.empty()) {
    if (extraIn.size() != atoms_.size()) {
      mprinterr("Error: Size of extra atom info (%zu) != # atoms (%zu)\n",
                 extraIn.size(), atoms_.size());
      return 1;
    }
    extra_ = extraIn;
  }
  return 0;
}

// Topology::SetNonbondInfo()
int Topology::SetNonbondInfo(NonbondParmType const& nonbondIn) {
  nonbond_ = nonbondIn;
  return 0;
}

// Topology::SetAtomBondInfo()
/** Set up bond information in the atoms array based on given BondArray.
  */
void Topology::SetAtomBondInfo(BondArray const& bonds) {
  // Add bonds based on array 
  for (BondArray::const_iterator bnd = bonds.begin(); bnd != bonds.end(); ++bnd) {
    atoms_[ bnd->A1() ].AddBond( bnd->A2() );
    atoms_[ bnd->A2() ].AddBond( bnd->A1() );
  }
}

// -----------------------------------------------------------------------------
// Topology::AddBondParam()
/** Create parameters for given bond based on element types. */
void Topology::AddBondParam(BondType& bnd, BP_mapType& bpMap)
{
  unsigned int bp_idx;
  Atom::AtomicElementType a1Elt = atoms_[bnd.A1()].Element();
  Atom::AtomicElementType a2Elt = atoms_[bnd.A2()].Element();
  std::set<Atom::AtomicElementType> Eset;
  Eset.insert( a1Elt );
  Eset.insert( a2Elt );
  // Has this bond parameter been defined?
  BP_mapType::iterator bp = std::find(bpMap.begin(), bpMap.end(), Eset);
  if (bp == bpMap.end()) { // Bond parameter Not defined
    bp_idx = bondparm_.size();
    bpMap.push_back( Eset );
    bondparm_.push_back( BondParmType(0.0, Atom::GetBondLength(a1Elt, a2Elt)) );
  } else
    bp_idx = bp - bpMap.begin();
  //mprintf("DEBUG:\t\t%i:[%s] -- %i:[%s] Cut=%f BPidx=%u\n",
  //        bnd.A1()+1, atoms_[bnd.A1()].c_str(), bnd.A2()+1, atoms_[bnd.A2()].c_str(),
  //        bondparm_[bp_idx].Req(), bp_idx);
  bnd.SetIdx( bp_idx );
}

// Topology::AssignBondParameters()
void Topology::AssignBondParameters() {
  mprintf("Warning: %s: Determining default bond distances from element types.\n", c_str());
  bondparm_.clear();
  // Hold indices into bondparm for unique element pairs
  BP_mapType bpMap;
  for (BondArray::iterator bnd = bondsh_.begin(); bnd != bondsh_.end(); ++bnd)
    AddBondParam( *bnd, bpMap ); 
  for (BondArray::iterator bnd = bonds_.begin(); bnd != bonds_.end(); ++bnd)
    AddBondParam( *bnd, bpMap );
} 

// Topology::GetBondsFromAtomCoords()
void Topology::GetBondsFromAtomCoords() {
  mprintf("\t%s: determining bond info from distances.\n",c_str());
  // ----- STEP 1: Determine bonds within residues
  for (std::vector<Residue>::iterator res = residues_.begin(); 
                                      res != residues_.end(); ++res) 
  {
    // Get residue start atom.
    int startatom = (*res).FirstAtom();
    // Get residue end atom.
    int stopatom = (*res).LastAtom();
    // Check for bonds between each atom in the residue.
    for (int atom1 = startatom; atom1 < stopatom - 1; ++atom1) {
      Atom::AtomicElementType a1Elt = atoms_[atom1].Element();
      // If this is a hydrogen and it already has a bond, move on.
      if (a1Elt==Atom::HYDROGEN && atoms_[atom1].Nbonds() > 0 )
        continue;
      for (int atom2 = atom1 + 1; atom2 < stopatom; ++atom2) {
        Atom::AtomicElementType a2Elt = atoms_[atom2].Element();
        double D2 = DIST2_NoImage(refCoords_.XYZ(atom1), refCoords_.XYZ(atom2) );
        double cutoff2 = Atom::GetBondLength(a1Elt, a2Elt) + offset_;
        cutoff2 *= cutoff2;
        if (D2 < cutoff2) {
          AddBond(atom1, atom2);
          // Once a bond has been made to hydrogen move on.
          if (a1Elt==Atom::HYDROGEN) break;
        }
      }
    }
  }
  // ----- STEP 2: Determine bonds between adjacent residues
  std::vector<Molecule>::iterator nextmol = molecules_.begin();
  if (!molecules_.empty())
    ++nextmol;
  for (std::vector<Residue>::iterator res = residues_.begin() + 1;
                                      res != residues_.end(); ++res)
  {
    // If molecule information is already present, check if first atom of 
    // this residue >= first atom of next molecule, which indicates this
    // residue and the previous residue are in different molecules.
    if ( (nextmol != molecules_.end()) && 
         ((*res).FirstAtom() >= (*nextmol).BeginAtom()) )
    {
      ++nextmol;
      continue;
    }
    // If this residue is recognized as solvent, no need to check previous or
    // next residue
    if ( (*res).NameIsSolvent() ) {
      ++res;
      if (res == residues_.end()) break;
      continue;
    }
    // Get previous residue
    std::vector<Residue>::iterator previous_res = res - 1;
    // If previous residue is recognized as solvent, no need to check previous.
    if ( (*previous_res).NameIsSolvent() ) continue;
    // Get previous residue start atom
    int startatom = (*previous_res).FirstAtom();
    // Previous residue stop atom, this residue start atom
    int midatom = (*res).FirstAtom();
    // This residue stop atom
    int stopatom = (*res).LastAtom();
    // Check for bonds between adjacent residues
    for (int atom1 = startatom; atom1 < midatom; atom1++) {
      Atom::AtomicElementType a1Elt = atoms_[atom1].Element();
      if (a1Elt==Atom::HYDROGEN) continue;
      for (int atom2 = midatom; atom2 < stopatom; atom2++) {
        Atom::AtomicElementType a2Elt = atoms_[atom2].Element();
        if (a2Elt==Atom::HYDROGEN) continue;
        double D2 = DIST2_NoImage(refCoords_.XYZ(atom1), refCoords_.XYZ(atom2) );
        double cutoff2 = Atom::GetBondLength(a1Elt, a2Elt) + offset_;
        cutoff2 *= cutoff2;
        if (D2 < cutoff2) 
          AddBond(atom1, atom2);
      }
    }
  }
  if (debug_>0)
    mprintf("\t%s: %zu bonds to hydrogen, %zu other bonds.\n",c_str(), 
            bondsh_.size(), bonds_.size());
}

// Topology::AddBond()
/** Create a bond between atom1 and atom2, update the atoms array.
  * For bonds to H always insert the H second.
  */
void Topology::AddBond(int atom1, int atom2) {
  // Check for duplicate bond
  for (Atom::bond_iterator ba = atoms_[atom1].bondbegin();
                           ba != atoms_[atom1].bondend(); ++ba)
    if ( *ba == atom2 ) {
      mprintf("Warning: Bond between atoms %i and %i already exists.\n", atom1+1, atom2+1);
      return;
    }
  bool a1H = (atoms_[atom1].Element() == Atom::HYDROGEN);
  bool a2H = (atoms_[atom2].Element() == Atom::HYDROGEN);
  //mprintf("\t\t\tAdding bond %i to %i (isH=%i)\n",atom1+1,atom2+1,(int)isH);
  // Update bonds arrays
  if (a1H || a2H) {
    if (a1H)
      bondsh_.push_back( BondType(atom2, atom1, -1) );
    else
      bondsh_.push_back( BondType(atom1, atom2, -1) );
  } else
    bonds_.push_back( BondType( atom1, atom2, -1 ) );
  // Update atoms
  atoms_[atom1].AddBond( atom2 );
  atoms_[atom2].AddBond( atom1 );
}

// Topology::VisitAtom()
// NOTE: Use iterator instead of atom num?
void Topology::VisitAtom(int atomnum, int mol) {
  // Return if this atom already has a molecule number
  if (!atoms_[atomnum].NoMol()) return;
  // Mark this atom as visited
  atoms_[atomnum].SetMol( mol );
  // Visit each atom bonded to this atom
  for (Atom::bond_iterator bondedatom = atoms_[atomnum].bondbegin();
                           bondedatom != atoms_[atomnum].bondend(); bondedatom++)
    VisitAtom(*bondedatom, mol);
}

// Topology::DetermineMolecules()
/** Determine individual molecules using bond information. Performs a 
  * recursive search over the bonds of each atom.
  */
int Topology::DetermineMolecules() {
  std::vector<Atom>::iterator atom;
  // Since this is always done only print when debugging
  if (debug_>0) mprintf("\t%s: determining molecule info from bonds.\n",c_str());
  // Reset molecule info for each atom
  for (atom = atoms_.begin(); atom != atoms_.end(); atom++)
    (*atom).SetMol( -1 );
  // Perform recursive search along bonds of each atom
  int mol = 0;
  int atomnum = 0;
  for (atom = atoms_.begin(); atom != atoms_.end(); atom++)
  {
    if ( (*atom).NoMol() ) {
      VisitAtom( atomnum, mol );
      ++mol;
    }
    ++atomnum;
  }
  if (debug_>0) {
    mprintf("\t%i molecules.\n",mol);
    if (debug_ > 1)
    for (atom = atoms_.begin(); atom != atoms_.end(); ++atom)
      mprintf("\t\tAtom %i assigned to molecule %i\n", atom - atoms_.begin(), (*atom).MolNum());
  }

  // Update molecule information
  molecules_.resize( mol );
  if (mol == 0) return 0;
  std::vector<Molecule>::iterator molecule = molecules_.begin();
  (*molecule).SetFirst(0);
  atom = atoms_.begin(); 
  int lastMol = (*atom).MolNum();
  int atomNum = 0;
  for (; atom != atoms_.end(); atom++)
  {
    if ( (*atom).MolNum() > lastMol ) {
      // Set last atom of molecule
      (*molecule).SetLast( atomNum );
      // Set first atom of next molecule
      ++molecule;
      (*molecule).SetFirst( atomNum );
      lastMol = (*atom).MolNum();
    } else if ( (*atom).MolNum()  < lastMol) {
      mprinterr("Error: Atom %u was assigned a lower molecule # than previous atom. This can\n"
                "Error:   happen if 1) bond information is incorrect or missing, or 2) if the\n"
                "Error:   atom numbering in molecules is not sequential. If topology did not\n"
                "Error:   originally contain bond info, 1) can potentially be fixed by\n"
                "Error:   increasing the bondsearch cutoff offset (currently %.3f). 2) can be\n"
                "Error:   fixed by either using the 'fixatomorder' command, or using\n"
                "Error:   the 'setMolecules' command in parmed.py.\n",
                atom - atoms_.begin() + 1, offset_);
      molecules_.clear();
      // Reset molecule info for each atom
      for (atom = atoms_.begin(); atom != atoms_.end(); atom++)
        (*atom).SetMol( -1 );
      return 1;
    }
    ++atomNum;
  }
  (*molecule).SetLast( atoms_.size() );
  return 0;
}

// Topology::AtomDistance()
void Topology::AtomDistance(int originalAtom, int atom, int dist, std::set<int> &excluded) const 
{
  // If this atom is already too far away return
  if (dist==4) return;
  // dist is less than 4 and this atom greater than original, add exclusion
  if (atom > originalAtom)
    excluded.insert( atom ); 
  // Visit each atom bonded to this atom
  for (Atom::bond_iterator bondedatom = atoms_[atom].bondbegin();
                           bondedatom != atoms_[atom].bondend();
                           bondedatom++)
    AtomDistance(originalAtom, *bondedatom, dist+1, excluded);
}

// Topology::DetermineExcludedAtoms()
/** For each atom, determine which atoms with greater atom# are within
  * 4 bonds (and therefore should be excluded from a non-bonded calc).
  */
void Topology::DetermineExcludedAtoms() {
  // A set is used since it automatically sorts itself and rejects duplicates.
  std::set<int> excluded_i;
  int natom = (int)atoms_.size();
  for (int atomi = 0; atomi < natom; atomi++) {
    excluded_i.clear();
    //mprintf("    Determining excluded atoms for atom %i\n",atomi+1);
    // AtomDistance recursively sets each atom bond distance from atomi
    AtomDistance(atomi, atomi, 0, excluded_i);
    atoms_[atomi].AddExclusionList( excluded_i );
    // DEBUG
    //mprintf("\tAtom %i Excluded:",atomi+1);
    //for (Atom::excluded_iterator ei = atoms_[atomi].excludedbegin(); 
    //                             ei != atoms_[atomi].excludedend(); ++ei)
    //  mprintf(" %i",*ei + 1);
    //mprintf("\n");
  } // END loop over atomi
}

// Topology::DetermineNumExtraPoints()
void Topology::DetermineNumExtraPoints() {
  n_extra_pts_ = 0;
  for (std::vector<Atom>::const_iterator atom = atoms_.begin();
                                         atom != atoms_.end(); ++atom)
    if ( (*atom).Element() == Atom::EXTRAPT ) ++n_extra_pts_;
}

// -----------------------------------------------------------------------------
// Topology::SetSolvent()
/** Set solvent information from atom mask. */
int Topology::SetSolvent(std::string const& maskexpr) {
  // Require molecule information
  if (molecules_.empty()) {
    mprinterr("Error: SetSolvent [%s]: No molecule information.\n", c_str());
    return 1;
  }
  // If maskexpr is empty this means remove all solvent information.
  if (maskexpr.empty()) {
    mprintf("Warning: Removing all solvent information from %s\n", c_str());
    for (std::vector<Molecule>::iterator mol = molecules_.begin(); 
                                         mol != molecules_.end(); ++mol)
      mol->SetNoSolvent();
    NsolventMolecules_ = 0;
    return 0;
  }
  // Setup mask
  AtomMask mask( maskexpr );
  SetupCharMask( mask );
  if (mask.None()) {
    mprinterr("Error: SetSolvent [%s]: Mask %s selects no atoms.\n", c_str(), maskexpr.c_str());
    return 1;
  }
  // Loop over all molecules
  NsolventMolecules_ = 0;
  int numSolvAtoms = 0;
  for (std::vector<Molecule>::iterator mol = molecules_.begin();
                                       mol != molecules_.end(); ++mol)
  {
    // Reset old solvent information.
    mol->SetNoSolvent();
    // If any atoms in this molecule are selected by mask, make entire
    // molecule solvent.
    for (int atom = mol->BeginAtom(); atom < mol->EndAtom(); ++atom) {
      if ( mask.AtomInCharMask( atom ) ) {
        mol->SetSolvent();
        ++NsolventMolecules_;
        numSolvAtoms += mol->NumAtoms();
        break;
      }
    }
  }

  mprintf("\tSolvent Mask [%s]: %i solvent molecules, %i solvent atoms\n",
          maskexpr.c_str(), NsolventMolecules_, numSolvAtoms);
  return 0;
}

// Topology::SetSolventInfo()
/** Determine which molecules are solvent based on residue name. */
int Topology::SetSolventInfo() {
  // Require molecule information
  if (molecules_.empty()) {
    mprinterr("Error: SetSolventInfo: No molecule information.\n");
    return 1;
  }
  // Loop over each molecule. Check if first residue of molecule is solvent.
  NsolventMolecules_ = 0;
  int numSolvAtoms = 0;
  for (std::vector<Molecule>::iterator mol = molecules_.begin();
                                       mol != molecules_.end(); mol++)
  {
    int firstRes = atoms_[ mol->BeginAtom() ].ResNum();
    if ( residues_[firstRes].NameIsSolvent() ) {
      mol->SetSolvent();
      ++NsolventMolecules_;
      numSolvAtoms += mol->NumAtoms();
    }
  }

  if (debug_>0) {
    if (NsolventMolecules_ == 0) 
      mprintf("    No solvent.\n");
    else
      mprintf("    %i solvent molecules, %i solvent atoms\n",NsolventMolecules_,numSolvAtoms);
  }
  return 0;
}

// -----------------------------------------------------------------------------
// Topology::SetupIntegerMask()
bool Topology::SetupIntegerMask(AtomMask &mask) const { 
  return ParseMask(refCoords_, mask, true);
}

// Topology::SetupCharMask()
bool Topology::SetupCharMask(AtomMask &mask) const {
  return ParseMask(refCoords_, mask, false);
}

// Topology::SetupIntegerMask()
bool Topology::SetupIntegerMask(AtomMask &mask, Frame const& frame) const {
  if (frame.empty()) return ParseMask(refCoords_, mask, true);
  return ParseMask( frame, mask, true );
}

// Topology::SetupCharMask()
bool Topology::SetupCharMask(AtomMask &mask, Frame const& frame) const {
  if (frame.empty()) return ParseMask(refCoords_, mask, false);
  return ParseMask( frame, mask, false );
}

// Topology::Mask_SelectDistance()
int Topology::Mask_SelectDistance( Frame const& REF, char *mask, bool within, 
                                    bool byAtom, double distance ) const 
{
  int endatom, resi;
  bool selectresidue;
  int atomi, idx, atomj;
  double d2;
  const double* i_crd;

  if (REF.empty()) {
    mprinterr("Error: No reference set for [%s], cannot select by distance.\n",c_str());
    return 1;
  }
  // Distance has been pre-squared.
  // Create temporary array of atom #s currently selected in mask. Also
  // reset mask, it will be the output mask.
  std::vector<unsigned int> selected;
  for (unsigned int i = 0; i < atoms_.size(); i++) {
    if (mask[i]=='T') {
      selected.push_back( i );
      mask[i] = 'F';
    }
  }
  if (selected.empty()) {
    mprinterr("Error: SelectAtomsWithin(%f): No atoms in prior selection.\n",distance);
    return 1;
  }
/*  if (debug_ > 1) {
    mprintf("\t\t\tDistance Op: Within=%i  byAtom=%i  distance^2=%lf\n",
            (int)within, (int)byAtom, distance);
    mprintf("\t\t\tInitial Mask=[");
    for (std::vector<unsigned int>::iterator at = selected.begin(); at != selected.end(); at++)
      mprintf(" %u",*at + 1);
    mprintf(" ]\n");
  }*/

  if (byAtom) { // Select by atom
    // Loop over all atoms
    int n_of_atoms = (int)atoms_.size();
#ifdef _OPENMP
#pragma omp parallel private(atomi, idx, atomj, d2, i_crd)
{
#pragma omp for
#endif
    for (atomi = 0; atomi < n_of_atoms; atomi++) {
      // No need to calculate if atomi already selected
      if (mask[atomi] == 'T') continue;
      // Loop over initially selected atoms
      i_crd = REF.XYZ( atomi );
      for (idx = 0; idx < (int)selected.size(); idx++) {
        atomj = selected[idx];
        d2 = DIST2_NoImage(i_crd, REF.XYZ(atomj));
        if (within) {
          if (d2 < distance) {
            mask[atomi] = 'T';
            break;
          }
        } else {
          if (d2 > distance) {
            mask[atomi] = 'T';
            break;
          }
        }
      }
    }
#ifdef _OPENMP
} // END pragma omp parallel
#endif
  } else { // Select by residue
    int n_of_res = (int)residues_.size();
#ifdef _OPENMP
#pragma omp parallel private(atomi, idx, atomj, d2, resi, selectresidue, endatom, i_crd)
{
#pragma omp for
#endif
    for (resi = 0; resi < n_of_res; resi++) {
      selectresidue = false;
      // Determine end atom for this residue
      endatom = residues_[resi].LastAtom();
      // Loop over mask atoms
      for (idx = 0; idx < (int)selected.size(); idx++) {
        atomj = selected[idx];
        i_crd = REF.XYZ( atomj );
        // Loop over residue atoms
        for (atomi = residues_[resi].FirstAtom(); atomi < endatom; atomi++) {
          d2 = DIST2_NoImage(REF.XYZ(atomi), i_crd);
          if (within) {
            if (d2 < distance) selectresidue = true;
          } else {
            if (d2 > distance) selectresidue = true;
          }
          if (selectresidue) break; 
        }
        if (selectresidue) break;
      }
      if (selectresidue) {
        for (atomi = residues_[resi].FirstAtom(); atomi < endatom; atomi++)
          mask[atomi] = 'T';
        continue;
      }
    }
#ifdef _OPENMP
} // END pragma omp parallel
#endif
  }
  return 0;
}

// Topology::Mask_AND()
void Topology::Mask_AND(char *mask1, char *mask2) const {
  //mprintf("\t\t\tPerforming AND on masks.\n");
  for (unsigned int i = 0; i < atoms_.size(); i++) {
    //mprintf(" [%c|%c]",mask1[i],mask2[i]);
    if (mask1[i]=='F' || mask2[i]=='F')
      mask1[i] = 'F';
    // Otherwise mask1 should already be T
  }
  //mprintf("\n");
}

// Topology::Mask_OR()
void Topology::Mask_OR(char *mask1, char *mask2) const {
  //mprintf("\t\t\tPerforming OR on masks.\n");
  for (unsigned int i = 0; i < atoms_.size(); i++) {
    if (mask1[i]=='T' || mask2[i]=='T')
      mask1[i] = 'T';
    else
      mask1[i] = 'F';
  }
}

// Topology::Mask_NEG()
void Topology::Mask_NEG(char *mask1) const {
  //mprintf("\t\t\tNegating mask.\n");
  for (unsigned int i = 0; i < atoms_.size(); i++) {
    if (mask1[i]=='T')
      mask1[i] = 'F';
    else
      mask1[i] = 'T';
  }
}

// Topology::MaskSelectResidues()
void Topology::MaskSelectResidues(NameType const& name, char *mask) const {
  //mprintf("\t\t\tSelecting residues named [%s]\n",*name);
  for (std::vector<Residue>::const_iterator res = residues_.begin();
                                            res != residues_.end(); res++)
  {
    if ( (*res).Name().Match( name ) ) {
      std::fill(mask + (*res).FirstAtom(), mask + (*res).LastAtom(), 'T');
    }
  }
}

// Topology::MaskSelectResidues()
// Mask args expected to start from 1
void Topology::MaskSelectResidues(int res1, int res2, char *mask) const {
  int endatom;
  int nres = (int) residues_.size();
  //mprintf("\t\t\tSelecting residues %i to %i\n",res1,res2);
  // Check start atom. res1 and res2 are checked by MaskToken
  if (res1 > nres) {
    if (debug_>0)
      mprintf("Warning: Select residues: res 1 out of range (%i)\n",res1);
    return;
  }
  // If last res > nres, make it nres
  if ( res2 >= nres )
    endatom = (int)atoms_.size();
  else
    endatom = residues_[res2-1].LastAtom();
  // Select atoms
  std::fill(mask + residues_[res1-1].FirstAtom(), mask + endatom, 'T');
}

// Topology::MaskSelectElements()
void Topology::MaskSelectElements( NameType const& element, char* mask ) const {
  unsigned int m = 0;
  for (std::vector<Atom>::const_iterator atom = atoms_.begin();
                                         atom != atoms_.end(); ++atom)
  {
    NameType atom_element( (*atom).ElementName() );
    if ( atom_element.Match( element ) )
      mask[m] = 'T';
    ++m;
  } 
}

// Topology::MaskSelectTypes()
void Topology::MaskSelectTypes( NameType const& type, char* mask ) const {
  unsigned int m = 0;
  for (std::vector<Atom>::const_iterator atom = atoms_.begin();
                                         atom != atoms_.end(); ++atom)
  {
    if ( (*atom).Type().Match( type ) )
      mask[m] = 'T';
    ++m;
  } 
}

// Topology::MaskSelectAtoms()
void Topology::MaskSelectAtoms( NameType const& name, char *mask) const {
  //mprintf("\t\t\tSelecting atoms named [%s]\n",*name);
  unsigned int m = 0;
  for (std::vector<Atom>::const_iterator atom = atoms_.begin();
                                         atom != atoms_.end(); atom++)
  {
    //mprintf("\t\t\t%u PARM[%s]  NAME[%s]",m,(*atom).c_str(),*name);
    if ( (*atom).Name().Match( name ) )
      mask[m] = 'T';
    //mprintf(" %c\n",mask[m]);
    ++m;
  } 
}

// Topology::MaskSelectAtoms()
// Mask args expected to start from 1
void Topology::MaskSelectAtoms(int atom1, int atom2, char *mask) const {
  int startatom, endatom;
  //mprintf("\t\t\tSelecting atoms %i to %i\n",atom1,atom2);
  if (atom1 > (int)atoms_.size()) {
    if (debug_>0) 
      mprintf("Warning: Select atoms: atom 1 out of range (%i)\n",atom1);
    return;
  }
  startatom = atom1 - 1;
  if (atom2 > (int)atoms_.size()) 
    //mprinterr("Error: Select atoms: atom 2 out of range (%i)\n",atom2)
    endatom = atoms_.size();
  else
    endatom = atom2;
  // Select atoms
  std::fill(mask + startatom, mask + endatom, 'T');
}

// Topology::ParseMask()
bool Topology::ParseMask(Frame const& REF, AtomMask &maskIn, bool intMask) const {
  std::stack<char*> Stack;
  char *pMask = 0; 
  char *pMask2 = 0;
  int err = 0;

  for (AtomMask::token_iterator token = maskIn.begintoken();
                                token != maskIn.endtoken(); ++token)
  {
    if (pMask==0) {
      // Create new blank mask
      pMask = new char[ atoms_.size() ];
      std::fill(pMask, pMask + atoms_.size(), 'F');
    }
    switch ( token->Type() ) {
      case MaskToken::ResNum : 
        MaskSelectResidues( token->Res1(), token->Res2(), pMask );
        break;
      case MaskToken::ResName :
        MaskSelectResidues( token->Name(), pMask );
        break;
      case MaskToken::AtomNum :
        MaskSelectAtoms( token->Res1(), token->Res2(), pMask );
        break;
      case MaskToken::AtomName :
        MaskSelectAtoms( token->Name(), pMask );
        break;
      case MaskToken::AtomType :
        MaskSelectTypes( token->Name(), pMask );
        break;
      case MaskToken::AtomElement :
        MaskSelectElements( token->Name(), pMask );
        break;
      case MaskToken::SelectAll :
        std::fill(pMask, pMask + atoms_.size(), 'T');
        break;
      case MaskToken::OP_AND :
        pMask2 = Stack.top();
        Stack.pop();
        Mask_AND( Stack.top(), pMask2 );
        delete[] pMask2;
        break;
      case MaskToken::OP_OR :
        pMask2 = Stack.top();
        Stack.pop();
        Mask_OR( Stack.top(), pMask2 );
        delete[] pMask2;
        break;
      case MaskToken::OP_NEG :
        Mask_NEG( Stack.top() );
        break;
      case MaskToken::OP_DIST :
        err = Mask_SelectDistance( REF, Stack.top(), token->Within(), 
                                   token->ByAtom(), token->Distance() );
        break;
      default:
        mprinterr("Error: Invalid mask token (Mask [%s], type [%s]).\n",
                  maskIn.MaskString(), token->TypeName() );
    }
    // If an error occurred, exit the loop.
    if (err != 0 ) break;
    // Check if this mask should now go on the stack
    if ( token->OnStack() ) {
      //mprintf("Pushing Mask on stack, last Token [%s]\n",token->TypeName());
      Stack.push( pMask );
      pMask = 0;
    }
  }
  // If pMask is not null it is probably a blank leftover
  if (pMask!=0) delete[] pMask;

  // If stack is empty here there was an error.
  if (Stack.empty()) {
    mprinterr("Error: Could not parse mask [%s].\n",maskIn.MaskString());
    return true;
  }

  // Top of the stack should point to the final mask
  pMask = Stack.top();
  Stack.pop();
  // Stack should be empty now
  if (!Stack.empty()) {
    mprinterr("Error: Mask stack is not empty.\n");
    while (!Stack.empty()) {
      delete[] Stack.top();
      Stack.pop();
    }
    delete[] pMask;
    return true;
  }
  if (err == 0) {
    if (intMask)
      maskIn.SetupIntMask( pMask, atoms_.size(), debug_ );
    else
      maskIn.SetupCharMask( pMask, atoms_.size(), debug_);
  }
  delete[] pMask;
  return (err != 0); // false if no error occurred
}

// -----------------------------------------------------------------------------
void Topology::ScaleDihedralK(double scale_factor) {
  for (DihedralParmArray::iterator dk = dihedralparm_.begin();
                                   dk != dihedralparm_.end(); ++dk)
    (*dk).Pk() *= scale_factor;
}

// Topology::ModifyByMap()
/** \return Pointer to new Topology based on this Topology, deleting atoms
  *         that are not in the given map (Map[newatom] = oldatom).
  */
Topology* Topology::ModifyByMap(std::vector<int> const& MapIn, bool setupFullParm) const {
  Topology *newParm = new Topology();

  newParm->parmName_ = parmName_;
  newParm->fileName_ = fileName_;
  // NOTE: Do NOT copy tag to avoid duplication.
  newParm->radius_set_ = radius_set_;
  newParm->debug_ = debug_;
  newParm->hasVelInfo_ = hasVelInfo_;
  newParm->nRepDim_ = nRepDim_;
  newParm->n_atom_types_ = n_atom_types_;

  // Reverse Atom map
  std::vector<int> atomMap( atoms_.size(),-1 );

  // Copy atoms from this parm that are in Mask to newParm.
  int oldres = -1;
  // TODO: Check the map size
  for (int newatom = 0; newatom < (int)MapIn.size(); newatom++) {
    int oldatom = MapIn[ newatom ];
    if (oldatom < 0) continue;
    // Store map of oldatom to newatom
    atomMap[oldatom] = newatom;
    // Copy oldatom 
    Atom newparmAtom = atoms_[oldatom];
    // Save oldatom residue number
    int curres = newparmAtom.ResNum();
    // Check if this old atom is in a different residue than the last. If so,
    // set new residue information.
    if ( curres != oldres ) {
      if (!newParm->residues_.empty())
        newParm->residues_.back().SetLastAtom( newatom );
      newParm->residues_.push_back( Residue(residues_[curres].OriginalResNum(),
                                            residues_[curres].Name(), newatom) );
      oldres = curres;
    }
    // Clear bond information from new atom
    newparmAtom.ClearBonds();
    // Set new atom num and residue num
    newparmAtom.SetResNum( newParm->residues_.size() - 1 );
    // Place new atom in newParm
    newParm->atoms_.push_back( newparmAtom );
  }
  if (newParm->atoms_.empty()) {
    mprintf("Warning: All atoms have been stripped.\n");
    return newParm;
  }

  // Set last residue last atom
  newParm->residues_.back().SetLastAtom( newParm->atoms_.size() );

  // Copy reference if present
  if (!refCoords_.empty()) {
    newParm->refCoords_.SetupFrameM( atoms_ );
    newParm->refCoords_.ModifyByMap( refCoords_, MapIn );
  }

  // NOTE: Since in the bond/angle/dihedral atom arrays the parm indices have 
  //       survived intact we can just include direct copies of all the 
  //       parameter arrays for now. May want to cull unused params later.

  // Set up new bond information
  newParm->bonds_ = StripBondArray( bonds_, atomMap );
  newParm->bondsh_ = StripBondArray( bondsh_, atomMap );
  newParm->SetAtomBondInfo( newParm->bonds_ );
  newParm->SetAtomBondInfo( newParm->bondsh_ );
  std::vector<int> parmMap( bondparm_.size(), -1 ); // Map[oldidx] = newidx
  StripBondParmArray( newParm->bonds_,  parmMap, newParm->bondparm_ );
  StripBondParmArray( newParm->bondsh_, parmMap, newParm->bondparm_ );
  //mprintf("DEBUG: Original bond parm array= %zu, new bond parm array = %zu\n",
  //        bondparm_.size(), newParm->bondparm_.size());
  // Give stripped parm the same pindex as original
  newParm->pindex_ = pindex_;
  newParm->nframes_ = nframes_;
  // Copy box information
  newParm->box_ = box_;
  // If we dont care about setting up full parm information, exit now.
  if (!setupFullParm) return newParm;

  // Set new molecule information based on new bonds
  if (newParm->DetermineMolecules()) {
    mprintf("Warning: Could not set up molecule information for stripped topology %s\n",
            newParm->c_str());
  }
  // Set new solvent information based on new molecules
  if (newParm->SetSolventInfo()) {
    mprintf("Warning: Could not set up solvent information for stripped topology %s\n",
            newParm->c_str());
  } 
  // Set up new angle info
  newParm->angles_ = StripAngleArray( angles_, atomMap );
  newParm->anglesh_ = StripAngleArray( anglesh_, atomMap );
  parmMap.assign( angleparm_.size(), -1 );
  StripAngleParmArray( newParm->angles_,  parmMap, newParm->angleparm_ );
  StripAngleParmArray( newParm->anglesh_, parmMap, newParm->angleparm_ );
  // Set up new dihedral info
  newParm->dihedrals_ = StripDihedralArray( dihedrals_, atomMap );
  newParm->dihedralsh_ = StripDihedralArray( dihedralsh_, atomMap );
  parmMap.assign( dihedralparm_.size(), -1 );
  StripDihedralParmArray( newParm->dihedrals_,  parmMap, newParm->dihedralparm_ );
  StripDihedralParmArray( newParm->dihedralsh_, parmMap, newParm->dihedralparm_ );
  // Set up nonbond info. First determine which atom types remain.
  if (nonbond_.HasNonbond()) {
    parmMap.clear();               // parmMap[oldtype]      = newtype
    std::vector<int> oldTypeArray; // oldTypeArray[newtype] = oldtype
    for (std::vector<Atom>::const_iterator atm = newParm->atoms_.begin();
                                           atm != newParm->atoms_.end(); ++atm)
    {
      int oldidx = atm->TypeIndex();
      if (oldidx >= (int)parmMap.size())
        parmMap.resize( oldidx+1, -1 );
      if (parmMap[oldidx] == -1) {
        parmMap[oldidx] = (int)oldTypeArray.size();
        oldTypeArray.push_back( oldidx );
      }
      //int newidx = parmMap[oldidx];
      //mprintf("DEBUG: '%s' Old type index=%i, new type index = %i\n", atm->c_str(), oldidx, newidx);
    }
    //mprintf("DEBUG: # new types %zu\n", oldTypeArray.size());
    // Set up new nonbond and nonbond index arrays.
    newParm->nonbond_.SetNtypes( oldTypeArray.size() );
    for (int a1idx = 0; a1idx != (int)oldTypeArray.size(); a1idx++)
    {
      int atm1 = oldTypeArray[a1idx];
      for (int a2idx = a1idx; a2idx != (int)oldTypeArray.size(); a2idx++)
      {
        int atm2 = oldTypeArray[a2idx];
        int oldnbidx = nonbond_.GetLJindex( atm1, atm2 );
        // NOTE: Certain routines in sander (like the 1-4 calcs) do NOT use
        //       the nonbond index array; instead they expect the nonbond
        //       arrays to be indexed like '(ibig*(ibig-1)/2+isml)', where
        //       ibig is the larger atom type index.
        int ibig = std::max(a1idx, a2idx) + 1;
        int isml = std::min(a1idx, a2idx) + 1;
        int testidx = (ibig*(ibig-1)/2+isml)-1;
        if (oldnbidx > -1) {
          // This is a traditional LJ 6-12 term. Because of the way the LJ 1-4
          // code is laid out in sander/pmemd the LJ matrix has to be laid out
          // indepdendent of the nonbond index array.
          newParm->nonbond_.AddLJterm( testidx, a1idx, a2idx, nonbond_.NBarray(oldnbidx) );
        } else {
          // This is an old LJ 10-12 hbond term. Add one to the LJ 6-12 matrix
          // and one to the hbond since that seems to be the convention.
          newParm->nonbond_.AddLJterm( testidx, a1idx, a2idx, NonbondType() );
          newParm->nonbond_.AddHBterm( a1idx, a2idx, nonbond_.HBarray((-oldnbidx)-1) );
        }
        //int newnbidx = newParm->nonbond_.GetLJindex( a1idx, a2idx );
        //mprintf("DEBUG: oldtypei=%i oldtypej=%i Old NB index=%i, newtypi=%i newtypej=%i new NB idx=%i testidx=%i\n", 
        //        atm1, atm2, oldnbidx, a1idx, a2idx, newnbidx, testidx);
      }
    }
    // Update atom type indices.
    for (std::vector<Atom>::iterator atm = newParm->atoms_.begin();
                                     atm != newParm->atoms_.end(); ++atm)
      atm->SetTypeIndex( parmMap[atm->TypeIndex()] );
  }
  // LES info - FIXME: Not sure if stripping this is valid so print a warning.
  if (lesparm_.HasLES()) {
    mprintf("Warning: LES info present. Stripped topology may not have correct LES info.\n");
    newParm->lesparm_.SetTypes( lesparm_.Ntypes(), lesparm_.FAC() );
    for (std::vector<int>::const_iterator old_it = MapIn.begin(); old_it != MapIn.end(); ++old_it)
    {
      if (*old_it >= 0)
        newParm->lesparm_.AddLES_Atom( lesparm_.Array()[*old_it] );
    }
  }
  // CAP info - dont support stripping such topologies right now
  if (cap_.HasWaterCap())
    mprintf("Warning: Stripping of CAP info not supported. Removing CAP info.\n");
  // CHAMBER info - Parameters remain intact
  if (chamber_.HasChamber()) {
    newParm->chamber_.SetChamber( chamber_.FF_Version(), chamber_.FF_Type() );
    newParm->chamber_.SetUB( StripBondArray(chamber_.UB(),atomMap), chamber_.UBparm() );
    newParm->chamber_.SetImproper( StripDihedralArray(chamber_.Impropers(),atomMap),
                                   chamber_.ImproperParm() );
    newParm->chamber_.SetLJ14( chamber_.LJ14() );
    if (chamber_.HasCmap()) {
      for (CmapArray::const_iterator cmap = chamber_.Cmap().begin();
                                     cmap != chamber_.Cmap().end(); ++cmap)
      {
        int newA1 = atomMap[ cmap->A1() ];
        if (newA1 != -1) {
          int newA2 = atomMap[ cmap->A2() ];
          if (newA2 != -1) {
            int newA3 = atomMap[ cmap->A3() ];
            if (newA3 != -1) {
              int newA4 = atomMap[ cmap->A4() ];
              if (newA4 != -1) {
                int newA5 = atomMap[ cmap->A5() ];
                if (newA5 != -1)
                  newParm->chamber_.AddCmapTerm( CmapType(newA1,newA2,newA3,
                                                          newA4,newA5,cmap->Idx()) );
              }
            }
          }
        }
      }
      // Only add CMAP grids if there are CMAP terms left.
      if (!newParm->chamber_.Cmap().empty()) {
        for (CmapGridArray::const_iterator g = chamber_.CmapGrid().begin();
                                           g != chamber_.CmapGrid().end(); ++g)
          newParm->chamber_.AddCmapGrid( *g );
      }
    }
  }
  // Amber extra info.
  if (!extra_.empty()) {
    for (std::vector<int>::const_iterator old_it = MapIn.begin(); old_it != MapIn.end(); ++old_it)
      if (*old_it >= 0)
        newParm->extra_.push_back( extra_[*old_it] );
  }
  
  // Setup excluded atoms list - Necessary?
  newParm->DetermineExcludedAtoms();

  // Determine number of extra points
  newParm->DetermineNumExtraPoints();

  return newParm;
}

/** \return BondArray with bonds for which both atoms are still present.
  * \param atomMap format Map[oldAtom]=newAtom
  */
BondArray Topology::StripBondArray(BondArray const& bondsIn, std::vector<int> const& atomMap) const {
  BondArray bondsOut;
  // Go through old array. Use atomMap to determine what goes into newArray.
  for (BondArray::const_iterator oldbond = bondsIn.begin(); oldbond != bondsIn.end(); ++oldbond) {
    int newA1 = atomMap[ oldbond->A1() ];
    if (newA1 != -1) {
      int newA2 = atomMap[ oldbond->A2() ];
      if (newA2 != -1)
        bondsOut.push_back( BondType(newA1, newA2, oldbond->Idx() ) );
    }
  }
  return bondsOut;
}

/** \return AngleArray with angles for which all atoms are still present.
  * \param atomMap format Map[oldAtom]=newAtom
  */
AngleArray Topology::StripAngleArray(AngleArray const& anglesIn, std::vector<int> const& atomMap) const {
  AngleArray anglesOut;
  for (AngleArray::const_iterator oldangle = anglesIn.begin(); oldangle != anglesIn.end(); ++oldangle) {
    int newA1 = atomMap[ oldangle->A1() ];
    if (newA1 != -1) {
      int newA2 = atomMap[ oldangle->A2() ];
      if (newA2 != -1) {
        int newA3 = atomMap[ oldangle->A3() ];
        if (newA3 != -1)
          anglesOut.push_back( AngleType(newA1, newA2, newA3, oldangle->Idx()) );
      }
    }
  }
  return anglesOut;
}

/** \return DihedralArray with dihedrals for which all atoms are still present.
  * \param atomMap format Map[oldAtom]=newAtom
  */
DihedralArray Topology::StripDihedralArray(DihedralArray const& dihIn, std::vector<int> const& atomMap) const {
  DihedralArray dihOut;
  for (DihedralArray::const_iterator olddih = dihIn.begin(); olddih != dihIn.end(); ++olddih) {
    int newA1 = atomMap[ olddih->A1() ];
    if (newA1 != -1) {
      int newA2 = atomMap[ olddih->A2() ];
      if (newA2 != -1) {
        int newA3 = atomMap[ olddih->A3() ];
        if (newA3 != -1) {
          int newA4 = atomMap[ olddih->A4() ];
          if (newA4 != -1) {
            // Since in Amber improper/end dihedrals are stored as negative #s,
            // atom index 0 cannot be in 3rd or 4th position. Reverse.
            if (olddih->Type() != DihedralType::NORMAL && (newA3 == 0 || newA4 == 0))
              dihOut.push_back( DihedralType( newA4, newA3, newA2, newA1, 
                                              olddih->Type(), olddih->Idx() ) );
            else
              dihOut.push_back( DihedralType( newA1, newA2, newA3, newA4, 
                                              olddih->Type(), olddih->Idx() ) );
          }
        }
      }
    }
  }
  return dihOut;
}

// Topology::StripBondParmArray()
void Topology::StripBondParmArray(BondArray& newBondArray, std::vector<int>& parmMap,
                                  BondParmArray& newBondParm) const
{
  for (BondArray::iterator bnd = newBondArray.begin();
                           bnd != newBondArray.end(); ++bnd)
  {
    int oldidx = bnd->Idx();
    int newidx = parmMap[bnd->Idx()];
    if (newidx == -1) { // This needs to be added to new parameter array.
      newidx = (int)newBondParm.size();
      parmMap[oldidx] = newidx;
      newBondParm.push_back( bondparm_[oldidx] );
    }
    //mprintf("DEBUG: Old bond parm index=%i, new bond parm index=%i\n", oldidx, newidx);
    bnd->SetIdx( newidx );
  }
}

// Topology::StripAngleParmArray()
void Topology::StripAngleParmArray(AngleArray& newAngleArray, std::vector<int>& parmMap,
                                   AngleParmArray& newAngleParm) const
{
  for (AngleArray::iterator ang = newAngleArray.begin();
                            ang != newAngleArray.end(); ++ang)
  {
    int oldidx = ang->Idx();
    int newidx = parmMap[ang->Idx()];
    if (newidx == -1) { // This needs to be added to new parameter array.
      newidx = (int)newAngleParm.size();
      parmMap[oldidx] = newidx;
      newAngleParm.push_back( angleparm_[oldidx] );
    }
    //mprintf("DEBUG: Old angle parm index=%i, new angle parm index=%i\n", oldidx, newidx);
    ang->SetIdx( newidx );
  }
}

// Topology::StripDihedralParmArray()
void Topology::StripDihedralParmArray(DihedralArray& newDihedralArray, std::vector<int>& parmMap,
                                      DihedralParmArray& newDihedralParm) const
{
  for (DihedralArray::iterator dih = newDihedralArray.begin();
                               dih != newDihedralArray.end(); ++dih)
  {
    int oldidx = dih->Idx();
    int newidx = parmMap[dih->Idx()];
    if (newidx == -1) { // This needs to be added to new parameter array.
      newidx = (int)newDihedralParm.size();
      parmMap[oldidx] = newidx;
      newDihedralParm.push_back( dihedralparm_[oldidx] );
    }
    //mprintf("DEBUG: Old dihedral parm index=%i, new dihedral parm index=%i\n", oldidx, newidx);
    dih->SetIdx( newidx );
  }
}

// Topology::AddBondArray()
void Topology::AddBondArray(BondArray const& barray, int atomOffset) {
  for (BondArray::const_iterator bond = barray.begin(); bond != barray.end(); ++bond)
    AddBond( bond->A1() + atomOffset, bond->A2() + atomOffset );
}

// Topology::AppendTop()
int Topology::AppendTop(Topology const& CurrentTop) {
  int atomOffset = (int)atoms_.size();
  // ATOMS
  for (atom_iterator atom = CurrentTop.begin(); atom != CurrentTop.end(); ++atom)
  {
    Atom CurrentAtom = *atom;
    Residue const& res = CurrentTop.Res( CurrentAtom.ResNum() );
    // Bonds need to be cleared and re-added.
    CurrentAtom.ClearBonds();
    AddTopAtom( CurrentAtom, res.OriginalResNum(), res.Name(), 0 );
  }
  // BONDS
  AddBondArray(CurrentTop.Bonds(),  atomOffset);
  AddBondArray(CurrentTop.BondsH(), atomOffset);
  // Re-set up this topology
  // TODO: Could get expensive for multiple appends.
  return CommonSetup(false);
}
