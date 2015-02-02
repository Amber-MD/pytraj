#include <algorithm> // sort
#include "AtomMap.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
AtomMap::AtomMap() : debug_(0) {}

/// Blank AtomMap for empty return of bracket operator
MapAtom AtomMap::EMPTYMAPATOM = MapAtom();

// AtomMap::operator[]()
/** Return a reference to the specifed MapAtom */
MapAtom& AtomMap::operator[](int idx) {
  if (idx < 0 || idx >= (int)mapatoms_.size()) {
    mprinterr("Error: AtomMap::operator[]: Index %i out of range.\n",idx);
    return EMPTYMAPATOM;
  }
  return mapatoms_[idx];
}

// AtomMap::InvalidElement()
bool AtomMap::InvalidElement() {
  if (mapatoms_.back().CharName() == 0) {
    mprinterr("Error: AtomMap: Mapping currently not supported for element %s\n",
              mapatoms_.back().ElementName());
    return true;
  }
  return false;
}

// AtomMap::Setup()
int AtomMap::Setup(Topology const& TopIn) {
  mapatoms_.clear();
  for (Topology::atom_iterator atom = TopIn.begin(); atom != TopIn.end(); atom++) {
    // This sets up 1 char atom name based on atom element
    mapatoms_.push_back( *atom );
    if (InvalidElement()) return 1;
  }
  return CheckBonds();
}

// AtomMap::SetupResidue()
int AtomMap::SetupResidue(Topology const& topIn, int resnum) {
  mapatoms_.clear();
  int firstAtom = topIn.Res(resnum).FirstAtom();
  int lastAtom = topIn.Res(resnum).LastAtom();
  //mprintf("DEBUG:\tResidue %i, atoms %i to %i\n", resnum + 1, firstAtom+1, lastAtom);
  for (int atom = firstAtom; atom < lastAtom; ++atom) {
    mapatoms_.push_back( topIn[atom] );
    if (InvalidElement()) return 1;
    // Add bonds for this residue 
    mapatoms_.back().ClearBonds();
    for (Atom::bond_iterator bndatm = topIn[atom].bondbegin();
                             bndatm != topIn[atom].bondend(); ++bndatm)
    {
      //mprintf("DEBUG:\t\tOriginal bond %u-%i", atom+1, *bndatm+1);
      if (*bndatm >= firstAtom && *bndatm < lastAtom) { 
        int newbndatm = *bndatm - firstAtom;
        mapatoms_.back().AddBond(newbndatm);
        //mprintf(", new bond %i-%i", mapatoms_.size(), newbndatm+1);
      }
      //mprintf("\n");
    }
  }
  return CheckBonds();
}

// AtomMap::ResetMapping()
void AtomMap::ResetMapping() {
  for (Marray::iterator matom = mapatoms_.begin();
                        matom != mapatoms_.end(); ++matom)
  {
    matom->SetNotMapped();
    matom->SetNotComplete();
  }
}

// AtomMap::BondIsRepeated()
/** Check if the atomID of the specified atom (bondedAtom) bonded to <atom> 
  * is the same as the atomID of any other non-mapped atom bonded to <atom>.
  */
bool AtomMap::BondIsRepeated(int atom, int bondedAtom) const {
  // If 1 or no bonds, atom cant possibly be repeated
  if (mapatoms_[atom].Nbonds() > 1) {
    for (Atom::bond_iterator bondedAtom2 = mapatoms_[atom].bondbegin();
                             bondedAtom2 != mapatoms_[atom].bondend(); bondedAtom2++)
    {
      if (mapatoms_[*bondedAtom2].IsMapped()) continue;
      if (mapatoms_[bondedAtom].AtomID() == mapatoms_[*bondedAtom2].AtomID())
        return true;
    }
  }
  return false;
}

// AtomMap::DetermineAtomIDs()
/** Give each atom an identifier (atomID) based on what atoms are bonded to 
  * it. The first part of the atomID is the atom itself, followed by an 
  * alphabetized list of bonded atoms. So C in O=C-H2 would be CHHO.
  * Then create a 'unique string' comprised of this atomID and the 
  * atomIDs of all bonded atoms (sorted). Once that is done determine
  * which unique strings are actually unique (i.e. they are not repeated
  * in this map). 
  */
void AtomMap::DetermineAtomIDs() {
  // Determine self IDs
  if (debug_>0) mprintf("ATOM IDs:\n");
  unsigned int anum = 1;
  for (Marray::iterator matom = mapatoms_.begin(); 
                        matom != mapatoms_.end(); ++matom)
  {
    std::string atomID;
    for (Atom::bond_iterator bondedAtom = matom->bondbegin();
                             bondedAtom != matom->bondend(); ++bondedAtom)
      atomID += mapatoms_[ *bondedAtom ].CharName();
    // Sort atom ID
    sort( atomID.begin(), atomID.end() );
    // Place current atom 1 char name at beginning
    atomID = matom->CharName() + atomID;
    matom->SetAtomID( atomID );
    if (debug_>0) mprintf("  Atom %u %4s : %s\n",anum, matom->c_str(), atomID.c_str());
    ++anum;
  }
  
  // Create a unique ID for each atom based on Atom IDs
  for (int ma1 = 0; ma1 < (int)mapatoms_.size(); ++ma1) {
    MapAtom& matom = mapatoms_[ma1];
    std::string unique = matom.AtomID();
    if (mapatoms_.size() > 9)
    {
      for (Atom::bond_iterator bondedAtom = matom.bondbegin();
                               bondedAtom != matom.bondend(); ++bondedAtom)
      {
        unique += mapatoms_[ *bondedAtom ].AtomID();
        // Go one more level through bonds for unique ID.
        // FIXME: This may only be optimal above a certain # of atoms
        MapAtom const& Batom = mapatoms_[ *bondedAtom ];
        for (Atom::bond_iterator ba2 = Batom.bondbegin(); ba2 != Batom.bondend(); ++ba2)
        {
          if (*ba2 != ma1) {
            unique += mapatoms_[ *ba2 ].AtomID();
            // For larger residues go one additional level.
            if (mapatoms_.size() > 19) {
              Atom const& Catom = mapatoms_[ *ba2 ];
              for (Atom::bond_iterator ca3 = Catom.bondbegin(); ca3 != Catom.bondend(); ++ca3)
              {
                if (ca3 != ba2 && *ca3 != ma1)
                  unique += mapatoms_[ *ca3 ].AtomID();
              }
            }
          }
        }
      }
    }
    // Do not sort first character (this atoms element ID char).
    sort( unique.begin() + 1, unique.end() );
    // NOTE: SetUnique also resets the duplicated counter.
    matom.SetUnique( unique );
  }

  // Determine which unique IDs are duplicated - set isUnique flag
  for (unsigned int i = 0; i < mapatoms_.size()-1; i++) {
    for (unsigned int j = i+1; j < mapatoms_.size(); j++) {
      if ( mapatoms_[i].Unique() == mapatoms_[j].Unique() ) {
        // This unique string is duplicated
        mapatoms_[i].IsDuplicated();
        mapatoms_[j].IsDuplicated();
      }
    }
  }

  // DEBUG
  if (debug_ > 0) {
    mprintf("UNIQUE IDs:\n");
    anum = 1;
    for (Marray::const_iterator matom = mapatoms_.begin();
                                matom != mapatoms_.end(); ++matom)
    {
      mprintf("  Atom %6u %4s [%3i]: %s", anum, matom->c_str(), matom->Nduplicated(),
              matom->Unique().c_str());
      if (matom->IsUnique()) mprintf(" UNIQUE!");
      mprintf("\n");
      ++anum;
    }
  }
}

// AtomMap::MarkAtomComplete()
/** If atom is mapped and all bonded atoms are mapped mark atom as completely 
  * mapped.
  * If printAtoms is true print isMapped value for this atom and all atoms
  * bonded to it.
  */
void AtomMap::MarkAtomComplete(int atom, bool printAtoms) {
  if (atom<0 || atom >= (int)mapatoms_.size()) return;
  if (!mapatoms_[atom].IsMapped() && !printAtoms) return;
  if ( mapatoms_[atom].Complete() && !printAtoms) return;
  int nunique = 0;
  for (Atom::bond_iterator bondedAtom = mapatoms_[atom].bondbegin();
                           bondedAtom != mapatoms_[atom].bondend(); bondedAtom++)
    if (mapatoms_[*bondedAtom].IsMapped())
      ++nunique;
  if (mapatoms_[atom].IsUnique() && nunique==mapatoms_[atom].Nbonds())
    mapatoms_[atom].SetComplete();
  if (printAtoms) {
    mprintf("  Atom %4i: [%s]-%1i |",atom+1,mapatoms_[atom].c_str(),
            (int)mapatoms_[atom].IsMapped());
    for (Atom::bond_iterator bondedAtom = mapatoms_[atom].bondbegin();
                           bondedAtom != mapatoms_[atom].bondend(); bondedAtom++)
    {
      mprintf(" %4i:[%s]-%1i",*bondedAtom+1,mapatoms_[*bondedAtom].c_str(),
              (int)mapatoms_[*bondedAtom].IsMapped());
    }
    if (mapatoms_[atom].Complete())
      mprintf(" Atom is completely mapped.");
    mprintf("\n");
  }
}

// AtomMap::CheckForCompleteAtoms()
/** Go through each atom in the map. If the atom is unique and all bonded
  * atoms are unique mark the atom as completely mapped.
  */
void AtomMap::CheckForCompleteAtoms() {
  bool printAtoms = (debug_ > 0);
  for (int atom = 0; atom < (int)mapatoms_.size(); atom++)
    MarkAtomComplete(atom,printAtoms);
}

// AtomMap::CheckBonds()
/** Checks that bonding information is present. Also checks for potential
  * chiral centers.
  */
int AtomMap::CheckBonds() {
  int total_bonds = 0;
  // Search for chiral centers by number of bonds
  for (Marray::iterator matom = mapatoms_.begin();
                        matom != mapatoms_.end(); ++matom)
  {
    // Sort the bonded atoms array by atom #
    matom->SortBonds();
    total_bonds += matom->Nbonds();
    if (matom->Nbonds() == 4) {
      // If >=3 bonds to single atoms, not chiral (e.g. -CH3)
      int N_single_atoms=0; // Count # bonds to single atoms
      for (Atom::bond_iterator bondedAtom = matom->bondbegin();
                               bondedAtom != matom->bondend(); ++bondedAtom)
      {
        if (mapatoms_[*bondedAtom].Nbonds() == 1)
          ++N_single_atoms;
      }
      if (N_single_atoms<3) {
        matom->SetChiral();
        for (Atom::bond_iterator bondedAtom = matom->bondbegin();
                               bondedAtom != matom->bondend(); ++bondedAtom)
          mapatoms_[*bondedAtom].SetBoundToChiral();
      }
    }
  }
  if (total_bonds == 0) {
    mprinterr("Error: No bond information present, required by AtomMap.\n");
    return 1;
  }

  // DEBUG
  if (debug_>0) {
    mprintf("AtomMap: Atom Bond information.\n");
    unsigned int anum = 1;
    for (Marray::const_iterator matom = mapatoms_.begin();
                                matom != mapatoms_.end(); ++matom)
    {
      mprintf("  Atom %s(%c)_%i has %i bonds.",matom->c_str(),matom->CharName(),
              anum, matom->Nbonds());
      if (matom->IsChiral()) mprintf(" CHIRAL");
      if (matom->BoundToChiral()) mprintf(" BOUND TO CHIRAL");
      mprintf("\n");
      for (Atom::bond_iterator bondedAtom = matom->bondbegin();
                               bondedAtom != matom->bondend(); ++bondedAtom)
      {
        mprintf("    to %s(%c)_%i\n",mapatoms_[*bondedAtom].c_str(),
                mapatoms_[*bondedAtom].CharName(), *bondedAtom+1);
      }
    }
  }
  return 0;
}

// -----------------------------------------------------------------------------
#ifdef DEBUG_ATOMMAP
static int recursionLevel_;
#endif
/** Recursive function to search for symmetric atoms.
  * \param at Atom to start search at.
  * \param Unique ID string of starting atom.
  * \param Selected 1 if at already visited, 0 if not.
  * \param symmatoms Output array of potentially symmetric atoms.
  */
void AtomMap::FindSymmetricAtoms(int at, std::string const& Unique,
                                 Iarray& Selected, Iarray& symmatoms) const
{
  // If this atom has already been selected, leave
  if (Selected[at]) return;
  Selected[at] = 1;
# ifdef DEBUG_ATOMMAP
  ++recursionLevel_;
  for (int i = 0; i < recursionLevel_; i++)
    mprintf("..");
  mprintf("Atom %i(%s)", at+1, mapatoms_[at].c_str());
# endif
  // Does this atom match the unique ID we are looking for?
  if (mapatoms_[at].Unique() == Unique) {
    symmatoms.push_back( at ); // NOTE: This index is relative to the residue 
#   ifdef DEBUG_ATOMMAP
    mprintf(" SYMM");
#   endif
  }
  // Recursively search through all atoms bonded to this atom unless they 
  // are a chiral center
# ifdef DEBUG_ATOMMAP
  mprintf("\n");
# endif
  for (Atom::bond_iterator bndatm = mapatoms_[at].bondbegin();
                           bndatm != mapatoms_[at].bondend();
                           ++bndatm)
  {
    if (!mapatoms_[*bndatm].IsChiral())
      FindSymmetricAtoms( *bndatm, Unique, Selected, symmatoms );
  }
}

/** Find groups of potentially symmetric atoms in residue.
  * \param topIn Topology.
  * \param SymmetricAtomIndices Groups of symmetric atoms will be added to this array.
  * \param res Residue # to search for symmetric atoms in.
  */
int AtomMap::SymmetricAtoms(Topology const& topIn,
                            AtomIndexArray& SymmetricAtomIndices,
//                          AtomMask const& tgtMask,
                            int res)
{
  enum atomStatusType { UNSELECTED = 0, NONSYMM, SYMM };
  int res_first_atom = topIn.Res(res).FirstAtom();
  // Are any of the residue atoms selected?
//  bool atomsAreSelected = false;
//  for (int ratom = res_first_atom; ratom != topIn.Res(res).LastAtom(); ++ratom)
//    if (SelectedIdx[ratom] != -1) {
//      atomsAreSelected = true;
//      break;
//    }
//  if (!atomsAreSelected) continue;
  if (debug_>0) mprintf("DEBUG: Residue %s\n", topIn.TruncResNameNum(res).c_str());
  // Create AtomMap of this residue to determine chiral centers, unique atom IDs etc
  if (SetupResidue(topIn, res) != 0) return 1;
  DetermineAtomIDs();
  // Potentially symmetric atom group; indices relative to this AtomMap.
  Iarray symmAtoms;
  // Symmetric atom group; indices relative to Topology.
  Iarray selectedSymmAtoms;
  // Current status of atoms in this residue.
  Iarray AtomStatus( Natom(), UNSELECTED );
  // Loop over all atoms in the residue
  for (int at = 0; at < Natom(); at++) {
    // If atom is unique in residue, mark non-symmetric 
    if (mapatoms_[at].IsUnique())
      AtomStatus[at] = NONSYMM;
    else if (AtomStatus[at] != SYMM) {
      Iarray Selected( Natom(), 0 );
      symmAtoms.clear();
      // Recursively search for other potentially symmetric atoms in residue.
      // The Selected array is used to keep track of which atoms have been
      // visited in this pass; this is used instead of AtomStatus so that
      // we can travel through atoms already marked as symmetric.
#     ifdef DEBUG_ATOMMAP
      recursionLevel_ = 0;
      mprintf("Starting recursive call for %i(%s)\n", at+1, mapatoms_[at].c_str());
#     endif
      FindSymmetricAtoms(at, mapatoms_[at].Unique(), Selected, symmAtoms);
#     ifdef DEBUG_ATOMMAP
      mprintf("Potentially symmetric:\n");
      for (Iarray::const_iterator sa = symmAtoms.begin(); sa != symmAtoms.end(); ++sa)
        mprintf("\t%8i %4s %8i\n", *sa + res_first_atom + 1,
                topIn[*sa + res_first_atom].c_str(),
                AtomStatus[ *sa ]);
#     endif
      // If only one atom, not symmetric.
      if (symmAtoms.size() == 1)
        AtomStatus[symmAtoms.front()] = NONSYMM;
      else if (symmAtoms.size() > 1) {
        // Store correct atom #s in selectedSymmAtoms
        selectedSymmAtoms.clear();
        for (Iarray::const_iterator sa = symmAtoms.begin();
                                    sa != symmAtoms.end(); ++sa)
        {
          selectedSymmAtoms.push_back( *sa + res_first_atom );
          AtomStatus[*sa] = SYMM;
        }
        // Add this group of symmetric atoms.
        SymmetricAtomIndices.push_back( selectedSymmAtoms );
      }
    }
  }
  if (debug_ > 0) {
    mprintf("DEBUG:\tResidue Atom Status:\n");
    for (int at = 0; at < Natom(); at++) {
      mprintf("\t%s", topIn.AtomMaskName(at + res_first_atom).c_str());
      switch (AtomStatus[at]) {
        case NONSYMM: mprintf(" Non-symmetric\n"); break;
        case SYMM   : mprintf(" Symmetric\n"); break;
        case UNSELECTED: mprintf(" Unselected\n");
      }
    }
  }
  return 0;
}
