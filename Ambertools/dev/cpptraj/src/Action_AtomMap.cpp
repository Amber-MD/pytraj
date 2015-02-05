// AtomMap
#include <list>
#include "Action_AtomMap.h"
#include "CpptrajStdio.h"
#include "TorsionRoutines.h"

// CONSTRUCTOR
Action_AtomMap::Action_AtomMap() :
  debug_(0),
  RefFrame_(0),
  TgtFrame_(0),
  maponly_(false),
  newFrame_(0),
  newParm_(0),
  stripParm_(0),
  rmsfit_(false),
  rmsdata_(0)
{}

void Action_AtomMap::Help() {
  mprintf("\t<target> <reference> [mapout <filename>] [maponly]\n"
          "\t[rmsfit [ rmsout <rmsout> ]]\n"
          "  Attempt to create a map from atoms in <target> to atoms in <reference> solely\n"
          "  based on how they are bonded; remap <target> so it matches <reference>.\n");
}

// DESTRUCTOR
Action_AtomMap::~Action_AtomMap() {
  if (newFrame_!=0) delete newFrame_;
  if (newParm_!=0) delete newParm_;
  if (stripParm_!=0) delete stripParm_;
}

// Action_AtomMap::mapBondsToUnique()
/** For each atom R in reference already mapped to unique atom T in 
  * target, try to match up non-mapped reference atoms r bonded to R to 
  * non-mapped target atoms t bonded to T. Repeat until no more reference 
  * atoms can be mapped in this way.
  * Checking is very strict in this routine; r and t must be the only 
  * possible match and the atomIDs must match.
  * \return the total number of atoms mapped.
  */
int Action_AtomMap::mapBondsToUnique(AtomMap& Ref, AtomMap& Tgt) {
  int numMappedAtoms=0;
  bool newSingle=true;

  while (newSingle) {

    // newSingle will be set back to true if any atoms are mapped
    newSingle=false;

    for (int ratom=0; ratom < Ref.Natom(); ratom++) {
      // Skip non-mapped atoms in Ref
      if (!Ref[ratom].IsMapped()) continue;
      int tatom = AMap_[ratom];
      // Check if map value is valid
      if (tatom<0) {
        mprintf("      Error: mapBondsToUnique: Ref %i:%s map value is invalid.\n",
                ratom+1, Ref[ratom].c_str());
        return -1;
      }
      // For each non-mapped atom bonded to Reference atom, try to 
      // find a matching non-mapped atom in Target by virtue of it being the
      // only possible match.
      for (Atom::bond_iterator r = Ref[ratom].bondbegin(); r != Ref[ratom].bondend(); r++)
      {
        // Check that bonded atom r is not already mapped
        if (Ref[*r].IsMapped()) continue;
        //mprintf("        Ref: Checking non-mapped %i:%s bonded to %i:%s\n",r,Ref->names[r],
        //        atom,Ref->names[atom]);
        // Check that non-mapped bonded ref atom r atomID is not the same as any 
        // other non-mapped bonded atomID.
        if ( Ref.BondIsRepeated(ratom, *r) ) continue;
        // At this point r is the only one of its kind bonded to atom.
        // Check if there is an analogous atom bonded to unique Target atom
        // tatom.
        for (Atom::bond_iterator t = Tgt[tatom].bondbegin(); t != Tgt[tatom].bondend(); t++)
        {
          // Check that bonded atom t is not already mapped 
          if (Tgt[*t].IsMapped()) continue;
          //mprintf("          Tgt: Checking non-mapped %i:%s bonded to %i:%s\n",t,Tgt->names[t],
          //        tatom,Tgt->names[tatom]);
          // Check that non-mapped bonded tgt atom t atomID is not the same as 
          // any other non-mapped bonded atomID.
          if ( Tgt.BondIsRepeated(tatom, *t) ) continue;
          // At this point t is the only one of its kind bonded to tatom.
          // Check if its atomID matches r. If so, map it.
          if ( Ref[*r].AtomID() == Tgt[*t].AtomID() ) {
            if (debug_>0) 
              mprintf("    Mapping tgt %i:%s to ref %i:%s based on single bond to unique.\n",
                      *t+1, Tgt[*t].c_str(), *r+1, Ref[*r].c_str());
            AMap_[*r] = *t;
            Ref[*r].SetMapped();
            Tgt[*t].SetMapped();
            newSingle = true;
            ++numMappedAtoms;
          }
        } // End loop over atoms bonded to tatom
      } // End loop over atoms bonded to ratom
      // Check if atom is completely mapped now
      Ref.MarkAtomComplete(ratom,false);
      Tgt.MarkAtomComplete(tatom,false);
    } // End loop over ref atoms

  } // End loop over newSingle
  return numMappedAtoms;
}        

// Action_AtomMap::mapChiral()
/** Given two atommaps and a map relating the two, find chiral centers for
  * which at least 3 of the atoms have been mapped. Assign the remaining
  * two atoms based on improper dihedrals. 
  * \return the total number of mapped atoms.
  * NOTE: ONLY WORKS FOR SP3
  */
int Action_AtomMap::mapChiral(AtomMap& Ref, AtomMap& Tgt) {
  int uR[5], uT[5], nR[4], nT[4];
  double dR[4], dT[4];
  int numMappedAtoms=0;

  for (int ratom=0; ratom < Ref.Natom(); ratom++) {
    // Skip non-mapped atoms
    if (!Ref[ratom].IsMapped()) continue;
    //mprintf("DBG: mapChiral: Ref atom %i:%s\n",atom,Ref->P->names[atom]);
    int tatom = AMap_[ratom];
    // Check that map value is valid
    if (tatom<0) {
      mprintf("      Error: mapChiral: Ref atom %i:%s map value is invalid.\n",
              ratom+1, Ref[ratom].c_str());
      return -1;
    }
    // If this Ref atom already completely mapped, skip
    if (Ref[ratom].Complete()) {
      // Sanity check - if Ref atom is completely mapped, target should be
      // unless # atoms in Tgt and Ref are different.
      if (!Tgt[tatom].Complete()) {
        mprintf("Warning: mapChiral: Ref atom %i:%s is complete but Tgt atom %i:%s is not.\n",
                ratom+1, Ref[ratom].c_str(), tatom+1, Tgt[tatom].c_str());
        //return 1;
      }
      continue;
    }
    // Check if this is a chiral center
    if (!Ref[ratom].IsChiral()) continue;
    // If target atom is not a chiral center (e.g. due to diff # atoms)
    // mapping by chirality is not important for this reference, let
    // mapByIndex handle it.
    if (!Tgt[tatom].IsChiral()) {
      mprintf("Warning: mapChiral: Ref atom %i:%s is chiral but Tgt atom %i:%s is not!\n",
              ratom+1, Ref[ratom].c_str(), tatom+1, Tgt[tatom].c_str());
      mprintf("         Marking Ref atom as non-chiral to try and map Tgt.\n");
      Ref[ratom].SetNotChiral();
      continue;
    }
    // Both atoms are chiral centers. Place bonded atoms (starting with 
    // central atom) in R and T.
    uR[0] = ratom;
    uT[0] = tatom;
    int nunique = 1;
    int notunique_r = 0;
    // Look for mapped bonded ref and target atoms, and nonmapped reference atoms
    for (Atom::bond_iterator r = Ref[ratom].bondbegin(); r != Ref[ratom].bondend(); r++)
    {
      int t = AMap_[*r];
      if (!Ref[*r].IsMapped()) {
        // Bonded atom r is not mapped 
        nR[notunique_r++] = *r;
      } else {
        // Bonded atom r is mapped. If a target was mapped to it
        // (i.e. it is the same atom) store it.
        if (t>=0) {
          if (Ref[*r].IsMapped() && Tgt[t].IsMapped()) {
            uR[nunique] = *r;
            uT[nunique] = t;
            ++nunique;
          }
        } 
      }
    }
    // Fill nT with nonmapped atoms from target
    int notunique_t = 0;
    for (Atom::bond_iterator tt = Tgt[tatom].bondbegin(); tt != Tgt[tatom].bondend(); tt++)
    {
      if (!Tgt[*tt].IsMapped()) 
        nT[notunique_t++] = *tt;
    }
    // notunique_r may not be the same as notunique_t if the # atoms is different
    if (notunique_r != notunique_t) 
      mprintf("Warning: Ref and Tgt do not have the same # of nonmapped atoms.\n");
    if (debug_>0) { 
      mprintf("  Potential Chiral center Ref=%i:%s Tgt=%i:%s  Mapped atoms=%i, non-Mapped=%i/%i\n",
              ratom+1, Ref[ratom].c_str(), tatom+1, Tgt[tatom].c_str(),
              nunique, notunique_r, notunique_t);
      for (int i=0; i<nunique; i++)
        mprintf("\t   Mapped\t%4i:%s %4i:%s\n", 
                uR[i]+1, Ref[uR[i]].c_str(), uT[i]+1, Tgt[uT[i]].c_str());
      for (int i=0; i<notunique_r; i++)
        mprintf("\tNotMappedRef\t%4i:%s\n", nR[i]+1, Ref[nR[i]].c_str());
      for (int i=0; i<notunique_t; i++)
        mprintf("\tNotMappedTgt\t         %4i:%4s\n", nT[i]+1, Tgt[nT[i]].c_str());
    }
    // If all atoms are unique no need to map
    // NOTE: Should be handled by complete check above.
    //if (nunique==5) continue;
    // Require at least 3 unique atoms for dihedral calc. 
    if (nunique<3) {
      if (debug_>0) 
        mprintf("    Warning: Center has < 3 mapped atoms, dihedral cannot be calcd.\n");
      continue;
    }
    // Calculate reference improper dihedrals
    for (int i=0; i<notunique_r; i++) {
      dR[i] = Torsion( RefFrame_->RefFrame().XYZ(uR[0]),
                       RefFrame_->RefFrame().XYZ(uR[1]), 
                       RefFrame_->RefFrame().XYZ(uR[2]),
                       RefFrame_->RefFrame().XYZ(nR[i]) );
      if (debug_>1) mprintf("    Ref Improper %i [%3i,%3i,%3i,%3i]= %lf\n",i,
                           uR[0]+1, uR[1]+1, uR[2]+1, nR[i]+1, dR[i]+1);
    }
    // Calculate target improper dihedrals
    for (int i=0; i<notunique_t; i++) {
      dT[i] = Torsion( TgtFrame_->RefFrame().XYZ(uT[0]),
                       TgtFrame_->RefFrame().XYZ(uT[1]),
                       TgtFrame_->RefFrame().XYZ(uT[2]),
                       TgtFrame_->RefFrame().XYZ(nT[i]) );
      if (debug_>1) mprintf("    Tgt Improper %i [%3i,%3i,%3i,%3i]= %lf\n",i,
                           uR[0]+1, uR[1]+1, uR[2]+1, nT[i]+1, dT[i]+1);
    }
    // Match impropers to each other using a cutoff. Note that all torsions
    // are in radians.
    // NOTE: 10.0 degrees seems reasonable? Also there is currently no 
    //       check for repeated deltas.
    for (int i=0; i<notunique_r; i++) {
      for (int j=0; j<notunique_t; j++) {
        double delta = dR[i] - dT[j];
        if (delta<0.0) delta=-delta;
        if (delta<0.17453292519943295769236907684886) {
          if (debug_>0)
            mprintf("    Mapping tgt atom %i:%s to ref atom %i:%s based on chirality.\n",
                    nT[j]+1, Tgt[nT[j]].c_str(), nR[i]+1, Ref[nR[i]].c_str() );
          AMap_[ nR[i] ] = nT[j];
          ++numMappedAtoms;
          // Once an atom has been mapped set its mapped flag
          Ref[nR[i]].SetMapped();
          Tgt[nT[j]].SetMapped();
        } else if (notunique_r == 1 && notunique_t == 1) {
          // This is the only non-mapped atom of the chiral center but for
          // some reason the improper dihedral doesnt match. Map it but warn
          // the user.
          mprintf("Warning: Ref %i:%s and Tgt %i:%s are the only unmapped atoms of chiral\n"
                  "Warning: centers %i:%s | %i:%s, but the improper dihedral angles do not\n"
                  "Warning: match (%.4f rad != %.4f rad). This can indicate structural problems\n"
                  "Warning: in either the target or reference. Mapping atoms, but it is\n"
                  "Warning: recommended the structures be visually inspected for problems.\n",
                  nR[i]+1, Ref[nR[i]].c_str(), nT[j]+1, Tgt[nT[j]].c_str(),
                  ratom+1, Ref[ratom].c_str(), tatom+1, Tgt[tatom].c_str(),
                  dR[i], dT[j]);
          AMap_[ nR[i] ] = nT[j];
          ++numMappedAtoms;
          Ref[nR[i]].SetMapped();
          Tgt[nT[j]].SetMapped();
        }
      }
    }
    // Check if ref atom or tgt atom is now completely mapped
    Ref.MarkAtomComplete(ratom,false);
    Tgt.MarkAtomComplete(tatom,false);
  } // End loop over ratom

  return numMappedAtoms;
}

// Action_AtomMap::mapUniqueRefToTgt()
/** If the number of atoms in Ref is different from Tgt, it is possible that
  * Tgt is missing atoms (or maybe vice versa). If the difference is not too
  * great it may be possible to look for an unmapped atom in Ref that has
  * same name and at least 1 matching bond (# bonds may be diff due to # atoms
  * so atomID cannot be used).
  * If only one name matches, probably safe to map it. 
  * \return 1 if the atom could be mapped, 0 otherwise.
  */
int Action_AtomMap::mapUniqueRefToTgt(AtomMap& Ref, AtomMap& Tgt, int ratom) {
  int match=-1;
  for (int tatom=0; tatom < Tgt.Natom(); tatom++) {
    //mprintf("DBG:        %i:%s %i:%s\n",t,Tgt->names[t],atom,Ref->names[atom]);
    // If atom #s are different Tgt atom could be unique but not mapped. Check
    // if tgt has already been mapped using Amap. 
    bool alreadyMapped=false;
    for (int at = 0; at < Ref.Natom(); at++) {
      if (AMap_[at] == tatom) {
        alreadyMapped = true;
        break;
      }
    }
    if (alreadyMapped) continue;
    // Check char name
    // NOTE: Should this be checking element instead?
    if ( Tgt[tatom].CharName() == Ref[ratom].CharName() ) {
      if (debug_>1) 
        mprintf("        Attempting match of Tgt %i:%s to Ref %i:%s\n",
                tatom+1,Tgt[tatom].c_str(), ratom+1, Ref[ratom].c_str());
      // Check that at least 1 bond is in common
      int commonBond=0;
      for (Atom::bond_iterator r = Ref[ratom].bondbegin(); r != Ref[ratom].bondend(); r++)
      {
        // Check Map for ref bonded atom
        int t = AMap_[ *r ];
        // If no mapping exists cant check it
        if (t<0) continue;
        if (debug_>1) 
          mprintf("          Ref %i:%s bonded to %i:%s (%i:%s in tgt)\n",
                  ratom+1, Ref[ratom].c_str(), 
                  *r+1,    Ref[*r].c_str(), 
                  t+1,     Tgt[t].c_str());
        // Loop over all bonds in current target atom. See if it matches index 
        // of a mapped ref bonded atom.
        for (Atom::bond_iterator tbond = Tgt[tatom].bondbegin(); 
                                 tbond != Tgt[tatom].bondend(); tbond++) 
        {
          if (debug_>1)
            mprintf("            Tgt %i:%s bonded to %i:%s\n",
                    tatom+1, Tgt[tatom].c_str(), *tbond+1, Tgt[*tbond].c_str());
          if (t == *tbond) 
            ++commonBond;
        }
      }
      if (commonBond==0) continue;
      // This Tgt Name matches and at least 1 bond in common with Ref atom
      // Check that a match has not yet been found for ref
      if (match!=-1) {
        mprintf("      Warning: mapUniqueRefToTgt: Ref %i:%s has multiple potential matches\n",
                ratom+1, Ref[ratom].c_str());
        mprintf("               among Tgt [%i:%s, %i:%s]\n",
                tatom+1, Tgt[tatom].c_str(), match+1, Tgt[match].c_str());
        return 0;
      }
      match = tatom;
    }
  }
  if (match==-1) return 0;
  // Only one match found - map it
  if (debug_>0) 
    mprintf("    Mapping target %i:%s to unique ref %i:%s\n",match+1,Tgt[match].c_str(),
            ratom+1, Ref[ratom].c_str());
  AMap_[ratom] = match;
  Ref[ratom].SetMapped();
  Tgt[match].SetMapped();
  return 1;
}

// Action_AtomMap::mapByIndex()
/** Given two atommaps and a map relating the two, attempt to map any remaining
  * incomplete atoms by assuming the atom indices in reference and target are
  * in similar orders. At this point all unique atoms should have been mapped.
  * First, for each reference atom R check if R is unique but not mapped and 
  * attempt to match it to a non-mapped target based on local bonding 
  * environment (mapUniqueRefToTgt). Lastly, for reference atom R mapped to 
  * target atom T, compare the non-mapped atoms bonded to R (r) to the 
  * non-mapped atoms bonded to T (t). If the unique IDs of r and t match, map 
  * them. Otherwise if there is only one potential match between r and t map 
  * them.
  * \return the number of atoms mapped this way. 
  */
int Action_AtomMap::mapByIndex(AtomMap& Ref, AtomMap& Tgt) {
  int match;
  int numAtomsMapped=0;

  for (int ratom=0; ratom<Ref.Natom(); ratom++) {
    int tatom = AMap_[ratom];
    // Check if no mapping exists for this atom 
    if (tatom<0) {
      // Check if reference atom is unique, but hasnt had a target mapped to it.
      // This can arise when the number of atoms in ref and tgt not equal.
      // If the difference in atoms is not too great (probably ~1) attempt
      // to look for a similar atom in tgt (name and index) and map it.
      if (Ref[ratom].IsUnique()) {
        mprintf("      Warning: mapByIndex: Atom %i:%s in reference is unique but not mapped!\n",
                ratom+1, Ref[ratom].c_str());
        if (mapUniqueRefToTgt(Ref,Tgt,ratom)) 
          ++numAtomsMapped;
      }
      continue;
    }
    // Skip over non-mapped atoms
    //if (!Ref->M[atom].isMapped) continue;

    // Check that num bonds match in Ref and target.
    // The # of bonds might not be equal if the # atoms in ref and tgt
    // not equal.
    if (Ref[ratom].Nbonds() != Tgt[tatom].Nbonds()) {
      mprintf(
        "\tWarning: mapByIndex: Ref atom %i:%s #bonds (%i) does not match Tgt atom %i:%s (%i)\n",
        ratom+1, Ref[ratom].c_str(), Ref[ratom].Nbonds(),
        tatom+1, Tgt[tatom].c_str(), Tgt[tatom].Nbonds()
      );
      //return 1;
    }
    // Skip completely mapped atoms - check that both Ref and Tgt are complete
    // NOTE: This is ok if #atoms in ref > #atoms in tgt but not the other way around.
    if (Ref[ratom].Complete()) {
      if (!Tgt[tatom].Complete()) {
        mprintf("Error: mapByIndex: Ref atom %i:%s is complete but Tgt atom %i:%s is not.\n",
                ratom+1,Ref[ratom].c_str(),tatom+1,Tgt[tatom].c_str() );
        //return 1;
      }
      continue;
    }
    // This atom is mapped, but bonded atoms are not completely mapped. Try
    // to map the unmapped reference atoms bonded to <atom> to the unmapped
    // target atoms bonded to <tatom>. 
    if (debug_>1)
      mprintf("DBG: Checking bonds of mapped Ref %i:%s (isChiral=%i)"
              " against mapped Tgt %i:%s (isChiral=%i)\n",
              ratom+1, Ref[ratom].c_str(), (int)Ref[ratom].IsChiral(),
              tatom+1, Tgt[tatom].c_str(), (int)Tgt[tatom].IsChiral());
    for (Atom::bond_iterator r = Ref[ratom].bondbegin(); r != Ref[ratom].bondend(); r++)
    {
      if (debug_>1) 
        mprintf("\t\tRefBond %i:%s [%1i]\n",*r+1,Ref[*r].c_str(),(int)Ref[*r].IsMapped());
      if (Ref[*r].IsMapped()) continue;
      // Dont map atoms that are single-bonded to chiral centers; let
      // mapChiral take care of them.
      if (Ref[ratom].IsChiral() && Ref[*r].Nbonds()==1) continue;
      match = -1;
      for (Atom::bond_iterator t = Tgt[tatom].bondbegin(); t != Tgt[tatom].bondend(); t++)
      {
        if (debug_>1) 
          mprintf("\t\t\tTgtBond %i:%s [%1i]\n",*t+1,Tgt[*t].c_str(),(int)Tgt[*t].IsMapped());
        if (Tgt[*t].IsMapped()) continue;
        // Atom r bonded to atom, atom t bonded to tatom. r and t are not
        // yet mapped. Check if names match
        // NOTE: Check elements instead?
        if ( Ref[*r].CharName() != Tgt[*t].CharName() ) continue;
        // If the uniqueIDs of bonded atom r and bonded atom t match, map them now
        // NOTE: Scan for repeats?
        if ( Ref[*r].Unique() == Tgt[*t].Unique() ) {
          match = *t;
          break;
        }
        // Store this atom t bonded to tatom as a potential match. If another
        // match has already been stored we cant tell these apart yet so ignore.
        if (match==-1) {
          match = *t;
        } else {
          mprintf("\tWarning: mapByIndex: Atom %i:%s bonded to Ref %i:%s has too many matches.\n",
                  *r+1,Ref[*r].c_str(), ratom+1, Ref[ratom].c_str());
          match = -1;
          break;
        }
      } // End loop tbond over bonds in target atom
      // If a match was found, Map it
      if (match!=-1) {
        if (debug_>0) mprintf("    Mapping Tgt %i:%s to Ref %i:%s based on name/bonding.\n",
                             match+1,Tgt[match].c_str(),*r+1,Ref[*r].c_str());
        AMap_[*r] = match;
        Ref[*r].SetMapped();
        Tgt[match].SetMapped();
        ++numAtomsMapped;
      }
    } // End loop over atoms bonded to Ref atom
    // Check if atom is completely mapped now
    Ref.MarkAtomComplete(ratom,false);
    Tgt.MarkAtomComplete(tatom,false);
  } // End loop over ratoms

  return numAtomsMapped;
}

// Action_AtomMap::MapUniqueAtoms()
/** Map unique atoms in reference to unique atoms in target. If no atoms
  * can be mapped in this way, attempt to guess a starting point based
  * first on uniqueID, then by chirality.
  * \return number of atoms mapped.
  */
int Action_AtomMap::MapUniqueAtoms(AtomMap& Ref, AtomMap& Tgt) {
  int numAtomsMapped=0;

  // Atoms have now been assigned IDs. Match up the unique strings in Ref with 
  // unique strings in target.
  for (int refatom=0; refatom < Ref.Natom(); refatom++) {
    AMap_[refatom] = -1;
    // If the ID of this reference atom is unique, look for same ID in target
    if (Ref[refatom].IsUnique()) {
      for (int targetatom=0; targetatom < Tgt.Natom(); targetatom++) {
        // If ID of thie target atom is unique, check if it matches reference atom ID
        if (Tgt[targetatom].IsUnique()) {
          if ( Tgt[targetatom].Unique() == Ref[refatom].Unique() ) {
            // Check that number of bonds is consistent
            if (Ref[refatom].Nbonds() != Tgt[targetatom].Nbonds()) {
              mprintf("\tWarning: AtomMap: Atoms R%i and T%i have same ID but different # bonds!\n",
                      refatom,targetatom);
            }
            AMap_[refatom] = targetatom;
            Ref[refatom].SetMapped();
            Tgt[targetatom].SetMapped();
            ++numAtomsMapped;
            if (debug_>0) {
              mprintf("    Mapping Tgt %i:%s to Ref %i:%s based on unique ID\n",
                      targetatom+1, Tgt[targetatom].c_str(),
                      refatom+1, Ref[refatom].c_str());
            }
          } // If unique strings match
        } // If target atom is unique
      } // Loop over target atoms
    } // If reference atom is unique
  } // Loop over reference atoms

  return numAtomsMapped;
}

// Action_AtomMap::MapWithNoUniqueAtoms()
/** If no unique atoms could be mapped it means the molecule is probably
  * very symmetric, so try to guess a good starting point. Map the first 
  * atoms that have a uniqueID duplicated only 1 time, preferably a chiral 
  * center. Try each pair of atoms and compute the resulting RMSD; use
  * the map with the lowest overall RMSD.
  * Note that the current implementation isn't very smart since it will 
  * try guess pairings that may have already been mapped in a previous
  * try.
  */
// NOTE: Also store the number of atoms mapped?
int Action_AtomMap::MapWithNoUniqueAtoms( AtomMap& Ref, AtomMap& Tgt ) {
  std::list<int> refGuess;
  std::list<int> tgtGuess;
  double lowestRMS = 0;
  std::vector<int> bestMap;
  int numAtomsMapped;

  mprintf("      Warning: No unique atoms found, usually indicates highly symmetric system.\n");
  mprintf("               Trying to guess starting point.\n");
  //mprintf("DEBUG: Ref has %i atoms, Tgt has %i\n",Ref->natom, Tgt->natom);
  // Get a list of atoms in ref duplicated only once, preferably chiral
  for (int refatom=0; refatom < Ref.Natom(); refatom++) {
    if (Ref[refatom].Nduplicated()==1) {
      if (Ref[refatom].IsChiral()) 
        refGuess.push_front(refatom);
      else 
        refGuess.push_back(refatom);
    }
  }
  if (refGuess.empty()) {
    mprintf("Error: AtomMap: Could not find starting point in reference.\n");
    return 1;
  }
  mprintf("Ref guess atoms:");
  for (std::list<int>::iterator r=refGuess.begin(); r!=refGuess.end(); r++)
    mprintf(" %i",(*r)+1);
  mprintf("\n");
  // Get a list of atoms in tgt duplicated only once, preferably chiral
  for (int targetatom=0; targetatom < Tgt.Natom(); targetatom++) {
    if (Tgt[targetatom].Nduplicated()==1) {
      if (Tgt[targetatom].IsChiral())
        tgtGuess.push_front(targetatom);
      else
        tgtGuess.push_back(targetatom);
    }
  }
  if (tgtGuess.empty()) {
    mprintf("Error: AtomMap: Could not find starting point in target.\n");
    return 1;
  }
  mprintf("Tgt guess atoms:");
  for (std::list<int>::iterator t=tgtGuess.begin(); t!=tgtGuess.end(); t++)
    mprintf(" %i",(*t)+1);
  mprintf("\n");
  // For each pair of atoms in refGuess and tgtGuess that have the same
  // ID string, guess that they are mapped and attempt to perform atom
  // mapping from there.
  for (std::list<int>::iterator r=refGuess.begin(); r!=refGuess.end(); r++) {
    for (std::list<int>::iterator t=tgtGuess.begin(); t!=tgtGuess.end(); t++) {
      if ( Ref[*r].Unique() == Tgt[*t].Unique() ) {
        // Reset any previous mapping
        for (int mapi=0; mapi < Ref.Natom(); mapi++) AMap_[mapi]=-1;
        Ref.ResetMapping();
        Tgt.ResetMapping();
        //mprintf("  Ref %i (%i) to Tgt %i (%i) MATCH!\n",*r,Ref->natom,*t,Tgt->natom); // DEBUG
        // Map this guess
        AMap_[(*r)] = (*t);
        Ref[(*r)].SetMapped();
        Tgt[(*t)].SetMapped();
        mprintf("    Mapping Tgt %i:%s to Ref %i:%s based on guess.\n",
                (*t)+1, Tgt[*t].c_str(),
                (*r)+1, Ref[*r].c_str());
        // Attempt to complete mapping based on the guess
        if ( MapAtoms(Ref,Tgt) ) return 1;
        // Count number of mapped atoms
        numAtomsMapped=0;
        for (int mapi=0; mapi < Ref.Natom(); mapi++) 
          if (AMap_[mapi]!=-1) ++numAtomsMapped;
        // If < 3 atoms mapped this will cause a problem with RMSD
        if (numAtomsMapped<3) continue;
        // Score this mapping with an RMSD ---------------------------------
        // Set up a reference/target frame containing only mapped atoms
        rmsRefFrame_.StripUnmappedAtoms(RefFrame_->RefFrame(), AMap_);
        rmsTgtFrame_.ModifyByMap(TgtFrame_->RefFrame(), AMap_);
        double RmsVal = rmsTgtFrame_.RMSD(rmsRefFrame_, false);
        mprintf("\tRMS fit (%i atoms) based on guess Tgt %i -> Ref %i, %lf\n",
                numAtomsMapped,(*t)+1, (*r)+1, RmsVal);
        // -----------------------------------------------------------------
        // If the current RmsVal is lower than the lowestRMS, store this map.
        if (bestMap.empty() || RmsVal < lowestRMS) {
          bestMap = AMap_;
          lowestRMS = RmsVal;
        }
      }
    } // End loop over tgt guesses
  } // End loop over ref guesses

  // If bestMap is null something went wrong. Otherwise set AMap to best map.
  if (bestMap.empty()) {
    mprinterr("Error: Could not guess starting point for atom mapping.\n");
    return 1;
  } else {
    AMap_ = bestMap;
  }
  return 0;
}

// Action_AtomMap::MapAtoms()
/** Map atoms in tgt to atoms in reference. Assumes that any uniquely 
  * identified atoms have already been mapped. First map unmapped atoms 
  * that are the only one of their kind bonded to a unique or already 
  * mapped atom (mapBondsToUnique). Then map atoms based on chirality; 
  * if any atoms are mapped in this way check to see if mapBondsToUnique 
  * finds new atoms. Last try to guess mapping based on bonds (mapByIndex), 
  * which will also attempt to map atoms in Ref that are unique but not 
  * mapped to atoms in Tgt (which can happen e.g. if Tgt is missing atoms).
  * Negative return values from mapXXX routines indicates error.
  * \return 0 on success, 1 on error.
  */
int Action_AtomMap::MapAtoms(AtomMap& Ref, AtomMap& Tgt) {
  bool mapatoms=true;
  int numAtomsMapped;
  int iterations=0;

  // DEBUG
  //char name[1024];
  //sprintf(name,"Ref.%i.mol2",iterations);
  //Ref->WriteMol2(name);
  //sprintf(name,"Tgt.%i.mol2",iterations);
  //Tgt->WriteMol2(name);
  // END DEBUG
  // Search for completely mapped atoms. If an atom and all atoms
  // it is bonded to are unique, mark the atom as completely mapped.
  RefMap_.CheckForCompleteAtoms();
  TgtMap_.CheckForCompleteAtoms();

  // Map remaining non-unique atoms
  while (mapatoms) {
    ++iterations;
    // First assign based on bonds to unique (already mapped) atoms.
    numAtomsMapped=mapBondsToUnique(Ref,Tgt);
    // DEBUG
    //sprintf(name,"Ref.%i.u.mol2",iterations);
    //Ref->WriteMol2(name);
    //sprintf(name,"Tgt.%i.u.mol2",iterations);
    //Tgt->WriteMol2(name);
    // END DEBUG
    if (debug_>0)
      mprintf("* [%3i] mapBondsToUnique: %i atoms mapped.\n",iterations,numAtomsMapped);
    if (numAtomsMapped<0) return 1;
    // Next assign based on chirality
    numAtomsMapped=mapChiral(Ref,Tgt);
    // DEBUG
    //sprintf(name,"Ref.%i.c.mol2",iterations);
    //Ref->WriteMol2(name);
    //sprintf(name,"Tgt.%i.c.mol2",iterations);
    //Tgt->WriteMol2(name);
    // END DEBUG
    if (debug_>0)
      mprintf("* [%3i]        mapChiral: %i atoms mapped.\n",iterations,numAtomsMapped);
    if (numAtomsMapped<0) return 1;
    if (numAtomsMapped>0) continue;
    // Last assign based on index/element
    numAtomsMapped=mapByIndex(Ref,Tgt);
    // DEBUG
    //sprintf(name,"Ref.%i.i.mol2",iterations);
    //Ref->WriteMol2(name);
    //sprintf(name,"Tgt.%i.i.mol2",iterations);
    //Tgt->WriteMol2(name);
    // END DEBUG
    if (debug_>0)
      mprintf("* [%3i]       mapByIndex: %i atoms mapped.\n",iterations,numAtomsMapped);
    if (numAtomsMapped<0) return 1;
    if (numAtomsMapped==0) mapatoms=false;
  }
  if (debug_>0) mprintf("* %i iterations.\n",iterations);
  return 0;
}

// Action_AtomMap::Init()
Action::RetType Action_AtomMap::Init(ArgList& actionArgs, TopologyList* PFL, DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  std::string rmsout;
  CpptrajFile outputfile;
  int refatom,targetatom;
  debug_ = debugIn; 
  RefMap_.SetDebug(debug_);
  TgtMap_.SetDebug(debug_);

  // Get Args
  std::string outputname = actionArgs.GetStringKey("mapout");
  maponly_ = actionArgs.hasKey("maponly");
  rmsfit_ = actionArgs.hasKey("rmsfit");
  if (rmsfit_)
    rmsout = actionArgs.GetStringKey("rmsout");
  std::string targetName = actionArgs.GetStringNext();
  std::string refName = actionArgs.GetStringNext();
  if (targetName.empty()) {
    mprinterr("Error: No target specified.\n");
    return Action::ERR;
  }
  if (refName.empty()) {
    mprinterr("Error: No reference specified.\n");
    return Action::ERR;
  }
  // Get Reference
  RefFrame_ = (DataSet_Coords_REF*)DSL->GetReferenceFrame( refName );
  if (RefFrame_ == 0) {
    mprinterr("Error: Could not get reference frame %s\n",refName.c_str());
    return Action::ERR;
  }
  // Get Target
  TgtFrame_ = (DataSet_Coords_REF*)DSL->GetReferenceFrame( targetName );
  if (TgtFrame_ == 0) {
    mprinterr("Error: Could not get target frame %s\n",targetName.c_str());
    return Action::ERR;
  }
 
  mprintf("    ATOMMAP: Atoms in trajectories associated with parm %s will be\n",
          TgtFrame_->Top().c_str());
  mprintf("             mapped according to parm %s.\n",RefFrame_->Top().c_str());
  if (!outputname.empty())
    mprintf("             Map will be written to %s\n",outputname.c_str());
  if (maponly_)
    mprintf("             maponly: Map will only be written, not used in trajectory read.\n");
  if (!maponly_ && rmsfit_) {
    mprintf("             rmsfit: Will rms fit mapped atoms in tgt to reference.\n");
    if (!rmsout.empty()) {
      rmsdata_ = DSL->AddSet(DataSet::DOUBLE, actionArgs.GetStringNext(), "RMSD");
      if (rmsdata_==0) return Action::ERR;
      DFL->AddSetToFile(rmsout, rmsdata_);
    }
  }

  // For each map, set up (get element for each atom, initialize map mem),
  // determine what atoms are bonded to each other via simple distance
  // cutoffs, the give each atom an ID based on what atoms are bonded to
  // it, noting which IDs are unique for that map. 

  if (RefMap_.Setup(RefFrame_->Top())!=0) return Action::ERR;
  //RefMap_.WriteMol2((char*)"RefMap.mol2\0"); // DEBUG
  RefMap_.DetermineAtomIDs();

  if (TgtMap_.Setup(TgtFrame_->Top())!=0) return Action::ERR;
  //TgtMap_.WriteMol2((char*)"TgtMap.mol2\0"); // DEBUG
  TgtMap_.DetermineAtomIDs();

  // Check if number of atoms in each map is equal
  if (RefMap_.Natom() != TgtMap_.Natom()) {
    mprintf("Warning: # atoms in reference (%i) not equal\n",
            RefMap_.Natom());
    mprintf("Warning:\tto # atoms in target (%i).\n",TgtMap_.Natom());
  }

  // Set up RMS frames to be able to hold max # of possible atoms 
  rmsRefFrame_.SetupFrame(RefMap_.Natom());
  rmsTgtFrame_.SetupFrame(RefMap_.Natom());

  // Allocate memory for atom map
  //   AMap_[reference]=target
  AMap_.resize( RefMap_.Natom(), -1); 
  // Map unique atoms
  int numMappedAtoms = MapUniqueAtoms(RefMap_, TgtMap_);
  if (debug_>0)
    mprintf("*         MapUniqueAtoms: %i atoms mapped.\n",numMappedAtoms);
  // If no unique atoms mapped system is highly symmetric and needs to be
  // iteratively mapped. Otherwise just map remaining atoms.
  if (numMappedAtoms==0) { 
    if (MapWithNoUniqueAtoms(RefMap_,TgtMap_)) return Action::ERR;
  } else {
    if (MapAtoms(RefMap_,TgtMap_)) return Action::ERR;
  }

  // Print atom map and count # mapped atoms
  numMappedAtoms = 0;
  outputfile.OpenEnsembleWrite(outputname, DSL->EnsembleNum());
  outputfile.Printf("%-6s %4s %6s %4s\n","#TgtAt","Tgt","RefAt","Ref");
  for (refatom=0; refatom < RefMap_.Natom(); refatom++) {
    targetatom = AMap_[refatom];
    if (targetatom < 0) 
      outputfile.Printf("%6s %4s %6i %4s\n","---","---",refatom+1,RefMap_[refatom].c_str());
    else
      outputfile.Printf("%6i %4s %6i %4s\n",targetatom+1,TgtMap_[targetatom].c_str(),
                            refatom+1, RefMap_[refatom].c_str());
    if (targetatom>=0) {
      //mprintf("* TargetAtom %6i(%4s) maps to RefAtom %6i(%4s)\n",
      //                targetatom,TgtMap_.P->names[targetatom],
      //                refatom,RefMap_.P->names[refatom]);
      ++numMappedAtoms;
    } //else {
    //  mprintf("* Could not map any TargetAtom to RefAtom %6i(%4s)\n",
    //                  refatom,RefMap_.P->names[refatom]);
    //}
  }
  outputfile.CloseFile();
  mprintf("      %i total atoms were mapped.\n",numMappedAtoms);
  if (maponly_) return Action::OK;

  // If rmsfit is specified, an rms fit of target to reference will be
  // performed using all atoms that were successfully mapped from 
  // target to reference.
  if (rmsfit_) {
    // Set up a reference frame containing only mapped reference atoms
    rmsRefFrame_.StripUnmappedAtoms(RefFrame_->RefFrame(), AMap_);
    mprintf("      rmsfit: Will rms fit %i atoms from target to reference.\n",numMappedAtoms);
    return Action::OK;
  }

  // Check if not all atoms could be mapped
  if (numMappedAtoms != RefMap_.Natom()) {
    // If the number of mapped atoms is less than the number of reference
    // atoms but equal to the number of target atoms, can modify the reference
    // frame to only include mapped atoms
    if (numMappedAtoms<RefMap_.Natom() && numMappedAtoms==TgtMap_.Natom()) {
      // Create mask that includes only reference atoms that could be mapped
      AtomMask M1;
      for (refatom = 0; refatom < RefMap_.Natom(); refatom++) {
        if (AMap_[refatom] != -1) M1.AddAtom(refatom);
      }
      // Strip reference parm
      mprintf("    Modifying reference '%s' topology and frame to match mapped atoms.\n",
              RefFrame_->FrameName().base());
      if (RefFrame_->StripRef( M1 )) return Action::ERR;
      // Since AMap[ ref ] = tgt but ref is now missing any stripped atoms,
      // the indices of AMap must be shifted to match
      int refIndex = 0; // The new index
      for (refatom = 0; refatom < RefMap_.Natom(); refatom++) {
        targetatom = AMap_[refatom];
        if (targetatom<0)
          continue;
        else
          AMap_[refIndex++]=targetatom;
      }
    } else {
      mprintf("Warning: AtomMap: Not all atoms were mapped. Frames will not be modified.\n");
      maponly_=true;
    }
  }

  if (!maponly_) {
    // Set up new Frame
    newFrame_ = new Frame();
    newFrame_->SetupFrameM( TgtFrame_->Top().Atoms() );

    // Set up new Parm
    newParm_ = TgtFrame_->Top().ModifyByMap(AMap_);
  }

  return Action::OK;
}

// Action_AtomMap::setup()
/** If the current parm does not match the target parm, deactivate. Otherwise
  * replace current parm with mapped parm.
  */
Action::RetType Action_AtomMap::Setup(Topology* currentParm, Topology** parmAddress) {
  if (maponly_) {
    mprintf("    ATOMMAP: maponly was specified, not using atom map during traj read.\n");
    return Action::OK;
  }
  if (currentParm->Pindex() != TgtFrame_->Top().Pindex() ||
      currentParm->Natom() != TgtFrame_->Top().Natom()) 
  {
    mprintf("    ATOMMAP: Map for parm %s -> %s (%i atom).\n",TgtFrame_->Top().c_str(),
            RefFrame_->Top().c_str(), TgtFrame_->Top().Natom());
    mprintf("             Current parm %s (%i atom).\n",currentParm->c_str(),
            currentParm->Natom());
    mprintf("             Not using map for this parm.\n");
    return Action::ERR;
  }
  if (rmsfit_) {
    mprintf("    ATOMMAP: rmsfit specified, %i atoms.\n",rmsRefFrame_.Natom());
    return Action::OK;
  }
  mprintf("    ATOMMAP: Map for parm %s -> %s (%i atom).\n",TgtFrame_->Top().c_str(),
          RefFrame_->Top().c_str(), TgtFrame_->Top().Natom());

  *parmAddress = newParm_;
  
  return Action::OK;
}

// Action_AtomMap::action()
/** Modify the current frame based on the atom map. 
  */
Action::RetType Action_AtomMap::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {
  if (maponly_) return Action::OK;

  // Perform RMS fit on mapped atoms only
  if (rmsfit_) {
    // Set target frame up according to atom map.
    rmsTgtFrame_.ModifyByMap(*currentFrame, AMap_);
    Matrix_3x3 Rot;
    Vec3 Trans, refTrans;
    double R = rmsTgtFrame_.RMSD(rmsRefFrame_, Rot, Trans, refTrans, false);
    currentFrame->Trans_Rot_Trans(Trans, Rot, refTrans);
    if (rmsdata_!=0)
      rmsdata_->Add(frameNum, &R);
    return Action::OK;
  }

  // Modify the current frame
  // TODO: Fix this since its probably busted for unmapped atoms
  newFrame_->SetCoordinatesByMap(*currentFrame, AMap_);
  *frameAddress = newFrame_;
  return Action::OK;
}

