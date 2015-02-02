#include "SymmetricRmsdCalc.h"
#include "DistRoutines.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
SymmetricRmsdCalc::SymmetricRmsdCalc() : debug_(0), fit_(true), useMass_(false) {}

// CONSTRUCTOR - For use when only RMSD is wanted.
SymmetricRmsdCalc::SymmetricRmsdCalc(AtomMask const& maskIn, bool fitIn, 
                                     bool useMassIn, Topology const& topIn, int debugIn) :
  debug_(debugIn), fit_(fitIn), useMass_(useMassIn)
{
  SetupSymmRMSD( topIn, maskIn, false ); // No remap warning
}

// SymmetricRmsdCalc::InitSymmRMSD()
int SymmetricRmsdCalc::InitSymmRMSD(bool fitIn, bool useMassIn, int debugIn)
{
  debug_ = debugIn;
  fit_ = fitIn;
  useMass_ = useMassIn;
  return 0;
}

// SymmetricRmsdCalc::SetupSymmRMSD()
/** Find potential symmetric atoms. All residues up to the last selected
  * residue are considered.
  */
int SymmetricRmsdCalc::SetupSymmRMSD(Topology const& topIn, AtomMask const& tgtMask, bool remapIn)
{
  // Allocate space for remapping selected atoms in target frame. This will
  // also put the correct masses in based on the mask.
  tgtRemap_.SetupFrameFromMask(tgtMask, topIn.Atoms());
  // Create map of original atom numbers to selected indices
  Iarray SelectedIdx( topIn.Natom(), -1 );
  int tgtIdx = 0;
  for (int originalAtom = 0; originalAtom != topIn.Natom(); ++originalAtom)
    if ( originalAtom == tgtMask[tgtIdx] )
      SelectedIdx[originalAtom] = tgtIdx++;
  if (debug_ > 0) {
    mprintf("DEBUG: Original atom -> Selected Index mapping:\n");
    for (int originalAtom = 0; originalAtom != topIn.Natom(); ++originalAtom)
      mprintf("\t%8i -> %8i\n", originalAtom + 1, SelectedIdx[originalAtom] + 1);
  }
  // Create initial 1 to 1 atom map for all selected atoms; indices in 
  // SymmetricAtomIndices will correspond to positions in AMap.
  AMap_.resize( tgtRemap_.Natom() );
  // Determine last selected residue.
  int last_res = topIn[tgtMask.back()].ResNum() + 1;
  mprintf("\tResidues up to %s will be considered for symmetry correction.\n",
          topIn.TruncResNameNum(last_res-1).c_str());
  // In each residue, determine which selected atoms are symmetric.
  SymmetricAtomIndices_.clear();
  AtomMap resmap;
  if (debug_ > 1) resmap.SetDebug(1);
  for (int res = 0; res < last_res; ++res) {
    AtomMap::AtomIndexArray residue_SymmetricGroups;
    if (resmap.SymmetricAtoms(topIn, residue_SymmetricGroups, res)) {
      mprinterr("Error: Finding symmetric atoms in residue '%s'\n",
                topIn.TruncResNameNum(res).c_str());
      return 1;
    }
    if (!residue_SymmetricGroups.empty()) {
      // Which atoms in symmetric groups are selected?
      bool resHasSelectedSymmAtoms = false;
      for (AtomMap::AtomIndexArray::const_iterator symmGroup = residue_SymmetricGroups.begin();
                                                   symmGroup != residue_SymmetricGroups.end();
                                                 ++symmGroup)
      {
        Iarray selectedAtomIndices;
        for (Iarray::const_iterator atnum = symmGroup->begin();
                                    atnum != symmGroup->end(); ++atnum)
        {
          if ( SelectedIdx[*atnum] != -1 )
            selectedAtomIndices.push_back( SelectedIdx[*atnum] ); // Store tgtMask indices
        }
        if (!selectedAtomIndices.empty()) {
          SymmetricAtomIndices_.push_back( selectedAtomIndices );
          resHasSelectedSymmAtoms = true;
        }
      }
      // If remapping and not all atoms in a residue are selected, warn user.
      // TODO: Should they just be considered even if not selected?
      if (remapIn && resHasSelectedSymmAtoms) {
        for (int atom = topIn.Res(res).FirstAtom(); atom != topIn.Res(res).LastAtom(); ++atom)
          if (SelectedIdx[atom] == -1) {
            mprintf("Warning: Not all atoms selected in residue '%s'. Re-mapped\n"
                    "Warning:   structures may appear distorted.\n", 
                    topIn.TruncResNameNum(res).c_str());
            break;
          }
      }
    }
  }
  if (debug_ > 0) {
    mprintf("DEBUG: Potential Symmetric Atom Groups:\n");
    for (AtomIndexArray::const_iterator symmatoms = SymmetricAtomIndices_.begin();
                                        symmatoms != SymmetricAtomIndices_.end();
                                        ++symmatoms)
    {
      mprintf("\t%8u) ", symmatoms - SymmetricAtomIndices_.begin());
      for (Iarray::const_iterator atom = symmatoms->begin();
                                  atom != symmatoms->end(); ++atom)
        mprintf(" %s(%i)", topIn.AtomMaskName(tgtMask[*atom]).c_str(), tgtMask[*atom] + 1);
      mprintf("\n");
    } 
  }
  return 0;
}

/** It is expected that TGT and REF correspond to each other. */
double SymmetricRmsdCalc::SymmRMSD(Frame const& selectedTGT, Frame& selectedREF) {
  selectedREF.CenterOnOrigin( useMass_ );
  return SymmRMSD_CenteredRef(selectedTGT, selectedREF);
}

// SymmetricRmsdCalc::SymmRMSD()
/** selectedTgt and centeredREF must correspond to each other. */
double SymmetricRmsdCalc::SymmRMSD_CenteredRef(Frame const& selectedTgt, Frame const& centeredREF)
{
  // Create initial 1 to 1 atom map for all atoms; indices in 
  // SymmetricAtomIndices will correspond to positions in AMap.
  for (int atom = 0; atom < (int)AMap_.size(); atom++)
    AMap_[atom] = atom;
  tgtRemap_.SetCoordinates(selectedTgt);
  // Calculate initial best fit RMSD if necessary
  if (fit_) {
    tgtRemap_.RMSD_CenteredRef(centeredREF, rotMatrix_, tgtTrans_, useMass_);
    // Since tgtRemap is moved to origin during RMSD calc and centeredREF
    // should already be at the origin, just rotate.
    tgtRemap_.Rotate( rotMatrix_ );
  }
  // Correct RMSD for symmetry
  for (AtomIndexArray::const_iterator symmatoms = SymmetricAtomIndices_.begin();
                                      symmatoms != SymmetricAtomIndices_.end(); ++symmatoms)
  {
    // For each array of symmetric atoms, determine the lowest distance score
#   ifdef DEBUGSYMMRMSD
    mprintf("    Symmetric atoms group %u starting with atom %i\n", 
            symmatoms - SymmetricAtomIndices_.begin(), tgtMask_[symmatoms->front()] + 1);
#   endif
    cost_matrix_.Initialize( symmatoms->size() );
    for (Iarray::const_iterator ta = symmatoms->begin(); ta != symmatoms->end(); ++ta)
    {
      for (Iarray::const_iterator ra = symmatoms->begin(); ra != symmatoms->end(); ++ra)
      { 
        double dist2 = DIST2_NoImage( centeredREF.XYZ(*ra), tgtRemap_.XYZ(*ta) );
#       ifdef DEBUGSYMMRMSD
        mprintf("\t\t%i to %i: %f\n", tgtMask_[*ta] + 1, tgtMask_[*ra] + 1, dist2);
#       endif
        cost_matrix_.AddElement( dist2 );
      }
    }
    Iarray resMap = cost_matrix_.Optimize();
#   ifdef DEBUGSYMMRMSD
    mprintf("\tMapping from Hungarian Algorithm:\n");
    for (Iarray::const_iterator ha = resMap.begin(); ha != resMap.end(); ++ha)
      mprintf("\t\tMap col=%u row=%i\n", ha - resMap.begin(), *ha);
#   endif
    // Fill in overall map
    Iarray::const_iterator rmap = resMap.begin();
    for (Iarray::const_iterator atmidx = symmatoms->begin();
                                atmidx != symmatoms->end(); ++atmidx, ++rmap)
    {
      AMap_[*atmidx] = (*symmatoms)[*rmap];
#     ifdef DEBUGSYMMRMSD
      mprintf("\tAssigned atom %i to atom %i\n", tgtMask_[*atmidx] + 1,
              tgtMask_[(*symmatoms)[*rmap]] + 1);
#     endif
    }
  }
# ifdef DEBUGSYMMRMSD
  mprintf("    Final Atom Mapping:\n");
  for (unsigned int ref = 0; ref < AMap_.size(); ++ref)
    mprintf("\t%u -> %i\n", tgtMask_[ref] + 1, tgtMask_[AMap_[ref]] + 1);
  mprintf("----------------------------------------\n");
# endif
  // Remap the target frame for symmetry, then calculate new RMSD.
  // TODO: Does the topology need to be remapped as well?
  double rmsdval;
  tgtRemap_.SetCoordinatesByMap(selectedTgt, AMap_);
  if (fit_)
    rmsdval = tgtRemap_.RMSD_CenteredRef( centeredREF, rotMatrix_, tgtTrans_, useMass_ );
  else
    rmsdval = tgtRemap_.RMSD_NoFit( centeredREF, useMass_ );
  return rmsdval;
}
