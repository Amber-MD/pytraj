#ifndef INC_ATOMMAP_H
#define INC_ATOMMAP_H
#include "MapAtom.h"
#include "Topology.h"
/// Used to set up mapping information for each atom.
class AtomMap {
  public:
    AtomMap();

    MapAtom& operator[](int);
    const MapAtom& operator[](int i) const { return mapatoms_[i]; } // FIXME: bounds
    /// \return the number of atoms in the AtomMap.
    int Natom() const { return (int)mapatoms_.size(); }
    /// Set the debug level of the AtomMap.
    void SetDebug(int d) { debug_ = d; }
    /// Setup AtomMap with all atoms from input Topology.
    int Setup(Topology const&);
    /// Setup AtomMap with just atoms from specified Residue.
    int SetupResidue(Topology const&,int);
    /// Reset any previously set mapping information.
    void ResetMapping();
    bool BondIsRepeated(int,int) const;
    void DetermineAtomIDs();
    void MarkAtomComplete(int,bool);
    void CheckForCompleteAtoms();
    typedef std::vector<int> Iarray;
    typedef std::vector<Iarray> AtomIndexArray;
    int SymmetricAtoms(Topology const&, AtomIndexArray&, int);
  private:
    /// Check if 1 char name set to 0, means unidentified element.
    bool InvalidElement();
    int CheckBonds();
    void FindSymmetricAtoms(int, std::string const&, Iarray&, Iarray&) const;

    static MapAtom EMPTYMAPATOM;
    typedef std::vector<MapAtom> Marray;
    Marray mapatoms_;
    int debug_;
};
#endif
