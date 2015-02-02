#ifndef INC_MAPATOM_H
#define INC_MAPATOM_H
#include "Atom.h"
/// Atom with extra information that can be used for mapping.
class MapAtom : public Atom {
  public :
    MapAtom();
    MapAtom(const MapAtom&);
    MapAtom(const Atom&);
    MapAtom& operator=(const MapAtom&);

    bool IsChiral()             const { return isChiral_; }
    bool BoundToChiral()        const { return boundToChiral_; }
    bool IsMapped()             const { return isMapped_; }
    bool Complete()             const { return complete_; }
    bool IsUnique()             const { return (Nduplicated_ == 0); }
    const std::string& AtomID() const { return atomID_; }
    const std::string& Unique() const { return unique_;}
    int Nduplicated()           const { return Nduplicated_; }
    char CharName()             const { return name_; }

    void IsDuplicated()                  { ++Nduplicated_; }
    void SetMapped()                     { isMapped_ = true; } 
    void SetComplete()                   { complete_ = true; }
    void SetChiral()                     { isChiral_ = true; }
    void SetBoundToChiral()              { boundToChiral_ = true; }
    void SetAtomID(std::string const& s) { atomID_ = s; }
    void SetUnique(std::string const& s) { unique_ = s; Nduplicated_ = 0; }
    void SetNotMapped()                  { isMapped_ = false; }
    void SetNotComplete()                { complete_ = false; }
    void SetNotChiral()                  { isChiral_ = false; }
  private:
    static const char AtomicElementChar[];
    bool isChiral_;      ///< true: Atom is a chiral center
    bool boundToChiral_; ///< true: Atom is bound to a chiral center.
    bool isMapped_;      ///< true: this atom has been mapped
    bool complete_;      ///< true: This atom an all bonded atoms have been mapped
    std::string atomID_; ///< ID created from this atom name, then bonded atom names
    std::string unique_; ///< ID created from this atomID, then bonded atomIDs
    int Nduplicated_;    ///< How many times is uniqueID duplicated, 0 = IsUnique
    char name_;          ///< 1 char name used to construct atomID/unique strings
};
#endif
