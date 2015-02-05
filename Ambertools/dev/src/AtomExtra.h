#ifndef INC_ATOMEXTRA_H
#define INC_ATOMEXTRA_H
#include "NameType.h"
/// Hold extra atom info.
class AtomExtra {
  public:
    AtomExtra() : join_(0), irotat_(0), atom_altloc_(' '), occupancy_(1.0), bfactor_(0.0) {}
    AtomExtra(float oc, float bf) :
      join_(0), irotat_(0), atom_altloc_(' '), occupancy_(oc), bfactor_(bf) {}
    AtomExtra(NameType const& it, int jo, int ir, char al) :
      itree_(it), join_(jo), irotat_(ir), atom_altloc_(al), occupancy_(1.0), bfactor_(0.0) {}
    NameType const& Itree() const { return itree_;       }
    int Join()              const { return join_;        }
    int Irotat()            const { return irotat_;      }
    char AtomAltLoc()       const { return atom_altloc_; }
    float Occupancy()       const { return occupancy_;   }
    float Bfactor()         const { return bfactor_;     }
  private:
    // Amber extra info.
    NameType itree_;
    int join_;
    int irotat_;
    // PDB extra info.
    char atom_altloc_; ///< Alternate location indicator.
    float occupancy_;  ///< Atom occupancy.
    float bfactor_;    ///< Atom B-factor.
};
#endif
