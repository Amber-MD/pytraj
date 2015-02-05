#ifndef INC_IMAGETYPES_H
#define INC_IMAGETYPES_H
/*! \file ImageTypes.h
    \brief Data types and enumerations used by imaging routines. 
 */
#include <vector>
namespace Image {
  typedef std::vector<int> PairType;
  enum Mode { BYMOL = 0, BYRES, BYATOM };
  inline const char* ModeString(Mode m) {
    if      (m == BYMOL) return "molecule";
    else if (m == BYRES) return "residue";
    else                 return "atom"; // BYATOM
  }
}
#endif
