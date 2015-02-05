#ifndef INC_IMAGEROUTINES_H
#define INC_IMAGEROUTINES_H
#include "Topology.h"
#include "ImageTypes.h"
namespace Image {
  PairType CreatePairList(Topology const&, Mode, AtomMask);
  Vec3 SetupTruncoct( Frame const&, AtomMask*, bool, bool);
  void Nonortho(Frame&, bool, Vec3 const&, Vec3 const&, Matrix_3x3 const&, Matrix_3x3 const&,
                bool, bool, bool, PairType const&);
  Vec3 Nonortho(Vec3 const&, bool, bool,
                Matrix_3x3 const&, Matrix_3x3 const&, Vec3 const&, double);
  int SetupOrtho(Box const&, Vec3&, Vec3&, bool);
  void Ortho(Frame&, Vec3 const&, Vec3 const&, Vec3 const&, bool, bool, PairType const&);
  Vec3 Ortho(Vec3 const&, Vec3 const&, Vec3 const&, Box const&);
  void UnwrapNonortho( Frame&, Frame&, PairType const&, 
                       Matrix_3x3 const&, Matrix_3x3 const&, bool, bool );
  void UnwrapOrtho( Frame&, Frame&, PairType const&, bool, bool );
}
#endif
