# distutil: language = c++

from pytraj.Topology cimport *
from pytraj.ImageTypes cimport *

cdef extern from "ImageRoutines.h" namespace "Image":
    PairType CreatePairList(const _Topology&, Mode, _AtomMask)
    _Vec3 SetupTruncoct(const _Frame&, _AtomMask*, bint, bint)
    void Nonortho(_Frame&, bint, const _Vec3&, const _Vec3&, const _Matrix_3x3&,
                  const _Matrix_3x3&,
                  bint, bint, bint, const PairType&)
    _Vec3 Nonortho(const _Vec3&, bint, bint,
                  const _Matrix_3x3&, const _Matrix_3x3&, const _Vec3&, double)
    int SetupOrtho(const _Box&, _Vec3&, _Vec3&, bint)
    void Ortho(_Frame&, const _Vec3&, const _Vec3&, const _Vec3&,
              bint, bint, const PairType&)
    _Vec3 Ortho(const _Vec3&, const _Vec3&, const _Vec3&, const _Box&)
    void UnwrapNonortho(_Frame&, _Frame&, const PairType&, const _Matrix_3x3&, const _Matrix_3x3&, bint, bint)
    void UnwrapOrtho(_Frame&, _Frame&, const PairType&, bint ,bint)
