# distutil: language = c++

cdef CreatePairList(const _Topology&, Mode mode, _AtomMask atm):
    return CreatePairList(top.thisptr[0], mode, atm.thisptr)

_Vec3 SetupTruncoct(const _Frame&, _AtomMask* , bint, bint)
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
