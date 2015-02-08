# distutils: language = c++


cdef class MapAtom(Atom):
    def __cinit__(self, *args):
        cdef Atom atom
        cdef MapAtom map_atom

        if not args:
            self.thisptr_ma = new _MapAtom()
        elif len(args) == 1:
            if isinstance(args[0], Atom):
                atom = args[0]
                self.thisptr_ma = new _MapAtom(atom.thisptr[0])
            elif isinstance(args[0], MapAtom):
                map_atom = args[0]
                self.thisptr_ma = new _MapAtom(map_atom.thisptr_ma[0])
            else:
                raise ValueError()
        else:
            raise ValueError()

        # recast pointers since MapAtom inherits from Atom
        self.thisptr = <_Atom*> self.thisptr_ma

    def copy(self):
        """Return a copy of instance
        >>>from pytraj.MapAtom import MapAtom
        >>>ma1 = MapAtom()
        >>>ma2 = ma1.copy()
        >>>print ma1 == ma2
        """
        cdef MapAtom matom = MapAtom()
        matom.thisptr_ma[0] = self.thisptr_ma[0]
        return matom

    def is_chiral(self):
        return self.thisptr_ma.IsChiral()

    def bound_to_chiral(self):
        return self.thisptr_ma.BoundToChiral()

    def is_mapped(self):
        return self.thisptr_ma.IsMapped()

    def complete(self):
        return self.thisptr_ma.Complete()

    def is_unique(self):
        return self.thisptr_ma.IsUnique()

    def unique(self):
        return self.thisptr_ma.Unique()

    def n_duplicated(self):
        return self.thisptr_ma.Nduplicated()

    def char_name(self):
        return self.thisptr_ma.CharName()

    def is_duplicated(self):
        self.thisptr_ma.IsDuplicated()

    def set_mapped(self):
        self.thisptr_ma.SetMapped()

    def set_complete(self):
        self.thisptr_ma.SetComplete()

    def set_chiral(self):
        self.thisptr_ma.SetChiral()

    def set_bound_to_chiral(self):
        self.thisptr_ma.SetBoundToChiral()

    def set_atom_id(self, string s):
        self.thisptr_ma.SetAtomID(s)

    def set_unique(self, string s):
        self.thisptr_ma.SetUnique(s)

    def set_not_mapped(self):
        self.thisptr_ma.SetNotMapped()

    def set_not_complete(self):
        self.thisptr_ma.SetNotComplete()

    def set_not_chiral(self):
        self.thisptr_ma.SetNotChiral()
