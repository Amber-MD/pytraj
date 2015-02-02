# distutils: language = c++
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as incr


cdef class TrajinList:
    def __cinit__(self):
        self.thisptr = new _TrajinList()
        self.py_free_mem = True

    def __dealloc__(self):
        if self.py_free_mem:
            del self.thisptr

    def clear(self):
        self.thisptr.Clear()

    def set_debug(self,int dIn):
        self.thisptr.SetDebug(dIn)

    def add_traj(self, string filename, ArgList arglist, TopologyList toplist):
        return self.thisptr.AddTrajin(filename, arglist.thisptr[0], toplist.thisptr[0])

    def add_ensemble(self, string s, ArgList arglist, TopologyList toplist):
        return self.thisptr.AddEnsemble(s, arglist.thisptr[0], toplist.thisptr[0])

    def __iter__(self):
        cdef Trajin trajin
        cdef cppvector[_Trajin*].const_iterator it
        it = self.thisptr.begin()
        while it != self.thisptr.end():
            trajin = Trajin()
            # use memoryview rather making instance copy
            trajin.baseptr_1 = deref(it)
            # recast trajin.baseptr0 too
            trajin.baseptr0 = <_TrajectoryFile*> trajin.baseptr_1
            yield trajin
            incr(it)

    @property
    def size(self):
        cdef cppvector[_Trajin*].const_iterator it
        it = self.thisptr.begin()
        s = 0
        while it != self.thisptr.end():
            s += 1
            incr(it)
        return s

    def __getitem__(self, int idx):
        """return Trajin instance
        TODO: return Trajin or Trajin_Single instance?
        """
        cdef Trajin trajin
        cdef cppvector[_Trajin*].const_iterator it
        it = self.thisptr.begin()

        if idx < 0 or idx >= self.size:
            raise ValueError("index is out of range")

        s = 0
        while it != self.thisptr.end():
            if idx == s:
                trajin = Trajin()
                # use memoryview rather making instance copy
                trajin.baseptr_1 = deref(it)
                # recast trajin.baseptr0 too
                trajin.baseptr0 = <_TrajectoryFile*> trajin.baseptr_1
                return trajin
            s += 1
            incr(it)
        return Trajin()

    def __setitem__(self, int idx, Trajin other):
        "TODO: validate"
        cdef cppvector[_Trajin*].const_iterator it
        cdef _Trajin* _trajinptr
        it = self.thisptr.begin()

        if idx < 0 or idx >= self.size:
            raise ValueError("index is out of range")

        s = 0
        while it != self.thisptr.end():
            if idx == s:
                _trajinptr = deref(it)
                _trajinptr[0] = other.baseptr_1[0]
            s += 1
            incr(it)

    def empty(self):
        return self.thisptr.empty()

    def mode(self, updatedmode=False):
        # Use "updatedmode" in case Dan Roe updates his TrajinList.Mode()
        
        TrajModeType_dict = {
                UNDEFINED : "UNDEFINED",
                NORMAL : "NORMAL",
                ENSEMBLE : "ENSEMBLE",
        }

        if not updatedmode:
            return TrajModeType_dict[self.thisptr.Mode()]
        else:
            raise NotImplementedError()

    def front(self):
        # TODO: add doc
        cdef Trajin trajin = Trajin()
        # create memoryview
        trajin.baseptr_1 = <_Trajin*> self.thisptr.front()
        # make sure two pointers pointing to the same address
        trajin.baseptr0 = <_TrajectoryFile*> trajin.baseptr_1
        return trajin

    def append(self, Trajin trajin):
        raise NotImplementedError("not yet")

    @property
    def max_frames(self):
        return self.thisptr.MaxFrames()

    def printlist(self):
        self.thisptr.List()
