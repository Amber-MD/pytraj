# distutils: language = c++
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as incr
from pytraj._utils cimport get_positive_idx

from pytraj.externals.six import string_types
from pytraj.cpptraj_dict import TrajModeDict, get_key



cdef class TrajinList:
    def __cinit__(self):
        self.thisptr = new _TrajinList()
        self.py_free_mem = True

    def __dealloc__(self):
        if self.py_free_mem:
            del self.thisptr

    def clear(self):
        self.thisptr.Clear()

    def add_traj(self, filename not None, top=None, arglist=None):
        cdef TopologyList tlist = TopologyList()
        cdef ArgList _arglist

        filename = filename.encode()
        if arglist is not None:
            _arglist = ArgList(arglist)
        else:
            _arglist = ArgList()
        if top is None:
            # use self.top
            tlist.add_parm(self.top)
        else:
            tlist.add_parm(top)
            # update self.top too
            self.top = top
        self.thisptr.AddTrajin(filename, _arglist.thisptr[0], tlist.thisptr[0])

    def add_ensemble(self, filename not None, top=None, arglist=None):
        cdef TopologyList tlist = TopologyList()
        cdef ArgList _arglist

        filename = filename.encode()
        if arglist is not None:
            if isinstance(arglist, ArgList):
                _arglist = <ArgList> arglist
            else:
                _arglist = ArgList(arglist)
        else:
            arglist = ArgList()

        if top is None:
            # use self.top
            tlist.add_parm(self.top)
        else:
            tlist.add_parm(top)
            # update self.top too
            self.top = top
        tlist.add_parm(top)
        self.thisptr.AddEnsemble(filename, _arglist.thisptr[0], tlist.thisptr[0])

    def __iter__(self):
        cdef Trajin trajin
        cdef cppvector[_Trajin*].const_iterator it
        it = self.thisptr.begin()

        while it != self.thisptr.end():
            trajin = Trajin()
            # FIXME: got segmentation fault if we set topology here
            # is this because we're using pointer?
            # we need to sub-class at Python level (not Cython level)
            #trajin.top = self.top
            #trajin.top.py_free_mem = False
            # use memoryview rather making instance copy
            trajin.baseptr_1 = deref(it)
            # recast trajin.baseptr0 too
            trajin.baseptr0 = <_TrajectoryFile*> trajin.baseptr_1
            yield trajin
            incr(it)

    def frame_iter(self):
        if self.top == None:
            raise ValueError("need to set top for TrajinList")
        for traj in self:
            for frame in traj:
                yield frame

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
        cdef int s = 0

        for traj in self:
            if s == idx:
                return traj
            s += 1

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

    def is_empty(self):
        return self.thisptr.empty()

    @property
    def mode(self):
        return get_key(self.thisptr.Mode(), TrajModeDict)

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
