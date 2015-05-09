# distutils: language = c++
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as incr
from pytraj._utils cimport get_positive_idx
from pytraj.trajs.Trajin_Single cimport _Trajin_Single, Trajin_Single
from pytraj.trajs.Trajin_Single cimport Trajin_Single as TrajectoryREMDIterator
from pytraj.TrajectoryREMDIterator cimport TrajectoryREMDIterator

from pytraj.TrajectoryIterator import TrajectoryIterator
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
        cdef Trajin_Single ts
        cdef cppvector[_Trajin*].const_iterator it
        it = self.thisptr.begin()

        while it != self.thisptr.end():
            ts = Trajin_Single()
            ts.baseptr_1 = deref(it) # baseptr_1 is from `Trajin`
            # need to cast other pointers too
            ts.baseptr0 = <_TrajectoryFile*> ts.baseptr_1
            # don't cast to `thisptr`, will get segfault
            # Why: because TrajinList might return Trajin_Multi (inherited from `Trajin`
            # but not `Trajin_Single`
            #ts.thisptr = <_Trajin_Single*> ts.baseptr_1
            ts.top = self.top.copy()
            yield ts
            incr(it)

    def _getitem_remd(self, idx):
        """return TrajectoryREMDIterator object
        """
        # Aim: to be used with `io.iterload_remd`
        cdef TrajectoryREMDIterator traj = TrajectoryREMDIterator()
        cdef Trajin_Single _traj
        cdef int s = 0

        for _traj in self:
            if s == idx:
                # casting
                traj.baseptr_1 = <_Trajin*> _traj.baseptr0
                # need to cast other pointers too
                traj.baseptr0 = <_TrajectoryFile*> traj.baseptr_1
                # dont' cast to _Trajin_Single* or will get segfault
                # we don't actually use _Trajin_Single here
                # ideally we can subclass `Trajin` but we can directly allocate class
                # having virtual method
                #traj.thisptr = <_Trajin_Single*> traj.baseptr_1
                #traj.py_free_mem = False
                traj.top = self.top.copy()
                return traj
            s += 1

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
