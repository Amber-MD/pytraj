# distutils: language = c++
from cython.operator cimport dereference as deref

from pytraj.externals.six import PY2, PY3, string_types

cdef class TopologyList:
    def __cinit__(self, py_free_mem=True):
        # py_free_mem is a flag to tell pytraj should free memory or let
        # cpptraj does
        # check ./CpptrajState.pyx

        cdef string filename
        self.thisptr = new _TopologyList()
        self.py_free_mem = py_free_mem

    def __dealloc__(self):
        if self.py_free_mem:
            del self.thisptr

    def clear(self):
        self.thisptr.Clear()

    def set_debug(self, int id):
        self.thisptr.SetDebug(id)

    def __getitem__(self, int idx):
        """return a copy of topology instance in TopologyList with index idx
        Input:
        =====
        idx : int
        """

        try:
           return (self.get_parm(idx)).copy()
        except:
            raise ValueError("index is out of range? do you have empty list?")

    def __setitem__(self, int idx, Topology other):
        cdef _Topology* topptr
        topptr = self.thisptr.GetParm(idx)
        topptr[0] = other.thisptr[0]

    #def __iter__(TopologyList self):
    #    cdef Topology top
    #    cdef int idx

    #    for i in range(self.size):
    #        top = self[idx]
    #        yield top

    #@property
    #def size(self):
    #    return self.thisptr.Size()

    def get_parm(self, arg):
        # TODO: checkbound
        """Return a Topology instance as a view to Topology instance in TopologyList
        If you made change to this topology, TopologyList would update this change too.
        """

        cdef Topology top = Topology()
        cdef int num
        cdef ArgList argIn
        # Explicitly del this pointer since we later let cpptraj frees memory
        # (we are use memoryview, not gettting a copy)
        del top.thisptr

        # since we pass C++ pointer to top.thisptr, we let cpptraj take care of 
        # freeing memory
        top.py_free_mem = False

        if isinstance(arg, (int, long)):
            num = arg
            #top.thisptr[0] = deref(self.thisptr.GetParm(num))
            # use memoryview instead making a copy
            top.thisptr = self.thisptr.GetParm(num)
            if not top.thisptr:
                raise IndexError("Out of bound indexing or empty list")
        if isinstance(arg, ArgList):
            argIn = arg
            # use memoryview instead making a copy
            top.thisptr = self.thisptr.GetParm(argIn.thisptr[0])
            #top.thisptr[0] = deref(self.thisptr.GetParm(argIn.thisptr[0]))
        return top

    def get_parm_from_pylist(self, list listin):
        cdef Topology top
        for top in listin:
            self.add_parm(top)

    def get_parm_by_index(self, ArgList argIn):
        """TODO: what is the difference between get_parm_by_index and  get_parm?"""
        cdef Topology top = Topology()
        top.py_free_mem = False
        top.thisptr = self.thisptr.GetParmByIndex(argIn.thisptr[0])
        return top

    def add_parm(self, *args):
        """Add parm file from file or from arglist or from Topology instance
        Input:
        =====
        filename :: str
        or
        arglist :: ArgList instance
        """
        # TODO: what's about adding Topology instance to TopologyList?
        cdef string filename
        cdef ArgList arglist
        cdef Topology top
        cdef Topology newtop

        if len(args) == 1:
            if isinstance(args[0], Topology):
                top = args[0]
                # let cpptraj frees memory since we pass a pointer
                newtop = top.copy()
                newtop.py_free_mem = False
                self.thisptr.AddParm(newtop.thisptr)
            else:
                filename = args[0].encode()
                self.thisptr.AddParmFile(filename)
        elif len(args) == 2:
            filename, arglist = args
            self.thisptr.AddParmFile(filename, arglist.thisptr[0])
        
    def info(self):
        self.thisptr.List()
