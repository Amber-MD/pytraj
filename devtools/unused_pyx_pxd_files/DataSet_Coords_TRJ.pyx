# distutils: language = c++
from ..decorators import for_testing
from ..externals.six import string_types


cdef class DataSet_Coords_TRJ(DataSet_Coords):
    def __cinit__(self, *args, **kwd):
        # TODO : do we really need to epoxe _DataSet and _DataSet_1D?
        # seriouly I need to use 4 pointers for class inheritance.
        # use pointer casting instead? (look ugly?)
        self.baseptr0 = <_DataSet*> new _DataSet_Coords_TRJ()
        # recast
        self.baseptr_1 = <_DataSet_Coords*> self.baseptr0
        self.thisptr = <_DataSet_Coords_TRJ*> self.baseptr0

        self.py_free_mem = True
        self.load(*args, **kwd)

    def __dealloc__(self):
        if self.py_free_mem:
            del self.thisptr
    
    @classmethod
    def alloc(cls):
        """return base class: DataSet"""
        cdef DataSet dset = DataSet()
        dset.baseptr0 = _DataSet_Coords_TRJ.Alloc()
        return dset

    def load(self, filename=None, top=None, arg=None):
        cdef Topology tmp_top
        cdef ArgList _arglist

        if filename is None and top is None and arg is None:
            pass
            # used for empty constructor in __init__
        else:
            if isinstance(top, string_types):
                self.top = Topology(top)
            elif isinstance(top, Topology):
                self.top = top.copy()

            if self.top.is_empty():
                raise ValueError("need to have non-empty topology file")

            # cast to Topology type so we can use cpptraj method
            tmp_top = <Topology> self.top

            if isinstance(filename, string_types):
                filename = filename.encode()
                if arg is None:
                    _arglist = ArgList()
                elif isinstance(arg, string_types):
                    _arglist = ArgList(arg)
                elif isinstance(arg, ArgList):
                    _arglist = arg
                else:
                    raise ValueError("arg must be None, string type or ArgList object")
                self.thisptr.AddSingleTrajin(filename, _arglist.thisptr[0], tmp_top.thisptr)
            elif isinstance(filename, (list, tuple)):
                # rename to avoid confusion
                filename_list = filename
                # recursive
                for fn in filename_list:
                    self.load(fn, top, arg)

    def add_trajin(self, Trajin trajin):
        """add memoryview for input trajin"""
        self.thisptr.AddInputTraj(trajin.baseptr_1)
