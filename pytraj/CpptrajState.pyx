# distutils: language = c++
from cython.operator cimport dereference as deref
from pytraj.externals.six import string_types


cdef class CpptrajState:
    """
    CpptrajState hold instances of:
    + TopologyList
    + DataSetList (having output data)
    + DataFileList

    TODO : add_trajin from Trajin or FrameArray instance too
    """
    def __cinit__(self):
        self.thisptr = new _CpptrajState()
        self.toplist = TopologyList(py_free_mem=False)
        self.datafilelist = DataFileList(py_free_mem=False)
        self.datasetlist = DataSetList(py_free_mem=False)

        # cpptraj will take care of memory deallocating from self.thisptr.PFL(FL, DSL, DFL)
        # We don't free memory again 
        # (example: self.toplist.thisptr and self.thisptr.PFL() point to the same address)
        # create memory view
        self.toplist.thisptr = self.thisptr.PFL()
        self.datasetlist.thisptr = self.thisptr.DSL()
        self.datafilelist.thisptr = self.thisptr.DFL()

    def __dealloc__(self):
        if self.thisptr is not NULL:
            del self.thisptr
    
    def is_empty(self):
        return self.thisptr.EmptyState()

    def add_trajin(self, arg, is_ensemble=False):
        # TODO: add trajector instance?
        cdef string filename
        cdef ArgList argIn
        
        if isinstance(arg, ArgList):
            argIn = arg
            self.thisptr.AddTrajin(argIn.thisptr[0], is_ensemble)
        elif isinstance(arg, string_types):
            filename = arg.encode()
            self.thisptr.AddTrajin(filename)
        else:
            raise NotImplementedError()

    def run_analyses(self):
        return self.thisptr.RunAnalyses()

    def get_trajinlist(self):
        """Return a copy of CpptrajState's TrajinList instance"""
        cdef TrajinList trajlist = TrajinList()

        # we need to let Cython know that cpptraj will free memory 
        # if I don't use this flag, the program will get segmentfault
        # Maybe Cython will try free memory of trajlist.thisptr twice? 
        # Reason? cpptraj does not have assignment operator for TrajinList class
        # --> ?

        # We used "get_trajlist" as method's name to remind that this will return an instance copy
        # Should we update other *.pyx files too?
        # so why do we need py_free_mem = False here?
        trajlist.py_free_mem = False
        
        if self.thisptr:
            trajlist.thisptr[0] = self.thisptr.InputTrajList()

            # does not work when using address (const)
            #trajlist.thisptr = &(self.thisptr.InputTrajList())
            return trajlist
        else:
            raise MemoryError("")

    def add_trajout(self, arg):
        """add trajout file
        
        Parameters
        ---------
        arg : str or ArgList object
        """
        cdef string filename
        cdef ArgList arglist

        if isinstance(arg, ArgList):
            arglist = arg
            return self.thisptr.AddTrajout(arglist.thisptr[0])
        elif isinstance(arg, string_types):
            filename = arg.encode()
            return self.thisptr.AddTrajout(filename)
        else:
            raise NotImplementedError()

    def add_reference(self, *args):
        """
        Parameters
        ---------
        filename : str
        arg : ArgList object, optional
        """
        cdef string name
        cdef ArgList arglist

        if len(args) == 1:
            if isinstance(args[0], string_types):
                name =  args[0].encode(0)
                self.thisptr.AddReference(name)
            else:
                raise NotImplementedError()
        elif len(args) == 2:
                name =  args[0].encode(0)
                if isinstance(args[1], string_types):
                    arglist = ArgList(args[1])
                else:
                    arglist = <ArgList> args[1]
                self.thisptr.AddReference(name, arglist.thisptr[0])
        else:
            raise NotImplementedError()

    def add_action(self, actobj, arglist):
        """
        Parameters
        ---------
        actobj : Action object or str
        arglist : ArgList object or str
        """
        # need to explicit casting to FunctPtr because self.thisptr.AddAction need to know type 
        # of variables
        cdef FunctPtr alloc_funct
        cdef ArgList _arglist 

        if isinstance(actobj, string_types):
            # if actobj is string, make Action object
            # then cast to FunctPtr
            from pytraj.action_dict import ADICT
            alloc_funct = ADICT[actobj]().alloc()
        else:
            alloc_funct = <FunctPtr> actobj.alloc()

        if isinstance(arglist, string_types):
            _arglist = ArgList(arglist)
        elif isinstance(arglist, ArgList):
            _arglist = arglist
        else:
            raise ValueError("must be string or ArgList object")

        return self.thisptr.AddAction(alloc_funct.ptr, _arglist.thisptr[0])

    def add_analysis(self, obj, ArgList arglist):
        """temp doc: add_analysis(self, obj, ArgList arglist)
        obj :: Action or Analysis instance
        """
        cdef ArgList _arglist 
        cdef FunctPtr alloc_funct = <FunctPtr> obj.alloc()

        if isinstance(arglist, string_types):
            _arglist = ArgList(arglist)
        elif isinstance(arglist, ArgList):
            _arglist = arglist
        else:
            raise ValueError("must be string or ArgList object")

        return self.thisptr.AddAnalysis(alloc_funct.ptr, _arglist.thisptr[0])

    def list_all(self, ArgList arglist):
        return self.thisptr.ListAll(arglist.thisptr[0])

    def clear_list(self, arglist='all'):
        return self.thisptr.ClearList(ArgList(arglist).thisptr[0])

    def remove_dataset(self, ArgList alist):
        return self.thisptr.RemoveDataSet(alist.thisptr[0])

    def run(self):
        return self.thisptr.Run()

    def write_all_datafiles(self):
        self.thisptr.MasterDataFileWrite()
