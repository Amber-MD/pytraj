# distutils: language = c++
from cython.operator cimport dereference as deref
from pytraj.externals.six import string_types


cdef class CpptrajState:
    """
    CpptrajState hold instances of:
    + TopologyList
    + FrameList (having reference frames)
    + DataSetList (having output data)
    + DataFileList

    TODO : add_trajin from Trajin or FrameArray instance too
    """
    def __cinit__(self):
        self.thisptr = new _CpptrajState()
        self.toplist = TopologyList(py_free_mem=False)
        #self.framelist = FrameList(py_free_mem=False)
        self.datafilelist = DataFileList(py_free_mem=False)
        self.datasetlist = DataSetList(py_free_mem=False)

        # cpptraj will take care of memory deallocating from self.thisptr.PFL(FL, DSL, DFL)
        # We don't free memory again (example: self.toplist.thisptr and self.thisptr.PFL() point to the same address)

        # create memory view
        self.toplist.thisptr = self.thisptr.PFL()
        #self.framelist.thisptr = self.thisptr.FL()
        self.datasetlist.thisptr = self.thisptr.DSL()
        self.datafilelist.thisptr = self.thisptr.DFL()

    def __dealloc__(self):
        if self.thisptr is not NULL:
            del self.thisptr
    
    def set_no_exit_on_error(self):
        self.thisptr.SetNoExitOnError()

    def set_no_progress(self):
        self.thisptr.SetNoProgress()

    def set_action_silence(self, bint b):
        self.thisptr.SetActionSilence(b)

    #def debug(self):
    #    return self.thisptr.Debug()

    def exit_on_error(self):
        return self.thisptr.ExitOnError()

    def is_empty(self):
        return self.thisptr.EmptyState()

    def add_trajin(self, arg, isEnsemble=False):
        # TODO: add trajector instance?
        cdef string filename
        cdef ArgList argIn
        
        if isinstance(arg, ArgList):
            argIn = arg
            self.thisptr.AddTrajin(argIn.thisptr[0], isEnsemble)
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
            #trajlist.thisptr = &(self.thisptr.InputTrajList())
            return trajlist
        else:
            raise MemoryError("")

    def add_trajout(self, arg):
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

    #def add_reference(self, arg):
    #    cdef string name
    #    cdef ArgList arglist

    #    if isinstance(arg, str):
    #        name = arg
    #        self.thisptr.AddReference(name)
    #    elif isinstance(arg, ArgList):
    #        arglist = arg
    #        self.thisptr.AddReference(arglist.thisptr[0])
    #    else:
    #        raise NotImplementedError()

    def add_action(self, obj, ArgList arglist):
        # need to explicit casting to FunctPtr because self.thisptr.AddAction need to know type 
        # of variables
        cdef FunctPtr alloc_funct = <FunctPtr> obj.alloc()
        return self.thisptr.AddAction(alloc_funct.ptr, arglist.thisptr[0])

    def add_analysis(self, obj, ArgList arglist):
        """temp doc: add_analysis(self, obj, ArgList arglist)
        obj :: Action or Analysis instance
        >>> obj = Action_Rmsd()
        """
        cdef FunctPtr alloc_funct = <FunctPtr> obj.alloc()
        return self.thisptr.AddAnalysis(alloc_funct.ptr, arglist.thisptr[0])

    # what is it?
    #def world_size(self):
    #    return self.thisptr.WorldSize()

    def list_all(self, ArgList arglist):
        return self.thisptr.ListAll(arglist.thisptr[0])

    def set_list_debug(self, ArgList arglist):
        return self.thisptr.SetListDebug(arglist.thisptr[0])

    def clear_list(self, ArgList arglist):
        return self.thisptr.ClearList(arglist.thisptr[0])

    def remove_data_set(self, ArgList alist):
        return self.thisptr.RemoveDataSet(alist.thisptr[0])

    def traj_length(self, string topname, vector[string] trajlist):
        return self.thisptr.TrajLength(topname, trajlist)

    def run(self):
        return self.thisptr.Run()

    def master_data_file_write(self):
        self.thisptr.MasterDataFileWrite()
