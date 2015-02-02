# distutils: language = c++
from pycpptraj.cpptraj_dict import TrajFormatDict, get_key

cdef class TrajectoryFile:
    """
    Base class that all input and output trajectories will inherit.

    Class hierarchy:
    TrajectoryFile --> Trajin --> Trajin_Single
    """

    def __cinit__(self):
        # ABC
        self._top = Topology()

        # let cpptraj free memory for self._top
        # Really?
        self._top.py_free_mem = False

        # we don't need to initialize self._top here
        # check "property" section
        # set memory view for self._top.thisptr
        #if self.baseptr0.TrajParm():
        #    self._top.thisptr = self.baseptr0.TrajParm()

    def __dealloc__(self):
        # ABC
        pass

    @classmethod
    def read_options(cls):
        _TrajectoryFile.ReadOptions()

    @classmethod
    def write_options(cls):
        _TrajectoryFile.WriteOptions()

    @classmethod
    def get_format(cls, arg):
        """
        Return format
        Parameters
        ---------
        arg : ArgList instance or string
        """
        cdef ArgList arglist
        cdef string s
        if isinstance(arg, ArgList):
            arglist = <ArgList> arg
            return _TrajectoryFile.GetFormatFromArg(arglist.thisptr[0])
        elif isinstance(arg, basestring):
            s = arg
            return _TrajectoryFile.GetFormatFromString(s)

    @classmethod
    def get_ext_for_type(cls, string key):
        """
        Parameters
        ----------
        value : str
        """
        return _TrajectoryFile.GetExtensionForType(TrajFormatDict[key])

    @classmethod
    def get_type_from_ext(cls, string e):
        """
        Parameters
        ----------
        value : str
        """
        ttype = _TrajectoryFile.GetTypeFromExtension(e)
        return get_key(ttype, TrajFormatDict)

    @classmethod
    def format_string(cls, string key):
        return _TrajectoryFile.FormatString(TrajFormatDict[key])

    def set_debug(self, int debug=0):
        self.baseptr0.SetDebug(debug)

    def set_trajfilename(self, string filename, bint is_read=True):
        self.baseptr0.SetTrajFileName(filename, is_read)

    # this decorator does not work
    # Python complains that 'property' is not callable
    # check the below solution
    #@property
    #def top(self):
    #    self._top.thisptr = self.baseptr0.TrajParm()
    #    return self._top

    #@top.setter
    #def top(self, Topology other):
    #    self.baseptr0.SetTrajParm(other.thisptr)

    property top:
        def __get__(self):
            if not self._top.is_empty():
                self._top.thisptr = self.baseptr0.TrajParm()
            return self._top

        def __set__(self, Topology other):
            # make a copy
            cdef Topology newtop = other.copy()
            # since we will pass a pointer to SetTrajParm method,
            # we let cpptraj frees memory
            newtop.py_free_mem = False

            self.baseptr0.SetTrajParm(newtop.thisptr)
            self._top.thisptr = self.baseptr0.TrajParm()

    def trajfilename(self):
        cdef FileName filename = FileName()
        if not filename.thisptr:
            raise MemoryError("Can not get Filename instance")
        filename.thisptr[0] = self.baseptr0.TrajFilename()
        return filename
