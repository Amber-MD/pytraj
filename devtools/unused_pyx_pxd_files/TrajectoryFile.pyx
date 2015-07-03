# distutils: language = c++
from pytraj.cpptraj_dict import TrajFormatDict, get_key
from pytraj.externals.six import string_types

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
    def _read_options(cls):
        _TrajectoryFile.ReadOptions()

    @classmethod
    def _write_options(cls):
        _TrajectoryFile.WriteOptions()

    @classmethod
    def _get_format(cls, arg):
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
        elif isinstance(arg, string_types):
            s = arg.encode()
            return _TrajectoryFile.GetFormatFromString(s)

    @classmethod
    def _get_ext_for_type(cls, key):
        """
        Parameters
        ----------
        value : str
        """
        return _TrajectoryFile.GetExtensionForType(TrajFormatDict[key]).decode()

    @classmethod
    def _get_type_from_ext(cls, e):
        """
        Parameters
        ----------
        value : str
        """
        e = e.encode()
        ttype = _TrajectoryFile.GetTypeFromExtension(e)
        return get_key(ttype, TrajFormatDict)

    @classmethod
    def _format_string(cls, key):
        return _TrajectoryFile.FormatString(TrajFormatDict[key]).decode()

    def _set_trajfilename(self, filename, bint is_read=True):
        filename = filename.encode()
        self.baseptr0.SetTrajFileName(filename, is_read)

    @property
    def filename(self):
        cdef FileName fname = FileName()
        fname.thisptr[0] = self.baseptr0.TrajFilename()
        return fname.__str__()

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
