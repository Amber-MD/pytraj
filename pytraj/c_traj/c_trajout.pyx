# distutils: language = c++
from pytraj.externals.six import string_types
from pytraj.c_dict import TrajFormatDict
from pytraj.utils.check_and_assert import file_exist


cdef class TrajectoryWriter:
    formats = TrajFormatDict.keys()
    """Writing output

    Parameters
    ----------
    filename: str
    format: str, optional, default='AMBERTRAJ'
        output format: %s
        if `format` is not provided, TrajectoryWriter will decide format based on extension.
        if not `format` and no extension, default format = AMBERTRAJ
    So the priority is format> extension > default

    """

    def __cinit__(self, *args, **kwd):
        self.thisptr = new _Trajout()
        self.count = 0

        if args or kwd:
            self.open(*args, **kwd)

    def __dealloc__(self):
        del self.thisptr

    def __enter__(self):
        return self

    def __exit__(self, arg1, arg2, arg3):
        self.close()

    @classmethod
    def help(cls):
        print "TrajFormat"
        print TrajFormatDict.keys()

    def open(self, filename='', top=Topology(), format='default',
             options=None, overwrite=False):

        cdef ArgList arglist
        cdef Topology top_

        filename = filename.encode("UTF-8")
        if not overwrite:
            # TODO : what's about 'append'?
            if file_exist(filename):
                err = "file exist and you're in overwrite=%s mode" % str(overwrite)
                raise RuntimeError(err)

        # check Topology
        if isinstance(top, string_types):
            top_ = Topology(top)
        elif isinstance(top, Topology):
            # assume this is Topology instance
            top_ = top

        local_dict = TrajFormatDict.copy()
        local_dict.get("", "")
        # make upper case in case user uses lower ones
        format= format.upper()

        if format == "PDB" or format == "MOL2":
            # add 'FILE' the end
            # 'PDBFILE' 'MOL2FILE'
            format += 'FILE'

        if options:
            if isinstance(options, string_types):
                inputstring = options
                arglist = <ArgList> ArgList(inputstring)
            elif isinstance(options, ArgList):
                arglist = <ArgList> options
            else:
                raise ValueError()
            self.thisptr.InitTrajWrite(filename, arglist.thisptr[0], top_.thisptr)
        else:
            self.thisptr.InitTrajWrite(filename, ArgList().thisptr[0], top_.thisptr)

        # real open
        self.thisptr.SetupTrajWrite(top_.thisptr, CoordinateInfo(), 0)

    def close(self):
        self.thisptr.EndTraj()

    def write(self, Frame frame):
        """

        Parameters
        ----------
        frame : Frame instance

        *args, **kwd: just dummy
        """
        self.thisptr.WriteFrame(self.count, frame.thisptr[0])
        self.count += 1
