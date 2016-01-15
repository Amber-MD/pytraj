# distutils: language = c++
from pytraj.externals.six import string_types
from pytraj.c_dict import TrajFormatDict
from pytraj.utils.check_and_assert import file_exist


cdef class TrajectoryWriter:
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

    def open(self, filename='', top=Topology(),
             format='infer',
             options='', overwrite=False):
        '''
        filename : str, output filename
        top : Topology
        format : str, default 'infer'
            if 'infer', determine file format based on extension.
            If can not detect extension, use AMBER mdcrd format
        options : str, additional keywords for writing file (good for pdb, mol2, ...)
        overwrite : bool, default False
        '''

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

        if format.lower() == 'infer':
            options += ''
        else:
            options = ' '.join((format.lower(), options))

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

    @classmethod
    def get_formats(cls):
        return list(TrajFormatDict.keys())
