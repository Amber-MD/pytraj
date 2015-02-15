# distutils: language = c++
from pytraj.externals.six import string_types
from pytraj.cpptraj_dict import TrajFormatDict
from pytraj.utils.check_and_assert import file_exist


cdef class Trajout:
    formats = TrajFormatDict.keys()
    """Writing output
    Parameters :
    filename: str

    fmt : str, optional, default='AMBERTRAJ'
        output format: %s 

    if `fmt` is not provided, Trajout will decide format based on extension.
    if not `fmt` and no extension, default format = AMBERTRAJ

    So the priority is fmt > extension > default
        
    """
    def __cinit__(self, *args, **kwd):
        self.thisptr = new _Trajout()

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

    #@property
    #def formats(self):
    #    """return a list of possible format to be used with self.open"""
    #    return TrajFormatDict.keys()
        
    def open(self, filename='', top=Topology(), fmt='UNKNOWN_TRAJ', 
             more_args=None, overwrite=False):
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
        fmt = fmt.upper()

        if fmt == "PDB" or fmt == "MOL2":
            # add 'FILE' the end
            # 'PDBFILE' 'MOL2FILE'
            fmt += 'FILE'

        #filename = filename.encode("UTF-8")
        if more_args:
            if isinstance(more_args, string_types):
                inputstring = more_args
                arglist = <ArgList> ArgList(inputstring) 
            elif isinstance(more_args, ArgList):
                arglist = <ArgList> more_args
            else:
                raise ValueError()
            #self.thisptr.InitTrajWrite(filename.encode("UTF-8"), arglist.thisptr[0], top_.thisptr, local_dict[fmt])
            self.thisptr.InitTrajWrite(filename, arglist.thisptr[0], top_.thisptr, local_dict[fmt])
        else:
            self.thisptr.InitTrajWrite(filename, top_.thisptr, local_dict[fmt])

    def close(self):
        self.thisptr.EndTraj()

    def writeframe(self, *args, **kwd):
        self.write_frame(*args, **kwd)

    def write_frame(self, int idx=0, Frame frame=Frame(), top=Topology()):
        """write trajout for Frame with given Topology
        Parameters:
        ----------
        frame : Frame instance
        top : Topology instance
        """
        cdef Topology _top
        # check Topology
        if isinstance(top, string_types):
            top_ = Topology(top)
        elif isinstance(top, Topology):
            # assume this is Topology instance
            top_ = <Topology> top
        if len(top) == 0:
            # we use `len` here since we don't know if this is string or 
            # Topology object
            raise ValueError("require non-empty topology")
        self.thisptr.WriteFrame(idx, top_.thisptr, frame.thisptr[0])

    def print_info(self, int idx):
        self.thisptr.PrintInfo(idx)

    def is_open(self):
        return self.thisptr.TrajIsOpen()

    def nframes_processed(self):
        return self.thisptr.NumFramesProcessed()
