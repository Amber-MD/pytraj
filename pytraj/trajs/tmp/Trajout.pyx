# distutils: language = c++
from pycpptraj.cpptraj_dict import TrajFormatDict


cdef class Trajout:
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

    @property
    def formats(cls):
        """return a list of possible format to be used with self.open"""
        return TrajFormatDict.keys()
        
    def open(self, string filename='', top=Topology(), fmt=None, more_args=None):
        cdef ArgList arglist
        cdef Topology top_ 

        # check Topology
        if isinstance(top, basestring):
            top_ = Topology(top)
        elif isinstance(top, Topology):
            # assume this is Topology instance
            top_ = top

        local_dict = TrajFormatDict.copy()
        local_dict.get("", "")

        if more_args:
            if isinstance(more_args, basestring):
                inputstring = more_args
                arglist = <ArgList> ArgList(inputstring) 
            elif isinstance(more_args, ArgList):
                arglist = <ArgList> more_args
            else:
                raise ValueError()
            self.thisptr.InitTrajWrite(filename, arglist.thisptr[0], top_.thisptr, local_dict[fmt])
        else:
            self.thisptr.InitTrajWrite(filename, top_.thisptr, local_dict[fmt])

    def close(self):
        self.thisptr.EndTraj()

    def writeframe(self, int idx=0, Frame frame=Frame(), top=Topology()):
        """write trajout for Frame with given Topology
        Parameters:
        ----------
        frame : Frame instance
        top : Topology instance
        """
        cdef Topology _top
        # check Topology
        if isinstance(top, basestring):
            top_ = Topology(top)
        elif isinstance(top, Topology):
            # assume this is Topology instance
            top_ = <Topology> top

        self.thisptr.WriteFrame(idx, top_.thisptr, frame.thisptr[0])

    def print_info(self, int idx):
        self.thisptr.PrintInfo(idx)

    def is_open(self):
        return self.thisptr.TrajIsOpen()

    def nframes_processed(self):
        return self.thisptr.NumFramesProcessed()
