# iistutils: languag]e = c++
from pycpptraj.cpptraj_dict import ParmFormatDict

cdef class ParmFile:
    def __cinit__(self):
        self.thisptr = new _ParmFile()

    def __dealloc__(self):
        del self.thisptr

    @classmethod
    def help(cls):
        print "read_options"
        cls.read_options()

        print
        print "write_options"
        cls.write_options()

    @classmethod
    def read_options(cls):
        _ParmFile.ReadOptions()

    @classmethod
    def write_options(cls):
        _ParmFile.WriteOptions()

    @classmethod
    def formats(cls):
        """return a list of supported parm formats"""
        return ParmFormatDict.keys()

    def readparm(self, Topology top=Topology(), string filename="AMBERPARM", *args):
        """readparm(Topology top=Topology(), string filename="", "*args)
        Return : None (update `top`)

        top : Topology instance
        filename : str, output filename
        arglist : ArgList instance, optional

        """
        cdef ArgList arglist
        cdef debug = 0

        if not args:
            self.thisptr.ReadTopology(top.thisptr[0], filename, debug)
        else:
            arglist = <ArgList> args[0]
            self.thisptr.ReadTopology(top.thisptr[0], filename, arglist.thisptr[0], debug)

    def write_prefix_topology(self, Topology top=Topology(), string prefix="default", fmt=""):
        """write_prefix_topology(Topology Top, string prefix, fmt="AMBER")
        TODO : automatically get ParmFormatDict for this doc
        ParmFormatDict = {
            "AMBER" : AMBERPARM,
            "PDB" : PDBFILE,
            "MOL2" : MOL2FILE,
            "CHARMMPSF" : CHARMMPSF,
            "CIF" : CIFFILE,
            "SDF" : SDFFILE,
            "UNKNOWN" : UNKNOWN_PARM,
                         }
        TODO : do we need this method?
        """
        cdef int debug = 0
        cdef int err  
        cdef ParmFormatType parmtype 
        
        if fmt.empty():
            parmtype = UNKNOWN_PARM
        else:
            try:
                parmtype = ParmFormatDict[fmt]
            except:
                print "supported keywords: ", self.formats
        # TODO : combine with write_topology
        err = self.thisptr.WritePrefixTopology(top.thisptr[0], prefix, parmtype, debug)
        if err == 1:
            print "Not supported or failed to write"

    def writeparm(self, Topology top=Topology(), string filename="default.top", 
                  ArgList arglist=ArgList(), string fmt=""):
        """writeparm(Topology top=top, string filename="default.top", 
                     ArgList arglist=ArgList(), string fmt="AMBER")")
        TODO : automatically get ParmFormatDict for this doc
        ParmFormatDict = {
            "AMBER" : AMBERPARM,
            "PDB" : PDBFILE,
            "MOL2" : MOL2FILE,
            "CHARMMPSF" : CHARMMPSF,
            "CIF" : CIFFILE,
            "SDF" : SDFFILE,
            "UNKNOWN" : UNKNOWN_PARM,
                         }
        """
        cdef int debug = 0
        cdef int err
        # change `fmt` to upper
        cdef ParmFormatType parmtype 
        
        if fmt.empty():
            parmtype = UNKNOWN_PARM
        else:
            try:
                fmt = fmt.upper()
                parmtype = ParmFormatDict[fmt]
            except:
                print "supported keywords: ", self.formats

        if top.is_empty():
            raise ValueError("empty topology")

        err = self.thisptr.WriteTopology(top.thisptr[0], filename, arglist.thisptr[0], parmtype, debug)
        if err == 1:
            print "Not supported or failed to write"

    def filename(self):
        cdef FileName filename = FileName()
        #del filename.thisptr
        filename.thisptr[0] = self.thisptr.ParmFilename()
        return filename
