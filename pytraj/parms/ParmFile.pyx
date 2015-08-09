# iistutils: languag]e = c++
from pytraj.cpptraj_dict import ParmFormatDict

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

    def readparm(self, filename="", top=Topology(), *args):
        """readparm(Topology top=Topology(), string filename="", "*args)
        Return : None (update `top`)

        top : Topology instance
        filename : str, output filename
        arglist : ArgList instance, optional

        """
        cdef ArgList arglist
        cdef debug = 0
        cdef Topology _top = <Topology> top

        filename = filename.encode()

        if not args:
            self.thisptr.ReadTopology(_top.thisptr[0], filename, debug)
        else:
            arglist = <ArgList> args[0]
            self.thisptr.ReadTopology(_top.thisptr[0], filename, arglist.thisptr[0], debug)

    def write_prefix_topology(self, Topology top=Topology(), prefix="default", format=""):
        """write_prefix_topology(Topology Top, string prefix, format="AMBER")
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

        prefix = prefix.encode()
        
        if format.empty():
            parmtype = UNKNOWN_PARM
        else:
            try:
                parmtype = ParmFormatDict[format]
            except:
                print "supported keywords: ", self.formats
        # TODO : combine with write_topology
        err = self.thisptr.WritePrefixTopology(top.thisptr[0], prefix, parmtype, debug)
        if err == 1:
            print "Not supported or failed to write"

    def writeparm(self, Topology top=Topology(), filename="default.top", 
                  ArgList arglist=ArgList(), format=""):
        """writeparm(Topology top=top, string filename="default.top", 
                     ArgList arglist=ArgList(), string format="AMBER")")
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
        # change `for` to upper
        cdef ParmFormatType parmtype 
        filename = filename.encode()
        
        if format== "":
            parmtype = UNKNOWN_PARM
        else:
            try:
                format= format.upper()
                parmtype = ParmFormatDict[format]
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
