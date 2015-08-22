# distutils: language = c++
from pytraj.cpptraj_dict cimport RetTypeAna, OKANALYSIS, ERRANALYSIS
from pytraj.externals.six import PY3
from pytraj.decorators import makesureABC
from pytraj.externals.six import string_types

cdef extern from "Analysis.h":
    ctypedef enum RetType "Analysis::RetType":
        OKANALYSIS "Analysis::OK"
        ERRANALYSIS "Analysis::ERR"

cdef class Analysis:
    """
    Original cpptraj doc:
    ====================
        The abstract base class that all other actions inherit. 
        By convention actions have 3 main phases: init, Setup, and DoAnalysis.
        Init is used to initialize the action, make sure that all arguments
        for the action are correct, and add any DataSets/DataFiles which will
        be used by the action. Setup will set up the action for a specific
        Topology file. DoAnalysis will perform the action on a given frame.
        A fourth function, Print, is for any additional calculations or output 
        the action may require once all frames are processed.

    pytraj doc:
    =============
    Add new action: add to pytraj/actions/ folder 
                    then update action in pytraj/actions/allactions
                    (TODO : allactions.py might be changed)
    """
    def __cinit__(self):
        # don't directly create instance of this ABC class.
        pass
        #self.baseptr = new _Analysis()

    def __dealloc__(self):
        # should I del pointer here or in subclass? 
        #del self.baseptr
        pass

    def __call__(self, *args, **kwd):
        return self._master(*args, **kwd)

    @makesureABC("Analysis")
    def read_input(self, command='', 
                   top=TopologyList(),
                   DataSetList dslist=DataSetList(), 
                   DataFileList dflist=DataFileList()):
        """
        Parameters
        ----------
        command : str
            Type of actions, mask, ... (Get help: Analysis_Clustering().help())
        top : Topology or TopologyList instance, default=TopologyList()
        dslist : DataSetList instance, default=DataSetList()
        dflist : DataFileList instance, default=DataFileList()
        debug : int, default=0
            debug option from cpptraj. (Do we need this?)
        """
        cdef ArgList arglist
        cdef TopologyList toplist
        cdef debug = 0
        cdef RetType STATUS

        if isinstance(top, Topology) or isinstance(top, string_types):
            toplist = TopologyList()
            toplist.add_parm(top)
        elif isinstance(top, TopologyList):
            toplist = <TopologyList> top
        else:
            raise ValueError('Topology must not be empty')

        if isinstance(command, string_types):
            arglist = ArgList(command)
        elif isinstance(command, ArgList):
            arglist = <ArgList> command
        else:
            raise ValueError()

        STATUS = self.baseptr.Setup(arglist.thisptr[0], 
                       dslist.thisptr, 
                       toplist.thisptr,
                       dflist.thisptr,
                       debug)
        if STATUS == ERRANALYSIS:
            raise ValueError()
        else:
            return STATUS

    @makesureABC("Analysis")
    def do_analysis(self):
        """
        Parameters: None
        """
        self.baseptr.Analyze()

    def _master(self, *args, **kwd):
        self.read_input(*args, **kwd)
        return self.do_analysis()
