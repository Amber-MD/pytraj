# distutils: language = c++
from pytraj.decorators import makesureABC
from pytraj.externals.six import string_types


cdef class Action:
    """
    Original cpptraj doc:
    ====================
        The abstract base class that all other actions inherit. 
        By convention actions have 3 main phases: init, Setup, and DoAction.
        Init is used to initialize the action, make sure that all arguments
        for the action are correct, and add any DataSets/DataFiles which will
        be used by the action. Setup will set up the action for a specific
        Topology file. DoAction will perform the action on a given frame.
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
        #self.baseptr = new _Action()

    def __dealloc__(self):
        # should I del pointer here or in subclass? 
        #del self.baseptr
        pass

    @makesureABC("Action")
    def read_input(self, command='', current_top=TopologyList(),
                   FrameList flist=FrameList(), 
                   DataSetList dslist=DataSetList(), 
                   DataFileList dflist=DataFileList(), 
                   int debug=0):
        """
        Parameters
        ----------
        command : str
            Type of actions, mask, ... (Get help: Action_Box().help())
        current_top : Topology or TopologyList instance, default=TopologyList()
        flist : FrameList instance, default=FrameList()
        dslist : DataSetList instance, default=DataSetList()
        dflist : DataFileList instance, default=DataFileList()
        debug : int, default=0
            debug option from cpptraj. (Do we need this?)
        """
        cdef ArgList arglist
        cdef TopologyList toplist

        if isinstance(current_top, Topology):
            toplist = TopologyList()
            toplist.add_parm(current_top)
        elif isinstance(current_top, TopologyList):
            toplist = <TopologyList> current_top

        if isinstance(command, string_types):
            #command = command.encode("UTF-8")
            arglist = ArgList(command)
        elif isinstance(command, ArgList):
            arglist = <ArgList> command

        return self.baseptr.Init(arglist.thisptr[0], toplist.thisptr, 
                       flist.thisptr, dslist.thisptr, dflist.thisptr,
                       debug)

    @makesureABC("Action")
    def process(self, Topology current_top=Topology(), Topology new_top=Topology()): 
        """
        Process input and do initial setup
        (TODO : add more doc)

        Parameters:
        ----------
        current_top : Topology instance, default (no default)
        new_top : new Topology instance, default=Topology()
            Need to provide this instance if you want to change topology
        """
        if "Strip" in self.__class__.__name__:
            # since `Action_Strip` will copy a modified version of `current_top` and 
            # store in new_top, then __dealloc__ (from cpptraj)
            # we need to see py_free_mem to False
            new_top.py_free_mem = False
        return self.baseptr.Setup(current_top.thisptr, &(new_top.thisptr))

    @makesureABC("Action")
    def do_action(self, int idx=0, current_frame=Frame(), Frame new_frame=Frame()):
        """
        Perform action on Frame. Depend on what action you want to perform, you might get
        new_frame or get data from dslist or dflist...
        TODO : add FrameArray
        Parameters:
        ----------
        idx : int, defaul=0 
            id of Frame
        current_frame : Frame instance need to be processed, default=Frame() 
        new_frame : Frame instance, defaul=Frame()
            if action change Frame, you need to have this
        >>> from pytraj._cast import cast_dataset
        >>> # dslist is DataSetList (list of DataSet instance)
        >>> d0 = cast_dataset(dslist[0], dtype='DataSet_double')
        """
        # debug
        cdef Frame frame
        new_frame.py_free_mem = False
        #current_frame.py_free_mem = False
        # got double-free memory error when not using above flag
        # end debug
        #return self.baseptr.DoAction(idx, current_frame.thisptr, &(new_frame.thisptr))
        if current_frame.is_empty():
            raise ValueError("require providing Frame or FrameArray")
        if isinstance(current_frame, Frame):
            frame = <Frame> current_frame
            frame.py_free_mem = False
            self.baseptr.DoAction(idx, frame.thisptr, &(new_frame.thisptr))
        else:
            # add check
            # assume Traj instance
            # TODO : check new_frame
            for frame in current_frame:
                frame.py_free_mem = False
                self.baseptr.DoAction(idx, frame.thisptr, &(new_frame.thisptr))

    @makesureABC("Action")
    def print_output(self):
        """Do we need this?"""
        self.baseptr.Print()

    # Do we really need this method?
    @classmethod
    def get_action_from_functptr(cls, FunctPtr funct):
        cdef Action act = Action()
        if funct.ptr() == NULL:
            raise ValueError("NULL pointer")
        act.baseptr = <_Action*> funct.ptr()
        return act

    def master(self, command='', int idx=0,
                   current_top=TopologyList(),current_frame=Frame(),
                   FrameList flist=FrameList(), 
                   DataSetList dslist=DataSetList(), 
                   DataFileList dflist=DataFileList(), 
                   new_top=Topology(),
                   new_frame=Frame(),
                   int debug=0):
        """combined all 3 steps
                master(command='', int idx=0,
                       current_top=TopologyList(),current_frame=Frame(),
                       FrameList flist=FrameList(), 
                       DataSetList dslist=DataSetList(), 
                       DataFileList dflist=DataFileList(), 
                       new_top=Topology(),
                       new_frame=Frame(),
                       int debug=0):
        """
        self.read_input(command=command, current_top=current_top, 
                        flist=flist, dslist=dslist,
                        dflist=dflist, debug=debug)
        self.process(current_top=current_top, new_top=new_top)
        self.do_action(idx=idx, current_frame=current_frame, new_frame=new_frame)
