# distutils: language = c++
from pytraj.decorators import makesureABC
from pytraj.externals.six import string_types
from pytraj.utils import is_generator

from pytraj.cast_dataset import cast_dataset
from pytraj import TrajinList
from pytraj.TrajReadOnly import TrajReadOnly


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
        self.n_frames = 0
        #self.baseptr = new _Action()

    def __dealloc__(self):
        # should I del pointer here or in subclass? 
        #del self.baseptr
        pass

    def __del__(self):
        del self.baseptr

    def __str__(self):
        txt = "< %s object >" % (self.__class__.__name__)
        return txt

    def __repr__(self):
        return self.__str__()

    def __call__(self, *args, **kwd):
        """
        >>> from pytraj import *
        >>> traj = io.load("../tz2.nc", "../tz2.parm7")
        >>> dslist = DataSetList.DataSetList()
        >>> adict['jcoupling']("outfile Jcoupling.dat kfile Karplus.txt", traj[0], traj.top, dslist=dslist)
        """
        return self.run(*args, **kwd)

    @makesureABC("Action")
    def read_input(self, command='', 
                   current_top=TopologyList(),
                   DataSetList dslist=DataSetList(), 
                   DataFileList dflist=DataFileList(), 
                   int debug=0):
        """
        Parameters
        ----------
        command : str
            Type of actions, mask, ... (Get help: Action_Box().help())
        current_top : Topology or TopologyList instance, default=TopologyList()
        #flist : FrameList instance, default=FrameList()
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
                       #flist.thisptr, dslist.thisptr, dflist.thisptr,
                       dslist.thisptr, dflist.thisptr,
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
    def do_action(self, current_frame=Frame(), Frame new_frame=Frame()):
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
        cdef int i
        cdef object traj, tmptraj

        new_frame.py_free_mem = False

        if isinstance(current_frame, Frame):
            frame = <Frame> current_frame
            frame.py_free_mem = False
            self.baseptr.DoAction(self.n_frames, frame.thisptr, &(new_frame.thisptr))
            self.n_frames += 1
        #elif hasattr(current_frame, 'n_frames'):
        elif isinstance(current_frame, (FrameArray, TrajReadOnly)):
            # Trajectory-like object
            traj = current_frame 
            for frame in traj:
                self.do_action(current_frame=frame, new_frame=new_frame)
        elif isinstance(current_frame, (list, tuple)):
            # creat alias to avoid con
            trajlist = current_frame
            # FIXME: correct `idx`
            # FIXME: ugly
            # make sure to check ActionList class to avoid duplication
            for tmptraj in trajlist:
                # recursive
                # why doesn't this work with chunk_iter?
                if hasattr(tmptraj, 'n_frames') or isinstance(tmptraj, Frame):
                    self.do_action(tmptraj, new_frame)
                else:
                    # chunk_iter
                    for tmptraj2 in tmptraj:
                        self.do_action(tmptraj2, new_frame)
        else:
            for frame in current_frame:
                self.do_action(current_frame, new_frame)

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

    def run(self, command='',
                  current_frame=Frame(),
                  current_top=Topology(),
                  dslist=DataSetList(), 
                  dflist=DataFileList(), 
                  new_top=Topology(),
                  new_frame=Frame(),
                  int debug=0,
                  update_frame=False,
                  quick_get=False):
        """
        TODO : (do we need this method?)
            + add doc
            + don't work with `chunk_iter`

        """
        if current_top.is_empty():
            _top = current_frame.top
        else:
            _top = current_top
        self.read_input(command=command, 
                        current_top=_top, 
                        dslist=dslist,
                        dflist=dflist, debug=debug)

        self.process(current_top=_top, new_top=new_top)
        self.do_action(current_frame, new_frame)

        # currently support only dtype = 'DOUBLE', 'MATRIX_DBL', 'STRING', 'FLOAT', 'INTEGER'
        # we get the last dataset from dslist
        # (if we call self.run() several times, the result will be dumped to dslist)
        # FIXME: add all dtype in cpptraj so we don't need to specify them
        if quick_get:
            idx = dslist.size - 1
            if hasattr(dslist[idx], 'dtype'):
                dtype = dslist[idx].dtype.upper()
                if dtype in ['DOUBLE', 'MATRIX_DBL', 'STRING', 'FLOAT', 'INTEGER']:
                    d0 = cast_dataset(dslist[idx], dtype=dtype)
                    return d0
                else:
                    # return what?
                    return None
            else:
                raise RuntimeError("don't know how to cast dataset")

    def reset_counter(self):
        self.new_frame = 0

    def master(self, *args, **kwd):
        """keep this method since some of examples uses them"""
        return self.run(*args, **kwd)
