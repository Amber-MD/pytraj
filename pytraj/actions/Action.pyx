# distutils: language = c++
from __future__ import print_function
from pytraj.decorators import makesureABC
from pytraj.externals.six import string_types
from pytraj.utils import is_generator
from pytraj.utils.check_and_assert import is_pytraj_trajectory
from pytraj.datasets.cast_dataset import cast_dataset
from pytraj._shared_methods import _frame_iter_master


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
        return self._master(*args, **kwd)

    @makesureABC("Action")
    def read_input(self, command='', 
                   top=TopologyList(),
                   DataSetList dslist=DataSetList(), 
                   DataFileList dflist=DataFileList(), 
                   int debug=0):
        """
        Parameters
        ----------
        command : str
            Type of actions, mask, ... (Get help: Action_Box().help())
        top : Topology or TopologyList instance, default=TopologyList()
        #flist : FrameList instance, default=FrameList()
        dslist : DataSetList instance, default=DataSetList()
        dflist : DataFileList instance, default=DataFileList()
        debug : int, default=0
            debug option from cpptraj. (Do we need this?)
        """
        cdef ArgList arglist
        cdef TopologyList toplist
        cdef RetType i_fail

        if isinstance(top, Topology):
            toplist = TopologyList()
            toplist.add_parm(top)
        elif isinstance(top, TopologyList):
            toplist = <TopologyList> top

        if isinstance(command, string_types):
            #command = command.encode("UTF-8")
            arglist = ArgList(command)
        elif isinstance(command, ArgList):
            arglist = <ArgList> command

        i_fail = self.baseptr.Init(arglist.thisptr[0], toplist.thisptr, 
                       dslist.thisptr, dflist.thisptr,
                       debug)

        if i_fail != OK:
            # check before do_action to avoid segfault
            raise ValueError("")
        else:
            return i_fail

    @makesureABC("Action")
    def process(self, Topology top=Topology(), Topology new_top=Topology()): 
        """
        Process input and do initial setup
        (TODO : add more doc)

        Parameters:
        ----------
        top : Topology instance, default (no default)
        new_top : new Topology instance, default=Topology()
            Need to provide this instance if you want to change topology
        """
        if "Strip" in self.__class__.__name__:
            # since `Action_Strip` will copy a modified version of `top` and 
            # store in new_top, then __dealloc__ (from cpptraj)
            # we need to see py_free_mem to False
            new_top.py_free_mem = False
        return self.baseptr.Setup(top.thisptr, &(new_top.thisptr))

    @makesureABC("Action")
    def do_action(self, current_frame=None, Frame new_frame=Frame(), 
            update_mass=True, Topology top=Topology()):
        """
        Perform action on Frame. 
        Parameters:
        ----------
        current_frame : Frame instance need to be processed, default=Frame() 

        new_frame : Frame instance, defaul=Frame()
            if action change Frame, you need to have this
        """
        # debug
        cdef Frame frame
        cdef int i
        cdef object traj, tmptraj, farray

        if self.__class__.__name__ == 'Action_Strip':
            # let cpptraj do its job for this special action
            new_frame.py_free_mem = False

        if isinstance(current_frame, Frame):
            frame = <Frame> current_frame
            # make sure to update frame mass
            if update_mass and not top.is_empty():
                frame.set_frame_mass(top)
            self.baseptr.DoAction(self.n_frames, frame.thisptr, &(new_frame.thisptr))
            self.n_frames += 1
        else:
            for frame in _frame_iter_master(current_frame):
                self.do_action(frame, new_frame)

    @makesureABC("Action")
    def print_output(self):
        """Do we need this?"""
        self.baseptr.Print()

    # Do we really need this method?
    @classmethod
    def _get_action_from_functptr(cls, FunctPtr funct):
        cdef Action act = Action()
        if funct.ptr() == NULL:
            raise ValueError("NULL pointer")
        act.baseptr = <_Action*> funct.ptr()
        return act

    def _master(self, command='',
                  current_frame=Frame(),
                  top=Topology(),
                  dslist=DataSetList(), 
                  dflist=DataFileList(), 
                  new_top=Topology(),
                  new_frame=Frame(),
                  update_mass=True,
                  int debug=0,
                  update_frame=False,
                  quick_get=False):
        """
        TODO : (do we need this method?)
            + add doc
            + don't work with `chunk_iter`

        """
        if top.is_empty():
            _top = current_frame.top
        else:
            _top = top
        self.read_input(command=command, 
                        top=_top, 
                        dslist=dslist,
                        dflist=dflist, debug=debug)

        self.process(top=_top, new_top=new_top)
        self.do_action(current_frame, new_frame, update_mass=update_mass, top=_top)

        # currently support only dtype = 'DOUBLE', 'MATRIX_DBL', 'STRING', 'FLOAT', 'INTEGER'
        # we get the last dataset from dslist
        # (if we call self.run() several times, the result will be dumped to dslist)
        # FIXME: add all dtype in cpptraj so we don't need to specify them
        if quick_get:
            idx = dslist.size - 1
            if hasattr(dslist[idx], 'dtype'):
                dtype = dslist[idx].dtype.upper()
                if dtype in ['DOUBLE', 'MATRIX_DBL', 'STRING', 'FLOAT', 'INTEGER',
                             'MATRIX_FLT', 'VECTOR']:
                    d0 = cast_dataset(dslist[idx], dtype=dtype)
                    return d0
                else:
                    # return what?
                    return None
            else:
                raise RuntimeError("don't know how to cast dataset")
        else:
            return dslist

    def reset_counter(self):
        self.n_frames = 0
