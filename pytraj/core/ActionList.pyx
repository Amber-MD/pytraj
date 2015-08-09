# distutils: language = c++
from cython.operator cimport dereference as deref

from . import TrajinList
from ..externals.six import string_types
from ..action_dict import ActionDict
from .._get_common_objects import _get_arglist
from .._shared_methods import _frame_iter_master

cdef class ActionList:
    def __cinit__(self):
        self.thisptr = new _ActionList()
        self.top_is_processed = False

    def __dealloc__(self):
        if self.thisptr:
            del self.thisptr

    def clear(self):
        self.thisptr.Clear()

    def add_action(self, action="", 
                         command="", 
                         top=None, 
                         DataSetList dslist=DataSetList(), 
                         DataFileList dflist=DataFileList(),
                         check_status=False):
        """
        Add action to ActionList

        Parameters:
        ==========
        action : str or Action object
        command : str or ArgList object
        top : str | Topology | TopologyList
        dslist : DataSetList 
        dflist : DataFileList
        check_status : bool, default=False
            return status of Action (0 or 1) if "True"
        """
        cdef object _action
        cdef int status

        if isinstance(action, string_types):
            # create action object from string
            _action = ActionDict()[action]
        else:
            _action = action

        cdef FunctPtr func = <FunctPtr> _action.alloc()
        cdef TopologyList toplist
        cdef ArgList _arglist

        if isinstance(top, Topology):
            toplist = TopologyList()
            toplist.add_parm(top)
        elif isinstance(top, TopologyList):
            toplist = top
        self.toplist = toplist
        # add function pointer: How?

        _arglist = _get_arglist(command)
        status = self.thisptr.AddAction(func.ptr, _arglist.thisptr[0], 
                                        toplist.thisptr,
                                        dslist.thisptr, dflist.thisptr)

        if check_status:
            # return "0" if sucess "1" if failed
            return status
        else:
            return None

    def process(self, Topology top):
        # let cpptraj free mem
        top.py_free_mem = False
        self.thisptr.SetupActions(&(top.thisptr))
        self.top_is_processed = True

    def do_actions(self, traj=Frame(), int idx=0, use_mass=True):
        cdef Frame frame
        cdef int i

        if not self.top_is_processed:
            self.process(self.toplist[0])

        if isinstance(traj, Frame):
            frame = <Frame> traj
            if use_mass:
                frame.set_frame_mass(self.toplist[0])
            self.thisptr.DoActions(&(frame.thisptr), idx)
        else:
            for i, frame in enumerate( _frame_iter_master(traj)):
                self.do_actions(frame, i) 

    def is_empty(self):
        return self.thisptr.Empty()

    @property
    def n_actions(self):
        return self.thisptr.Naction()

    def _action_alloc(self, int i):
        cdef FunctPtr func = FunctPtr()
        if i >= self.n_actions:
            raise IndexError("index must be < " + str(self.n_actions)) 
        func.ptr = self.thisptr.ActionAlloc(i)
        return func
