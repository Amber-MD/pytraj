# distutils: language = c++
from cython.operator cimport dereference as deref

# TODO : double-check C++ code
from pytraj import TrajinList
from pytraj.externals.six import string_types
from pytraj.action_dict import ActionDict

cdef class ActionList:
    def __cinit__(self):
        self.thisptr = new _ActionList()

    def __dealloc__(self):
        if self.thisptr:
            del self.thisptr

    def clear(self):
        self.thisptr.Clear()

    def add_action(self, action="", 
                         command="", 
                         top=None, 
                         DataSetList dlist=DataSetList(), 
                         DataFileList dflist=DataFileList()):
        """
        Add action to ActionList

        Parameters:
        ==========
        actionobj :: Action object
        arglist :: ArgList instance
        toplist :: TopologyList instance
        #flist :: FrameList instance
        dlist :: DataSetList 
        dflist :: DataFileList
        """
        cdef object _action

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
        # add function pointer: How?

        if isinstance(command, ArgList):
            _arglist = command
        else:
            # try creating arglist
            _arglist = ArgList(command)

        return self.thisptr.AddAction(func.ptr, _arglist.thisptr[0], 
                                      toplist.thisptr,
                                      dlist.thisptr, dflist.thisptr)

    def process(self, Topology top):
        # let cpptraj free mem
        top.py_free_mem = False
        return self.thisptr.SetupActions(&(top.thisptr))

    def do_actions(self, traj=Frame(), int idx=0):
        # TODO : read cpptraj code to check memory stuff
        # set py_free_mem = False to let cpptraj does its job

        cdef Frame frame
        cdef int i

        if len(traj) == 0:
            raise ValueError("empty Frame/Traj/List, what can I do with this?")
        if isinstance(traj, Frame):
            frame = <Frame> traj
            frame.py_free_mem = False
            self.thisptr.DoActions(&(frame.thisptr), idx)
        elif hasattr(traj, 'n_frames'):
            for i, frame in enumerate(traj):
                self.do_actions(frame, i)
        elif isinstance(traj, (list, tuple, TrajinList)):
            for tmtraj in traj:
                self.do_actions(tmtraj)

    def listinfo(self):
        self.thisptr.List()

    def is_empty(self):
        return self.thisptr.Empty()

    @property
    def n_actions(self):
        return self.thisptr.Naction()

    def cmd_string(self, int i):
        return self.thisptr.CmdString(i)

    def action_alloc(self, int i):
        # TODO : do we need to expose this method here?
        # return func_ptr
        cdef FunctPtr func = FunctPtr()
        if i >= self.n_actions:
            raise IndexError("index must be < " + str(self.n_actions)) 
        func.ptr = self.thisptr.ActionAlloc(i)
        return func
