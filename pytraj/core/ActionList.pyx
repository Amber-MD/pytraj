# distutils: language = c++
from cython.operator cimport dereference as deref

from ..externals.six import string_types
from ..action_dict import ActionDict
from .._shared_methods import iterframe_master

def _get_arglist(arg):
    if isinstance(arg, ArgList):
        return arg
    else:
        return ArgList(arg)

cdef class ActionList:
    def __cinit__(self):
        self.thisptr = new _ActionList()
        self.top_is_processed = False

    property data:
        def __get__(self):
            return self._dslist

    def __init__(self, commands=None, Topology top=None,
                 DatasetList dslist=DatasetList(), DataFileList dflist=DataFileList()):
        """not done yet

        Parameters
        ----------
        actionlist : a list of tuple
        top : Topology
        dslist : DatasetList, optional
            hold data for actions
        dflist : DataFileList, optional
            hold datafiles

        Examples
        --------
        >>> import pytraj as pt
        >>> from pytraj import ActionList
        >>> list_of_commands = ['autoimage',
                                'rmsd first @CA',
                                'hbond :3,8,10']
        >>> alist = ActionList(list_of_commands, traj.top, dslist=dslist)
        >>> for frame in traj:
        >>>     alist.do_actions(frame)
        """
        self._dslist = dslist
        self._dflist = dflist

        if commands is not None and top is not None:
            for command in commands:
                command = command.rstrip().lstrip()
                try:
                    action, cm = command.split(" ", 1)
                except ValueError:
                    action = command.split(" ", 1)[0]
                    cm = ''
                action = action.rstrip().lstrip()
                self.add_action(action, command=cm,
                                top=top, dslist=dslist, dflist=dflist)

    def __dealloc__(self):
        if self.thisptr:
            del self.thisptr

    def clear(self):
        self.thisptr.Clear()

    def add_action(self, action="", 
                         command="", 
                         top=None, 
                         DatasetList dslist=DatasetList(), 
                         DataFileList dflist=DataFileList(),
                         check_status=False):
        """Add action to ActionList

        Parameters
        ----------
        action : str or Action object
        command : str or ArgList object
        top : str | Topology | TopologyList
        dslist : DatasetList 
        dflist : DataFileList
        check_status : bool, default=False
            return status of Action (0 or 1) if "True"
        """
        cdef object _action
        cdef int status
        cdef _ActionInit actioninit_
        actioninit_ = _ActionInit(dslist.thisptr[0], dflist.thisptr[0])

        if isinstance(action, string_types):
            # create action object from string
            _action = ActionDict()[action]
        else:
            _action = action

        cdef FunctPtr func = <FunctPtr> _action.alloc()
        cdef ArgList _arglist

        self.top = top
        # add function pointer: How?

        _arglist = _get_arglist(command)
        status = self.thisptr.AddAction(func.ptr, _arglist.thisptr[0], 
                                        actioninit_)

        if check_status:
            # return "0" if sucess "1" if failed
            return status
        else:
            return None

    def process(self, Topology top, crdinfo={}, n_frames_t=0):
        # let cpptraj free mem
        cdef _ActionSetup actionsetup_
        cdef CoordinateInfo crdinfo_
        cdef Box box
        cdef bint has_velocity, has_time, has_force

        box = crdinfo.get('box', top.box)
        has_velocity = crdinfo.get('has_velocity', False)
        has_time = crdinfo.get('has_time', False)
        has_force = crdinfo.get('has_force', False)

        crdinfo_ = CoordinateInfo(box.thisptr[0], has_velocity, has_time, has_force)

        #top._own_memory = False

        actionsetup_ = _ActionSetup(top.thisptr, crdinfo_, n_frames_t)
        self.thisptr.SetupActions(actionsetup_)

    def do_actions(self, traj=Frame(), int idx=0, use_mass=True):
        cdef _ActionFrame actionframe_
        cdef Frame frame
        cdef int i

        if not self.top_is_processed:
            self.process(self.top)

        if isinstance(traj, Frame):
            frame = <Frame> traj
            if use_mass:
                frame.set_frame_mass(self.top)
            actionframe_ = _ActionFrame(frame.thisptr)
            self.thisptr.DoActions(idx, actionframe_)
        else:
            for i, frame in enumerate( iterframe_master(traj)):
                self.do_actions(frame, i) 

    def is_empty(self):
        return self.thisptr.Empty()

    @property
    def n_actions(self):
        return self.thisptr.Naction()
