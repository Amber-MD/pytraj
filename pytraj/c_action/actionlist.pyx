# distutils: language = c++
from cython.operator cimport dereference as deref

from ..externals.six import string_types
from .c_action import ActionDict
from ..shared_methods import iterframe_master


def _get_arglist(arg):
    if isinstance(arg, ArgList):
        return arg
    else:
        return ArgList(arg)


def pipe(traj, commands, DatasetList dslist=DatasetList(), frame_indices=None):
    '''create frame iterator from cpptraj's commands.

    This method is useful if you want cpptraj pre-processing your Trajectory before
    throwing it to your own method.

    Parameters
    ----------
    commands : a list of strings of cpptraj's Action commands
    traj : Trajectory or any iterable that produces Frame
    dslist : CpptrajDatasetList, optional

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> for frame in pt.pipe(traj, ['autoimage', 'rms', 'center :1']): pass

    Above example is similiar to cpptraj's command::

         cpptraj -i EOF<<
         parm tz2.ortho.parm
         trajin tz2.ortho.nc
         autoimage
         rms
         center :1
         EOF

    You can design your own method::

        def new_method(traj, ...):
            for frame in traj:
                do_some_thing_fun_with(frame)

        fi = pt.pipe(traj, ['autoimage', 'rms', 'center :1'])

        # perform action with pre-processed frames (already autoimaged, then rms fit to
        # 1st frame, then center at box center.
        data = new_method(fi, ...)
    '''
    cdef Frame frame
    cdef ActionList actlist

    if frame_indices is None:
        fi = traj
    else:
        fi = traj.iterframe(frame_indices=frame_indices)

    if isinstance(commands, (list, tuple)):
        commands = commands
    elif isinstance(commands, string_types):
        commands = [line.lstrip().rstrip()
                    for line in commands.split('\n') if line.strip() != '']

    actlist = ActionList(commands, top=traj.top, dslist=dslist)
    for frame in iterframe_master(fi):
        actlist.compute(frame)
        yield frame


def compute(lines, traj, *args, **kwd):
    """perorm a series of cpptraj's actions on trajectory

    Parameters
    ----------
    lines : {str, list of string}
    traj : Trajectory-like
    process : {None, int}, default None
        if not None, there will have progess bar if user is using jupyter notebook.
        The value of ``progess`` shows how often the bar is displayed. Make sure to choose large ``progess`` number to avoid slowing down
        your calculation.
    *args, **kwd : more arguments

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2()
    >>> # cpptraj command style
    >>> data = pt.compute('''
    ...              rms
    ...              radgyr
    ...              molsurf
    ...              ''', traj)

    >>> # a list of commands
    >>> data = pt.compute([
    ...              'rms',
    ...              'radgyr',
    ...              'molsurf',], traj)

    >>> # a list of commands with reference
    >>> # need to explicitly add ``reference`` keyword and specify `ref=`
    >>> # calculate rms and use traj[3] as reference
    >>> data = pt.compute([
    ...              'rms myrms reference @CA',
    ...              'radgyr myrg @CA nomax',
    ...              'molsurf',], traj, ref=traj[3])
    >>> data['myrms']
    array([ 3.46476336,  3.4343108 ,  2.94523273, ...,  4.21848857,
            4.4566457 ,  3.94477017])
    >>> data['myrms'][3]
    0.0

    >>> # a list of commands with reference
    >>> # can also specify 'refindex'
    >>> # calculate rms and use traj[3] as reference
    >>> data = pt.compute([
    ...              'rms myrms refindex 0 @CA',
    ...              'radgyr myrg @CA nomax',
    ...              'molsurf',], traj, ref=traj[3])
    >>> data['myrms']
    array([ 3.46476336,  3.4343108 ,  2.94523273, ...,  4.21848857,
            4.4566457 ,  3.94477017])
    >>> data['myrms'][3]
    0.0
    """
    cdef DatasetList dslist

    # frequency to make the bar
    # None or an int
    freq =  kwd.get("progress")
    color =  kwd.get("color")

    if color is None:
        color = '#0080FF'
    else:
        kwd.pop('color')

    if isinstance(lines, (list, tuple, string_types)):
        ref = kwd.get('ref')
        if ref is not None:
            if isinstance(ref, Frame):
                reflist = [ref, ]
            else:
                # list/tuplex
                reflist = ref
        else:
            reflist = []

        dslist = DatasetList()

        if reflist:
            for ref_ in reflist:
                ref_dset = dslist.add_new('reference')
                ref_dset.top = traj.top
                ref_dset.add_frame(ref_)

        # create Frame generator
        fi = pipe(traj, commands=lines, dslist=dslist)

        # just iterate Frame to trigger calculation.
        # this code is for fun.
        if freq is not None:
           from pytraj.utils.progress import make_bar, init_display
           if hasattr(fi, 'n_frames'):
               max_frames = fi.n_frames
           elif hasattr(traj, 'n_frames'):
               max_frames = traj.n_frames
           else:
               # inaccurate max_frames
               max_frames = 1000000
           init_display(color)

        for idx, _ in enumerate(fi):
            if freq is not None:
                if idx % freq == 0:
                    make_bar(idx, max_frames)
                if idx == max_frames - 1:
                    make_bar(max_frames, max_frames)
            pass

        # remove ref
        return dslist[len(reflist):].to_dict()

    elif callable(lines):
        return lines(traj, *args, **kwd)


cdef class ActionList:
    def __cinit__(self):
        self.thisptr = new _ActionList()
        self.is_setup = False
        self.n_frames = 0

    property data:
        '''Store data (CpptrajDatasetList). This is for internal use.
        '''
        def __get__(self):
            return self._dslist

    def __init__(self, commands=None, Topology top=None,
                 DatasetList dslist=DatasetList(), DataFileList dflist=DataFileList(),
                 crdinfo={}):
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
        >>> from pytraj.datasets import CpptrajDatasetList
        >>> dslist = CpptrajDatasetList()
        >>> traj = pt.datafiles.load_tz2_ortho()
        >>> list_of_commands = ['autoimage',
        ...                     'rmsd first @CA',
        ...                     'hbond :3,8,10']
        >>> alist = ActionList(list_of_commands, traj.top, dslist=dslist)
        >>> for frame in traj:
        ...     alist.compute(frame)
        """
        self._dslist = dslist
        self._dflist = dflist
        self._crdinfo = crdinfo
        self.top = top

        if commands is not None and top is not None:
            for command in commands:
                command = command.rstrip().lstrip()
                try:
                    action, cm = command.split(" ", 1)
                except ValueError:
                    action = command.split(" ", 1)[0]
                    cm = ''
                action = action.rstrip().lstrip()
                self.add(action, command=cm,
                         top=top, dslist=dslist, dflist=dflist)

    def __dealloc__(self):
        if self.thisptr:
            del self.thisptr

    def add(self, action="",
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

        Examples
        --------
        >>> actlist = ActionList()
        >>> actlist.add('radgyr', '@CA', top=traj.top, dslist=dslist) # doctest: +SKIP
        """
        cdef Action action_
        cdef int status
        cdef _ActionInit actioninit_
        actioninit_ = _ActionInit(dslist.thisptr[0], dflist.thisptr[0])

        if isinstance(action, string_types):
            # create action object from string
            action_ = ActionDict()[action]
        else:
            action_ = action

        cdef ArgList _arglist

        self.top = top if top is not None else self.top

        _arglist = _get_arglist(command)
        # let ActionList free memory
        action_.own_memory = False
        status = self.thisptr.AddAction(action_.baseptr, _arglist.thisptr[0],
                                        actioninit_)

        if check_status:
            # return "0" if sucess "1" if failed
            return status
        else:
            return None

    def setup(self, Topology top, crdinfo={}, n_frames_t=0, bint exit_on_error=True):
        '''perform Topology checking and some stuff
        '''
        # let cpptraj free mem
        cdef _ActionSetup actionsetup_
        cdef CoordinateInfo crdinfo_
        cdef Box box
        cdef bint has_velocity, has_time, has_force

        if not crdinfo:
            crdinfo = self._crdinfo
        else:
            crdinfo = crdinfo

        crdinfo2 = dict((k, v) for k, v in crdinfo.items())
        if 'box' not in crdinfo2:
            crdinfo2['box'] = top.box

        crdinfo_ = CoordinateInfo(crdinfo2)

        actionsetup_ = _ActionSetup(top.thisptr, crdinfo_.thisptr[0], n_frames_t)
        self.thisptr.SetupActions(actionsetup_, exit_on_error)

    def compute(self, traj=Frame()):
        '''perform a series of Actions on Frame or Trajectory
        '''
        cdef _ActionFrame actionframe_
        cdef Frame frame

        if not self.is_setup:
            self.setup(self.top)
            # make sure to make is_setup True after processing
            # if not, pytraj will try to setup for every Frame
            self.is_setup = True

        if isinstance(traj, Frame):
            frame = <Frame> traj
            actionframe_ = _ActionFrame(frame.thisptr, self.n_frames)
            self.thisptr.DoActions(self.n_frames, actionframe_)
            self.n_frames += 1
        else:
            for frame in iterframe_master(traj):
                self.compute(frame)

    def post_process(self):
        self.thisptr.PrintActions()
