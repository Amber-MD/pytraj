"""this file has commonly used actions such as rmsd calculation, 
randomizeions, strip atoms, ..."""

from __future__ import print_function, absolute_import
from pytraj.Topology import Topology
from .TopologyList import TopologyList
from .ArgList import ArgList
from pytraj.Frame import Frame
#from pytraj.Trajin_Single import Trajin_Single
from pytraj.FrameArray import FrameArray
from pytraj.actions import allactions
from pytraj import adict
from pytraj.DataSetList import DataSetList

# external
from pytraj.externals.six import string_types

__all__ = ['strip', 'fit', 'get_subframe', 'randomize_ions']

def strip(arg, mask):
    """
    TODO: can not modify oldtop, Python to pass it as value
    TODO: validate
    Modify top and farray

    Strip atoms
    Parameters
    ---------
    Topology instance, or Frame instance, or Trajin_Single instance

    Out:
    ----
    Return : None
    """
    if isinstance(arg, Topology):
        top = arg.copy()
    elif isinstance(arg, FrameArray):
        top = arg.top.copy()

    toplist = TopologyList()
    toplist.add_parm(top)

    stripact = allactions.Action_Strip()
    stripact.read_input(ArgList("strip " + mask), toplist)
    stripact.process(toplist[0], top)

    if isinstance(arg, FrameArray):
        for i in range(arg.size):
            tmp = arg[i]
            stripact.do_action(tmp, i)
            arg[i] = tmp
        # need to update arg.top
        arg.top = top.copy()

    if isinstance(arg, Topology):
        return top.copy()

def fit(frame, ref=None):
    """fit Frame intance to reference Frame

    Parameters
    ---------
    frame : Frame instance
    ref : reference Frame (default = None)
    """
    if not ref:
        raise ValueError("missing reference Frame")
    else:
        # TODO : fitting
        pass

def get_subframe(frame=None, mask=None, top=None, atommask=None):
    # TODO : move to `io.py`
    return frame.get_subframe(mask=mask, top=top, atommask=atommask)


def randomize_ions(frame=Frame(), top=Topology(), command=""):
    """randomize_ions for given Frame with Topology
    Return : None
    Parameters
    ---------
    frame : Frame instance, default=Frame()
        frame coords will be modified

    top : Topology instance, default=Topology()

    >>> from pytraj.misc import randomize_ions
    >>> randomize_ions(frame, top, command="randomizeions @Na+ around :1-16 by 5.0 overlap 3.0")
    """
    act = allactions.Action_RandomizeIons()
    act.master(command=command,
               current_top=top,
               current_frame=frame,
               )

def action_help(action=None):
    # where should we put this method? putting here seems not really reasonable
    from pytraj import adict

    if action is None:
        print ("give the name of Action to get help")
        print ("Example: action_help(): all action keywords")
        print ("Example: action_help('rmsd'): help for Action_Rmsd")
        print (adict.keys())
    else:
        adict[action].help()

def get_action_dict():
    actdict = {}
    for key in allactions.__dict__.keys():
        if "Action_" in key:
            act = key.split("Action_")[1]
            # add Action classes
            actdict[act] = allactions.__dict__["Action_" + act]
    return actdict

# add action_dict
action_dict = get_action_dict()

def show_code(func, get_txt=False):
    """show code of func or module"""
    import inspect
    txt = inspect.getsource(func)
    if not get_txt:
        print (txt)
    else:
        return txt

def calculate(action=None, command=None, traj=None, top=None, **kwd):
    # TODO : should write universal help's method
    """
    quick way to get data
    Parameters:
    action : Action object or str, default=None
    command : str, default=None
    traj : Trajectory object (FrameArray, TrajReadOnly, ...) or list, tuple of traj object
    top : topology

    Use `calculate(ahelp=True)` or `calculate(ahelp='action name')` for help

    """
    from pytraj import adict
    if action is None and command is None and traj is None and top is None:
        if not kwd:
            #
            #print (calculate.__doc__)
            print ()
            print (adict.keys())
            print ()
            print ("use calculate(key=action_name) for help")
        else:
            adict[kwd['key'].lower()].help()
    else:
        if top is None:
            try:
               top = traj.top
            except:
                # list, tuple of traj objects
                top = traj[0].top
        if traj is None:
            raise ValueError("must have trajectory object")
        if isinstance(action, string_types):
            # convert to action
            act = adict[action]
        else:
            act = action
        return act(command, traj, top, quick_get=True)

def to_string_ss(arr0):
    """
    arr0 : ndarray
    """
    ss = ['None', 'Para', 'Anti', '3-10', 'Alpha', 'Pi', 'Turn', 'Bend']
    len_ss = len(ss)
    ssdict = dict(zip(range(len_ss), ss))
    return map(lambda idx: ssdict[idx], arr0)

def calc_dssp(command="", traj=None, dtype='int'):
    dslist = DataSetList()
    adict['dssp'](command, 
                  current_frame=traj, current_top=traj.top, 
                  dslist=dslist)
    dtype = dtype.upper()
    arr0 = dslist.get_dataset(dtype="integer")
    if dtype in ['INT', 'INTERGER']:
        return arr0
    elif dtype in ['STRING', 'STR']:
        shape = arr0.shape
        tmplist = [[x for x in to_string_ss(arr)] for arr in arr0]
        return tmplist
    else:
        raise NotImplementedError("dtype = integer, int, string, str")

def simple_plot(d0, *args, **kwd):
    # TODO : return object so we can update axis, label, ..
    from pytraj import _import

    has_plot, plt = _import('matplotlib.pyplot')
    if not has_plot:
        raise RuntimeError("require matplotlib installed")
    fig = plt.pyplot.plot(range(d0.size), d0[:], *args, **kwd)
    plt.pyplot.show()


def frame_iter(self, start=0, stop=-1, stride=1):
    """iterately get Frames with start, stop, stride 
    Parameters
    ---------
    start : int (default = 0)
    chunk : int (default = 1)
    stop : int (default = max_frames - 1)
    """
    frame = Frame(self.top.n_atoms)
    if stop == -1 or stop >= self.n_frames:
        stop = self.n_frames - 1

    i = start
    # use `with self` in case needed to open/close file
    with self:
        while i <= stop:
            if hasattr(self, 'read_traj_frame'):
                # cpptraj Traj-like object
                self.read_traj_frame(i, frame)
            else:
                # FrameArray object
                frame = self[i]
            yield frame
            i += stride
