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
from pytraj import adict, analdict
from pytraj.DataSetList import DataSetList

from pytraj._utils import set_world_silent

# external
from pytraj.externals.six import string_types

def info(obj):
    """get `help` for obj
    Useful for Actions and Analyses
    
    Since we use `set_worl_silent` to turn-off cpptraj' stdout, we need 
    to turn on to use cpptraj's help methods
    """
    if isinstance(obj, string_types):
        if obj in adict.keys():
            # make Action object
            _obj = adict[obj]
        elif obj in analdict.keys():
            # make Analysis object
            _obj = analdict[obj]
        else:
            raise ValueError("keyword must be an Action or Analysis")
    else:
        # assume `obj` hasattr `help`
        _obj = obj

    if hasattr(_obj, 'help'):
        set_world_silent(False)
        _obj.help()
        set_world_silent(True)
    else:
        raise ValueError("object does not have `help` method")

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

def get_atts(obj):
    """get methods and atts from obj but excluding special methods __"""
    atts_dict = dir(obj)
    return [a for a in atts_dict if not a.startswith("__")]
