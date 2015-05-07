"""This is a thin wrapper of Trajin_Single
We need to sub-class Trajin_Single to use Trajectory
(we called Trajin_Single from Trajectory, so we can not call Trajectory back from 
Trajin_Single)
"""
from __future__ import absolute_import
from pytraj.trajs.Trajin_Single import Trajin_Single
from pytraj._action_in_traj import ActionInTraj
from pytraj.action_dict import ActionDict
from pytraj.Frame import Frame
from pytraj.AtomMask import AtomMask
from pytraj.externals.six import string_types
from pytraj.exceptions import PytrajMemviewError


class TrajectoryIterator(Trajin_Single, ActionInTraj):
    def __init__(self, *args, **kwd):
        pass

    @property
    def topology(self):
        """traditional name for Topology file"""
        return self.top

    @topology.setter
    def topology(self, newtop):
        self.top = newtop

    def frame_iter(self, start=0, stop=-1, stride=1, mask=None, autoimage=False):
        if autoimage:
            act = ActionDict()['autoimage']
        for frame in super(TrajectoryIterator, self).frame_iter(start, stop, stride):
            if autoimage:
                act(current_frame=frame, top=self.top)
            if mask is not None:
                if isinstance(mask, string_types):
                    atm = self.top(mask)
                else:
                    try:
                        atm = AtomMask()
                        atm.add_selected_indices(mask)
                    except TypeError:
                        raise PytrajMemviewError()
                frame2 = Frame(atm.n_atoms)
                frame2.set_coords(frame, atm)
                yield frame2
            else:
                yield frame
