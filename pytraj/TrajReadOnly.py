"""This is a thin wrapper of Trajin_Single
We need to sub-class Trajin_Single to use FrameArray
(we called Trajin_Single from FrameArray, so we can not call FrameArray back from 
Trajin_Single)
"""
from __future__ import absolute_import
from pytraj.trajs.Trajin_Single import Trajin_Single
from pytraj._action_in_traj import ActionInTraj

class TrajReadOnly(Trajin_Single, ActionInTraj):
    def __init__(self, *args, **kwd):
        pass

    @property
    def topology(self):
        """traditional name for Topology file"""
        return self.top

    @topology.setter
    def topology(self, newtop):
        self.top = newtop
