from __future__ import absolute_import
from pytraj.utils import has_, require
from pytraj.FrameArray import FrameArray
from pytraj.exceptions import PytrajRequireObject
from ._load_pseudo_parm import load_pseudo_parm

def load_MDAnalysis(its_obj):
    """load MDAnalysis' Universe object to pytra's traj object"""
    if not has_("MDAnalysis"):
        require("MDAnalysis")
    else:
        from MDAnalysis import Universe 
        if not isinstance(its_obj, Universe):
            raise PytrajRequireObject("Universe")

        # creat pseudotop
        pseudotop = load_pseudo_parm(its_obj)

        # creat atom group
        ag = its_obj.atoms

        # load xyz coords
        farray = FrameArray()
        farray.top = pseudotop

        for _ in its_obj.trajectory:
            farray.load_xyz(ag.positions.flatten())
        return farray
