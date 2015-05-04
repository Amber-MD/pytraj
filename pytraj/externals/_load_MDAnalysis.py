from __future__ import absolute_import
from pytraj.utils import has_, require, _import_numpy
from pytraj.FrameArray import FrameArray
from pytraj.exceptions import PytrajRequireObject
from ._load_pseudo_parm import load_pseudo_parm
from ..Frame import Frame

# MDAnalysis needs numpy. So we always have numpy when using this
_, np = _import_numpy()

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

        farray = FrameArray()
        farray.top = pseudotop
        for _ in its_obj.trajectory:
            frame = Frame(farray.top.n_atoms)
            # set box for each Frame
            frame.boxview[:] = farray.top.box[:]
            # load xyz coords
            frame.buffer2d[:] = ag.positions
            farray.append(frame)
        return farray
