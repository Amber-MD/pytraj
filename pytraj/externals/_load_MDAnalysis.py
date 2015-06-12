from __future__ import absolute_import

def load_MDAnalysis(its_obj, top=None):
    """load MDAnalysis' Universe object to pytra's traj object

    Notes
    -----
    All coords will be loaded

    See Also
    --------
    load_MDAnalysisIterator
    """
    from pytraj.utils import has_, require, _import_numpy
    from pytraj.Trajectory import Trajectory
    from pytraj.exceptions import PytrajRequireObject
    from ._load_pseudo_parm import load_pseudo_parm
    from ..Frame import Frame
    
    # don't import here since we import load_pseudo_parm in 
    # TrajectoryMDAnalysisIterator
    #from ..trajs.TrajectoryMDAnalysisIterator import TrajectoryMDAnalysisIterator
    
    # MDAnalysis needs numpy. So we always have numpy when using this
    _, np = _import_numpy()

    if not has_("MDAnalysis"):
        require("MDAnalysis")
    else:
        from MDAnalysis import Universe 
        if not isinstance(its_obj, Universe):
            raise PytrajRequireObject("Universe")

        # creat pseudotop
        if top is None:
            pseudotop = load_pseudo_parm(its_obj)
        else:
            pseudotop = top

        # creat atom group
        ag = its_obj.atoms

        farray = Trajectory()
        farray.top = pseudotop
        for _ in its_obj.trajectory:
            frame = Frame(farray.top.n_atoms)
            # set box for each Frame
            frame.boxview[:] = farray.top.box[:]
            # load xyz coords, let numpy do automatically casting
            frame.xyz[:] = ag.positions
            # we don't need to make copy=True since we already created
            # frame and `farray` can 'keep' it
            farray.append(frame, copy=False)
        return farray
