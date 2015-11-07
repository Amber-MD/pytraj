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

    from MDAnalysis import Universe
    from pytraj.Trajectory import Trajectory
    from ..Frame import Frame

    # don't import here since we import load_pseudo_parm in
    # TrajectoryMDAnalysisIterator
    #from ..trajs.TrajectoryMDAnalysisIterator import TrajectoryMDAnalysisIterator

    # MDAnalysis needs numpy. So we always have numpy when using this
    if not isinstance(its_obj, Universe):
        raise ValueError("must be a Universe")

    # creat pseudotop
    if top is None:
        raise ValueError("need a Topology or pdb/mol2/...")
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
