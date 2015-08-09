from pytraj.utils import has_
from pytraj.warnings import PytrajWarningMissing
from ..Trajectory import Trajectory
from ..Frame import Frame
from ..Topology import Topology
from ..utils.context import goto_temp_folder
from .six import string_types


def load_ParmEd(parmed_obj, as_traj=False, **kwd):
    """return pytraj's Topology or Trajectory objects

    Parameters
    ---------
    parmed_obj : ParmEd's Structure object
    save_and_reload: bool, default True
        if True, save `parmed_obj` to mol2 file and reload
        if False, internal convert. Might have bug
    as_traj: bool, default False
        if True, return pytraj.api.Trajectory
        if False, return Topology

    >>> import parmed as pmd
    >>> import pytraj as pt
    >>> p = pmd.download_PDB("1l2y")
    >>> top = pt.load_ParmEd(p, save_and_reload=True) 
    """
    import parmed as pmd
    from parmed.amber import AmberParm

    if isinstance(parmed_obj, string_types):
        import parmed as pmd
        # reserve **kwd for `structure=True`
        parmed_obj = pmd.load_file(parmed_obj, **kwd)
    # faster
    with goto_temp_folder():
        if isinstance(parmed_obj, AmberParm):
            fname = 'tmp.parm7'
        else:
            fname = 'tmppdb.pdb'
        parmed_obj.save(fname)
        top = Topology(fname)
    if as_traj:
        from pytraj import Trajectory
        coords = parmed_obj.coordinates
        coords = coords.reshape(1, *coords.shape)
        return Trajectory(xyz=coords, top=top)
    else:
        return top

def _load_parmed(parm_name):
    has_parmed = has_("parmed")
    if has_parmed:
        from parmed import load_file
        return load_file(parm_name)
    else:
        if not has_parmed:
            PytrajWarningMissing("`parmed`")
        return None


def to_ParmEd(pytraj_top):
    # TODO: exten to gromacs, charmm too
    # need to change extension
    """convert to ParmEd object"""
    from pytraj.utils.context import goto_temp_folder
    from pytraj.parms.ParmFile import ParmFile
    import parmed as chem

    # I am not a fan of saving/loading again but this might be best choice
    with goto_temp_folder():
        fname = "tmp_pytrajtop.prmtop"
        ParmFile().writeparm(pytraj_top, fname, format="")
        return chem.load_file(fname)
