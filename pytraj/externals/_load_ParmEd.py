from pytraj.utils import has_
from pytraj.warnings import PytrajWarningMissing
from ._load_pseudo_parm import load_pseudo_parm
from ..Trajectory import Trajectory
from ..Frame import Frame
from ..Topology import Topology
from ..utils.context import goto_temp_folder


def load_ParmEd(parmed_obj, save_and_reload=True):
    """return pytraj's Topology or Trajectory objects

    Parameters
    ---------
    parmed_obj : ParmEd's Structure object
    save_and_reload: bool, default True
        if True, save `parmed_obj` to mol2 file and reload
        if False, internal convert. Might have bug

    >>> import parmed as pmd
    >>> import pytraj as pt
    >>> p = pmd.download_PDB("1l2y")
    >>> top = pt.load_ParmEd(p, save_and_reload=True) 
    """
    if save_and_reload:
        # faster
        with goto_temp_folder():
            fname = 'tmppdb.mol2'
            parmed_obj.save(fname)
            return Topology(fname)
    else:
        return load_pseudo_parm(parmed_obj)


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
