from pytraj.utils import has_
from pytraj.warnings import PytrajWarningMissing
from ._load_pseudo_parm import load_pseudo_parm
from ..Trajectory import Trajectory
from ..Frame import Frame

def load_ParmEd(parmed_obj, restype="top"):
    """return pytraj's Topology or Trajectory objects

    Parameters
    ---------
    parmed_obj : ParmEd's Structure object
    restype : str {'top', 'traj'}
       return type
    """
    ptop = load_pseudo_parm(parmed_obj)
    if restype.lower() == 'top':
        return ptop
    elif restype.lower() == 'traj':
        if parmed_obj.coords is None:
            raise ValueError("can not convert to Traj with None-coords")
        else:
            fa = Trajectory()
            fa.top = ptop
            frame = Frame()
            frame.set_from_crd(parmed_obj.coords)
            fa.append(frame)
            return fa
    else:
        raise ValueError("only support `top` or `traj` keyword")

def _load_chem(parm_name):
    has_parmed = has_("chemistry")
    if has_parmed:
        from chemistry import load_file
        return load_file(parm_name)
    else:
        if not has_parmed:
            PytrajWarningMissing("`chemistry`")
        return None

def to_ParmEd(pytraj_top):
    # TODO: exten to gromacs, charmm too
    # need to change extension
    """convert to ParmEd object"""
    from pytraj.utils.context import goto_temp_folder
    from pytraj.parms.ParmFile import ParmFile
    import chemistry as chem

    # I am not a fan of saving/loading again but this might be best choice
    with goto_temp_folder():
        fname = "tmp_pytrajtop.prmtop"
        ParmFile().writeparm(pytraj_top, fname, fmt="")
        return chem.load_file(fname)
