from pytraj.utils import has_
from pytraj.warnings import PytrajWarningMissing
from ._load_pseudo_parm import load_pseudo_parm
from ..FrameArray import FrameArray
from ..Frame import Frame

def load_ParmEd(parmed_obj, restype="top"):
    """return pytraj's Topology or FrameArray objects

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
            fa = FrameArray()
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
