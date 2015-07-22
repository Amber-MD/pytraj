from pytraj.utils import has_
from pytraj.warnings import PytrajWarningMissing
from ._load_pseudo_parm import load_pseudo_parm
from ..Trajectory import Trajectory
from ..Frame import Frame
from ..Topology import Topology
from ..utils.context import goto_temp_folder


def load_ParmEd(parmed_obj, restype="top", save_and_reload=True):
    """return pytraj's Topology or Trajectory objects

    Parameters
    ---------
    parmed_obj : ParmEd's Structure object
    restype : str {'top', 'traj'}
       return type
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
            if parmed_obj.coordinates is None:
                none_coords = True
            else:
                none_coords = False

            if hasattr(parmed_obj, 'write_parm') and none_coords:
                fname = 'tmpparm.parm7'
                parmed_obj.write_parm(fname)
                return Topology(fname)
            elif hasattr(parmed_obj, 'coordinates') and none_coords:
                n_atoms = len(parmed_obj.atoms)
                need_to_reset_coords = True
                # make fake pdb
                import numpy as np
                parmed_obj.coordinates = np.array([0. for _ in range(3 * n_atoms)]
                                                 ).reshape(n_atoms, 3)
                parmed_obj.write_pdb(fname)
                return Topology(fname)
            else:
                need_to_reset_coords = False

            parmed_obj.write_pdb(fname)

            if restype == 'top':
                result = Topology(fname)
            elif restype == 'traj':
                if none_coords:
                    raise ValueError("need coordinates to load to traj")
                result= Trajectory(fname, fname)
            else:
                raise ValueError()

            if need_to_reset_coords:
                parmed_obj.coordinates = None
        return result

    else:
        ptop = load_pseudo_parm(parmed_obj)
        if restype.lower() == 'top':
            return ptop
        elif restype.lower() == 'traj':
            if hasattr(parmed_obj, 'coordinates'):
                coords = parmed_obj.coordinates
            elif hasattr(parmed_obj, 'coords'):
                coords = parmed_obj.coords
            if coords is None:
                raise ValueError("can not convert to Traj with None-coords")
            else:
                fa = Trajectory()
                fa.top = ptop
                try:
                    shape = coords.shape
                except AttributeError:
                    import numpy as np
                    coords = np.asarray(coords)
                shape = coords.shape
                if len(shape) in [1, 2]:
                    coords = coords.reshape(1, fa.top.n_atoms, 3)
                    shape = coords.shape
                print(shape)
                fa._allocate(shape[0], shape[1])
                fa.update_coordinates(coords)
                return fa
        else:
            raise ValueError("only support `top` or `traj` keyword")


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
