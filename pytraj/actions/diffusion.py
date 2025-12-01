"""
Diffusion and transport analysis functions
"""
from .base import *

__all__ = [
    'diffusion', 'tordiff', 'toroidal_diffusion'
]


@super_dispatch()
def diffusion(traj,
              mask='',
              dtype='dataset',
              top=None,
              frame_indices=None):
    """compute diffusion

    Parameters
    ----------
    traj : Trajectory-like
    mask : str, optional
        atom mask
    dtype : str, default 'dataset'
        return data type
    top : Topology, optional
    frame_indices : array-like, optional

    Returns
    -------
    out : DatasetList or ndarray

    Notes
    -----
    This method is equal to `cpptraj` with command
    cpptraj.run('diffusion {}'.format(mask))
    """
    action = c_action.Action_Diffusion()
    c_dslist = CpptrajDatasetList()
    action.read_input(mask, top=traj.top, dslist=c_dslist)
    action.setup(traj.top)

    for frame in traj:
        action.compute(frame)

    action.post_process()
    return get_data_from_dtype(c_dslist, dtype=dtype)


@super_dispatch()
def tordiff(traj=None, mask="", mass=False, out=None, diffout=None, time=1.0, extra_options="", dtype='dict', top=None, frame_indices=None):
    """Calculate diffusion using the toroidal-view-preserving scheme.

    Parameters
    ----------
    traj : Trajectory-like
    mask : str, atom mask
        Mask to select molecules for the calculation.
    mass : bool, default False
        Use center of mass if True, else geometric center.
    out : str, optional
        Output filename for average diffusion sets.
    diffout : str, optional
        Output filename for diffusion results.
    time : float, default 1.0
        Time between frames (ps).
    extra_options : str, optional
        Additional cpptraj options for future extensibility.
    dtype : str, default 'dataset'
        Output data type.
    top : Topology, optional
    frame_indices : array-like, optional

    Returns
    -------
    dict or DatasetList depending on dtype
    """
    command = (CommandBuilder()
               .add("TOR") # dataset name
               .add(mask, condition=bool(mask))
               .add("out", out, condition=out is not None)
               .add("diffout", diffout, condition=diffout is not None)
               .add("mass", condition=mass)
               .add("time", str(time), condition=time != 1.0)
               .add(extra_options, condition=bool(extra_options))
               .build())

    action_datasets, _ = do_action(traj, command, c_action.Action_ToroidalDiffusion)
    return get_data_from_dtype(action_datasets, dtype=dtype)


toroidal_diffusion = tordiff