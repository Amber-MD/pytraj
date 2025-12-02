"""
Diffusion and transport analysis functions
"""
from .base import *

__all__ = [
    'diffusion', 'tordiff', 'toroidal_diffusion'
]


@super_dispatch()
def diffusion(traj,
              mask="",
              tstep=1.0,
              individual=False,
              top=None,
              dtype='dataset',
              frame_indices=None):
    '''compute diffusion

    Parameters
    ----------
    traj : Trajectory-like or iterator
    mask : str, default '' (all atoms)
    tstep : float, time step between frames, default 1.0 ps
    individual : bool, default False
    top : Topology, optional
    dtype : str, default 'dataset'
    frame_indices : array or None

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> data = pt.diffusion(traj, dtype='dict')
    >>> data['X']
    array([ 0.        ,  0.87027302,  1.64626022,  2.26262651,  2.98068114,
            3.57075535,  4.07030655,  4.71894512,  5.42302306,  6.01310377])
    '''
    time_step = 'time ' + str(tstep)
    individual_option = 'individual' if individual else ''

    # add 'df' as label
    label = 'df'
    command = ' '.join((mask, label, time_step, individual_option))

    action_datasets, _ = do_action(traj, command, c_action.Action_Diffusion)

    # make the label nicer
    for dataset in action_datasets:
        dataset.key = dataset.key.replace('[', '').replace(']', '').replace(label, '')
    return get_data_from_dtype(action_datasets, dtype=dtype)


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