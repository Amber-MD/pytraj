"""
Energy analysis functions
"""
from .base import *

__all__ = ['ene_decomp']


@super_dispatch()
def ene_decomp(traj=None,
               mask="",
               savecomponents=False,
               out=None,
               extra_options="",
               dtype='dataset',
               top=None,
               frame_indices=None):
    """Perform energy decomposition analysis.

    Parameters
    ----------
    traj : Trajectory-like
    mask : str, atom mask
        Mask to select atoms for energy decomposition.
    savecomponents : bool, default False
        Save individual energy components if True.
    out : str, optional
        Output filename for the energy decomposition results.
    extra_options : str, optional
        Additional cpptraj options for future extensibility.
    dtype : str, default 'dataset'
        Output data type.
    top : Topology, optional
    frame_indices : array-like, optional

    Returns
    -------
    DatasetList or ndarray
    """
    command = (CommandBuilder().add(mask).add(
        "savecomponents", condition=savecomponents).add(
            "out", out, condition=out
            is not None).add(extra_options,
                             condition=bool(extra_options)).build())
    action_datasets, _ = do_action(traj, command, c_action.Action_EneDecomp)
    return get_data_from_dtype(action_datasets, dtype=dtype)
