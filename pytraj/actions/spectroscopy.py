"""
Spectroscopy analysis functions
"""
from .base import *

__all__ = [
    'infraredspec'
]


@super_dispatch()
def infraredspec(traj=None, mask="", out=None, maxlag=None, tstep=None, rawout=None, extra_options=None, dtype='dict', top=None, frame_indices=None):
    """Compute the infrared spectrum.

    Parameters
    ----------
    traj : Trajectory-like
    mask : str, atom mask
    out : str, optional
        Output filename for the IR spectrum.
    maxlag : int, optional
        Maximum lag for the calculation.
    tstep : float, optional
        Time step between frames.
    rawout : str, optional
        Output filename for raw data.
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
    # Build the command string
    command = (CommandBuilder()
               .add(mask)
               .add("out", out, condition=out is not None)
               .add("maxlag", str(maxlag), condition=maxlag is not None)
               .add("tstep", str(tstep), condition=tstep is not None)
               .add("rawout", rawout, condition=rawout is not None)
               .add(extra_options, condition=extra_options is not None)
               .build())

    # Initialize the action
    action = c_action.Action_InfraredSpectrum()
    action_datasets = CpptrajDatasetList()
    top = traj.top
    action.read_input(command, top=top, dslist=action_datasets)

    # Set up the action
    crdinfo = dict(has_velocity=True)
    action.setup(top, crdinfo=crdinfo)

    # Execute the action for each frame
    for frame in traj:
        action.compute(frame)

    # Finalize the action
    action.post_process()

    return get_data_from_dtype(action_datasets, dtype=dtype)