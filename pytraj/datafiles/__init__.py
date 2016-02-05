from __future__ import absolute_import
import os
from ..trajectory_iterator import TrajectoryIterator
from ..utils.context import tempfolder
from .datafiles import DataFile, DataFileList
from .load_samples import *
from .load_cpptraj_file import load_cpptraj_file

__all__ = ['load_cpptraj_state', 'load_cpptraj_output', 'Ala3_crd',
           'Ala3_crd_top', 'tz2_ortho_nc', 'tz2_ortho_parm7',
           'load_cpptraj_file']

mydir = os.path.dirname(os.path.abspath(__file__))

Ala3_crd = os.path.join(mydir, "Ala3", "Ala3.crd")
Ala3_crd_top = os.path.join(mydir, "Ala3", "Ala3.top")
tz2_ortho_nc = os.path.join(mydir, "tz2", "tz2.ortho.nc")
tz2_ortho_parm7 = os.path.join(mydir, "tz2", "tz2.ortho.parm7")

trajin_template = """
parm %s
trajin %s
"""

pairlist = {
    'tc5b_trajin': ["./data/Tc5b.x", "./data/Tc5b.top"],
    'DPDP_trajin': ["./data/DPDP.nc", "./data/DPDP.parm7"],
    'tz2_ortho_trajin': ["./data/tz2.ortho.nc", "./data/tz2.ortho.parm7"],
}


def _get_trajin_text(alist, template=trajin_template):
    fname = os.path.abspath(alist[0])
    topname = os.path.abspath(alist[1])
    return trajin_template % (topname, fname)


tc5b_trajin = _get_trajin_text(pairlist['tc5b_trajin'])
DPDP_trajin = _get_trajin_text(pairlist['DPDP_trajin'])
tz2_ortho_trajin = _get_trajin_text(pairlist['tz2_ortho_trajin'])


def load_cpptraj_output(txt, dtype=None):
    """load output from cpptraj

    Parameters
    ----------
    txt : str
        cpptraj's trajin
    dtype : str, return data type

    Returns
    -------
    if dtype is 'state', return CpptrajState
    if dtype is 'ndarray', return ndarray and so on

    """
    from pytraj.core.c_core import _load_batch
    from pytraj.datasetlist import DatasetList
    from pytraj import ArgList

    commands = list(filter(lambda x: x, txt.split("\n")))

    for idx, line in enumerate(commands):
        if 'parm' in line:
            arglist = ArgList(line)
            # use absolute path
            fname = os.path.abspath(arglist.get_string_key('parm'))
            commands[idx] = " ".join(('parm', fname))

        if 'trajin' in line:
            arglist = ArgList(line)
            # use absolute path
            relative_fname = arglist.get_string_key('trajin')
            the_rest_of_line = ' '.join(line.split(relative_fname)[1:])
            fname = os.path.abspath(relative_fname)
            commands[idx] = " ".join(('trajin', fname, the_rest_of_line))

    txt = "\n".join([line for line in commands])
    state = _load_batch(txt, traj=None)

    state.run()

    if dtype == 'state':
        out = state
    else:
        out = DatasetList(state.datasetlist)
    return out


def load_cpptraj_state(txt, traj=None):
    """load text to CpptrajState

    Examples
    --------
    >>> # load from string
    >>> import pytraj as pt
    >>> state  = pt.datafiles.load_cpptraj_state('''
    ...                        parm data/tz2.ortho.parm7
    ...                        trajin data/tz2.ortho.nc
    ...                        autoimage
    ...                        center origin
    ...                        rms''')
    >>> state.run()
    CpptrajState, include:
    <datasetlist: 2 datasets>

    >>> # load from string with given TrajectoryIterator
    >>> traj = pt.datafiles.load_tz2_ortho()
    >>> state  = pt.datafiles.load_cpptraj_state('''
    ...                        autoimage
    ...                        center origin
    ...                        rms''', traj=traj)
    >>> state.run()
    CpptrajState, include:
    <datasetlist: 2 datasets>
    """
    if os.path.exists(txt):
        # txt is a file
        return load_cpptraj_file(txt)
    else:
        from pytraj.core.c_core import _load_batch
        return _load_batch(txt, traj=traj)
