from __future__ import absolute_import
import os
from ..TrajectoryIterator import TrajectoryIterator
from ..utils.context import goto_temp_folder

__all__ = ['Ala3_crd', 'Ala3_crd_top', 'tz2_ortho_nc', 'tz2_ortho_parm7']

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
    'tc5b_trajin': ["./data/md1_prod.Tc5b.x", "./data/Tc5b.top"],
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


def load_result_from_cpptraj_state(txt, dtype=None):
    """load output from cpptraj

    Parameters
    ----------
    txt : str
        cpptraj's trajin
    dtype : str, return data type

    Returns
    -------
    if with_state=False: return DatasetList
    if with_state=True: return [CpptrajState, DatasetList]

    """
    from pytraj.io import load_cpptraj_file
    from pytraj.datasetlist import DatasetList
    from pytraj import ArgList

    command_list = list(filter(lambda x: x, txt.split("\n")))

    for idx, line in enumerate(command_list):
        if 'parm' in line:
            arglist = ArgList(line)
            # use absolute path
            fname = os.path.abspath(arglist.get_string_key('parm'))
            command_list[idx] = " ".join(('parm', fname))

        if 'trajin' in line:
            arglist = ArgList(line)
            # use absolute path
            fname = os.path.abspath(arglist.get_string_key('trajin'))
            command_list[idx] = " ".join(('trajin', fname))

    txt = "\n".join([line for line in command_list])

    with goto_temp_folder():
        with open("tmp.in", 'w') as fh:
            fh.write(txt)
        state = load_cpptraj_file("tmp.in")
        state.run()

        if dtype == 'cpptraj_dataset':
            out = state.datasetlist
        else:
            out = DatasetList(state.datasetlist)
        return out

def cpptraj_dry_run(txt):
    '''for speed comparison
    '''
    from pytraj.io import load_cpptraj_file
    from pytraj import ArgList
    command_list = list(filter(lambda x: x, txt.split("\n")))

    for idx, line in enumerate(command_list):
        if 'parm' in line:
            arglist = ArgList(line)
            # use absolute path
            fname = os.path.abspath(arglist.get_string_key('parm'))
            command_list[idx] = " ".join(('parm', fname))

        if 'trajin' in line:
            arglist = ArgList(line)
            # use absolute path
            fname = os.path.abspath(arglist.get_string_key('trajin'))
            command_list[idx] = " ".join(('trajin', fname))

    _txt = "\n".join([line for line in command_list])

    with goto_temp_folder():
        with open("tmp.in", 'w') as fh:
            fh.write(_txt)
        state = load_cpptraj_file("tmp.in")
        state.run()


def load_outtraj(txt, top=None):
    '''
    Examples
    --------
    >>> import pytraj as pt
    >>> pt.datafiles.load_outtraj(tc5b_trajin + """
    trajout test.nc
    """, top=traj.top)
    '''
    from pytraj.io import load_cpptraj_file
    from pytraj import ArgList
    import pytraj as pt

    command_list = list(filter(lambda x: x, txt.split("\n")))
    arglist = ArgList(txt)

    filelist = []
    _top = None
    for line in command_list:
        if 'trajout' in line:
            arg = ArgList(line)
            filelist.append(arg.get_string_key('trajout'))

        if 'parm' in line:
            _top = ArgList(line).get_string_key('parm')

    if top is None:
        top = _top

    with goto_temp_folder():
        with open("tmp.in", 'w') as fh:
            fh.write(txt)
        state = load_cpptraj_file("tmp.in")
        state.run()
        traj = pt.iterload(filelist[0], top)
        return traj


load_cpptraj_output = load_result_from_cpptraj_state
