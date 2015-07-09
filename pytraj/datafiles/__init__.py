from __future__ import absolute_import
import os
from .. TrajectoryIterator import TrajectoryIterator

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
        'tc5b_trajin' : ["./data/md1_prod.Tc5b.x", "./data/Tc5b.top"],
        'DPDP_trajin' : ["./data/DPDP.nc", "./data/DPDP.parm7"],
        'tz2_ortho_trajin' : ["./data/tz2.ortho.nc", "./data/tz2.ortho.parm7"],
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
    """
    from pytraj.utils.context import goto_temp_folder
    from pytraj.io import load_cpptraj_file
    from pytraj.datasetlist import DatasetList

    with goto_temp_folder():
        with open("tmp.in", 'w') as fh:
            fh.write(txt)
        state = load_cpptraj_file("tmp.in")
        state.run()
        if dtype == 'cpptraj_dataset':
            return state.datasetlist
        return DatasetList(state.datasetlist)

load_cpptraj_output = load_result_from_cpptraj_state
