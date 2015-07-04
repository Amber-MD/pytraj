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

fname = os.path.abspath("./data/md1_prod.Tc5b.x")
topname = os.path.abspath("./data/Tc5b.top")
tc5b_trajin = trajin_template % (topname, fname)

fname = os.path.abspath("./data/DPDP.nc")
topname = os.path.abspath("./data/DPDP.parm7")
DPDP_trajin = trajin_template % (topname, fname)

def load_result_from_cpptraj_state(txt):
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
        return DatasetList(state.datasetlist)

load_cpptraj_output = load_result_from_cpptraj_state
