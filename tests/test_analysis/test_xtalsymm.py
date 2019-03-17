import os
import pytraj as pt
from pytraj.testing import cpptraj_test_dir, aa_eq


def test_xtalsymm(tmpdir):
    with tmpdir.as_cwd():
        x6dky_parm = os.path.join(cpptraj_test_dir, 'x6dky.parm7')
        x6dky_tname = os.path.join(cpptraj_test_dir, 'mdXtal.nc')
        x6dky_refname = os.path.join(cpptraj_test_dir, 'mdXtal.inpcrd')
        state = pt.load_cpptraj_state("""
        parm {parm}
        trajin {trajin}
        reference {reference}
        xtalsymm :1-16 reference group P22(1)2(1) na 2 nb 1 nc 1
        trajout mdAsuOnly.pdb
        """.format(
            parm=x6dky_parm, trajin=x6dky_tname, reference=x6dky_refname))
        state.run()

        traj = pt.load(x6dky_tname, x6dky_parm)
        ref = pt.load(x6dky_refname, traj.top)
        pt.xtalsymm(
            traj,
            mask=':1-16',
            ref=ref,
            options="group P22(1)2(1) na 2 nb 1 nc 1")
        c_traj = pt.load('mdAsuOnly.pdb')
        aa_eq(traj.xyz, c_traj.xyz, decimal=2)
