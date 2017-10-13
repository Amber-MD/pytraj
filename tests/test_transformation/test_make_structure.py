from __future__ import print_function
import unittest
import pytraj as pt
from utils import fn
from pytraj.testing import aa_eq
from pytraj.testing import cpptraj_test_dir


class TestMakeStructure(unittest.TestCase):
    def test_makestructure(self):
        # FIXME: What does this test?
        # https://github.com/Amber-MD/cpptraj/issues/27
        # load only 1st frame
        traj = pt.iterload(fn('Tc5b.x'), fn('Tc5b.top'))

        # polyproline II dihedral to residues 1-13
        t0 = traj[:1].copy()
        pt.make_structure(t0, 'pp2:1-13')
        pt.write_traj(
            'dummy_test0.pdb', t0, options='model', overwrite=True)

        t1 = traj[:1].copy()
        pt.make_structure(t1, "chi1:8:N:CA:CB:CG:35")
        pt.write_traj(
            'dummy_test1.pdb', t1, options='model', overwrite=True)

    def test_makestructure_with_reference(self):
        tz2_parm7 = cpptraj_test_dir + '/tz2.parm7'
        ref_rst7 = cpptraj_test_dir + '/tz2.rst7'
        trajin_rst7 = cpptraj_test_dir + '/Test_MakeStructure/pp2.rst7.save'

        # cpptraj
        command = """
        parm {}
        reference {}
        trajin {}
        
        makestructure "ref:1-13:tz2.rst7"
        rmsd reference 
        createcrd mycrd
        """.format(tz2_parm7, ref_rst7, trajin_rst7)
        state = pt.load_cpptraj_state(command)
        state.run()
        cpp_traj = state.data['mycrd']

        # pytraj
        traj = pt.load(trajin_rst7, top=tz2_parm7)
        ref = pt.load(ref_rst7, top=tz2_parm7)
        pt.make_structure(traj, "ref:1-13:1-13", ref=ref)
        pt.rmsd(traj, ref=ref)
        aa_eq(traj.xyz, cpp_traj.xyz)


if __name__ == "__main__":
    unittest.main()
