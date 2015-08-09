from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
import pytraj.common_actions as pyca


class Test(unittest.TestCase):

    def test_1(self):
        import parmed as pmd
        import sander

        traj = pt.iterload("./data/Tc5b.crd", "./data/Tc5b.top")
        p = pmd.load_file(traj.top.filename)
        inp = sander.gas_input(8)
        coords = traj[0].coords

        fname = "tmp.parm7"

        with pt.utils.context.goto_temp_folder():
            for k in range(20, 100):
                p.bonds[3].type.k = k
                p.remake_parm()
                with sander.setup(p, coords, None, inp):
                    ene, frc = sander.energy_forces()
                    print (ene.bond, ene.dihedral)

if __name__ == "__main__":
    unittest.main()
