from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq

try:
    import sander
    has_sander = True
except ImportError:
    has_sander = False

try:
    import parmed as pmd
    has_parmed = True
except ImportError:
    has_parmed = False


@unittest.skipIf(not has_sander, 'skip if not having sander')
class TestUpdateDihedral(unittest.TestCase):

    def test_update_dihedral_parm(self):
        traj = pt.iterload("./data/Tc5b.crd", "./data/Tc5b.top")
        p = pmd.load_file(traj.top.filename)
        inp = sander.gas_input(8)
        coords = traj[0].xyz

        fname = "tmp.parm7"

        with pt.utils.context.tempfolder():
            for k in range(20, 100):
                p.bonds[3].type.k = k
                p.remake_parm()
                with sander.setup(p, coords, None, inp):
                    ene, frc = sander.energy_forces()


if __name__ == "__main__":
    unittest.main()
