from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
import pytraj.common_actions as pyca


class Test(unittest.TestCase):
    def test_0(self):
        traj = pt.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")

        # convert to mutable traj
        t0 = traj[:1]

        for deg in range(-170, 170, 10):
            pt.rotate_dihedral(t0, "custom:3:phi:" + str(deg))
            _deg = pt.calc_phi(t0, '3', dtype='ndarray')[0]
            dih = pt.dihedral(t0, ':2@C :3@N :3@CA :3@C')[0]
            aa_eq(deg, _deg)
            aa_eq(deg, dih)

            aa_eq(pt.calc_psi(traj[:1]).values, pt.calc_psi(t0))

    def test_1(self):
        # different from test_0 a bit in `mask`
        traj = pt.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")

        # convert to mutable traj
        t0 = traj[:1]

        for deg in range(-170, 170, 10):
            pt.rotate_dihedral(t0, "3:phi:" + str(deg))
            _deg = pt.calc_phi(t0, '3', dtype='ndarray')[0]
            dih = pt.dihedral(t0, ':2@C :3@N :3@CA :3@C')[0]
            print(deg, _deg, dih)
            aa_eq(deg, _deg)
            aa_eq(deg, dih)

            aa_eq(pt.calc_psi(traj[:1]).values, pt.calc_psi(t0))


if __name__ == "__main__":
    unittest.main()
