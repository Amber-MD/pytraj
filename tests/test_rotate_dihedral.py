from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq


class Test(unittest.TestCase):

    def test_0(self):
        traj = pt.iterload("./data/Tc5b.x", "./data/Tc5b.top")

        # convert to mutable traj
        t0 = traj[:1]

        for deg in range(-170, 170, 10):
            pt.rotate_dihedral(t0, "custom:3:phi:" + str(deg))
            _deg = pt.calc_phi(t0, '3', dtype='ndarray')[0]
            dih = pt.dihedral(t0, ':2@C :3@N :3@CA :3@C')[0]
            aa_eq(deg, _deg)
            aa_eq(deg, dih)

            aa_eq(pt.calc_psi(traj[:1]).values, pt.calc_psi(t0).values)

    def test_1(self):
        # different from test_0 a bit in `mask`
        traj = pt.iterload("./data/Tc5b.x", "./data/Tc5b.top")

        # convert to mutable traj
        t0 = traj[:1]

        for deg in range(-170, 170, 10):
            pt.rotate_dihedral(t0, "3:phi:" + str(deg))
            _deg = pt.calc_phi(t0, '3', dtype='ndarray')[0]
            dih = pt.dihedral(t0, ':2@C :3@N :3@CA :3@C')[0]
            aa_eq(deg, _deg)
            aa_eq(deg, dih)

            aa_eq(pt.calc_psi(traj[:1]).values, pt.calc_psi(t0).values)

    def test_2(self):
        traj = pt.iterload("./data/Tc5b.nat.crd", "./data/Tc5b.top")
        t0 = traj[:1]
        pt.set_dihedral(t0, resid='4', dihedral_type='phi', deg=120)
        dih = pt.calc_phi(t0, resrange='4').values[0]
        assert abs(dih - 120) < 1E-3
        t0.save('test.pdb', options='model', overwrite=True)


if __name__ == "__main__":
    unittest.main()
