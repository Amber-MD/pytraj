from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq


class Test(unittest.TestCase):

    def test_0(self):
        traj = pt.iterload("./data/Tc5b.x", "./data/Tc5b.top")

        traj(0, 8, 2, mask='@CA').save('output/test0.nc', overwrite=True)
        pt.write_traj('./output/test1.nc',
                      traj(0,
                           8,
                           2,
                           mask='@CA'),
                      overwrite=True)

        new_top = traj.top._get_new_from_mask('@CA')
        t0 = pt.iterload('./output/test0.nc', new_top)
        t1 = pt.iterload('./output/test1.nc', new_top)

        aa_eq(t0.xyz, t1.xyz)


if __name__ == "__main__":
    unittest.main()
