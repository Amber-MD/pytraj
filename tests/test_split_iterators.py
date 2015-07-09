from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
import pytraj.common_actions as pyca


class Test(unittest.TestCase):

    def test_0(self):
        import numpy as np
        from glob import glob
        fname = "./data/md1_prod.Tc5b.x"
        ftop = "./data/Tc5b.top"
        traj = pt.iterload(fname, ftop)

        # naive
        assert np.all(pt.calc_center_of_mass(list(traj.split_iterators(4))) ==
                      pt.calc_center_of_mass(traj))

        # with mask
        assert np.all(pt.calc_center_of_mass(list(traj.split_iterators(4, mask='@CA'))) ==
                      pt.calc_center_of_mass(traj, '@CA'))

        # with mask and rmsfit
        ilist = list(
            traj.split_iterators(n_chunks=4, mask='!@H=', rmsfit=(traj[0], '@CA')))
        arr0 = pt.calc_center_of_mass(ilist)
        arr1 = pt.calc_center_of_mass(traj(rmsfit=(traj[0], '@CA')), '!@H=')
        aa_eq(arr0, arr1)

if __name__ == "__main__":
    unittest.main()
