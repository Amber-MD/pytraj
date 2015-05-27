from __future__ import print_function
import unittest
from pytraj import io as mdio

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        dslist = traj.search_hbonds()
        print (dslist.groupby("SER").to_dict())

    def test_1(self):
        import pytraj.common_actions as pyca
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        dslist = pyca.calculate("multidihedral", traj)

        # groupby legend
        dlg = dslist.groupby("phi", mode='legend')
        dlg2 = dslist("phi", mode='legend')
        assert sorted(dlg.keys()) == sorted(dlg2.keys())
        for d0 in dlg:
            assert ('phi' in d0.legend) == True

        # groupby aspect
        key = "chip"
        mode = 'aspect'
        dnew = dslist.groupby(key, mode)
        for d0 in dnew:
            assert (key in d0.aspect) == True

        key = "omega"
        mode = 'aspect'
        dnew = dslist.groupby(key, mode)
        for d0 in dnew:
            assert (key in d0.aspect) == True

        # m_torsion
        key = "torsion"
        mode = "scalar_mode"
        dnew = dslist.groupby(key, mode)
        for d0 in dnew:
            assert (key in d0.scalar_mode) == True

    def test_2(self):
        import pytraj.common_actions as pyca
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        d = traj.calc_multidihedral(dtype='dataset')
        d2 = d.groupby("phi:2+")
        assert sorted(d2.keys()) == sorted(['phi:2', 'phi:20'])

if __name__ == "__main__":
    unittest.main()
