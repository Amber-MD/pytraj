from __future__ import print_function
import unittest
from pytraj import io as mdio

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        dslist = traj.search_hbonds()
        print (dslist.groupby("SER").to_dict())

    def test_0(self):
        import pytraj.common_actions as pyca
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        dslist = pyca.calculate("multidihedral", "", traj) 

        # groupby legend
        dlg = dslist.groupby("phi", mode='legend')
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

if __name__ == "__main__":
    unittest.main()
