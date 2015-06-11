import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.datasets import *
from pytraj.testing import aa_eq

class Test(unittest.TestCase):
    def test_0(self):
        from pytraj import datasets
        ddict = datasets.__dict__
        keys = [key for key in ddict if 'DataSet' in key]
        print (keys)

        # remove base classes
        useless_keys = ['DataSet', 'DataSet_1D', 'DataSet_2D', 'DataSet_3D', 
                        'DataSet_Coords', 'DataSet_Modes']
        for _key in useless_keys:
            if _key in keys:
                keys.remove(_key)

        print (keys)

        # test create
        for key in keys:
            print (key)
            d0 = ddict[key]()
            print (d0.name, d0.dtype)

        # test resize
        d_double = DataSet_double()
        d_float = DataSet_float()
        d_int = DataSet_integer()
        d_v = DataSet_Vector()
        d_s = DataSet_string()
        
        N = 100
        d_double.resize(N)
        d_float.resize(N)
        d_int.resize(N)
        d_v.resize(N)
        d_s.resize(N)

        assert d_double.size == d_float.size == d_int.size == N
        assert d_v.size == d_s.size == N

        # test copy
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        d = traj.calc_dssp(dtype='dataset')
        print (d.keys())
        print (d.get_dtypes())

        ds_new = d.groupby('integer', mode='dtype')
        d_new_cp = ds_new[0].copy()
        aa_eq(d_new_cp.data, ds_new[0].data)

        ds_new = d.groupby('float', mode='dtype')
        d_new_cp = ds_new[0].copy()
        aa_eq(d_new_cp.data, ds_new[0].data)

        ds_v = traj.calc_vector("@CA @CB")
        dcp = ds_v[0].copy()

        # shape
        d0 = d[0]
        assert d0.shape == (d0.size,)
        import numpy as np
        assert np.abs((np.mean(d0.values) == d0.avg())) < 1E-4
        print (np.mean(d0.values))
        print (d0.avg())
        print (np.sum(d0))
        print (np.array_split(d0, 3))

if __name__ == "__main__":
    unittest.main()
