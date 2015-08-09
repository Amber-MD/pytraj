import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.datasets import *
from pytraj.datasets.DataSetList import DataSetList
from pytraj.testing import aa_eq


class Test(unittest.TestCase):
    def test_0(self):
        from pytraj import datasets
        ddict = datasets.__dict__
        keys = [key for key in ddict if 'DataSet' in key]
        print(keys)

        # remove base classes
        useless_keys = ['DataSet', 'DataSet_1D', 'DataSet_2D', 'DataSet_3D',
                        'DataSet_Coords', 'DataSet_Modes', 'DataSetList']
        for _key in useless_keys:
            if _key in keys:
                keys.remove(_key)

        print(keys)

        # test create
        for key in keys:
            print(key)
            d0 = ddict[key]()
            print(d0.name, d0.dtype)

        # test resize
        d_double = DatasetDouble()
        d_float = DatasetFloat()
        d_int = DatasetInteger()
        d_v = DatasetVector()
        d_s = DatasetString()

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
        print(d.keys())
        print(d.get_dtypes())

        ds_v = traj.calc_vector("@CA @CB")
        dcp = ds_v[0].copy()
        d0 = ds_v[0]
        assert dcp.name == d0.name
        assert dcp.legend == d0.legend
        assert dcp.aspect == d0.aspect

        # shape
        d0 = d[0]
        assert d0.shape == (d0.size, )
        import numpy as np
        assert np.abs((np.mean(d0.values) - d0.avg())) < 1E-4
        print(np.mean(d0.values))
        print(d0.avg())
        print(np.sum(d0))
        print(np.array_split(d0, 3))

        # from_array
        d_double = DatasetDouble()
        d_float = DatasetFloat()
        d_int = DatasetInteger()

        d_double.from_array_like([1, 2, 3])
        d_float.from_array_like([1, 2, 3])
        d_int.from_array_like([1, 2, 3])
        assert d_double.size == d_float.size == d_int.size
        arr = [1, 2, 3]
        # values
        aa_eq(d_double, arr)
        aa_eq(d_float, arr)
        aa_eq(d_int, arr)
        print(d_double.values, d_float.values, d_int.values)

        # mean_with_error
        d2_double = d_double.copy()
        d2_double.values[:] *= 2.
        print(d2_double.values)
        assert (d2_double.mean_with_error(d_double) == (3., 1.))

        # chunk_average
        d3_double = DatasetDouble()
        d3_double.from_array_like(range(10))
        aa_eq(d3_double.chunk_average(5), np.array([0.5, 2.5, 4.5, 6.5, 8.5]))


if __name__ == "__main__":
    unittest.main()
