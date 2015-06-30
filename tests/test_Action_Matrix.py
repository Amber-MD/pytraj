from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import allactions
from pytraj.datasets.cast_dataset import cast_dataset
from pytraj import adict
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.decorators import no_test
from pytraj.datasets.DataSetList import DataSetList

farray = TrajectoryIterator(top=Topology("./data/Tc5b.top"),
                            filename='data/md1_prod.Tc5b.x',
                            )


class TestRadgyr(unittest.TestCase):
    #@no_test

    def test_0(self):
        dslist = DataSetList()
        act = adict['matrix']
        act(command="byres @CA", current_frame=farray,
            top=farray.top, dslist=dslist)

        d1 = cast_dataset(dslist[0], dtype="matrix double")
        print(d1.size)
        print(dir(d1))
        print(d1.n_cols, d1.n_rows)
        print(d1.dtype)
        print(d1.ndim)
        print(d1.mkind)
        print(d1.format)
        # TODO : add assert to make sure reproducing cpptraj output

        for i in range(d1.size):
            print(d1[i])

        # another way
        d0 = adict['matrix']("byres @CA", farray, farray.top, quick_get=True)
        print(d0.size)
        print(dir(d0))
        print(d0.n_cols, d0.n_rows)
        print(d0.dtype)
        print(d0.ndim)
        print(d0.mkind)
        print(d0.format)

        assert_almost_equal(d0, d1)
        d2 = adict['distance'](
            ":2@CA :10@CA", farray, farray.top, quick_get=True)
        print(d2.dtype)
        #assert d2[:] != d0[:]

    def test_1(self):
        dslist = DataSetList()
        act = adict['matrix']
        act(command="byres @CA", current_frame=farray,
            top=farray.top, dslist=dslist)
        act.print_output()
        d0 = dslist[0]
        print(d0.dtype)
        print(cast_dataset(d0, dtype=d0.dtype))
        print(dslist.get_dataset(0))

        for i in range(d0.size):
            print(d0[i])

        #print (d0.scalar_type)
        #print (d0)
        arr0 = []
        for _d in d0:
            arr0.append(_d)

        arr1 = []
        for i in range(d0.size):
            arr1.append(d0[i])

        print(arr0[:10])
        print(arr1[:10])
        print(len(arr0))
        assert_almost_equal(arr0, arr1)

        print(d0.get_element(10, 10))
        for i in range(d0.n_rows):
            for j in range(d0.n_cols):
                d0.get_element(i, j)

        fullmat = d0.get_full_matrix()
        print(type(fullmat))
        print(len(fullmat))

        assert_almost_equal(arr1[:20], fullmat[:20])
        try:
            from pytraj.plottingplot_matrix import plot_matrix
            from pytraj.plottingbase import plt
            ax0 = plot_matrix(d0)
            print(ax0)
            # plt.show()
            plt.savefig("./output/test_saveplot.png")
        except:
            print("don't have numpy, matplotlib. Ignore")

if __name__ == "__main__":
    unittest.main()
