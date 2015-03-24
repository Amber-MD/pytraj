import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal

class Test(unittest.TestCase):
    def test_0(self):
        # TODO: add assert
        traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        dslist = DataSetList()
        act = adict['pairwise']
        act("@CA", traj, dslist=dslist)
        print (dslist.size)

        for ds in dslist:
            if hasattr(ds, 'mkind'):
                print (ds.name, ds.dtype, ds.mkind)
            else:
                print (ds.name, ds.dtype)

        print (dslist.get_legends())
        d0 = dslist['PW_00000[EMAP]'][0]
        print (d0)
        print (d0.size)

        try:
            from pytraj.plots import plot_matrix
            ax0 = plot_matrix(d0)
            from matplotlib import pyplot
            pyplot.show()
        except:
            raise ImportError("need matplotlib")

if __name__ == "__main__":
    unittest.main()
