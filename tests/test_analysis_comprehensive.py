from __future__ import print_function
import unittest
from pytraj import DatasetList
from pytraj import adict
from pytraj import calculate
from pytraj import analdict
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir
import pytraj.common_actions as pyca
from pytraj.datasets.DataSetList import DataSetList


class Test(unittest.TestCase):
    def test_0(self):
        # histogram
        # TODO : assert
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        dslist = DataSetList()
        calculate("distance", traj, "@CA @CB mylovelyname", dslist=dslist)
        d0 = dslist[0]
        #print(d0.to_pyarray())
        dmax = d0.max()
        dmin = d0.min()
        #print(dmax, dmin)
        act = analdict['hist']
        command = "%s min %s max %s bins %s step 50 out ./output/test.out" % (
            d0.name, dmin, dmax, 100)
        #print(command)
        act(command, dslist=dslist)
        #print(dslist.get_legends())
        #print(dslist[1])
        #print(dslist[1].tolist())

    def test_1(self):
        # Corr
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        dslist = calculate("distance", traj, "@CA @CB mylovelyname")
        act = analdict['corr']
        command = "mylovelyname out mydummyoutput.txt"
        act(command, dslist=dslist)
        #print(dslist[-1].to_pyarray())


if __name__ == "__main__":
    unittest.main()
