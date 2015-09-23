from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.datasets.DatasetList import DatasetList


class TestOrderParameter(unittest.TestCase):
    def test_0(self):
        # TODO : add assert
        # results seem wrong

        traj = mdio.iterload('./data/DOPC.rst7', './data/DOPC.parm7')

        command = '''lipiorder out ./output/sn2_dir.dat z scd
        :OL2@C12 :OL2@H2R :OL2@H2S :OL2@C13 :OL2@H3R :OL2@H3S \
        :OL2@C14 :OL2@H4R :OL2@H4S :OL2@C15 :OL2@H5R :OL2@H5S \
        :OL2@C16 :OL2@H6R :OL2@H6S :OL2@C17 :OL2@H7R :OL2@H7S \
        :OL2@C18 :OL2@H8R :OL2@H8S :OL2@C19 :OL2@H9R :OL2@H9R \
        :OL2@C110 :OL2@H10R :OL2@H10R :OL2@C111 :OL2@H11R :OL2@H11S \
        :OL2@C112 :OL2@H12R :OL2@H12S :OL2@C113 :OL2@H13R :OL2@H13S \
        :OL2@C114 :OL2@H14R :OL2@H14S :OL2@C115 :OL2@H15R :OL2@H15S \
        :OL2@C116 :OL2@H16R :OL2@H16S :OL2@C117 :OL2@H17R :OL2@H17S \
        :OL2@C118 :OL2@H18R :OL2@H18S
        '''

        act = adict['orderparameter']
        dslist = DatasetList()
        act(command, traj, dslist=dslist)


if __name__ == '__main__':
    unittest.main()
