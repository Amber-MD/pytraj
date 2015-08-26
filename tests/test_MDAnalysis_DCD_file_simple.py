from __future__ import print_function
import unittest
from pytraj import io as mdio
from pytraj.testing import test_if_having, no_test


class Test(unittest.TestCase):
    @test_if_having("MDAnalysis")
    def test_0(self):
        # PSF, DCD
        from MDAnalysis import Universe
        from MDAnalysisTests.datafiles import PSF, DCD
        u = Universe(PSF, DCD)

        t = mdio.load_MDAnalysisIterator(u)

        for idx, f in enumerate(t):
            #print(idx, f)

        # hanged out if use iteration again
        for idx, f in enumerate(t):
            #print(idx, f)

    @test_if_having("MDAnalysis")
    def test_1(self):
        # GRO, TRR
        from MDAnalysis import Universe
        from MDAnalysisTests.datafiles import GRO, TRR
        u = Universe(GRO, TRR)

        t = mdio.load_MDAnalysisIterator(u)

        for idx, f in enumerate(t):
            #print(idx, f)

        for idx, f in enumerate(t):
            #print(idx, f)


if __name__ == "__main__":
    unittest.main()
