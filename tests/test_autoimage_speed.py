import unittest  # no test
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
import pytraj.io as io
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.decorators import test_if_having, no_test
from pytraj.compat import izip as zip
from pytraj.utils import Timer
from pytraj.externals.six.moves import range


class Test(unittest.TestCase):
    def test_1(self):
        traj = mdio.iterload(
            "./data/tz2.truncoct.nc", "./data/tz2.truncoct.parm7")

        fa = traj[:]
        for i in range(100):
            fa += traj[:]

        fa2 = fa.copy()

        @Timer()
        def normal():
            fa.autoimage()

        @Timer()
        def supposed_faster():
            fa2._autoimage_faster()

        assert_almost_equal(fa.xyz, fa2.xyz)

        normal()
        supposed_faster()

        print(fa2, fa)
        # Conclusion: not much faster


if __name__ == "__main__":
    pass
    # unittest.main()
