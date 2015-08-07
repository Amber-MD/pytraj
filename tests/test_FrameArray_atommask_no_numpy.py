import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.utils import _import, flatten_list


class Test(unittest.TestCase):
    def test_0(self):
        has_np, np = _import('numpy')
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        arr0 = traj['@CA'].xyz

        if not has_np:
            assert isinstance(arr0, list) == True
        else:
            assert hasattr(np, 'ndarray')
            assert isinstance(arr0[:], np.ndarray) == True

        # make sure that we did the right thing
        topCA = traj.top.strip_atoms("!@CA", copy=True)
        naked_traj = mdio.iterload(("./data/stripAllButCA.Tc5b.x"), topCA)

        for idx, frame in enumerate(naked_traj):
            assert_almost_equal(flatten_list(arr0[idx]), frame.coords)

        # it's time to see the output :D
        print(arr0)
        print(type(arr0))


if __name__ == "__main__":
    unittest.main()
