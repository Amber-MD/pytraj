from __future__ import print_function
import unittest
from pytraj import io as mdio
from pytraj.decorators import no_test, test_if_having, test_if_path_exists


class Test(unittest.TestCase):
    def test_0(self):
        from pytraj._common_actions import _get_top
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        saved_top = traj.top.copy()

        # get top from traj
        _top = _get_top(traj, None)
        assert _top.n_atoms == saved_top.n_atoms

        # get top from top
        _top = _get_top(None, traj.top)
        assert _top.n_atoms == saved_top.n_atoms

        # get top from file
        _top = _get_top(None, "./data/Tc5b.top")
        assert _top.n_atoms == saved_top.n_atoms

        # get top from list of trajs
        _top = _get_top([traj, 1], None)
        assert _top.n_atoms == saved_top.n_atoms

        # get top from list of trajs 2
        _top = _get_top([1, traj, 0], None)
        assert _top.n_atoms == saved_top.n_atoms

        # get top from frame_iter
        _top = _get_top(traj(), None)
        assert _top is not None

        # get top from chunk_iter
        _top = _get_top(traj.iterchunk(), None)
        assert _top is None


if __name__ == "__main__":
    unittest.main()
