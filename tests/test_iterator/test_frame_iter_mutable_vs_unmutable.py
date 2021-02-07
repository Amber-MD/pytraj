import pytraj as pt
from pytraj.testing import aa_eq
from pytraj import Trajectory, TrajectoryIterator

from utils import fn


def test():
    # TrajectoryIterator object
    traj0 = pt.iterload(fn('tz2.ortho.nc'), fn('tz2.ortho.parm7'))
    assert isinstance(traj0, TrajectoryIterator)

    # Trajectory object
    traj1 = traj0[:]
    assert isinstance(traj1, Trajectory)

    # check same coords
    ref = traj0[0]

    for f0, f1 in zip(
            traj0(autoimage=True, rmsfit=(ref, '@CA,C,N')),
            traj1(autoimage=True, rmsfit=(ref, '@CA,C,N'))):
        aa_eq(f0.xyz, f1.xyz)
        assert f0.rmsd_nofit(f1) == 0.
