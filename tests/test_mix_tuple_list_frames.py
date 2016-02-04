import unittest
import pytraj as pt
from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal

count = 0


def count_frames(traj):
    global count
    if isinstance(traj, Frame):
        count += 1
    elif isinstance(traj, Trajectory):
        for i, frame in enumerate(traj):
            count_frames(frame)
    elif isinstance(traj, (list, tuple)):
        for tmtraj in traj:
            count_frames(tmtraj)
    else:
        for i, frame in enumerate(traj):
            count_frames(frame)


def main():
    global count
    traj = mdio.iterload("./data/Tc5b.x", "./data/Tc5b.top")
    boring_list = [traj[0], traj[1], traj, traj(1, 6, 2),
                   traj.iterchunk(chunksize=4)]
    count_frames(boring_list)


if __name__ == "__main__":
    main()
