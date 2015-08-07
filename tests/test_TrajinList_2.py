from pytraj import *
from pytraj.trajs.Trajout import Trajout
from pytraj.TrajinList import TrajinList as _TrajinList


class TrajinList2(list):
    def write(self, filename="", format='UNKNOWN_TRAJ', top=None):
        with Trajout(filename=filename, format=format, top=top) as trajout:
            for traj in self:
                for frame in traj:
                    trajout.write(0, frame, top)


class TrajinList:
    def __init__(self):
        self._tlist = TrajinList()
    # overwrite __iter__ to assign Topology for frame

    def __iter__(self):
        pass


if __name__ == "__main__":
    traj = io.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
    tlist = TrajinList2()
    tlist.append(traj)
    tlist.append(traj)
    print(tlist)

    tlist.write("./output/test_TrajList2.x", top=traj.top)
