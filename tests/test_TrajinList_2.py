from pytraj import *
from pytraj.trajs.Trajout import Trajout

class TrajinList2(list):
    def write(self, filename="", fmt='UNKNOWN_TRAJ', top=None):
        with Trajout(filename=filename, fmt=fmt, top=top) as trajout:
            for traj in self:
                for frame in traj:
                    trajout.writeframe(0, frame, top)

if __name__ == "__main__":
    traj = io.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
    tlist = TrajinList2()
    tlist.append(traj)
    tlist.append(traj)
    print (tlist)

    tlist.write("./output/test_TrajList2.x", top=traj.top)
