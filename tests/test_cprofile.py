import cProfile, pstats
from pytraj import io as mdio
from pytraj import *

def load_test():
    traj = mdio.load("./data/NuG2/test.x.000", "./data/NuG2/NuG2.top")
    act = adict['surf']
    act("@CA", traj)
    farray = traj[:]

cProfile.runctx("load_test()", globals(), locals(), "Profile.prof")
s = pstats.Stats("Profile.prof")
s.strip_dirs().sort_stats("time").print_stats()
