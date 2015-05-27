# no unittest (and no travis test)
from __future__ import print_function
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq
from pytraj.utils import Timer
from pytraj._shared_methods import _frame_iter

def test_0():
    traj = mdio.iterload("./data/nogit/tip3p/md.trj", "./data/nogit/tip3p/tc5bwat.top")
    traj2 = mdio.iterload("./data/nogit/tip3p/md.trj", "./data/nogit/tip3p/tc5bwat.top")
    fa = traj[:]
    print (traj)

    stop = 100
    chunk = 17

    @Timer()
    def normal(chunk, stop):
        for chunk in traj.chunk_iter(chunk=chunk, stop=stop):
            pass

    @Timer()
    def iterate(stop):
        i = 0
        for frame in traj:
            i += 1
            if i == stop:
                break

    @Timer()
    def frame_iter(stop):
        for frame in traj.frame_iter(stop=stop):
            pass

    @Timer()
    def _shared_frame_iter(stop):
        for frame in _frame_iter(traj, stop=stop):
            pass

    @Timer()
    def fa_frame_iter(stop):
        for frame in fa.frame_iter(stop=stop):
            pass

    @Timer()
    def fa_normal_iterating(stop):
        i = 0
        for frame in fa:
            i += 1
            if i == stop:
                break
            pass

    def check_eq(stop, chunk):
        for c0, c1 in zip(traj.chunk_iter(stop=stop, chunk=chunk), traj2.chunk_iter(stop=stop, chunk=chunk)):
            aa_eq(c0.xyz, c1.xyz)

    print ("chunk: trajiter")
    normal(chunk, stop)

    print ("frame_iter: trajiter")
    frame_iter(stop)

    print ("_shared_frame_iter: trajiter")
    _shared_frame_iter(stop)

    print ("normal iterating with stop: trajiter")
    iterate(stop)

    print ("normal iterate with stop: fa")
    fa_normal_iterating(stop)

    print ("frame_iter fa")
    fa_frame_iter(stop)

    #check_eq(stop, chunk) # YES

if __name__ == "__main__":
    test_0()
