from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.utils import has_
from pytraj.decorators import test_if_having
from pytraj.utils import Timer

print ("pytraj version = 0.1.2.dev0")
 
if has_("mdtraj"):
    import mdtraj as md
    print ("mtrajd version = %s" % md.version.full_version)
    traj = mdio.load("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")[:]
    m_top = md.load_prmtop("./data/Tc5b.top") 
    m_traj = md.load_mdcrd("./data/md1_prod.Tc5b.x", m_top)

def Run(func, msg):
    print (msg)
    with Timer() as t0:
        func(traj)

    time_0 = t0.time_gap()

    with Timer() as t1:
        func(m_traj)
    time_1 = t1.time_gap()
    print ("pytraj (speed up) vs mdtraj = %s" % (float(time_1)/time_0))

class Test(unittest.TestCase):
    @test_if_having("mdtraj")
    def test_0(self):
        def iter_traj(traj):
            for frame in traj:
                pass
        Run(iter_traj, "iter_traj")

    @test_if_having("mdtraj")
    def test_1(self):
        def get_xyz(traj):
            traj.xyz
        Run(get_xyz, "get_xyz")

    @test_if_having("mdtraj")
    def test_2(self):
        def tolist(traj):
            traj.xyz.tolist()
        Run(tolist, "tolist")

    @test_if_having("mdtraj")
    def test_3(self):
        def save_nc(traj):
            traj.save("./output/x_speed.nc")
        Run(save_nc, "save .nc")

        def save_xtc(traj):
            traj.save("./output/x_speed.xtc")
        Run(save_xtc, "save .xtc")

        def save_xyz(traj):
            traj.save("./output/x_speed.xyz")
        Run(save_xyz, "save .xyz")

        def save_dcd(traj):
            traj.save("./output/x_speed.dcd")
        Run(save_dcd, "save .dcd")

        def save_binpos(traj):
            traj.save("./output/x_speed.binpos")
        Run(save_binpos, "save .binpos")

    @test_if_having("mdtraj")
    def test_4(self):
        def n_frames(traj):
            traj.n_frames
        Run(n_frames, 'n_frames')

if __name__ == "__main__":
    unittest.main()
