from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from numpy.testing import assert_almost_equal
from pytraj.utils import has_
from pytraj.decorators import test_if_having
from pytraj.utils import Timer

print ("pytraj version = 0.1.2.dev0")
 
if has_("mdtraj"):
    import mdtraj as md
    print ("mtrajd version = %s" % md.version.full_version)
    # use Trajectory (in memory for comparison)
    traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")[:]
    #traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
    m_top = md.load_prmtop("./data/Tc5b.top") 
    m_traj = md.load_mdcrd("./data/md1_prod.Tc5b.x", m_top)

def Run(func, msg, n_times=50, test_load=False):
    print (msg)
    my_ratio = 0
    for _ in range(n_times):
        if not test_load:
            with Timer() as t0:
                func(traj)

            time_0 = t0.time_gap()

            with Timer() as t1:
                func(m_traj)
            time_1 = t1.time_gap()
            my_ratio += time_1 / time_0
        else:
            txt = "./data/md1_prod.fit_to_first.Tc5b.x"
            fa = Trajectory()
            fa.top = traj.top
            with Timer() as t0:
                fa.load(txt)
            time_0 = t0.time_gap()

            with Timer() as t1:
                f_mdtraj = md.load_mdcrd(txt)
            time_1 = t1.time_gap()
            assert_almost_equal(fa.xyz, f_mdtraj.xyz)
            my_ratio += time_1 / time_0
    print ("pytraj (speed up) vs mdtraj = %s" % (my_ratio/n_times))
    print ("mdtraj (speed up) vs pytraj = %s" % (1./(my_ratio/n_times)))
    print ()

class Test(unittest.TestCase):
    @test_if_having("mdtraj")
    def test_load(self):
        def load(test_load=True):
            pass
        Run(load, "test_load_traj: text file")

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

    @test_if_having("mdtraj")
    def test_5(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")[:]
        m_top = md.load_prmtop("./data/Tc5b.top") 
        m_traj = md.load_mdcrd("./data/md1_prod.Tc5b.x", m_top)

        print ("get xyz for each frame")
        m_f = m_traj[0]
        p_f = traj[0]

        my_ratio = 0
        n_times = 10

        for _idx in range(n_times):
            with Timer() as t0:
                #a_p = p_f.xyz
                a_p = p_f[:]

            with Timer() as t0:
                a_p = p_f.xyz

            with Timer() as t1:
                # mdtraj * 10
                a_m = m_f.xyz * 10

            from numpy.testing import assert_almost_equal
            assert_almost_equal(a_p, a_m[0], decimal=4)

            ratio = t0.time_gap() / t1.time_gap()
            my_ratio += ratio
        print ("pytraj (speed up) vs mdtraj = %s" % (my_ratio/n_times))
        print ("mdtraj (speed up) vs pytraj = %s" % (1./(my_ratio/n_times)))


if __name__ == "__main__":
    unittest.main()
