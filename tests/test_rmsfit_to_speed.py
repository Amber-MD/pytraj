import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal as aa_eq
from pytraj.utils import Timer
from pytraj.decorators import test_if_path_exists

tip3pdir = "./data/nogit/tip3p/"
class Test(unittest.TestCase):
    @test_if_path_exists(tip3pdir)
    def test_frame_fit(self):
        traj = mdio.iterload(tip3pdir + "/md.trj", tip3pdir + "/tc5bwat.top")[:200]
        traj.autoimage()
        f00 = traj[0].copy()
        f01 = traj[0].copy()
        fa0 = traj.copy()
        fa1 = traj.copy()
        fa2 = traj.copy()

        _f00 = f00.copy() # make copy for reference
        @Timer()
        def mode_pytraj_rms():
            fa0.rmsfit(_f00, "@CA,C,N,O", mode='pytraj')

        _f00 = f00.copy()
        @Timer()
        def mode_cpptraj_rms():
            fa2.rmsfit(_f00, "@CA,C,N,O", mode='cpptraj')

        act = adict['rmsd']

        _f00 = f00.copy()
        @Timer()
        def cpptraj_rms():
            # cpptraj use 1st frame as default reference
            act("@CA,C,N,O", [_f00, fa1], top=fa1.top)

        print ("mode_pytraj_rms")
        mode_pytraj_rms()
        print ("mode_cpptraj_rms")
        mode_cpptraj_rms()
        print ("direct cpptraj")
        cpptraj_rms()

        aa_eq(fa0.xyz, fa1.xyz)
        aa_eq(fa0.xyz, fa2.xyz)

if __name__ == "__main__":
    unittest.main()
