from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils.context import tempfolder
from pytraj.utils import eq, aa_eq
from pytraj.utils import Timer
from pytraj.compat import zip


class Test_iterchunk_autoimage(unittest.TestCase):

    def setUp(self):
        self.traj = pt.iterload("data/tz2.ortho.nc", "data/tz2.ortho.parm7")
        self.mask = '@C,N,CA,O'

    def test_only_autoimage(self):
        traj = self.traj.copy()
        fa = traj[:]
        ref00 = traj[0]
        ref01 = ref00.copy()

        # test autoimage
        for chunk in traj.iterchunk(chunksize=2, autoimage=True):
            pass
        fa0 = traj[-2:]
        fa0.autoimage()
        aa_eq(chunk.xyz, fa0.xyz)

    def test_only_rmsfit(self):
        traj = self.traj.copy()
        fa = traj[:]
        ref00 = traj[0]

        # test rmsfit
        ref00 = traj[0]
        for chunk in traj.iterchunk(chunksize=2, rmsfit=(ref00, self.mask)):
            pass

        fa0 = traj[-2:]
        ref00 = traj[0]
        fa0.rmsfit(ref=ref00, mask=self.mask)
        aa_eq(chunk.xyz, fa0.xyz)

    def test_rmsfit_with_autoimage(self):
        traj = self.traj.copy()

        # test rmsfit and autoimage
        ref0 = traj[0]
        # need to autoimage reference frame first
        pt.autoimage(ref0, top=traj.top)
        for chunk in traj.iterchunk(chunksize=2,
                                    rmsfit=(ref0, self.mask),
                                    autoimage=True):
            pass

        fa0 = traj[-2:]
        fa0.autoimage()
        ref01 = traj[0]
        pt.autoimage(ref01, top=traj.top)
        fa0.rmsfit(ref=ref01, mask=self.mask)
        aa_eq(chunk.xyz, fa0.xyz[-2:])

    def test_rmsfit_with_autoimage_compared_to_cpptraj(self):
        # assert to cpptraj: need to set mass
        traj = self.traj.copy()

        txt = '''
        parm {0}
        trajin {1}
        autoimage
        rms first {2} mass
        trajout tmp.nc
        '''.format(traj.top.filename, traj.filename, self.mask)

        with tempfolder():
            state = pt.datafiles.load_cpptraj_output(txt, dtype='state')
            state.run()
            # need to load to memory (not iterload)
            saved_traj = pt.load('tmp.nc', traj.top)

        fa1 = traj[:]
        fa1.autoimage()
        pt.superpose(fa1, ref=0, mask=self.mask, mass=True)

        aa_eq(fa1.xyz, saved_traj.xyz)

        fa_saved_nowat = saved_traj['!:WAT']
        fa1_nowat = fa1['!:WAT']
        aa_eq(fa_saved_nowat.xyz, fa1_nowat.xyz)


if __name__ == "__main__":
    unittest.main()
