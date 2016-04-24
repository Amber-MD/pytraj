from __future__ import print_function
import os
import unittest
import pytraj as pt
from pytraj.testing import amberhome, aa_eq


# adapted test from $AMBERHOME/test/sanderapi/test.py
def assert_close(x, y, tol=1E-4):
    print("computed = %s, expected = %s" % (x, y))
    assert abs(x - y) < tol


try:
    import sander
    has_sander = True
except ImportError:
    has_sander = False


@unittest.skipIf(not has_sander, 'skip if not having sander')
class TestSander(unittest.TestCase):

    @unittest.skipIf(not amberhome, 'skip since there is no AMBERHOME')
    def test_GB(self):
        # compare to saved test: GB
        topfile = os.path.join(amberhome, "test/gb7_trx/prmtop_an")
        rstfile = os.path.join(amberhome, "test/gb7_trx/trxox.2.4ns.x")
        traj = pt.load(rstfile, topfile)
        options = sander.gas_input(7)
        options.cut = 9999.0
        options.saltcon = 0.2
        options.gbsa = 1
        edict = pt.esander(traj=traj,
                                        mm_options=options,
                                        prmtop=topfile)
        assert_close(edict['bond'][0], 631.8993, tol=3E-4)
        assert_close(edict['angle'][0], 898.2543, tol=3E-4)
        assert_close(edict['surf'][0], 33.8338, tol=3E-4)
        assert_close(edict['gb'][0], -1943.0838, tol=3E-4)

        # dummy test to make sure `esander` can work with list
        edict2 = pt.esander(traj=[traj, ],
                                         mm_options=options,
                                         prmtop=topfile,
                                         top=traj.top)
        edict3 = pt.esander(traj=traj(),
                                         mm_options=options,
                                         prmtop=topfile,
                                         top=traj.top)
        edict4 = pt.esander(traj=[traj[:5], traj[5:]],
                                         mm_options=options,
                                         prmtop=topfile,
                                         top=traj.top)
        edict5 = pt.esander(traj=[traj[:5], traj(start=5)],
                                         mm_options=options,
                                         prmtop=topfile,
                                         top=traj.top)
        # test dtype
        dslist = pt.esander(traj=[traj, ],
                                         mm_options=options,
                                         prmtop=topfile,
                                         top=traj.top,
                                         dtype='dataset')
        assert edict == edict2
        assert edict == edict3
        assert edict == edict4
        assert sorted(dslist.to_dict()) == sorted(edict)

    @unittest.skipIf(not amberhome, 'skip since there is no AMBERHOME')
    def test_PME(self):
        # compare to saved test: PME
        topfile = os.path.join(amberhome, "test/4096wat/prmtop")
        rstfile = os.path.join(amberhome, "test/4096wat/eq1.x")
        traj = pt.iterload(rstfile, topfile)
        options = sander.pme_input()
        options.cut = 8.0
        edict = pt.esander(traj=traj, mm_options=options)
        assert_close(edict['bond'][0], 0., tol=3E-4)
        assert_close(edict['vdw'][0], 6028.9517, tol=3E-4)
        assert_close(edict['elec'][0], -45371.5995, tol=3E-4)

    @unittest.skipIf(not amberhome, 'skip since there is no AMBERHOME')
    def test_PME_with_energy_decompositionr(self):
        # compare to saved test: PME
        topfile = os.path.join(amberhome, "test/4096wat/prmtop")
        rstfile = os.path.join(amberhome, "test/4096wat/eq1.x")
        traj = pt.iterload(rstfile, topfile)
        options = sander.pme_input()
        options.cut = 8.0
        edict = pt.energy_decomposition(traj=traj, mm_options=options)
        assert_close(edict['bond'][0], 0., tol=3E-4)
        assert_close(edict['vdw'][0], 6028.9517, tol=3E-4)
        assert_close(edict['elec'][0], -45371.5995, tol=3E-4)

    @unittest.skipIf(not amberhome, 'skip since there is no AMBERHOME')
    def test_GB_QMMM(self):
        # compare to saved test: GB + QMMM
        topfile = os.path.join(amberhome, "test/qmmm2/lysine_PM3_qmgb2/prmtop")
        rstfile = os.path.join(amberhome,
                               "test/qmmm2/lysine_PM3_qmgb2/lysine.crd")
        traj = pt.load(rstfile, topfile)

        options = sander.gas_input(1)
        options.cut = 99.0
        options.ifqnt = 1
        qm_options = sander.qm_input()
        qm_options.iqmatoms[:3] = [8, 9, 10]
        qm_options.qm_theory = "PM3"
        qm_options.qmcharge = 0
        qm_options.qmgb = 2
        qm_options.adjust_q = 0

        edict = pt.esander(traj=traj,
                                        mm_options=options,
                                        qm_options=qm_options)
        assert_close(edict['bond'][0], 0.0016, tol=3E-4)
        assert_close(edict['vdw'][0], 0.1908, tol=3E-4)
        assert_close(edict['vdw_14'][0], 3.7051, tol=3E-4)
        assert_close(edict['elec'][0], -4.1241, tol=3E-4)
        assert_close(edict['elec_14'][0], 65.9137, tol=3E-4)
        assert_close(edict['gb'][0], -80.1406, tol=3E-4)
        assert_close(edict['scf'][0], -11.9100, tol=3E-4)

    @unittest.skipIf(not amberhome, 'skip since there is no AMBERHOME')
    def test_PME_QMMM(self):
        # compare to saved test: PME + QMMM
        topfile = os.path.join(amberhome,
                               "test/qmmm2/MechEm_nma-spcfwbox/prmtop")
        rstfile = os.path.join(amberhome,
                               "test/qmmm2/MechEm_nma-spcfwbox/inpcrd")
        traj = pt.load(rstfile, topfile)

        options = sander.pme_input()
        options.cut = 8.0
        options.ifqnt = 1
        options.jfastw = 4

        qm_options = sander.QmInputOptions()
        qm_options.qm_theory = "PDDG-PM3"
        qm_options.qmmask = ":1-2"
        qm_options.qmcharge = 0
        qm_options.scfconv = 1e-10
        qmmm_tight_p_conv = 1
        qm_options.qmmm_int = 5

        edict = pt.esander(traj=traj,
                                        mm_options=options,
                                        qm_options=qm_options)
        assert_close(edict['bond'][0], 605.7349, tol=3E-4)
        assert_close(edict['vdw_14'][0], 0.0000, tol=3E-4)
        assert_close(edict['elec_14'][0], 0.0000, tol=3E-4)
        assert_close(edict['elec'][0], -7409.7167, tol=3E-1)
        assert_close(edict['scf'][0], -37.1277, tol=3E-4)

    @unittest.skipIf(not amberhome, 'skip since there is no AMBERHOME')
    def test_gbneck2nu(self):
        # compare to saved test: GBneck2nu
        topfile = os.path.join(amberhome, "test/gbneck2nu/1hji/prmtop")
        rstfile = os.path.join(amberhome, "test/gbneck2nu/1hji/min.r")
        traj = pt.load(rstfile, topfile)
        options = sander.gas_input(8)
        options.cut = 9999.0
        edict = pt.esander(traj=traj,
                                        mm_options=options,
                                        prmtop=topfile)
        assert_close(edict['gb'][0], -2287.6880, tol=3E-4)
        assert_close(edict['gb'][0], -2287.6880, tol=3E-4)
        assert_close(edict['elec'][0], -1659.5740, tol=3E-4)
        assert_close(edict['vdw'][0], 384.2512, tol=3E-4)

    def test_frame_indices(self):
        traj = pt.iterload('data/tz2.nc', 'data/tz2.parm7')
        frame_indices = [0, 6, 7, 4, 5]

        data_without_frame_indices = pt.esander(traj, igb=8)
        data_with_frame_indices = pt.esander(
            traj,
            igb=8,
            frame_indices=frame_indices)
        data_with_frame_indices_2 = pt.esander(
            traj[frame_indices],
            igb=8)

        for key in data_without_frame_indices:
            aa_eq(data_without_frame_indices[key][frame_indices],
                  data_with_frame_indices[key])
            aa_eq(data_without_frame_indices[key][frame_indices],
                  data_with_frame_indices_2[key])


if __name__ == "__main__":
    unittest.main()
