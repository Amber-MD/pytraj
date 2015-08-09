from __future__ import print_function
import unittest
from pytraj import io as mdio
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
import pytraj.common_actions as pyca


# adapted test from $AMBERHOME/test/sanderapi/test.py
def assert_close(x, y, tol=1E-4):
    print("computed = %s, expected = %s" % (x, y))
    assert abs(x - y) < tol


class Test(unittest.TestCase):

    @test_if_having("sander")
    @test_if_having("parmed")
    def test_2(self):
        # compare to saved test: GB
        import os
        from pytraj.testing import amberhome
        import sander

        if amberhome:
            topfile = os.path.join(amberhome, "test/gb7_trx/prmtop_an")
            rstfile = os.path.join(amberhome, "test/gb7_trx/trxox.2.4ns.x")
            traj = mdio.load(rstfile, topfile)
            options = sander.gas_input(7)
            options.cut = 9999.0
            options.saltcon = 0.2
            options.gbsa = 1
            edict = pyca.energy_decomposition(
                traj=traj, input_options=options, parm=topfile)
            print(edict)
            assert_close(edict['bond'][0], 631.8993, tol=3E-4)
            assert_close(edict['angle'][0], 898.2543, tol=3E-4)
            assert_close(edict['surf'][0], 33.8338, tol=3E-4)
            assert_close(edict['gb'][0], -1943.0838, tol=3E-4)

            # dummy test to make sure `energy_decomposition` can work with list
            edict2 = pyca.energy_decomposition(
                traj=[traj, ], input_options=options, parm=topfile, top=traj.top)
            edict3 = pyca.energy_decomposition(
                traj=traj(), input_options=options, parm=topfile, top=traj.top)
            edict4 = pyca.energy_decomposition(
                traj=[traj[:5], traj[5:]], input_options=options, parm=topfile, top=traj.top)
            edict5 = pyca.energy_decomposition(
                traj=[traj[:5], traj(start=5)], input_options=options, parm=topfile, top=traj.top)
            # test dtype
            dslist = pyca.energy_decomposition(traj=[traj, ], input_options=options,
                                               parm=topfile, top=traj.top, dtype='dataset')
            assert edict == edict2
            assert edict == edict3
            assert edict == edict4
            assert sorted(dslist.to_dict()) == sorted(edict)
        else:
            print("has not set AMBERHOME or can not find test folder. skip")

    @test_if_having("sander")
    @test_if_having("parmed")
    def test_3(self):
        # compare to saved test: PME
        import os
        from pytraj.testing import amberhome
        import sander

        if amberhome:
            topfile = os.path.join(amberhome, "test/4096wat/prmtop")
            rstfile = os.path.join(amberhome, "test/4096wat/eq1.x")
            traj = mdio.iterload(rstfile, topfile)
            options = sander.pme_input()
            options.cut = 8.0
            edict = pyca.energy_decomposition(traj=traj, input_options=options)
            assert_close(edict['bond'][0], 0., tol=3E-4)
            assert_close(edict['vdw'][0], 6028.9517, tol=3E-4)
            assert_close(edict['elec'][0], -45371.5995, tol=3E-4)
        else:
            print("has not set AMBERHOME or can not find test folder. skip")

    @test_if_having("sander")
    @test_if_having("parmed")
    def test_4(self):
        # compare to saved test: GB + QMMM
        import os
        from pytraj.testing import amberhome
        import sander

        if amberhome:
            topfile = os.path.join(
                amberhome, "test/qmmm2/lysine_PM3_qmgb2/prmtop")
            rstfile = os.path.join(
                amberhome, "test/qmmm2/lysine_PM3_qmgb2/lysine.crd")
            traj = mdio.load(rstfile, topfile)

            options = sander.gas_input(1)
            options.cut = 99.0
            options.ifqnt = 1
            qmmm_options = sander.qm_input()
            qmmm_options.iqmatoms[:3] = [8, 9, 10]
            qmmm_options.qm_theory = "PM3"
            qmmm_options.qmcharge = 0
            qmmm_options.qmgb = 2
            qmmm_options.adjust_q = 0

            edict = pyca.energy_decomposition(
                traj=traj, input_options=options, qmmm_options=qmmm_options)
            assert_close(edict['bond'][0], 0.0016, tol=3E-4)
            assert_close(edict['vdw'][0], 0.1908, tol=3E-4)
            assert_close(edict['vdw_14'][0], 3.7051, tol=3E-4)
            assert_close(edict['elec'][0], -4.1241, tol=3E-4)
            assert_close(edict['elec_14'][0], 65.9137, tol=3E-4)
            assert_close(edict['gb'][0], -80.1406, tol=3E-4)
            assert_close(edict['scf'][0], -11.9100, tol=3E-4)
        else:
            print("has not set AMBERHOME or can not find test folder. skip")

    @test_if_having("sander")
    @test_if_having("parmed")
    def test_5(self):
        # compare to saved test: PME + QMMM
        import os
        from pytraj.testing import amberhome
        import sander

        if amberhome:
            topfile = os.path.join(
                amberhome, "test/qmmm2/MechEm_nma-spcfwbox/prmtop")
            rstfile = os.path.join(
                amberhome, "test/qmmm2/MechEm_nma-spcfwbox/inpcrd")
            traj = mdio.load(rstfile, topfile)

            options = sander.pme_input()
            options.cut = 8.0
            options.ifqnt = 1
            options.jfastw = 4

            qmmm_options = sander.QmInputOptions()
            qmmm_options.qm_theory = "PDDG-PM3"
            qmmm_options.qmmask = ":1-2"
            qmmm_options.qmcharge = 0
            qmmm_options.scfconv = 1e-10
            qmmm_tight_p_conv = 1
            qmmm_options.qmmm_int = 5

            edict = pyca.energy_decomposition(
                traj=traj, input_options=options, qmmm_options=qmmm_options)
            assert_close(edict['bond'][0], 605.7349, tol=3E-4)
            assert_close(edict['vdw_14'][0], 0.0000, tol=3E-4)
            assert_close(edict['elec_14'][0], 0.0000, tol=3E-4)
            assert_close(edict['elec'][0], -7409.7167, tol=3E-1)
            assert_close(edict['scf'][0], -37.1277, tol=3E-4)
        else:
            print("has not set AMBERHOME or can not find test folder. skip")

    @test_if_having("sander")
    @test_if_having("parmed")
    def test_6(self):
        # compare to saved test: GBneck2nu
        import os
        from pytraj.testing import amberhome
        import sander

        if amberhome:
            topfile = os.path.join(amberhome, "test/gbneck2nu_test/prmtop")
            rstfile = os.path.join(amberhome, "test/gbneck2nu_test/min.r")
            traj = mdio.load(rstfile, topfile)
            options = sander.gas_input(8)
            options.cut = 9999.0
            edict = pyca.energy_decomposition(
                traj=traj, input_options=options, parm=topfile)
            assert_close(edict['gb'][0], -2287.6880, tol=3E-4)
            assert_close(edict['gb'][0], -2287.6880, tol=3E-4)
            assert_close(edict['elec'][0], -1659.5740, tol=3E-4)
            assert_close(edict['vdw'][0], 384.2512, tol=3E-4)

if __name__ == "__main__":
    unittest.main()
