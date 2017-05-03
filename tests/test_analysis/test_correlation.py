from __future__ import print_function
import unittest
import pytraj as pt
from utils import fn
from pytraj.testing import aa_eq, tempfolder
from utils import tc5b_trajin, tc5b_top


def test_corr():
    cm = """
    parm {}
    trajin {}
    distance d0 @2 @3
    distance d1 @4 @7
    corr d0 d1 out test.out
    corr d0 out test2.out
    """.format(tc5b_top, tc5b_trajin)

    with tempfolder():
        # exclude DatasetTopology (1st)
        state = pt.datafiles.load_cpptraj_state(cm).run()
        cout = state.data[1:]

        traj = pt.iterload(tc5b_trajin, tc5b_top)
        dslist = pt.calc_distance(traj, ['@2 @3', '@4, @7'])

        pout = pt.xcorr(dslist[0], dslist[1])
        # corr d0, d1
        aa_eq(pout, cout[2])

        # corr d0, d0
        pout = pt.xcorr(dslist[0], dslist[0])
        aa_eq(pout, cout[4])

        # autocorr d0, d0
        pout = pt.acorr(dslist[0])
        aa_eq(pout, cout[4])
