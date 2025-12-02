from __future__ import print_function
import unittest
import pytraj as pt
from utils import fn
from pytraj.testing import aa_eq, tempfolder, load_cpptraj_reference_data
from utils import tc5b_trajin, tc5b_top, fn


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


def test_correlation_cpptraj_reference():
    """Test correlation against cpptraj reference using tz2 data"""
    # Load tz2 trajectory - same as cpptraj Test_Corr
    traj = pt.load(fn('tz2.nc'), fn('tz2.parm7'))

    # Calculate distance :2 :12 (cpptraj: distance d1 :2 :12)
    distances = pt.calc_distance(traj, ':2 :12')

    # Calculate auto-correlation (cpptraj: corr d1 d1 out corr.dat)
    corr_result = pt.acorr(distances)

    # Load cpptraj reference data
    ref_data = load_cpptraj_reference_data('Test_Corr', 'corr.dat.save')

    # Compare correlation results
    aa_eq(corr_result, ref_data, decimal=4)
