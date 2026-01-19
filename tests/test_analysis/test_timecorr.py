from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.testing import aa_eq, tempfolder, load_cpptraj_reference_data

# local
from utils import fn

tc5b_trajin = fn('Tc5b.x')
tc5b_top = fn('Tc5b.top')


def test_timecorr():

    with tempfolder():
        # center of mass
        trajin = """
        parm {}
        trajin {}
        vector center v0
        timecorr vec1 v0
        """.format(tc5b_top, tc5b_trajin)

        state = pt.load_cpptraj_state(trajin)
        state.run()
        cpptraj_output = state.data

        traj = pt.iterload(tc5b_trajin, tc5b_top)
        dslist0 = pt.center_of_mass(traj)
        data = pt.timecorr(dslist0, dslist0)
        aa_eq(data, cpptraj_output[-1].values)

        # 2 vectors
        cm = """
        parm {}
        trajin {}
        vector v0 :2 :5
        vector v1 :3 :7
        timecorr vec1 v0 vec2 v1
        """.format(tc5b_top, tc5b_trajin)
        state = pt.load_cpptraj_state(cm)
        state.run()
        cpptraj_output = state.data

        dslist0 = pt.vector.vector(traj, [':2 :5', ':3 :7'])
        data = pt.timecorr(dslist0[0], dslist0[1])
        aa_eq(data, cpptraj_output[-1].values)

        # corrplane
        cm = """
        parm {}
        trajin {}
        vector v0 @2,@5,@9 corrplane
        vector v1 @3,@7,@20 corrplane
        timecorr vec1 v0 vec2 v1
        """.format(tc5b_top, tc5b_trajin)

        state = pt.load_cpptraj_state(cm)
        state.run()
        cpptraj_output = state.data

        dslist0 = pt.vector.vector(
            traj, ['@2,@5,@9 corrplane', '@3,@7,@20 corrplane'])
        dslist1 = pt.vector.corrplane(traj, ['@2,@5,@9', '@3,@7,@20'])
        data0 = pt.timecorr(dslist0[0], dslist0[1])
        data1 = pt.timecorr(dslist1[0], dslist1[1])
        aa_eq(data0, cpptraj_output[-1].values)
        aa_eq(data1, cpptraj_output[-1].values)


def test_timecorr_cpptraj_reference():
    """Test timecorr against cpptraj reference using tz2 data"""
    # Load tz2 trajectory - same as cpptraj Test_Timecorr
    traj = pt.load(fn('tz2.nc'), fn('tz2.parm7'))

    # Calculate vector :2 :4 (cpptraj: vector v1 :2 :4)
    vectors = pt.vector.vector(traj, ':2 :4')

    # Calculate auto time correlation (cpptraj: analyze timecorr vec1 v1)
    timecorr_result = pt.timecorr(vectors, vectors)

    # Load cpptraj reference data - v1.auto.dat.save
    ref_data = load_cpptraj_reference_data('Test_Timecorr', 'v1.auto.dat.save')

    assert ref_data is not None, "cpptraj reference data Test_Timecorr/v1.auto.dat.save must be available"

    # Compare time correlation results
    aa_eq(timecorr_result, ref_data, decimal=4)