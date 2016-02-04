from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq


class Test(unittest.TestCase):

    def test_0(self):

        # center of mass
        trajin = """
        parm data/Tc5b.top
        trajin data/Tc5b.x
        vector center v0
        timecorr vec1 v0
        """

        state = pt.load_cpptraj_state(trajin)
        state.run()
        cpptraj_output = state.data

        traj = pt.iterload("./data/Tc5b.x", "./data/Tc5b.top")
        dslist0 = pt.center_of_mass(traj)
        data = pt.timecorr(dslist0, dslist0)
        aa_eq(data, cpptraj_output[-1].values)

        # 2 vectors
        cm = """
        parm data/Tc5b.top
        trajin data/Tc5b.x
        vector v0 :2 :5
        vector v1 :3 :7
        timecorr vec1 v0 vec2 v1
        """
        state = pt.load_cpptraj_state(cm)
        state.run()
        cpptraj_output = state.data

        dslist0 = pt.calc_vector(traj, [':2 :5', ':3 :7'])
        data = pt.timecorr(dslist0[0], dslist0[1])
        aa_eq(data, cpptraj_output[-1].values)

        # corrplane
        cm = """
        parm data/Tc5b.top
        trajin data/Tc5b.x
        vector v0 @2,@5,@9 corrplane
        vector v1 @3,@7,@20 corrplane
        timecorr vec1 v0 vec2 v1
        """

        state = pt.load_cpptraj_state(cm)
        state.run()
        cpptraj_output = state.data

        dslist0 = pt.calc_vector(traj,
                                 ['@2,@5,@9 corrplane', '@3,@7,@20 corrplane'])
        dslist1 = pt.vector.corrplane(traj, ['@2,@5,@9', '@3,@7,@20'])
        data0 = pt.timecorr(dslist0[0], dslist0[1])
        data1 = pt.timecorr(dslist1[0], dslist1[1])
        aa_eq(data0, cpptraj_output[-1].values)
        aa_eq(data1, cpptraj_output[-1].values)


if __name__ == "__main__":
    unittest.main()
