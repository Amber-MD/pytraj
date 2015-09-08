from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
import pytraj.common_actions as pyca


class TestActionListRMSD(unittest.TestCase):
    def test_0(self):
        traj = pt.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        standard_rmsd = pt.rmsd(traj, mask='@CA')

        def test_rmsd(input_traj):
            from pytraj.actions.CpptrajActions import Action_Rmsd
            from pytraj.datasets import DataSetList
            dslist = DataSetList()
            act = Action_Rmsd()
            act.read_input('first @CA', top=input_traj.top, dslist=dslist)
            act.process(input_traj.top)

            for frame in input_traj:
                act.do_action(frame)
            return (dslist.values)

        def test_rmsd_actlist(input_traj):
            from pytraj.actions.CpptrajActions import Action_Rmsd
            from pytraj.core.ActionList import ActionList
            from pytraj.datasets import DataSetList

            alist = ActionList()
            dslist = DataSetList()
            act = Action_Rmsd()
            alist.add_action(act, 'first @CA',
                             top=input_traj.top,
                             dslist=dslist)

            for frame in input_traj:
                alist.do_actions(frame)
            return (dslist.values)

        rmsd0 = test_rmsd(traj)
        rmsd1 = test_rmsd(traj[:])
        rmsd2 = test_rmsd_actlist(traj)
        rmsd3 = test_rmsd_actlist(traj[:])
        t0 = traj[:]
        rmsd4 = test_rmsd_actlist(t0)
        aa_eq(standard_rmsd, rmsd0)
        aa_eq(standard_rmsd, rmsd1)
        aa_eq(standard_rmsd, rmsd2)
        aa_eq(standard_rmsd, rmsd3)
        aa_eq(standard_rmsd, rmsd4)


if __name__ == "__main__":
    unittest.main()
