from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
import pytraj.common_actions as pyca


class Test(unittest.TestCase):

    def test_0(self):
        traj = pt.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        print(pt.rmsd(traj, mask='@CA'))

        def test_rmsd(input_traj):
            from pytraj.actions.CpptrajActions import Action_Rmsd
            from pytraj.datasets import DataSetList
            dslist = DataSetList()
            act = Action_Rmsd()
            act.read_input('first @CA', top=input_traj.top, dslist=dslist)
            act.process(input_traj.top)

            for frame in input_traj:
                act.do_action(frame)
            print(dslist.values)

        def test_rmsd_actlist(input_traj):
            from pytraj.actions.CpptrajActions import Action_Rmsd
            from pytraj.core.ActionList import ActionList
            from pytraj.datasets import DataSetList

            alist = ActionList()
            dslist = DataSetList()
            act = Action_Rmsd()
            alist.add_action(act, 'first @CA', top=input_traj.top, dslist=dslist)

            for frame in input_traj:
                alist.do_actions(frame)
            print(dslist.values)


        test_rmsd(traj)
        test_rmsd(traj[:])
        test_rmsd_actlist(traj)
        test_rmsd_actlist(traj[:])
        #t0 = traj[:]
        #test_rmsd_actlist(t0)
        #print(t0.xyz)

if __name__ == "__main__":
    unittest.main()
