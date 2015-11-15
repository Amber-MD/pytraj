#!/usr/bin/env python
from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
'''for using reference frame internally (without loading from file in Action)
'''


class TestReferenceFrame(unittest.TestCase):

    def test_reference(self):
        traj = pt.iterload("./data/tz2.nc", "./data/tz2.parm7")

        text = '''
        parm {0}
        trajin {1}
        reference {1} 3 3 [ref0]
        reference {1} 2 2 [ref1]
        rms ref [ref0]  @CA
        rms ref [ref1]  @CA

        # can not use name now
        #reference {1} 3 3 name ref0
        #reference {1} 2 2 name ref1
        #rms ref ref0  @CA
        #rms ref ref1  @CA
        '''.format(traj.top.filename, traj.filename)

        state = pt.load_cpptraj_state(text)
        state.run()
        print(state.data.keys())

        rmsd_0 = pt.rmsd(traj, ref=2, mask='@CA')
        print(rmsd_0)
        print(state.data[-1])

    @unittest.skip('not finished yet')
    def test_reference_2(self):
        traj = pt.iterload("./data/tz2.nc", "./data/tz2.parm7")
        from pytraj.actions.CpptrajActions import Action_Rmsd
        from pytraj.datasets.DatasetList import DatasetList
        act = Action_Rmsd()
        ref = traj[2]
        dslist = DatasetList()
        dslist.add_set('ref_frame', 'myref')
        dslist[-1].add_frame(ref)
        act('myrmsd ref myref @CA', traj, top=traj.top, dslist=dslist)
        print(dslist)


if __name__ == "__main__":
    unittest.main()
