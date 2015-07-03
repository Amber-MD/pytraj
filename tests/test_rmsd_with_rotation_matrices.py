from __future__ import print_function
import os
import pytraj as pt
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
import pytraj.common_actions as pyca

txt = """
parm %s
trajin %s
reference %s
rms reference @CA savematrices
"""


class Test(unittest.TestCase):

    def test_0(self):
        from pytraj.utils.context import goto_temp_folder
        fname = os.path.abspath("data/md1_prod.Tc5b.x")
        topname = os.path.abspath("data/Tc5b.top")
        refname = os.path.abspath("data/Tc5b.crd")
        trajin = txt % (topname, fname, refname)
        traj = pt.iterload(fname, topname)
        ref = pt.iterload(refname, topname)[0]

        with goto_temp_folder():
            with open("test.in", 'w') as fh:
                fh.write(trajin)
            state = pt.io.load_cpptraj_file("test.in")
            state.run()

            # cpptraj output
            dslist0 = pt.datasetlist.DatasetList(state.datasetlist[1:])

            # pytraj output
            dslist1 = pt.rmsd_with_rotation_matrices(traj, '@CA', ref=ref)
            dslist2 = pt.rmsd_with_rotation_matrices(traj, '@CA', ref=ref, dtype='dict')

            print(pt.tools.rmsd_1darray(dslist0[0].values.flatten(),
                                        dslist1[0].values.flatten()))
            print(dslist0[1].shape)
            print(dslist1[1].shape)
            print(pt.tools.rmsd_1darray(dslist0[1].values.flatten(),
                                        dslist1[1].values.flatten()))

if __name__ == "__main__":
    unittest.main()
