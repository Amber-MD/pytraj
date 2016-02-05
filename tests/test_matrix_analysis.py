#!/usr/bin/env python
from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj.testing import cpptraj_test_dir
from pytraj.utils import tempfolder

not_byres = '''
matrix correl @N @C out mtest.4.dat
matrix correl @CA out mtest.9.dat
matrix covar @N @C out mtest.5.dat
matrix covar @CA out mtest.10.dat
matrix mwcovar @N @C out mtest.6.dat
matrix mwcovar @CA out mtest.11.dat
matrix idea @CA out mtest.8.dat
matrix distcovar :1-4@CA out mtest.13.dat
matrix dist @N @C out mtest.12.dat
matrix dist @CA out mtest.7.dat
matrix dist @N @CA out mtest.3.dat bymask
matrix dist @CA out mtest.2.dat bymask
'''

byres_cm = '''
matrix dist @CA out mtest.0.dat byres
matrix dist @N @CA out mtest.1.dat byres
'''

all_commands = not_byres + byres_cm

# return a list of non-blank lines
command_list = list(filter(lambda x: x, all_commands.split("\n")))


class TestMatrixConprehensive(unittest.TestCase):

    def test_matrix(self):
        import numpy as np
        from pytraj import ArgList
        from pytraj import matrix as ma
        from pytraj.externals.six import iteritems

        matrix_test_dir = cpptraj_test_dir + "/Test_Matrix/"
        top_file = matrix_test_dir + "/1rrb_vac.prmtop"
        crd_file = matrix_test_dir + "/1rrb_vac.mdcrd"
        traj = pt.iterload(crd_file, top_file)

        with tempfolder():
            state = pt.load_cpptraj_state(all_commands, traj)
            state.run()

            state_byres = pt.load_batch(traj, byres_cm)
            state_byres.run()

            byres_matlist = []

            # no byres keyword
            for idx, line in enumerate(command_list):
                arg = ArgList(line)
                # get function
                act_key = arg.get_string_key("matrix")
                slist = arg.get_string_key('out').split(".")
                mask = arg.get_next_mask()
                fname = ".".join((slist[0], slist[-1], slist[1]))
                # get correct name
                func = ma.__dict__[act_key]

                # get command
                command = line.split(act_key)[1]
                matout = func(traj, command, dtype='ndarray')

                # cpptraj output has only 3 digits after decimal
                if 'byres' not in command:
                    aa_eq(matout.flatten(), state.data[idx + 1].values)
                else:
                    # save data for asserting later
                    byres_matlist.append(matout)

            byres_arr = np.array(byres_matlist, dtype='f8')
            # only take byres datasets
            saved_matbyres = state_byres.data[[1, 3]].values
            aa_eq(byres_arr, saved_matbyres)


if __name__ == "__main__":
    unittest.main()
