from __future__ import print_function
import unittest
from pytraj.base import *
from pytraj import adict
from pytraj import io as mdio
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists
from pytraj.testing import cpptraj_test_dir, duplicate_traj
import pytraj.common_actions as pyca

cpptraj_trajin = """
matrix dist @CA out mtest.0.dat byres
matrix dist @N @CA out mtest.1.dat byres
matrix dist @CA out mtest.2.dat bymask
matrix dist @N @CA out mtest.3.dat bymask
matrix correl @N @C out mtest.4.dat
matrix covar @N @C out mtest.5.dat
matrix mwcovar @N @C out mtest.6.dat
matrix dist @CA out mtest.7.dat
matrix idea @CA out mtest.8.dat
matrix correl @CA out mtest.9.dat
matrix covar @CA out mtest.10.dat
matrix mwcovar @CA out mtest.11.dat
matrix dist @N @C out mtest.12.dat
matrix distcovar :1-4@CA out mtest.13.dat
"""

# return a list of non-blank lines
command_list = list(filter(lambda x: x, cpptraj_trajin.split("\n")))
print(command_list)


class Test(unittest.TestCase):
    @test_if_path_exists(cpptraj_test_dir)
    def test_0(self):
        import numpy as np
        from pytraj import ArgList
        from pytraj import matrix_analysis as ma
        from pytraj.externals.six import iteritems

        new_dict = dict((v, k) for k, v in iteritems(ma.default_key_dict))

        matrix_test_dir = cpptraj_test_dir + "/Test_Matrix/"
        top_file = matrix_test_dir + "/1rrb_vac.prmtop"
        crd_file = matrix_test_dir + "/1rrb_vac.mdcrd"
        traj = mdio.iterload(crd_file, top_file)

        for line in command_list:
            arg = ArgList(line)
            # get function
            act_key = arg.get_string_key("matrix")
            slist = arg.get_string_key('out').split(".")
            mask = arg.get_mask_next()
            fname = ".".join((slist[0], slist[-1], slist[1]))
            # get correct name
            saved_file_name = matrix_test_dir + fname + ".save"
            saved_mat = np.loadtxt(saved_file_name).transpose()
            func = ma.__dict__[new_dict[act_key]]
            # get command
            command = line.split(act_key)[1]
            print(line)
            print("command = %s, func = %s" % (command, func))
            print("saved file dir = '%s'" % saved_file_name)
            mat_out = func(traj, command, dtype='ndarray')

            if 'byres' in command:
                # TODO: simplify
                mat_out = mat_out[0]

            if 'bymask' in command:
                pass
            else:
                aa_eq(mat_out.flatten(), saved_mat.flatten())


if __name__ == "__main__":
    unittest.main()
