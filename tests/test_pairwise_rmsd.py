from __future__ import print_function
import unittest
import pytraj as pt
from pytraj.utils import eq, aa_eq
from pytraj.decorators import no_test, test_if_having, test_if_path_exists

txt = '''
parm ./data/Tc5b.top
trajin ./data/md1_prod.Tc5b.x
rms2d @CA metric_holder rmsout tmp.out
'''

class TestPairWiseRMSD(unittest.TestCase):
    def testTwoTrajTypes(self, txt=txt):
        '''test different metrics with different traj objects
        '''
        funclist = [pt.iterload, pt.load]

        for func in funclist: 
            traj = func("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
            for metric in ['rms', 'nofit', 'dme']:
                d0 = pt.pairwise_rmsd(traj(mask='@CA'), metric=metric)
                d1 = pt.pairwise_rmsd(traj, mask='@CA', metric=metric)
                d2 = pt.pairwise_rmsd(traj(), mask='@CA', metric=metric)

                txt0 = txt.replace('metric_holder', metric)
                state = pt.datafiles.load_cpptraj_output(txt0, dtype='state')
                d3 = state.datasetlist[-1].values

                aa_eq(d0, d1)
                aa_eq(d0, d2)
                aa_eq(d0, d3)


if __name__ == "__main__":
    unittest.main()
