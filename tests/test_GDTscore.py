import unittest
from array import array
from pytraj.base import *
from pytraj import io as mdio
from pytraj.decorators import no_test
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.common_actions import calc_score

class Test(unittest.TestCase):
    def test_0(self):
        traj = mdio.iterload("./data/md1_prod.Tc5b.x", "./data/Tc5b.top")
        # write pdb files for TMalign program so we can compare our result to TMalign
        # ./TMalign -A test_gdt_0.pdb -B test_gdt_1.pdb

        # do our calculation
        # score = 'gdtscore', 'tmscore' or 'maxsub'
        # need to add assert
        score = 'tmscore'
        ref = traj[9]
        tmscore = calc_score(traj, ref=ref, mask="@CA", top=traj.top, score=score)
        print (tmscore)
        print (tmscore[8])
        assert_almost_equal([tmscore[8],], [0.38941,], decimal=2) # 0.38941: from TMalign

if __name__ == "__main__":
    unittest.main()
    #pass
