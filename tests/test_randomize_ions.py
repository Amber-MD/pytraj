import unittest
from pytraj.common_actions import randomize_ions
from pytraj.io import load
from pytraj import adict


class TestRandomizeIons(unittest.TestCase):

    def test_0(self):
        # get `traj` instance (Trajectory)
        traj = load(filename="./Test_RandomizeIons/adh206.tip3p.rst7.gz",
                    top="./Test_RandomizeIons/adh206.ff10.tip3p.parm7.gz")
        # get 1st frame from `traj`
        frame0 = traj[0]
        frame0_1 = traj[0].copy()
        fsaved = frame0.copy()

        # randomize ions for frame0
        randomize_ions(traj=frame0,
                       top=traj.top.copy(),
                       command="@Na+ around :1-16 by 5.0 overlap 3.0",)

        # make sure to reproduce cpptraj output
        savedframe = load(filename="./Test_RandomizeIons/random.crd.save",
                          top="./Test_RandomizeIons/adh206.ff10.tip3p.parm7.gz")[0]

        # another way
        adict['randomizeions'](
            "@Na+ around :1-16 by 5.0 overlap 3.0", frame0_1, traj.top)
        print(frame0.rmsd(savedframe))
        print(frame0_1.rmsd(savedframe) < 1E-3)
        assert frame0.rmsd(savedframe) < 1E-3

        # all atoms
        _rmsd = frame0.rmsd(fsaved)
        print(_rmsd)
        _rmsd = frame0.rmsd(fsaved, use_mass=True)
        print(_rmsd)

        # only P
        _rmsd = frame0.rmsd(fsaved, traj.top("@P"))
        print(_rmsd)

        # another way
        _rmsd = frame0.rmsd(fsaved, top=traj.top, mask="@P")
        print(_rmsd)

if __name__ == "__main__":
    unittest.main()
