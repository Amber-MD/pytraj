import unittest
from pytraj.common_actions import randomize_ions
from pytraj.io import load


def randions():
    # get `traj` instance (Trajectory)
    traj = load(filename="../tests/Test_RandomizeIons/adh206.tip3p.rst7.gz",
                top="../tests/Test_RandomizeIons/adh206.ff10.tip3p.parm7.gz")
    # get 1st frame from `traj`
    frame0 = traj[0]

    # randomize ions for frame0
    randomize_ions(
        traj=frame0,
        top=traj.top,
        command="randomizeions @Na+ around :1-16 by 5.0 overlap 3.0", )

    # make sure to reproduce cpptraj output
    savedframe = load(
        filename="../tests/Test_RandomizeIons/random.crd.save",
        top="../tests/Test_RandomizeIons/adh206.ff10.tip3p.parm7.gz")[0]

    assert frame0.rmsd(savedframe) < 1E-3


if __name__ == "__main__":
    randions()
