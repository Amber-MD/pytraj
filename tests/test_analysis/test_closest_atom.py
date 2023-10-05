import unittest
from pytraj import *
import pytraj as pt
from utils import fn
from pytraj.testing import aa_eq
import numpy as np


class Test(unittest.TestCase):
    # tests to see if function can find an atom if we search with a coordinate
    # very close to the atom's xyz coordinates
    def test_0(self):
        traj = pt.iterload(fn('Tc5b.x'), fn('Tc5b.top'))
        frame = traj[0]

        atom_to_find = 1
        coords_x, coords_y, coords_z = frame.atom(atom_to_find)
        close_idx = pt.closest_atom(
            traj.top, traj[0],
            (coords_x + 0.01, coords_y + 0.01, coords_z + 0.01))
        aa_eq(close_idx, atom_to_find)

    # ensuring that the function computes min distance atom accurately, finds closest point
    # to origin both directly and using the method.
    def test_1(self):
        traj = pt.iterload(fn('Tc5b.x'), fn('Tc5b.top'))
        origin = (0, 0, 0)
        for frame in traj:
            distances = []
            for atm in traj.top.atom_indices(""):
                atm_crd = frame.atom(atm)
                distance = np.sqrt(atm_crd[0]**2 + atm_crd[1]**2 +
                                   atm_crd[2]**2)
                distances.append(distance)
            direct_answer = distances.index(min(distances))
            func_answer = pt.closest_atom(traj.top, frame, origin)
            aa_eq(direct_answer, func_answer)

    # same as test 1 but for trp cage
    def test_2(self):
        traj = pt.datafiles.load_trpcage()
        origin = (0, 0, 0)
        for frame in traj:
            distances = []
            for atm in traj.top.atom_indices(""):
                atm_crd = frame.atom(atm)
                distance = np.sqrt(atm_crd[0]**2 + atm_crd[1]**2 +
                                   atm_crd[2]**2)
                distances.append(distance)
            direct_answer = distances.index(min(distances))
            func_answer = pt.closest_atom(traj.top, frame, origin)
            aa_eq(direct_answer, func_answer)

    # same as test 1 but for rna
    def test_3(self):
        traj = pt.datafiles.load_rna()
        origin = (0, 0, 0)
        for frame in traj:
            distances = []
            for atm in traj.top.atom_indices(""):
                atm_crd = frame.atom(atm)
                distance = np.sqrt(atm_crd[0]**2 + atm_crd[1]**2 +
                                   atm_crd[2]**2)
                distances.append(distance)
            direct_answer = distances.index(min(distances))
            func_answer = pt.closest_atom(traj.top, frame, origin)
            aa_eq(direct_answer, func_answer)

    # makes sure exception is raised for empty frame
    def test_4(self):
        traj = pt.datafiles.load_rna()
        emptyTop = pt.Topology()
        with self.assertRaises(ValueError):
            pt.closest_atom(emptyTop, traj[0], (0, 0, 0))

    # make sure exception raised for argument without 3 coordinate tuple
    def test_5(self):
        traj = pt.datafiles.load_rna()
        with self.assertRaises(ValueError):
            pt.closest_atom(traj.top, traj[0], (0, 0, 0, 0))
        with self.assertRaises(ValueError):
            pt.closest_atom(traj.top, traj[0], "(0, 0, 0)")

    # makes sure that mask argument indeed returns an atom with that mask
    def test_6(self):
        traj = pt.datafiles.load_trpcage()
        closest_idx = pt.closest_atom(traj.top, traj[0], (0, 0, 0), "@CA")
        assert (traj.top.atom(closest_idx).name == "CA")

    def test_7(self):
        # find closest atom to origin in a given trajectory
        traj = pt.iterload(fn('Tc5b.x'), fn('Tc5b.top'))
        frame = traj[0]
        print(pt.closest_atom(traj.top, traj[0], (0, 0, 0)))


if __name__ == "__main__":
    unittest.main()
