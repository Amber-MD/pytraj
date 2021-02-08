import pytraj as pt
from utils import fn
from array import array
from pytraj import *
import numpy as np
from pytraj.testing import aa_eq


class TestCHARMM:
    def test_psf(self):
        top = pt.load_topology(fn('ala3.psf'))
        list(top.residues)
        atm = AtomMask("@CA")
        top._set_integer_mask(atm)

        atm.invert_mask()
        frame = Frame(atm.n_atoms)
        frame[:10] = np.asarray(array('d', list(range(30)))).reshape(10, 3)

    def test_write(self, tmpdir):
        traj = pt.iterload(fn("ala3.dcd"), fn('ala3.psf'))
        with tmpdir.as_cwd():
            output = "dummy_save_charmm_to_amber.x"
            traj.save(output, overwrite=True)
            # test loading
            trajamber = pt.iterload(output,
                                      fn('ala3.psf'))
            for i in range(traj.n_frames):
                aa_eq(trajamber[i].xyz, traj[i].xyz, decimal=3)
