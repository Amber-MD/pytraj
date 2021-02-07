from __future__ import print_function
import pytraj as pt
from utils import fn
import numpy as np
import unittest
from pytraj.testing import aa_eq
import pytest


class TestSearchHbonds(unittest.TestCase):
    def test_hbond_general(self):
        traj = pt.iterload(fn('DPDP.nc'), fn('DPDP.parm7'))
        dslist = pt.search_hbonds(traj, dtype='dataset')
        for key in dslist.keys():
            if 'UU' not in key:
                assert len(dslist[key].values) == traj.n_frames
        mydict = dslist.to_dict()
        mydict_np = dslist.to_dict()
        assert len(mydict.keys()) == dslist.size
        assert len(mydict_np.keys()) == dslist.size

        for key in mydict.keys():
            mydict[key] = np.asarray(mydict[key])
            aa_eq(mydict[key], mydict_np[key])

        # raise if dtype='hbond' and series=False
        with pytest.raises(ValueError):
            pt.hbond(traj, series=False, dtype='hbond')

    def test_hbonds_with_image(self):
        traj = pt.iterload(fn('tz2.ortho.nc'), fn('tz2.ortho.parm7'))

        hbonds_0 = pt.search_hbonds(traj(autoimage=True))
        hbonds_1 = pt.search_hbonds(traj, image=True)
        aa_eq(hbonds_0.values, hbonds_1.values)

    def test_hbonds_solvent_bridge(self):
        traj = pt.iterload(fn('tz2.ortho.nc'), fn('tz2.ortho.parm7'))
        # find solvent bridge between residue 10 and 11.
        hb = pt.hbond(
            traj, ':10,11', solvent_donor=':WAT', solvent_acceptor=':WAT')
        assert hb._get_bridge() ==  \
         [[('WAT208', 'THR10', 'TRP11')],
          [('WAT208', 'THR10', 'TRP11')],
          [('WAT208', 'THR10', 'TRP11')],
          [('WAT208', 'THR10', 'TRP11')],
          [],
          [],
          [],
          [],
          [],
          [('WAT266', 'THR10', 'TRP11')]]

    def test_hbonds_from_pdb(self):
        traj = pt.load(fn('1L2Y.pdb'))
        hb = pt.search_hbonds(traj)

        state = pt.load_cpptraj_state('''
        parm {pdb}
        trajin {pdb}
        hbond series
        '''.format(pdb=fn('1L2Y.pdb')))
        state.run()

        for data_p, data_cpp in zip(hb.data, state.data[1:]):
            assert len(
                data_p) == traj.n_frames == 38, 'size of dataset must be 38'
            aa_eq(data_p, data_cpp.values)

        # make sure distances are smaller than cutoff
        distance_cutoff = 2.5
        angle_cutoff = 135.
        hb = pt.search_hbonds(traj)
        distances = pt.distance(traj, hb.get_amber_mask()[0])
        angles = pt.angles(traj, hb.get_amber_mask()[1])
        dist_indices = np.where(distances > distance_cutoff)
        angle_indices = np.where(angles < angle_cutoff)

        saved_donor_acceptors = [
            'ASP9_OD2-ARG16_NH1-HH12', 'ASP9_OD2-ARG16_NH2-HH22',
            'ASP9_OD2-ARG16_NE-HE', 'ASP9_OD2-ARG16_NH2-HH21',
            'ASP9_OD2-ARG16_NH1-HH11'
        ]

        donor_acceptors = pt.search_hbonds(traj, ':9,16').donor_acceptor
        assert saved_donor_acceptors == donor_acceptors, 'saved_donor_acceptors'

        aa_eq(hb.total_solute_hbonds(), hb.data['total_solute_hbonds'])


if __name__ == "__main__":
    unittest.main()
