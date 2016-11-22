import pytraj as pt
from pytraj.testing import aa_eq
from utils import fn
import pytest


def test_analyze_modes():

    for mask in ['@CA', '@N,CA', '@CA,H']:
        command = '''
        matrix mwcovar name tz2 {}
        diagmatrix tz2 name my_modes vecs 20
        # modes fluct name my_modes
        '''.format(mask)
        
        traj = pt.iterload(fn('tz2.nc'), fn('tz2.parm7'))

        # cpptraj
        state = pt.load_cpptraj_state(command, traj)
        state.run()
        cpp_dict = state.data[1:].to_dict()
        c_dslist = state.data[1:]
        cpp_modes = c_dslist[1]

        # pytraj
        mat =  pt.matrix.mwcovar(traj, mask)
        indices = traj.top.select(mask)
        dslist = pt.matrix.diagonalize(mat, n_vecs=20, scalar_type='mwcovar',
                mass=traj.top.mass[indices],
                dtype='dataset')
        evecs, evals = dslist[0].eigenvectors, dslist[0].eigenvalues
        aa_eq(cpp_dict['tz2'], mat)
        aa_eq(evals, cpp_modes.eigenvalues)
        aa_eq(evecs, cpp_modes.eigenvectors)

        with pytest.raises(AssertionError):
            # wrong mass
            dslist = pt.matrix.diagonalize(mat, n_vecs=20, scalar_type='mwcovar',
                    mass=traj.top.mass[indices].tolist() + [1],
                    dtype='dataset')
        with pytest.raises(AssertionError):
            # mass is None
            dslist = pt.matrix.diagonalize(mat, n_vecs=20, scalar_type='mwcovar',
                    mass=None,
                    dtype='dataset')
