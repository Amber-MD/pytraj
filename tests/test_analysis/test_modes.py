import pytraj as pt
from pytraj.testing import aa_eq, cpptraj_test_dir
from utils import fn
import pytest
import numpy as np


def aa_eq_abs(arr0, arr1):
    aa_eq(np.abs(arr0), np.abs(arr1))


def test_analyze_modes():
    for mask in ['@CA', '@N,CA', '@CA,H']:
        command = '''
        matrix mwcovar name tz2 {}
        diagmatrix tz2 name my_modes vecs 20
        modes fluct name my_modes
        '''.format(mask)

        traj = pt.iterload(fn('tz2.nc'), fn('tz2.parm7'))

        # cpptraj
        state = pt.load_cpptraj_state(command, traj)
        state.run()
        state.run()
        cpp_dict = state.data[1:].to_dict()
        c_dslist = state.data[1:]
        cpp_modes = c_dslist[1]

        # pytraj
        mat = pt.matrix.mwcovar(traj, mask)
        indices = traj.top.select(mask)
        dslist = pt.matrix.diagonalize(
            mat,
            n_vecs=20,
            scalar_type='mwcovar',
            mass=traj.top.mass[indices],
            dtype='dataset')
        evecs, evals = dslist[0].eigenvectors, dslist[0].eigenvalues
        aa_eq_abs(cpp_dict['tz2'], mat)
        aa_eq_abs(evals, cpp_modes.eigenvalues)
        aa_eq_abs(evecs, cpp_modes.eigenvectors)

        with pytest.raises(AssertionError):
            # wrong mass
            dslist = pt.matrix.diagonalize(
                mat,
                n_vecs=20,
                scalar_type='mwcovar',
                mass=traj.top.mass[indices].tolist() + [1],
                dtype='dataset')
        with pytest.raises(AssertionError):
            # mass is None
            dslist = pt.matrix.diagonalize(
                mat,
                n_vecs=20,
                scalar_type='mwcovar',
                mass=None,
                dtype='dataset')

        fluct = pt.analyze_modes(
            'fluct', evecs, evals, scalar_type='mwcovar', dtype='dataset')
        p_rms = fluct['FLUCT_00001[rms]'].values
        c_rms = c_dslist['FLUCT_00003[rms]'].values
        aa_eq_abs(p_rms, c_rms)


def test_mode_disp():
    tz2_evecs = cpptraj_test_dir + '/Test_Analyze_Modes/tz2.evecs.dat'
    c_dslist = pt.io.read_data(tz2_evecs)
    c_modes = c_dslist[0]
    disp_dict = pt.analyze_modes(
        'displ', c_modes.eigenvectors, c_modes.eigenvalues, dtype='dict')

    command = """
    readdata {} name tz2modes
    analyze modes displ name tz2modes
    """.format(tz2_evecs)
    state = pt.load_cpptraj_state(command)
    state.run()
    c_disp_dict = state.data[1:].to_dict()

    for key in disp_dict.keys():
        aa_eq_abs(disp_dict[key], c_disp_dict[key])
