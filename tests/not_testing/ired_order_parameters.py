#!/usr/bin/env python
from __future__ import print_function
import numpy as np
import pytraj as pt
from pytraj.all_actions import _ired


def NH_order_paramters(traj, vector_pairs, order=2, tstep=1., tcorr=10000.):
    '''the results seem ok.
    '''

    # compute N-H vectors and ired matrix
    vecs_and_mat = pt.ired_vector_and_matrix(traj, vector_pairs, order=order)
    state_vecs = vecs_and_mat[:-1].values
    mat_ired = vecs_and_mat[-1]
    mat_ired /= mat_ired[0, 0]

    # get eigenvalues and eigenvectors
    modes = pt.matrix.diagonalize(mat_ired, n_vecs=len(state_vecs))[0]
    evals, evecs = modes.eigenvalues, modes.eigenvectors

    data = pt._ired(state_vecs, modes=(evals, evecs), tcorr=tcorr, tstep=tstep)
    order = [d.values.copy() for d in data if 'S2' in d.key][0]
    return (order)


def main():
    # add filenames
    parmfile = '../cpptraj/test/Test_IRED/1IEE_A_prot.prmtop'
    trajfile = '../cpptraj/test/Test_IRED/1IEE_A_test.mdcrd'

    # load to TrajectoryIterator, first 4000 frames
    traj = pt.iterload(trajfile, parmfile)

    # create N-H pairs
    h_indices = pt.select_atoms(traj.top, '@H')
    n_indices = pt.select_atoms(traj.top, '@H') - 1
    nh_indices = list(zip(n_indices, h_indices))
    print(NH_order_paramters(traj, nh_indices, tcorr=2000.))


if __name__ == "__main__":
    main()
