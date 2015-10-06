#!/usr/bin/env python
from __future__ import print_function
import numpy as np
import pytraj as pt
from pytraj.common_actions import _ired

def main():
    # add filenames
    parmfile =  '../cpptraj/test/Test_IRED/1IEE_A_prot.prmtop'
    trajfile = '../cpptraj/test/Test_IRED/1IEE_A_test.mdcrd'

    #parmfile = '6LYT.HIP.99SB.JSC.mb3.top'
    #trajfile = 'output_md1_prod/md1_prod.x'

    # load to TrajectoryIterator, first 4000 frames
    #traj = pt.iterload(trajfile, parmfile, frame_slice=(0, 4000))
    traj = pt.iterload(trajfile, parmfile)
    print(traj)

    # create N-H pairs
    h_indices = pt.select_atoms(traj.top, '@H')
    n_indices = pt.select_atoms(traj.top, '@H') - 1
    nh_indices = list(zip(n_indices, h_indices))

    # compute N-H vectors and ired matrix
    vecs_and_mat = pt.ired_vector_and_matrix(traj, mask=nh_indices, order=2)
    state_vecs = vecs_and_mat[:-1].values
    mat_ired = vecs_and_mat[-1]
    mat_ired /= mat_ired[0, 0]

    # get eigenvalues and eigvenvectors
    evals, evecs = np.linalg.eigh(mat_ired)

    # need to sort a numpy array bit to match to cpptraj's order
    evals = evals[::-1]
    evecs = evecs[:, ::-1].T
    print(evals)

    data = _ired(state_vecs, modes=(evals, evecs), tcorr=4000, tstep=1.)
    print(data.keys())
    order = data['IRED_00127[S2]']
    print(order)


if __name__ == "__main__":
    main()
