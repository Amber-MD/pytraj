#!/usr/bin/env python
import pytraj as pt
from sklearn.cluster import k_means
from pytraj._fast_iterframe import _fast_iterptr
from memory_profiler import profile


@profile
def do_clustering():
    traj = pt.iterload('GAAC3.5000frames.nc', 'GAAC.topo')
    n_frames = 50

    # as an example, we extract only first 50 frames and do clustering for @P atoms
    t0 = traj[:n_frames, '@P']

    # need to convert to 2D array, shape=(n_frames, n_features)
    X = pt.tools.as_2darray(t0)

    # specify the total number of clusters
    n_clusters = 10

    # do it
    data = k_means(X, n_clusters,
                   init='k-means++',
                   precompute_distances='auto',
                   n_init=10,
                   max_iter=300,
                   verbose=False,
                   tol=0.0001,
                   random_state=None,
                   copy_x=False,
                   n_jobs=-1,
                   return_n_iter=False)

    x0 = data[0]

    # write output, dump all frames to a single pdb file, each snapshot is seperated by
    # 'MODEL' keyword. good for visualizing in VMD.

    with pt.Trajout('tes.pdb',
                    mode='model',
                    top=t0.top,
                    overwrite=True) as trajout:
        x = x0.reshape(n_clusters, t0.n_atoms, 3)
        # to avoid extra copying, we use _fast_iterptr (for advanced users)
        for idx, frame in enumerate(_fast_iterptr(x, t0.n_atoms,
                                                  range(n_clusters))):
            trajout.write_frame(idx, frame)

do_clustering()
