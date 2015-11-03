#!/usr/bin/env python
import pytraj as pt
from sklearn.cluster import k_means
from pytraj import iterframe_from_array
import numpy as np


def do_clustering():
    # set random seed for reproducibility
    np.random.seed(1)

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

    # convert 2D to 3D array
    x0 = pt.tools.as_3darray(data[0])
    print('cluster index: ', data[1].tolist())

    # write output, dump all frames to a single pdb file, each snapshot is seperated by
    # 'MODEL' keyword. good for visualizing in VMD.
    # first, need to create a frame iterator for ``pytraj.write_traj``
    frame_iter = pt.iterframe_from_array(x0, t0.n_atoms, range(n_clusters))
    pt.write_traj('output.pdb', frame_iter,
                  top=t0.top,
                  options='model',
                  overwrite=True)


do_clustering()
