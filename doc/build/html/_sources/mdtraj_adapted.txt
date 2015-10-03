.. _mdtraj_adapted

Tutorials adapted from ``mdtraj``
================================

.. contents:: Table

Clustering with scipy
---------------------

This example is adapted from `mdtraj <http://mdtraj.org/latest/examples/clustering.html>`_

.. ipython:: python

    from __future__ import print_function
    %matplotlib inline
    import matplotlib
    import pytraj as pt

    import scipy
    import matplotlib.pyplot as plt
    import scipy.cluster.hierarchy
    
    # load data

    traj = pt.iterload('data/tz2.ortho.nc', 'data/tz2.ortho.parm7')
    traj
    
    # calculate pairwise rmsd with `autoimage=True`
    # generally we only need to cluster a subset of atoms.
    # in this example, we chose 'CA' atoms
    
    distances = pt.pairwise_rmsd(traj(mask='@CA', autoimage=True))

    # use `scipy` to perform clustering
    linkage = scipy.cluster.hierarchy.ward(distances)

    scipy.cluster.hierarchy.dendrogram(linkage, no_labels=True, count_sort='descendent')
