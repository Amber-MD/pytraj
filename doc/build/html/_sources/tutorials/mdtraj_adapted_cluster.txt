.. _mdtraj_adapted

Tutorials adapted from ``mdtraj``
================================

.. content:: Table

Clustering with scipy
---------------------

original link: `mdtraj <http://mdtraj.org/latest/examples/clustering.html>`_

.. ipython:: python

    from __future__ import print_function
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
    # cluster for 'CA' atoms
    
    distances = pt.pairwise_rmsd(traj(autoimage=True), mask='@CA')

    # use `scipy` to perform clustering
    linkage = scipy.cluster.hierarchy.ward(distances)

    scipy.cluster.hierarchy.dendrogram(linkage, no_labels=True, count_sort='descendent')

    # cluster for all atoms but H
    
    distances = pt.pairwise_rmsd(traj(autoimage=True), mask='!@H=')
    linkage = scipy.cluster.hierarchy.ward(distances)
    scipy.cluster.hierarchy.dendrogram(linkage, no_labels=True, count_sort='descendent')
