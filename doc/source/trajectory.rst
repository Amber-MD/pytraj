Trajectory
==========

try ``pytraj`` online:

.. image:: http://mybinder.org/badge.svg
   :target: http://mybinder.org/repo/hainm/notebook-pytraj

.. contents::

.. currentmodule:: pytraj

.. ipython:: python
    :suppress:

    import numpy as np
    np.set_printoptions(precision=4, suppress=True)


Overview
========

There are two types of trajecoties in ``pytraj``.

- an immutable :class:`pytraj.TrajectoryIterator` (work-horse of pytraj).
This class offers out-of-core data store with easiness to load data to memory. 

- a mutable :class:`pytraj.Trajectory`. This class hold in-memory data. It is suitable for loading a chunk of snapshots.

Why are there two trajectory classes?
-------------------------------------

See our rationale here :ref:`design_trajectory`

Try it
------

.. ipython:: python

    import pytraj as pt
    traj = pt.iterload('data/tz2.ortho.nc', 'data/tz2.ortho.parm7')
    traj
    pt.rmsd(traj, ref=0, mask='@CA')

pytraj is able to detect single file (mol2, pdb) to load as TrajectoryIterator.

::

    >>> pt.iterload('my_pdb.pdb') 
    >>> pt.iterload('your_mol2.mol2') 


TrajectoryIterator is something like `range(start, stop[, step]) <https://docs.python.org/3/library/stdtypes.html#range>`_ in python3 or
`xrange(start, stop[, step]) <https://docs.python.org/2/library/functions.html#xrange>`_ in python2

.. ipython:: python

    for i in range(0, 8, 2): print(i)
    for f in traj(0, 8, 2): print(f)

However, TrajectoryIterator is much more than that, you can slice atoms:

.. ipython:: python

    for f in traj(0, 8, 2, mask='@CA'): print(f, f.xyz[0])

To load all frames to memory at once, use ``[]`` notation:

.. ipython:: python
    
    traj[0:8:2, '@CA']

How to perform analysis with TrajectoryIterator? It's very simple. For example, if you want to calculate
rmsd to 3rd frame (index starts from 0) for all atoms, just:

.. ipython:: python

    pt.rmsd(traj, ref=3)

You have TB of data and want to speed up your calculation, just add more cpus

.. ipython:: python

    pt.pmap(pt.rmsd, traj, ref=3, n_cores=4)

How to get raw coordinates?

.. code-block:: python

    >>> traj.xyz
    >>> traj[[1, 3, 5]].xyz

See also
--------

`Fancy indexing <trajectory_slice>`_
