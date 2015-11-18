.. _design_trajectory:

Trajectory design
=================

try ``pytraj`` online:

.. image:: http://mybinder.org/badge.svg
   :target: http://mybinder.org/repo/hainm/notebook-pytraj

.. ipython:: python
    :suppress:

    import pytraj as pt
    traj = pt.iterload('data/tz2.nc', 'data/tz2.parm7')

.. currentmodule:: pytraj

Why does pytraj have two Trajectory classes?
--------------------------------------------

Regular MD simulation usually produce GB to TB data. So loading all the data
to memory is not really a good choice. So we first design `TrajectoryIterator`
as an out-of-core trajetory holder (same as `cpptraj`).

But sometimes we need to load a small chunk of data to memory for fast
calculation, for structure editting (rotate dihedral, ...), we need 2nd
`Trajectory` class to hold the in-memory data.

.. code-block:: python 

    >>> import pytraj as pt
    >>> traj = pt.iterload('./tz2.ortho.nc', './tz2.ortho.parm7') # out-of-core
    >>> traj
    <pytraj.TrajectoryIterator, 10 frames, include:
    <Topology: 5293 atoms, 1704 residues, 1692 mols, PBC with box type = ortho>

    >>> traj[[1, 3, 5], '@CA'] # load frame 1, 3, 5 with only CA atoms to memory
    <pytraj.Trajectory, 3 frames, include:
    <Topology: 12 atoms, 12 residues, 12 mols, PBC with box type = ortho>>

You can perform most of analysis with those two. Note that any actions modifying
coordinates don't apply to `TrajectoryIterator`

.. code-block:: python
    
    pt.rmsd(traj, ref=0, mask='@CA')
    mem_traj = traj[[1, 3, 5], '@CA']
    pt.rmsd(mem_traj, ref=0) 

Advantages of TrajectoryIterator
--------------------------------

* immutable (like python's tuple)
* use very little RAM
* still have random access to frame/frames
* easy to `parallelize <parallel>`_

.. ipython:: python
    
    traj
    traj[0]
    traj[[0, 3, 88]]
    traj[0:40:2]

(to be continued)
