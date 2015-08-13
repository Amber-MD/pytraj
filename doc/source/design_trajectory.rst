Trajectory design
=================

.. currentmodule:: pytraj

**Why does pytraj have two Trajectory classes?**

Regular MD simulation usually produce GB to TB data. So loading all the data
to memory is not really a good choice. So we first design :class:`TrajectoryIterator`
as an out-of-core trajetory holder (same as `cpptraj`).

But sometimes we need to load a small chunk of data to memory for fast
calculation, for structure editting (rotate dihedral, ...), we need 2nd
:class:`Trajectory` class to hold the in-memory data.

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
coordinates don't apply to :class:`TrajectoryIterator`

.. code-block:: python
    
    pt.rmsd(traj, ref=0, mask='@CA')
    mem_traj = traj[[1, 3, 5], '@CA']
    pt.rmsd(mem_traj, ref=0) 

(to be continued)
