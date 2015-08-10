Atom mask selection
===================

check `parmed's doc <http://parmed.github.io/ParmEd/html/amber.html#amber-mask-syntax>`_

**Examples adapted from parmed's website**
------------

Load trajectory
.. code-block:: python

    >>> import pytraj as pt
    >>> traj = pt.iterload('./tz2.ortho.nc', './tz2.ortho.parm7')
    >>> traj
    <pytraj.TrajectoryIterator, 10 frames, include:
    <Topology: 5293 atoms, 1704 residues, 1692 mols, 5300 bonds, PBC with box type = ortho>>

*Residue selections*
.. code-block:: python

    >>> traj[':1,3,6-10'] # Select residues 1, 3, 6, 7, 8, 9, and 10
    <pytraj.Trajectory, 10 frames, include:
    <Topology: 108 atoms, 7 residues, 3 mols, 107 bonds, PBC with box type = ortho>>

    >>> traj[':ALA,LYS,10']  # Select all residues with names ALA or LYS, and residue 10
    <pytraj.Trajectory, 10 frames, include:
    <Topology: 58 atoms, 3 residues, 3 mols, 55 bonds, PBC with box type = ortho>>

