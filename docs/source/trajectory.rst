Trajectory
============


TrajectoryIterator is work-horse of pytraj. This class offers out-of-core data store with easiness to load data to memory. 


.. toctree::
   :maxdepth: 1

   atom_selection

.. code-block:: python

    >>> import pytraj as pt
    >>> traj = pt.iterload('./tz2.nc', './tz2.parm7')
    >>> traj
    <pytraj.TrajectoryIterator with 101 frames: 
    <Topology with 223 atoms, 13 residues, 1 mols, 230 bonds, non-PBC>>
    >>> pt.rmsd(traj, ref=0, mask='@CA')
    array([  1.94667955e-07,   2.54596866e+00,   4.22333034e+00, ...,
             4.97189564e+00,   5.53947712e+00,   4.83201237e+00])


eytraj is able to detect single file (mol2, pdb) to load as TrajectoryIterator too.

.. code-block:: python

    >>> pt.iterload('my_pdb.pdb') 
    >>> pt.iterload('your_mol2.mol2') 


TrajectoryIterator is something like `range <https://docs.python.org/3/library/stdtypes.html#range>`_ in python3 or
`xrange <https://docs.python.org/2/library/functions.html#xrange>`_ in python2

.. code-block:: python

    >>> for i in range(0, 8, 2): print(i)
    ...
    0
    2
    4
    6
    >>> for f in traj(0, 8, 2): print(f)
    ...
    <Frame with 223 atoms>
    <Frame with 223 atoms>
    <Frame with 223 atoms>
    <Frame with 223 atoms>

However, TrajectoryIterator is much more than ``range`` or ``xrange``, you can slice atoms:

.. code-block:: python

    >>> for f in traj(0, 8, 2, mask='@CA'): print(f)
    ...
    <Frame with 12 atoms>
    <Frame with 12 atoms>
    <Frame with 12 atoms>
    <Frame with 12 atoms>

You can load a chunk of data into memory:

.. code-block:: python

    >>> traj[1:3, '!@H=']
    <pytraj.Trajectory, 2 frames, include:
    <Topology: 117 atoms, 13 residues, 1 mols, 124 bonds, non-PBC>>

How to perform analysis with TrajectoryIterator? It's very simple. For example, if you want to calculate
rmsd to 3rd frame (index starts from 0) for all atoms, just:

.. code-block:: python

    >>> pt.rmsd(traj, ref=3)
    array([ 5.87272014,  5.57581904,  4.35580747, ...,  8.17707783,
            8.69956761,  7.78286185])

You have TB of data and want to speed up your calculation, just add more cpus

.. code-block:: python

    >>> pt.pmap(n_cores=4, func=pt.radgyr, traj=traj)

How to get raw coordinates?

.. code-block:: python

    >>> traj.xyz
    >>> traj[[1, 3, 5]].xyz

