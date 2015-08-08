Trajectory
============


``TrajectoryIterator`` is work-horse of ``pytraj``. This class offers out-of-core data store with easiness to load data to memory. 

.. code-block:: python

    >>> import pytraj as pt
    >>> traj = pt.iterload('./tz2.nc', './tz2.parm7')
    >>> traj
    <pytraj.TrajectoryIterator with 101 frames: 
    <Topology with 223 atoms, 13 residues, 1 mols, 230 bonds, non-PBC>>
    >>> pt.rmsd(traj, ref=0, mask='@CA')
    array([  1.94667955e-07,   2.54596866e+00,   4.22333034e+00, ...,
             4.97189564e+00,   5.53947712e+00,   4.83201237e+00])


pytraj is able to detect single file (mol2, pdb) to load as `TrajectoryIterator` too. ::
    >>> pt.iterload('my_pdb.pdb') 
    >>> pt.iterload('your_mol2.mol2') 
