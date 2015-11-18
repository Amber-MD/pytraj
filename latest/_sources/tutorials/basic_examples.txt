.. _basic_examples:

Basic examples
==============

try ``pytraj`` online:

.. image:: http://mybinder.org/badge.svg
   :target: http://mybinder.org/repo/hainm/notebook-pytraj

.. contents::

Load a Topology and Trajectory
------------------------------

.. ipython:: python

    import pytraj as pt
    traj = pt.load('tz2.nc', 'tz2.parm7')
    traj

Select atoms
------------

.. ipython:: python
    
    import pytraj as pt
    traj = pt.load('tz2.nc', 'tz2.parm7')
    # get indices for backbone H atoms
    h_indices = pt.select_atoms(traj.top, '@H')
    h_indices
    # get indices for backbone N atoms
    n_indices = pt.select_atoms(traj.top, '@N')
    n_indices

DSSP analysis
-------------

.. ipython:: python

    import pytraj as pt
    pdb = pt.load_pdb_rcsb("1l2y")
    result = pt.dssp(pdb)
    print(result)

Hbond analysis
--------------

.. ipython:: python
 
   # search hbonds for all residues
   hbonds = pt.search_hbonds(pdb)
   hbonds
   # print first few hbond
   hbonds.donor_aceptor[:5]
   hbonds.values

   # search hbonds between residue 9 and 16
   h = pt.search_hbonds(pdb, ':9,16')
   h.donor_aceptor

.. include:: load_pdb_rcsb.rst

More examples
-------------
`check our github page <https://github.com/Amber-MD/pytraj/tree/master/examples>`_
