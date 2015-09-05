Atom mask selection
===================

.. contents::

.. note:: check `parmed's doc <http://parmed.github.io/ParmEd/html/amber.html#amber-mask-syntax>`_ for detail explanation.

.. note:: Examples adapted from parmed's website and the comments are from `parmed's doc`

Atom mask selection for ``Trajectory``
--------------------------------------

.. ipython:: python

    import pytraj as pt
    traj = pt.iterload('data/tz2.ortho.nc', 'data/tz2.ortho.parm7')
    print(traj)

    traj.top.residue_names
    list(traj.top.atom_names)[:20]

    # Residue selections
    traj[':1,3,6-10']  # Select residues 1, 3, 6, 7, 8, 9, and 10
    traj[':ALA,LYS,10']  # Select all residues with names ALA or LYS, and residue 10

    # Atom selections
    traj['@1-100,150-160'] #  Select the first 100 atoms as well as the 11 atoms 150-160
    traj['@CA,C,10-15'] # Select all atoms whose names are CA or C, and 6 atoms 10-15

    # Logical operators
    traj['(:1-100)&(@CA)'] # Selects atoms named CA in the first 100 residues

Atom mask selection for ``pytraj.Topology``
-------------------------------------------

Get atom_indices
~~~~~~~~~~~~~~~~

.. ipython:: python
    
    top = traj.top
    # Topology can be loaded alone by:
    #top = pt.load_topology('data/tz2.ortho.parm7')
    top

    # Residue selections, return atom_indices (0-based)
    top.select(':1,3,6-10')
    top.select(':ALA,LYS,10')
    top.select('@1-100,150-160')
    top.select('(:1-100)&(@CA)')
    # select based on distance, need to set reference frame
    top.set_reference_frame(traj[0])
    # pick all atoms around residue 1, distance cutoff < 5.0 Angstrom
    top.select(':1 <:5.0')

Get new ``Topology``
~~~~~~~~~~~~~~~~~~~~

.. ipython:: python

    top[':1,3,6-10']
    top[':ALA,LYS,10']
    top['@1-100,150-160']
    top['(:1-100)&(@CA)']
    top[':1 <:5.0']
    
