Atom mask selection
===================

.. contents::

check `parmed's doc <http://parmed.github.io/ParmEd/html/amber.html#amber-mask-syntax>`_
for detail explanation.

**Examples adapted from parmed's website**
------------------------------------------
.. note:: the comments are from `parmed's doc`

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
