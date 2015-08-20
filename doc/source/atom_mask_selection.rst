Atom mask selection and trajectory slicing
==========================================

check `parmed's doc <http://parmed.github.io/ParmEd/html/amber.html#amber-mask-syntax>`_
for detail explanation.

**Examples adapted from parmed's website**
------------------------------------------
.. note:: the comments are from `parmed's doc`

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

**Trajecotry slicinge**
----------------------=

.. note:: `[ ]` slicing notation will load frames into memory

.. ipython:: python
    
    # get 1st frame
    traj[0]

    # get last frame
    traj[-1]

    # get frames 1, 3, 7 with only CA atoms
    traj[[1, 3, 7], '@CA']

    # skip every 2 frames
    traj[::2]

.. note:: `( )` notation to create frame iterator

.. ipython:: python

    traj(0, 8, 2)
    for frame in traj(0, 8, 2): print(frame)
    
