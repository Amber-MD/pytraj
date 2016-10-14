Atom mask selection
===================

try ``pytraj`` online:

.. include:: mybinder.rst

.. contents::

.. note:: Section describing ``Amber mask syntax`` is direct excerpt from `parmed's website <http://parmed.github.io/ParmEd/html/amber.html#amber-mask-syntax>`_

Amber mask syntax
-----------------

The Amber programs use a selection syntax very similar to Chimera's selection
syntax, referred to as an *Amber mask*. The character ``:`` is used for residue
selections and ``@`` is used for atom selections.  Both numbers and names are
supported in the atom and residue selections (where the serial numbers range
from 1 to *N*, where *N* is either the total number of atoms or residues,
depending on what is being selected).

The atom/residue name or number list can be specified as a comma-delimited list
of entries (*inclusive* number ranges are also supported).  You can also use
logical operators ``&`` (and), ``|`` (or), and ``!`` (not) and group expressions
with parenthesis to denote order-of-operations.

Residue selections
~~~~~~~~~~~~~~~~~~

Examples of residue selections, accompanied by comments annotating what is being
selected, is shown below::

    :1,3,6-10    Select residues 1, 3, 6, 7, 8, 9, and 10
    :ALA,LYS,10  Select all residues with names ALA or LYS, and residue 10
    :1-3,6-10    Select residues 1, 2, 3, 6, 7, 8, 9, and 10

Atom selections
~~~~~~~~~~~~~~~

Examples of residue selections, accompanied by comments annotating what is being
selected, is shown below::

    @1-100,150-160   Select the first 100 atoms as well as the 11 atoms 150-160
    @CA,C,10-15      Select all atoms whose names are CA or C, and 6 atoms 10-15

Atom type selections
~~~~~~~~~~~~~~~~~~~~

In some cases, it is more succinct to select by atom *type* name, rather than
atom name. Be careful though! Atom type names can change from force field to
force field, and are considered an implementation detail.  While atom *names*
should never change (and are based on the PDB standard, typically), atom *type*
names can be set and/or changed to whatever the force field developer wants.

You have been warned. If you wish to push on, though, you can use the ``@%``
characters to select from type names rather than atom names. Note, no indices
(or index ranges) can be used, as atom types have no defined order.  Examples
are shown below::

    @%CT,CX         Select all atoms whose types are CT or CX

Logical operators
~~~~~~~~~~~~~~~~~

You often want to select certain atoms from certain residues. You can use the
``&`` binary operator to do this, as shown in examples below::

    (:1-100)&(@CA)      Selects atoms named CA in the first 100 residues
    (:1-100)&(@CA,CB)   Selects atoms named CA or CB in the first 100 residues
    (:ASP,GLU)&(@CA)    Selects atoms named CA in all ASP or GLU residues

This is so common, you can use a simple shortcut, demonstrated below, giving
equivalent expressions to the examples above::

    :1-100@CA
    :1-100@CA,CB
    :ASP,GLU@CA

You can also expand your selection using the ``|`` (or) operator. For example::

    :1-100|@CA      Selects all atoms in the first 100 residues AND all atoms named CA

Don't get confused between what the ``&`` (and) and ``|`` (or) operators do.
The above example selects all atoms if it satisfies the criteria of being in the
first 100 residues *or* if it satisfies the criteria of having the name ``CA``,
meaning it selects all atoms in the first 100 residues *and* all atoms named
``CA``.

You can also select all atoms *except* a particular selection::

    !(:LIG)         Selects all atoms that are NOT in a residue named LIG

Wild-cards
~~~~~~~~~~

Let's suppose you want to select all hydrogen atoms. These atoms usually have
many names, but they almost always start with "H", so we would like a way to
select all atoms that start with the letter "H".  We can do this with
*wild-cards*, which is the ``=`` character in Amber masks, and works for both
atoms and residues. For example::

    :AS=@H=   Select all atoms whose residue starts with AS AND atom name starts with H

Distance-based selections
~~~~~~~~~~~~~~~~~~~~~~~~~

Sometimes you also want to be able to select all atoms within a certain distance
(or all atoms within a *residue* that is a certain distance) from a specific set
of atoms.  For this, the ``<`` and ``>`` operators are used to denote distances
that are less than or greater than a particular value, respectively.

This is probably the most complex part of the Amber mask syntax, but it is also
quite powerful.  In addition to the ``<`` and ``>`` operators, you must use
either ``@`` or ``:`` to indicate whether you want to select individual *atoms*
or whole *residues* satisfying the distance.  When selecting residues, the whole
residue is selected if *any atom* satisfies the distance criteria.  The general
syntax is ``SELECTION<:DISTANCE``, where ``SELECTION`` is an Amber mask (``:``
can also be ``@`` to select only atoms), and ``DISTANCE`` is a floating point
distance that defines the criteria.  Examples are shown below with annotations::

    :8<@5.0     All atoms that are LESS THAN 5 Angstroms away from any atom in residue 8
    :1>:10.0    All atoms in any residue in which one atom is MORE THAN 10 angstroms away from residue 1

Examples: Atom mask selection for ``Trajectory``
------------------------------------------------

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

Example: Atom mask selection for ``pytraj.Topology``
----------------------------------------------------

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
    top.set_distance_mask_reference(traj[0])
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
    
