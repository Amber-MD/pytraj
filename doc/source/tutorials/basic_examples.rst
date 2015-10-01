.. _basic_examples:

Basic examples
==============

try ``pytraj`` online:

.. image:: http://mybinder.org/badge.svg
   :target: http://mybinder.org/repo/hainm/notebook-pytraj

.. contents::

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
 
   hbonds = pt.search_hbonds(pdb)
   hbonds

More examples
-------------
`check our github page <https://github.com/Amber-MD/pytraj/tree/master/examples>`_
