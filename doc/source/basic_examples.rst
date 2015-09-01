.. _basic_examples:

Basic examples
==============

.. contents::

DSSP analysis
-------------

.. ipython:: python

    import pytraj as pt
    pdb = pt.load_pdb_rcsb("1l2y")
    result = pt.dssp(pdb)
    print(result)

hbond analysis
--------------

.. ipython:: python
 
   hbonds = pt.search_hbonds(pdb)
   hbonds

More examples
-------------
`check our github page <https://github.com/Amber-MD/pytraj/tree/master/examples>`_
