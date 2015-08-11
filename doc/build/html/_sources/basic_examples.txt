Basic examples
==============

DSSP analysis
-------------

.. ipython:: python

    import pytraj as pt
    pdb = pt.load_pdb_rcsb("1l2y")
    result = pt.dssp(pdb)
    print(result)

`More examples <https://github.com/Amber-MD/pytraj/tree/master/examples>`_
--------------------------------------------------------------------------
