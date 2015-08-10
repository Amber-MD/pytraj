Parallel calculation
--------------------

.. code-block:: python

    In [2]: pt.pmap(n_cores=4, func=pt.radgyr, traj=traj)
    Out[2]:
    [(0, array([ 18.91114428,  18.93654996])),
     (1, array([ 18.84969884,  18.90449256])),
     (2, array([ 18.8568644 ,  18.88917208])),
     (3, array([ 18.9430491 ,  18.88878079,  18.91669565,  18.87069722]))]
