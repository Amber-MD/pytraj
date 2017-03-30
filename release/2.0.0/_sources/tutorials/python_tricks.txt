.. _python_tricks:


.. contents::

Python tricks that make your life easier
========================================

.. ipython:: python
    
    import pytraj as pt
    traj = pt.load_sample_data('tz2')[:3]
    traj

Unpacking
---------

.. ipython:: python

    a, b, c = traj[0], traj[1], traj[2]
    a, b, c

Unpacking dictionary in python function
---------------------------------------

For new users, it's confusing when you see ``f(*args, **kwd)``. This is Python sytax to
unpack list/tuple and dictionary respectively

.. ipython:: python

     def f(x, y, z=0, a=100):
         print('x = {}'.format(x))
         print('y = {}'.format(y))
         print('z = {}'.format(z))
         print('a = {}'.format(a))

         return x + y + z + a

     my_list = [4, 6]
     my_dict = {'z': 8, 'a': 9}

     # call your function by unpacking list and dictionary
     f(*my_list, **my_dict)

     # which is equal to
     f(4, 6, z=8, a=9)

Iterating over Trajectory index and value pairs (enumerate)
-----------------------------------------------------------

.. ipython:: python

    for i, frame in enumerate(traj):
        print(i, frame)

Zipping
-------

.. ipython:: python
 
    h_indices = pt.select_atoms(traj.top, '@H')
    n_indices = h_indices - 1
    list(zip(h_indices, n_indices))[:-5]

Sliding windows
---------------

.. ipython:: python
 
    for value in pt.tools.n_grams([0, 2, 4, 5, 7, 9], 2):
        print(value)

Flattening lists
----------------

.. ipython:: python
 
    a = [[0, 54, 6], [[6, 7, 8], [9, 7, 5]]]
    print(pt.tools.flatten(a))

List comprehensions
-------------------

.. ipython:: python

    [x**2 for x in [3., 5., 7.]]

Dictionary comprehensions
-------------------------

.. ipython:: python

    a = ['one', 'two', 'three']
    b = [1, 2, 3]
    print(dict(zip(a, b)))

Get value from Dictionary if key does not exist
-----------------------------------------------

.. ipython:: python

    a = dict(x=3, y=4)
    a
    a.get('x')
    a.get('z', 100)

Set
---

.. ipython:: python

    set(res.name for res in traj.top.residues)

Filter data
-----------

.. ipython:: python

    a = [3, 8, 2, 9]
    list(filter(lambda x: x < 5, a))

Combinations 
------------

.. ipython:: python

    from itertools import combinations
    a = [3, 8, 2, 9]
    list(combinations(a, 3))

Chain several iterators
-----------------------

.. ipython:: python

    from itertools import chain
    traj = pt.iterload('tz2.nc', 'tz2.parm7')
    fi_0 = traj(0, 3)
    fi_0
    fi_1 = traj(5, 8)
    fi_1
    for frame in chain(fi_0, fi_1): print(frame)
