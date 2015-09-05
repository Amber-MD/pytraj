.. _faq:

Frequently Asked Questions
==========================

.. contents::

Why do most of examples use ``pytraj.iterload`` instead of ``pytraj.load``
--------------------------------------------------------------------------

Because we encourage user to use out-of-core calculation for memory saving. Most of
the time we don't need change trajecotry's coordinate, so it's OK to just use
immutable trajectory.

So does not ``pytraj`` support much for in-memory calculation?
--------------------------------------------------------------

Yes, ``pytraj`` does support in-memory calculation very well. Just consider
``pytraj.Trajectory`` as a mutable Python's list and ``pytraj.TrajectoryIterator`` as
immutable Python's tuple.

The difference between ``traj[1:11:2]`` and ``traj(1, 11, 2)``?
---------------------------------------------------------------

``traj[1:11:2]`` returns a new ``pytraj.Trajectory`` while ``traj(1, 11, 2)`` returns a
FrameIterator for lazy loading. You can use both for analysis::

    pt.radgyr(traj[1:11:2])
    pt.radgyr(traj(1, 11, 2))

It's too slow to convert Frame.xyz to numpy array, about 3 us for only 10 atom-Frame
------------------------------------------------------------------------------------

It's constant time to cast to numpy array. It takes similar time to cast from 1
million-atom Frame too. This is not problem since you only need cast once without making
extra data copying.

What's the idiom to load data to memory for specific atoms?
-----------------------------------------------------------

.. ipython:: python
    
    import pytraj as pt
    traj = pt.iterload('tz2.nc', 'tz2.parm7')
    traj
    traj[:8, '@CA']

List comprehensions give identical frames when iterating trajecotry
-------------------------------------------------------------------

This is our intention for memory saving. Use ``copy`` method.
(But what do you need list comprehensions for in this case?)

.. code-block:: python

    # same frames
    [frame for frame in traj]

    # different frames
    [frame.copy() for frame in traj]

Get undefined symbol error when import ``pytraj``
-------------------------------------------------

Probably you use old ``libcpptraj``. You can export LD_LIBRARY_PATH

.. code-block:: bash

    # supposed you have libcpptraj in /home/lib/cpptraj/, do
    $ export LD_LIBRARY_PATH=/home/lib/cpptraj/:$LD_LIBRARY_PATH

Error: Could not determine trajectory format ".nc"
--------------------------------------------------

``libcpptraj`` needs to be installed with NetCDF library. If you encouter this, send us
email or open issue on github

ImportError: libhdf5_hl.so.9: cannot open shared object file
------------------------------------------------------------

Something is wrong with ``conda``. Try ``conda install libnetcdf``.

ValueError: Big-endian buffer not supported on little-endian compiler
---------------------------------------------------------------------

When you are using memoryview, make sure to use correct type. Just google this error.
