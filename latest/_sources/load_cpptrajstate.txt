.. _load_cpptrajstate:

.. contents::

Run cpptraj's batch
===================

This tutorial is for those who love running cpptraj in batch mode. Try it online:

.. image:: http://mybinder.org/badge.svg
   :target: http://mybinder.org/repo/hainm/notebook-pytraj

.. note:: This tutorial is for advanced users and the syntax might be changed in future.

.. ipython:: python
    :suppress:

    import numpy as np
    np.set_printoptions(precision=4, suppress=True)


if you want to run cpptraj's batch mode like below::

    parm tz2.ortho.parm7
    trajin tz2.ortho.nc
    autoimage
    rms first @CA
    radgyr @H,C,N,O
    molsurf 

you can create a 'CpptrajState'

Create state with existing ``TrajectoryIterator`` class
-------------------------------------------------------

.. ipython:: python

    # load traj first
    import pytraj as pt
    traj = pt.iterload('./tz2.ortho.nc', 'tz2.ortho.parm7')
    traj
    # create text (do not need to have ``parm`` and ``trajin`` info)
    text = '''
    autoimage
    rms first @CA
    radgyr @H,C,N,O
    molsurf @CA'''
    # load CpptrajState
    state = pt.load_batch(traj, text)
    # performa calculation
    state.run()
    state
    # get some data
    state.data
    state.data.keys()
    state.data[0]
    state.data[-1]
    # convert data to regular numpy array
    state.data.values


Create state without ``TrajectoryIterator`` class
-------------------------------------------------

.. ipython:: python
    
    # suppose you have 
    text = '''
    parm tz2.parm7
    trajin tz2.nc
    dihedral phi :1@C  :2@N  :2@CA :2@C
    dihedral psi   :1@N  :1@CA :1@C  :2@N
    dihedral omega   :1@CA :1@C  :2@N  :2@CA
    distance end_to_end :1@N :7@N
    rms first :1-1584
    strip !@H=
    createcrd mycrd
    '''
    import pytraj as pt
    state = pt.load_cpptraj_state(text)
    # need to explicit call run
    state.run()
    # All datasets are stored in ``state.data``
    state.data
    # if you already label your Dataset, you can access to them by using dict-like acessing
    # get Dataset with label `end_to_end` (distance)
    state.data['end_to_end']
    # get raw values (usually numpy array)
    state.data['end_to_end'].values

    # get `mycrd`
    state.data['mycrd']

    for dataset in state.data:
        print(dataset)
