.. _load_cpptrajstate:

Run cpptraj's batch
===================

if you want to run cpptraj's batch mode like below::

    parm tz2.ortho.parm7
    trajin tz2.ortho.nc
    autoimage
    rms first @CA
    radgyr @H,C,N,O
    molsurf 

you can create a 'CpptrajState'

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
    text
    # load CpptrajState
    state = pt.load_batch(traj, text)
    # performa calculation
    state.run()
    state
    # get some data
    state.datasetlist
    state.datasetlist.keys()
    state.datasetlist[0]
    state.datasetlist[-1]
