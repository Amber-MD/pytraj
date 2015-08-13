Read and write
==============

(coming soon)

.. ipython:: python
    
    # read amber format
    traj = pt.iterload('data/tz2.ortho.nc', 'data/tz2.ortho.parm7')
    traj

    # write to CHARMM dcd format
    traj.save('test.dcd', overwrite=True)

    # write more than two trajs
    pt.write_traj('test2.dcd', [traj[:3], traj[5:]], top=traj.top)
