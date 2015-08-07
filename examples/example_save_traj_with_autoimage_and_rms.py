import pytraj as pt

traj = pt.iterload('../tests//data/tz2.ortho.nc',
                   '../tests/data/tz2.ortho.parm7')
pt.write_traj("output/test.nc", traj(autoimage=True,
                                     rmsfit=(0, '@CA,N,O')),
              overwrite=True)
