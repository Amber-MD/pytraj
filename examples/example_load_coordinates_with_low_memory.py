import unittest # for travis (you don't need to add this)

# load TrajectoryIterator for saving memory
import pytraj as pt
trajiter = pt.iterload('../tests/data/tz2.ortho.nc', '../tests/data/tz2.ortho.parm7')
# load coordinates to memory
xyz = pt.get_coordinates(trajiter(mask='@CA', autoimage=True, rmsfit=(0, '!@H=')))
