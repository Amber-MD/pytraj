import os
import pytraj as pt
from pytraj.utils.context import goto_temp_folder
from pytraj._get_common_objects import _get_top

MIN_IN = """
min.in
&cntrl
  imin = 1,
  maxcyc = 500,
  ncyc = 250,
  ntb = 0,
  igb = 8,
  cut = 999.0
/
"""

def minimize(traj, engine='sander', input=None, top=None):
    from pytraj import Trajectory

    if not isinstance(traj, Trajectory):
        raise ValueError("support only mutable Trajectory")

    _top = _get_top(traj, top)

    if input is not None:
        min_in = input
    else:
        min_in = MIN_IN

    if engine == 'sander' or engine == 'pmemd':
        _engine = "$AMBERHOME/bin/" + engine
    else:
        _engine = engine

    with goto_temp_folder():
        with open("min.in", 'w') as fh:
            fh.write(min_in)

        pt.write_parm("tmp.prmtop", _top)

        for frame in traj:
            pt.write_traj("tmp_frame.rst7", frame, top=_top, overwrite=True)
            os.system("%s -O -p tmp.prmtop -c tmp_frame.rst7.1 -r min.r -i min.in" % _engine)
            f0 = pt.load("min.r", traj.top)[0]
            # update coords
            frame.xyz[:] = f0.xyz
