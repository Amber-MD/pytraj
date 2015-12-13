import os
import pytraj as pt
from pytraj.utils.context import goto_temp_folder
from pytraj.get_common_objects import get_topology

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
    """
    >>> from pytraj.amber_wrap import minimize
    >>> minimize(traj)

    >>> minimize(traj, engine='pmemd')
    """

    from pytraj import Trajectory

    if not isinstance(traj, Trajectory):
        raise ValueError("support only mutable Trajectory")

    _top = get_topology(traj, top)

    if input is not None:
        min_in = input
    else:
        min_in = MIN_IN

    if engine in ['sander', 'pmemd']:
        _engine = "$AMBERHOME/bin/" + engine
    else:
        _engine = engine

    with goto_temp_folder():
        with open("min.in", 'w') as fh:
            fh.write(min_in)

        pt.write_parm("tmp.prmtop", _top)

        for frame in traj:
            pt.write_traj("tmp_frame.rst7", frame, top=_top, overwrite=True)
            os.system(
                "%s -O -p tmp.prmtop -c tmp_frame.rst7.1 -r min.r -i min.in" %
                (_engine))
            f0 = pt.load("min.r", traj.top)[0]
            # update coords
            frame.xyz[:] = f0.xyz


leapin = """
source leaprc.ff14SB
set default PBradii mbondi3
x = loadpdb %s
saveamberparm x tmp.top tmp.crd
quit
"""


def prmtop_from_tleap(fname, leapin=leapin, verbose=False):
    import os
    import subprocess
    import pytraj as pt

    try:
        amberhome = os.environ['AMBERHOME']
    except KeyError:
        raise KeyError("must set AMBERHOME")

    tleap = amberhome + '/bin/tleap'

    fname = os.path.abspath(fname)

    with goto_temp_folder():
        leapin = leapin % fname

        with open("_leap.in", 'w') as f:
            f.write(leapin)

        with open(os.devnull, 'wb') as devnull:
            if not verbose:
                subprocess.check_call(
                    [tleap, ' -f _leap.in'],
                    stdout=devnull,
                    stderr=subprocess.STDOUT)
            else:
                subprocess.check_call([tleap, ' -f _leap.in'])
        return pt.load_topology("tmp.top")
