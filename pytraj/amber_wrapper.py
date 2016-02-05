import os
import pytraj as pt
from pytraj.utils.context import tempfolder
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

    with tempfolder():
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


leap_template = """
source leaprc.ff14SB
set default PBradii mbondi3
x = loadpdb {pdbfile}
saveamberparm x {pdbfile}.prmtop {pdbfile}.rst7
quit
"""


def prmtop_from_tleap(filename, leapin=None, verbose=False):
    '''make prmtop file from pdb

    Parameters
    ----------
    filename : str, pdb filename
    leapin : str, optional, default None
        leap input
        if None, use::

            source leaprc.ff14SB
            set default PBradii mbondi3
            x = loadpdb {pdbfile}
            saveamberparm x {pdbfile}.prmtop {pdbfile}.rst7
            quit
    verbose : bool, default False
        if False, suppress tleap output
    '''
    import os
    import subprocess
    import pytraj as pt

    if leapin is None:
        leapin = leap_template

    amberhome = os.environ.get('AMBERHOME')
    if amberhome is None:
        raise RuntimeError('must set AMBERHOME')

    tleap = amberhome + '/bin/tleap'

    filename = os.path.abspath(filename)

    with tempfolder():
        leapin = leapin.format(pdbfile=filename)

        with open("tmp_leap.in", 'w') as f:
            f.write(leapin)

        with open(os.devnull, 'wb') as devnull:
            if not verbose:
                subprocess.check_call(
                    [tleap, ' -f tmp_leap.in'],
                    stdout=devnull,
                    stderr=subprocess.STDOUT)
            else:
                x = subprocess.check_call([tleap, ' -f tmp_leap.in'])
        if verbose:
            print(x)
        return pt.load_topology("tmp.top")
