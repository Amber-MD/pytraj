import pytraj as pt
from pytraj.analysis import topology_analysis as ta

# local
from utils import fn


def test_info():
    traj = pt.load(fn('tz2.ortho.nc'), fn('tz2.ortho.parm7'))

    print(ta.atominfo(traj.top, '@10', ref=traj[0]))
    print(ta.resinfo(traj.top, ':10', ref=traj[0]))
    print(ta.bondinfo(traj.top, '@%CT', ref=traj[0]))
    print(ta.angleinfo(traj.top, '@2', ref=traj[0]))
    print(ta.dihedralinfo(traj.top, '@N @CA @CB @%H1', ref=traj[0]))
