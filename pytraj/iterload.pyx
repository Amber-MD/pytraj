# distutils: language=c++
from libcpp.string cimport string
from pytraj.Topology cimport Topology
from pytraj.Trajin_Single cimport Trajin_Single

def _iterload(Topology top, traj, int start=0, chunk=None):
    """Iterately loading trajectory from file with provided topology file
    We use cpptraj "Trajin_Single" class to fast reading
    """
    #traj = traj.encode("UTF=8")
    ts = Trajin_Single()
    ts.load(traj, top)
    if chunk is None:
        return ts.__iter__()
    else:
        return ts.frame_iter(start, chunk)
