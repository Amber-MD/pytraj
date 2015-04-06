"""Load mdtraj traj object

"""
from pytraj.utils import has_
from pytraj import FrameArray, Topology, Atom, Frame

# not sure if we need this `load_mdtraj` since cpptraj can do anything :D
# might need to move to Cython level for faster loading

if not has_("mdtraj") and not has_("numpy"):
    # mdtraj need numpy
    print ("need to have mdtraj and numpy")

def load_mdtraj(m_traj):
    """load_mdtraj traj object

    Parameters
    ---------
    m_traj : Trajectory object from mdtraj 
    """
    farray = FrameArray()

    # convert to pseudo-topology
    # to fully use Topology object in pytraj, we can do:
    # >>> farray.top = Topology(top_name) # or
    # >>> from pytraj import io
    # >>> farray.top = io.load(top_name) 

    pseudotop = Topology()
    for m_atom in m_traj.top.atoms:
        atom = Atom(m_atom.name, m_atom.name)
        m_res = m_atom.residue
        pseudotop.add_atom(atom=atom, resid=m_res.index, resname=m_res.name)

    # load coords
    farray = FrameArray(m_traj.xyz, pseudotop)
    return farray
