import numpy as np
from pytraj import CpptrajState
from pytraj.cast_dataset import cast_dataset
from pytraj import ArgList

# load prmtop file
def load_parm(fname=None):
    """Loading Amber toplogy file
    Parameters:
    ==========
    fname :: str (Amber topology filename)

    """
    top = Topology(top)
    return top
    
def create_state(top=None, trajin=None, ref=None):
    """Return CpptrajState instance

    Parameters:
    ========
    top :: topology filename (or topology instance or either?)
    trajin :: trajectory filename (or trajectory instance or either?)
    ref :: reference filename (or ReferenceFrame instance or either?)
    """
    state = CpptrajState()
    state.toplist.add_parm(top)
    state.add_trajin(trajin)

    if ref:
        state.add_reference(ref)
    return state

def do_calculation(action=None, command=None, state=None):
    """Perform specific actions with cpptraj state
    Return data list (should convert to Python array or numpy array?)

    Parameters:
    ========
    action :: instance of Action sublcass (e.g Action_Rmsd())
    input :: str
    state :: CpptrajState intance
    >>> state = create_state(top="./data/Tc5b.top", trajin="./data/md1_prod.Tc5b.x", ref=None)
    >>> distance = do_calculation(action=Action_Distance(), input="distance :2@CA :10@CA", state=state)
    """
    state.add_action(action, ArgList(command))
    state.run()
    d1 = cast_dataset(state.datasetlist[0], dtype="general")
    return np.asarray(d1.data)
