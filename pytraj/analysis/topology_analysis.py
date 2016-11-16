from ..core.c_core import CpptrajState, Command
from ..utils.context import capture_stdout

def compare_topology(top0, top1):
    ''' top0, top1 are :class:`pytraj.Topology`

    Notes
    -----
    Experiment 
    '''

    # TODO : make copies of topologies?
    state = CpptrajState()
    t0 = state.data.add('topology', name='top0')
    t0.data = top0

    t1 = state.data.add('topology', name='top1')
    t1.data = top1

    with capture_stdout() as (out, _):
        with Command() as command:
            command.dispatch(state, 'comparetop top0 top1')
    return out.read()
