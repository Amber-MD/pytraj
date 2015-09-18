from __future__ import print_function, absolute_import
import numpy as np
from .externals.six import iteritems

supported_types = [
    x for x in
    'minimage dipole center corrplane box boxcenter ucellx ucelly ucellz principal'.split()]

def _2darray_to_atommask_groups(seq):
    '''[[0, 3], [4, 7]] turns to ['@1 @4', '@5 @8']
    '''
    for arr in seq:
        # example: arr = [0, 3]; turns ot '@1 @4'
        yield '@' + str(arr[0]+1) + ' @' +  str(arr[1]+1)


def vector_mask(traj=None, mask="", frame_indices=None, dtype='ndarray', top=None):
    """compute vector between two maskes

    Parameters
    ----------
    traj : Trajectory-like or iterable that produces Frame
    mask: str or array of string or array of intergers, shape (n_vectors, 2)
        vector maskes
    frame_indices : array-like or iterable that produces integer number
        frame indices
    dtype : str, default 'ndarray'
        output dtype 
    top : Topology, optional, default None

    Rerturns
    --------
    if mask is a string, return 2D ndarray, shape (n_frames, 3)
    if mask is a list of strings or a 2D ndarray, return 3D ndarray, shape (n_vectors, n_frames, 3)

    Examples
    --------
    >>> # calcualte N-H vector
    >>> from pytraj import vector_analysis as va
    >>> import numpy as np
    >>> n_indices = pt.select_atoms(traj.top, '@N')
    >>> h_indices = n_indices + 1

    >>> # create n-h pair for vector calculation
    >>> n_h_pairs = np.array(list(zip(n_indices, h_indices)))
    >>> data_vec = va.vector_mask(traj, n_h_pairs, dtype='ndarray')

    >>> # compute vectors for specific frame indices (0, 4)
    >>> data_vec = va.vector_mask(traj, n_h_pairs, frame_indices=[0, 4], dtype='ndarray')
    """
    from ._get_common_objects import _get_top, _get_data_from_dtype, _get_fiterator
    from ._get_common_objects import _get_list_of_commands
    from .datasets.DatasetList import DatasetList as CpptrajDatasetList
    from .actions.CpptrajActions import Action_Vector
    from .core.ActionList import ActionList

    fi = _get_fiterator(traj, frame_indices)
    _top = _get_top(fi, top)
    dslist = CpptrajDatasetList()
    template_command = ' mask '

    cm_arr = np.asarray(mask)
    if cm_arr.dtype.kind != 'i':
        list_of_commands = _get_list_of_commands(mask)
    else:
        if cm_arr.ndim != 2:
            raise ValueError('if mask is a numpy.ndarray, it must have ndim = 2')
        list_of_commands = _2darray_to_atommask_groups(cm_arr)

    actlist = ActionList()

    for command in list_of_commands:
        act = Action_Vector()
        _command = command + template_command
        actlist.add_action(act, _command, _top, dslist=dslist)
    actlist.do_actions(fi)
    return _get_data_from_dtype(dslist, dtype=dtype)

template = '''
def vector_%s(traj=None, command="", frame_indices=None, dtype='ndarray', top=None):
    """
    Parameters
    ----------
    traj : Trajectory-like
    command : cpptraj command
    top : {str, Topology}, optional, default None
    *args, **kwd: more arguments
    """
    from ._get_common_objects import _get_top, _get_data_from_dtype, _get_fiterator
    from ._get_common_objects import _get_list_of_commands
    from .datasets.DatasetList import DatasetList as CpptrajDatasetList
    from .actions.CpptrajActions import Action_Vector
    from .core.ActionList import ActionList

    fi = _get_fiterator(traj, frame_indices)
    _top = _get_top(fi, top)
    dslist = CpptrajDatasetList()
    template_command = ' %s '

    list_of_commands = _get_list_of_commands(command)
    actlist = ActionList()

    for command in list_of_commands:
        act = Action_Vector()
        _command = command + template_command
        actlist.add_action(act, _command, _top, dslist=dslist)
    actlist.do_actions(fi)
    return _get_data_from_dtype(dslist, dtype=dtype)
'''

for key in supported_types:
    my_func_str = template % (key, key)
    g_dict = globals()
    exec(my_func_str)

del key
