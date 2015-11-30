from __future__ import print_function, absolute_import

__all__ = []

supported_dihedral_types = [
    x
    for x in
    'multidihedral phi psi chip omega alpha beta gamma delta epsilon zeta nu1 nu2 chin'.split(
    )
]

chinu = {
    'A': ("O4'", "C1'", "N9", "C4"),
    'G': ("O4'", "C1'", "N9", "C4"),
    'C': ("O4'", "C1'", "N1", "C2"),
    'T': ("O4'", "C1'", "N1", "C2"),
    'U': ("O4'", "C1'", "N1", "C2"),
}

template = '''
from .decorators import _register_pmap
from ._get_common_objects import _super_dispatch

@_register_pmap
@_super_dispatch()
def calc_{0}(traj=None, resrange="",
            range360=False,
            top=None, dtype='dataset',
            frame_indices=None):
    """
    Parameters
    ----------
    traj : Trajectory-like or anything that makes iterframe_master(traj) return Frame
    resrange : string or iterable, cpptraj resrange, default ""
        if resrange is a string, use index of 1-based
        if resrange is a python sequence (range, list, array, ...),
           use index of 0-based
    range360 : bool, default False
        use -180/180 or 0/360
    top : {str, Topology}, optional, default None
    frame_indices : {None, 1D array-like}, default None
        if not None, only compute for given indices

    Examples
    --------
    >>> import pytraj.dihedral_analysis as da
    >>> da.calc_phi(traj)
    >>> da.calc_psi(traj, resrange="3-10")
    >>> da.calc_chip(traj, resrange="3-10")
    >>> da.calc_chip(traj, resrange="3-10", range360=True)
    >>> da.calc_chip(traj, resrange="3-10", dtype='dict')
    >>> da.calc_multidihedral(traj, resrange="3-10")
    >>> da.calc_multidihedral(traj, resrange=range(5, 10))
    >>> # assert
    >>> from pytraj import common_actions as pyca
    >>> phi0 = pyca.calc_multidihedral(traj, "phi", dtype='dataset')
    >>> phi1 = da.calc_phi(traj, dtype='dataset')
    >>> from pytraj.testing import aa_eq
    >>> for key in phi0.keys():
    >>>     aa_eq(phi0[key], phi1[key])

    See Also
    --------
    """

    from .datasets.DatasetList import DatasetList
    from .actions.CpptrajActions import Action_MultiDihedral
    from ._get_common_objects import _get_data_from_dtype
    from .compat import string_types
    from .utils import is_int

    if range360:
        _range360 = 'range360'
    else:
        _range360 = ''

    if resrange:
        if is_int(resrange):
            resrange = [resrange,]

        if isinstance(resrange, string_types):
            _resrange = "resrange " + str(resrange)
        else:
            from pytraj.utils import convert as cv
            _resrange = "resrange " + cv.array_to_cpptraj_range(resrange)
    else:
        _resrange = ""

    dslist = DatasetList()
    template_command = '{1} '
    template_command = " ".join((template_command, _resrange, _range360))

    act = Action_MultiDihedral()
    # we dont need to set mass here
    act(template_command, traj, top=top, dslist=dslist)
    act.post_process()
    return _get_data_from_dtype(dslist, dtype=dtype)
'''

for key in supported_dihedral_types:
    if key != 'multidihedral':
        my_func_str = template.format(key, key)
    else:
        my_func_str = template.format(key, " ")
    g_dict = globals()
    exec(my_func_str)

del key
