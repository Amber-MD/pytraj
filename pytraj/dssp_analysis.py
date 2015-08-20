from __future__ import absolute_import
from ._base_result_class import BaseAnalysisResult
from ._get_common_objects import _get_data_from_dtype, _get_top
from .utils import _import_numpy
from .utils.convert import array_to_cpptraj_atommask as to_cpptraj_mask
from pytraj.compat import string_types

_, np = _import_numpy()


class DSSPAnalysisResult(BaseAnalysisResult):
    """
    Notes
    -----
    class's name might be changed
    """

    def to_dict(self, dtype='int'):
        """
        Return a dict of numpy.ndarray

        Parameters
        ----------
        dtype : str, {'int', 'string'}
        """
        if dtype == 'string':
            return to_string_ss(self.dslist.filter(
                lambda x: 'int' in x.dtype.name).to_dict())
        elif dtype == 'int':
            return self.dslist.filter(
                lambda x: 'int' in x.dtype.name).to_dict()
        else:
            raise NotImplementedError()

    def to_ndarray(self, dtype='string'):
        """
        Return a numpy.ndarray

        Parameters:
        dtype : str, {'string', 'int'}
        """
        if dtype == 'string':
            return to_string_ss(self.dslist.filter(
                lambda x: 'int' in x.dtype.name).to_ndarray())
        elif dtype == 'int':
            return self.dslist.filter(
                lambda x: 'int' in x.dtype.name).to_ndarray()
        else:
            raise NotImplementedError()

    def to_ndarray_per_frame(self, dtype='string'):
        return self.to_ndarray(dtype).T

    def average(self):
        """
        Return a `pytraj.datasetlist.DatasetList` object having average value
        for each frame for each type of secondary structure
        """
        return self.dslist.grep("avg")

    @property
    def residues(self):
        return np.array(self.dslist.grep('res', mode='aspect').keys())

    def values_per_frame(self, restype='string'):
        return np.vstack((self.residues, self.to_ndarray(restype).T))

    def values_per_residue(self, restype='string'):
        return np.vstack((self.residues, self.to_ndarray(restype).T)).T


def calc_dssp(traj=None, mask="", top=None, dtype='ndarray', *args, **kwd):
    """return dssp profile for frame/traj

    Parameters
    ----------
    mask: str
    traj : {Trajectory, Frame, mix of them}
    dtype : str {'ndarray', 'dataset', 'dict', 'dataframe', '_dssp_class'}, default 'ndarray'
        return data type, for regular user, just use default one (ndarray)

    Returns
    -------
    numpy.ndarray, shape=(n_residues + 1,) (array of residue names and DSSP characters)

    Examples
    --------
    >>> import pytraj as pt
    >>> d = pt.dssp(traj, ":2-10", dtype='ndarray')
    >>> print(d)
    >>> pt.dssp(traj, ":2-10", dtype='dict')
    >>> pt.dssp(traj, ":2-10", dtype='dataframe')
    >>> pt.dssp(traj, ":2-10", dtype='dataset')

    Notes
    -----
    ========= ======= ========= =======================
    Character Integer DSSP_Char SS_type
    ========= ======= ========= =======================
    0         0       ' '       None
    b         1       'E'       Parallel Beta-sheet
    B         2       'B'       Anti-parallel Beta-sheet
    G         3       'G'       3-10 helix
    H         4       'H'       Alpha helix
    I         5       'I'       Pi (3-14) helix
    T         6       'T'       Turn
    S         7       'S'       Bend
    ========= ======= ========= =======================
    """
    from pytraj.datasets.DataSetList import DataSetList as CpptrajDatasetList
    from pytraj.actions.CpptrajActions import Action_DSSP
    _, np = _import_numpy()

    if not isinstance(mask, string_types):
        mask = to_cpptraj_mask(mask)

    command = mask

    _top = _get_top(traj, top)
    dslist = CpptrajDatasetList()

    Action_DSSP()(command,
                  current_frame=traj,
                  top=_top,
                  dslist=dslist, *args, **kwd)

    # replace legend to something nicer
    for legend, dset in dslist.iteritems():
        if 'DSSP' in legend:
            legend = legend.replace("DSSP_00000[", "")
            legend = legend.replace("]", "_avg")
            dset.legend = legend.lower()
    dtype = dtype.lower()

    if dtype == 'ndarray':
        # get all dataset from DatSetList if dtype == integer
        arr0 = dslist.grep("integer", mode='dtype').values
        keys = dslist.grep("integer", mode='dtype').keys()
        return np.array([keys, [to_string_ss(arr) for arr in arr0]],
                        dtype='object')
    if dtype == '_dssp_class':
        return DSSPAnalysisResult(_get_data_from_dtype(dslist,
                                                       dtype='dataset'))
    else:
        return _get_data_from_dtype(dslist, dtype=dtype)


def to_string_ss(arr0):
    """
    arr0 : {ndarray, dict of ndarray}
    """
    _, np = _import_numpy()
    #ss = ['None', 'Para', 'Anti', '3-10', 'Alpha', 'Pi', 'Turn', 'Bend']
    ss = ["0", "b", "B", "G", "H", "I", "T", "S"]
    len_ss = len(ss)
    ssdict = dict(zip(range(len_ss), ss))

    if np:

        def myfunc(key):
            return ssdict[key]

        if not isinstance(arr0, dict):
            return np.vectorize(myfunc)(arr0)
        else:
            new_dict = {}
            for key in arr0.keys():
                new_dict[key] = to_string_ss(arr0[key])
            return new_dict
    else:
        raise ImportError("require numpy")
