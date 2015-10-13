from __future__ import absolute_import
import numpy as np
from ._base_result_class import BaseAnalysisResult
from ._get_common_objects import _get_data_from_dtype, _get_topology, _get_fiterator
from .utils.convert import array_to_cpptraj_atommask as to_cpptraj_mask
from pytraj.compat import string_types
from pytraj import DatasetList


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
            return _to_string_secondary_structure(self.data.filter(
                lambda x: 'int' in x.dtype.name).to_dict())
        elif dtype == 'int':
            return self.data.filter(
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
            return _to_string_secondary_structure(self.data.filter(
                lambda x: 'int' in x.dtype.name).to_ndarray())
        elif dtype == 'int':
            return self.data.filter(
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
        return self.data.grep("avg")

    @property
    def residues(self):
        return np.array(self.data.grep('res', mode='aspect').keys())

    def values_per_frame(self, restype='string'):
        return np.vstack((self.residues, self.to_ndarray(restype).T))

    def values_per_residue(self, restype='string'):
        return np.vstack((self.residues, self.to_ndarray(restype).T)).T


def calc_dssp(traj=None, mask="", frame_indices=None, dtype='ndarray', top=None):
    """return dssp profile for frame/traj

    Parameters
    ----------
    traj : Trajectory-like
    mask: str
        atom mask
    frame_indices : {None, array-like}, default None, optional
        specify frame numbers for calculation.
        if None, do all frames
    dtype : str, default 'ndarray'
        return data type, for regular user, just use default one (ndarray).
        use dtype='dataset' if wanting to get secondary structure in integer format

    Returns
    -------
    out_0: ndarray, shape=(n_residues,)
        residue names
    out_1: ndarray, shape=(n_residues, n_frames)
        DSSP for each residue
    out_2 : pytraj.DatasetList
        average value for each secondary structure type

    Examples
    --------
    >>> import pytraj as pt
    >>> residues, ss, _ = pt.dssp(traj, ":2-10")
    >>> residues
    array(['LEU:2', 'TYR:3', 'ILE:4', 'GLN:5', 'TRP:6', 'LEU:7', 'LYS:8',
           'ASP:9', 'GLY:10'],
          dtype='<U6')
    >>> ss
    array([['0', '0', '0', ..., '0', '0', '0'],
           ['H', 'H', 'H', ..., 'H', 'H', 'H'],
           ['H', 'H', 'H', ..., 'H', 'H', 'H'],
           ...,
           ['H', 'H', 'H', ..., 'H', 'H', 'H'],
           ['T', 'T', 'T', ..., 'T', 'H', 'T'],
           ['0', '0', '0', ..., '0', '0', '0']],
          dtype='<U1')

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
    from pytraj.datasets.DatasetList import DatasetList as CpptrajDatasetList
    from pytraj.actions.CpptrajActions import Action_DSSP

    if not isinstance(mask, string_types):
        mask = to_cpptraj_mask(mask)

    command = mask

    _top = _get_topology(traj, top)
    fi = _get_fiterator(traj, frame_indices)
    dslist = CpptrajDatasetList()

    Action_DSSP()(command,
                  fi,
                  top=_top,
                  dslist=dslist)

    # replace key to something nicer
    for key, dset in dslist.iteritems():
        if 'DSSP' in key:
            key = key.replace("DSSP_00000[", "")
            key = key.replace("]", "_avg")
            dset.key = key.lower()
    dtype = dtype.lower()

    if dtype == 'ndarray':
        # get all dataset from DatSetList if dtype == integer
        arr0 = dslist.grep("integer", mode='dtype').values
        keys = dslist.grep("integer", mode='dtype').keys()
        avg_dict = DatasetList(dslist.grep('_avg'))
        return np.asarray(keys), np.asarray([_to_string_secondary_structure(arr) for arr
            in arr0]), avg_dict
    if dtype == '_dssp_class':
        return DSSPAnalysisResult(_get_data_from_dtype(dslist,
                                                       dtype='dataset'))
    else:
        return _get_data_from_dtype(dslist, dtype=dtype)


def _to_string_secondary_structure(arr0):
    """
    arr0 : {ndarray, dict of ndarray}
    """
    #ss = ['None', 'Para', 'Anti', '3-10', 'Alpha', 'Pi', 'Turn', 'Bend']
    ss = ["0", "b", "B", "G", "H", "I", "T", "S"]
    len_ss = len(ss)
    ssdict = dict(zip(range(len_ss), ss))

    myfunc = lambda key: ssdict[key]

    if not isinstance(arr0, dict):
        return np.vectorize(myfunc)(arr0)
    else:
        new_dict = {}
        for key in arr0.keys():
            new_dict[key] = _to_string_secondary_structure(arr0[key])
        return new_dict
