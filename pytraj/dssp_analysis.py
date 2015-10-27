from __future__ import absolute_import
import numpy as np
from ._base_result_class import BaseAnalysisResult
from ._get_common_objects import _get_data_from_dtype, _get_topology, _get_fiterator
from .utils.convert import array_to_cpptraj_atommask as to_cpptraj_mask
from pytraj.compat import string_types
from pytraj import DatasetList
from .decorators import _register_openmp


@_register_openmp
def calc_dssp(traj=None,
              mask="",
              frame_indices=None,
              dtype='ndarray',
              simplified=False,
              top=None):
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
    simplified : bool, default False
        if True, use simplified codes, only has 'H', 'E' and 'C'
        if False, use all DSSP codes

    Returns
    -------
    out_0: ndarray, shape=(n_residues,)
        residue names
    out_1: ndarray, shape=(n_frames, n_residues)
        DSSP for each residue
    out_2 : pytraj.DatasetList
        average value for each secondary structure type

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.load_pdb_rcsb('1l2y')
    >>> residues, ss, _ = pt.dssp(traj, ":2-10")
    >>> residues
    array(['LEU:2', 'TYR:3', 'ILE:4', 'GLN:5', 'TRP:6', 'LEU:7', 'LYS:8',
           'ASP:9', 'GLY:10'], 
          dtype='<U6')
    >>> ss # doctest: +SKIP
    array([['0', 'H', 'H', ..., 'H', 'T', '0'],
           ['0', 'H', 'H', ..., 'H', 'T', '0'],
           ['0', 'H', 'H', ..., 'H', 'T', '0'],
           ...,
           ['0', 'H', 'H', ..., 'H', 'T', '0'],
           ['0', 'H', 'H', ..., 'H', 'H', '0'],
           ['0', 'H', 'H', ..., 'H', 'T', '0']],
          dtype='<U1')

    >>> residues, ss, _ = pt.dssp(traj, mask=range(100))

    >>> traj = pt.fetch_pdb('1l2y')
    >>> residues, ss, _ = pt.dssp(traj, simplified=True)
    >>> ss[0].tolist() # first frame
    ['C', 'H', 'H', 'H', 'H', 'H', 'H', 'H', 'C', 'C', 'H', 'H', 'H', 'H', 'C', 'C', 'C', 'C', 'C', 'C']


    Notes
    -----
    ========= ======= ========= =======================
    Character Integer DSSP_Char Seconday structure type
    ========= ======= ========= =======================
    0         0       '0'       None
    b         1       'E'       Parallel Beta-sheet
    B         2       'B'       Anti-parallel Beta-sheet
    G         3       'G'       3-10 helix
    H         4       'H'       Alpha helix
    I         5       'I'       Pi (3-14) helix
    T         6       'T'       Turn
    S         7       'S'       Bend
    ========= ======= ========= =======================

    Simlified codes::

        - 'H': include 'H', 'G', 'I' (helix)
        - 'E': include 'E', 'B' (strand)
        - 'C': include 'T', 'S' or '0' (coil)

    Simlified codes will be mostly used for visulization in other packages.
    """
    from pytraj.datasets.DatasetList import DatasetList as CpptrajDatasetList
    from pytraj.actions.CpptrajActions import Action_DSSP

    if not isinstance(mask, string_types):
        mask = to_cpptraj_mask(mask)

    command = mask

    _top = _get_topology(traj, top)
    fi = _get_fiterator(traj, frame_indices)
    dslist = CpptrajDatasetList()

    Action_DSSP()(command, fi, top=_top, dslist=dslist)

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
        ss_array = np.asarray([_to_string_secondary_structure(arr, simplified=simplified)
            for arr in arr0]).T
        return np.asarray(keys), ss_array, avg_dict
    else:
        return _get_data_from_dtype(dslist, dtype=dtype)


# _s0 = ['None', 'Para', 'Anti', '3-10', 'Alpha', 'Pi', 'Turn', 'Bend']
_s1 = ["0", "b", "B", "G", "H", "I", "T", "S"]

def _to_string_secondary_structure(arr0, simplified=False):
    """
    arr0 : ndarray
    """
    if not simplified:
        ssdict = {0: '0', 1: 'b', 2: 'B', 3: 'G', 4: 'H', 5: 'I', 6: 'T', 7: 'S'}
    else:
        ssdict = {0: 'C', 1: 'E', 2: 'E', 3: 'H', 4: 'H', 5: 'H', 6: 'C', 7: 'C'}

    return np.vectorize(lambda key: ssdict[key])(arr0)

def dssp_full_residues(traj, *args, **kwd):
    '''mostly for visulization. 
    Status: not finished yet
    '''
    residues, ss, _ = calc_dssp(traj, *args, **kwd)
    pass
