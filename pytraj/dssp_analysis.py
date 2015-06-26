from . _base_result_class import BaseAnalysisResult
from ._get_common_objects import _get_data_from_dtype, _get_top
from .utils import _import_numpy


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
            return to_string_ss(self.dslist.grep("int", mode='dtype').to_dict())
        elif dtype == 'int':
            return self.dslist.grep("int", mode='dtype').to_dict()
        else:
            raise NotImplementedError()

    def to_ndarray(self, dtype='string'):
        """
        Return a numpy.ndarray

        Parameters:
        dtype : str, {'string', 'int'}
        """
        if dtype == 'string':
            return to_string_ss(self.dslist.grep("int", mode='dtype').to_ndarray())
        else:
            return self.dslist.grep("int", mode='dtype').to_ndarray()

    def average(self):
        """
        Return a `pytraj.datasetlist.DatasetList` object having average value
        for each frame for each type of secondary structure
        """
        return self.dslist.grep("avg")


def calc_dssp(traj=None, command="", top=None, dtype='ndarray', *args, **kwd):
    """return dssp profile for frame/traj

    Parameters
    ----------
    command : str
    traj : {Trajectory, Frame, mix of them}
    dtype : str {'dataset', 'ndarray', 'dict', 'dataframe'}, default 'ndarray'

    Returns
    -------
    (try it)

    Examples
    --------
        calc_dssp(traj, ":2-10", dtype='ndarray')
        calc_dssp(traj, ":2-10", dtype='dict')
        calc_dssp(traj, ":2-10", dtype='dataframe')
        calc_dssp(traj, ":2-10", dtype='dataset')

    Notes
    -----
    Character Integer DSSP_Char SS_type
    0         0       ' '       None
    b         1       'E'       Parallel Beta-sheet
    B         2       'B'       Anti-parallel Beta-sheet
    G         3       'G'       3-10 helix
    H         4       'H'       Alpha helix
    I         5       'I'       Pi (3-14) helix
    T         6       'T'       Turn
    S         7       'S'       Bend

    See Also
    --------
    Amber15 manual: http://ambermd.org/doc12/Amber15.pdf (page 588)
    """
    from pytraj.datasets.DataSetList import DataSetList as CpptrajDatasetList
    from pytraj.actions.CpptrajActions import Action_DSSP
    _, np = _import_numpy()

    _top = _get_top(traj, top)
    dslist = CpptrajDatasetList()

    Action_DSSP()(command,
                  current_frame=traj, 
                  top=_top,
                  dslist=dslist,
                  *args, **kwd)

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
        return np.array([to_string_ss(arr) for arr in arr0])
    if dtype == 'dssp_class':
        return DSSPAnalysisResult(_get_data_from_dtype(dslist, dtype='dataset'))
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
