from __future__ import absolute_import, print_function, division

from .action_dict import ActionDict
from .externals.six import string_types
from .datasets.DataSetList import DataSetList
from ._get_common_objects import _get_data_from_dtype
from ._base_result_class import BaseAnalysisResult

adict = ActionDict()

__all__ = ['HbondAnalysisResult', 'search_hbonds', 'search_nointramol_hbonds',
           'search_hbonds_noseries']


class HbondAnalysisResult(BaseAnalysisResult):
    """Hold data for HbondAnalysisResult

    Examples
    --------
    >>> import pytraj as pt
    >>> traj = pt.load_pdb_rcsb("1l2y")
    >>> h = pt.hbonds.HbondAnalysisResult(traj.search_hbonds())
    >>> h
    <pytraj.hbonds.HbondAnalysisResult
    donor_aceptor pairs : 31>
    >>> 
    >>> h.lifetime(cut=(0.5, 1.0))
    {'ARG16_O-TRP6_NE1-HE1': 0.9473684210526315,
     'ILE4_O-LYS8_N-H': 0.7631578947368421,
     'LEU7_O-GLY10_N-H': 0.8947368421052632,
     'TYR3_O-LEU7_N-H': 1.0}
    >>>
    >>> h.grep(['ARG', 'TYR']).dslist.to_dict()
    """

    def __str__(self):
        root_msg = "<pytraj.hbonds.HbondAnalysisResult"
        more_info = "donor_aceptor pairs : %s>" % len(self.donor_aceptor)
        return root_msg + "\n" + more_info

    def __repr__(self):
        return self.__str__()

    @property
    def donor_aceptor(self):
        return self.dslist.grep(["solventhb", "solutehb"],
                                mode='aspect').keys()

    def lifetime(self, cut=None):
        """return a dict with keys as donor_aceptor pairs and
        values as fraction of frames (vs total) having hbond. This fraction
        within `cut`. If `cut` is a single number, cut=(cut, 1.0). If `cut`
        is a tuple, cut=(min, max)
        """
        c = self.dslist.count(1)
        n_frames = self.dslist[0].size

        result_dict = dict((key, c[key] / n_frames)
                           for key in self.donor_aceptor)

        if cut is None:
            return result_dict
        else:

            def func(result_dict, cut=cut):
                d = {}
                for k in result_dict.keys():
                    if isinstance(cut, tuple):
                        if cut[0] <= result_dict[k] <= cut[1]:
                            d[k] = result_dict[k]
                    else:
                        if result_dict[k] >= cut:
                            d[k] = result_dict[k]
                return d

            return func(result_dict, cut=cut)

    def to_amber_mask(self):
        """convert donor_aceptor pair mask to amber mask to calculate
        distance (for example: 'ARG16_O-TRP6_NE1-HE1' will be ':16@O :6@HE1')
        """
        raise NotImplementedError("not yet")


def _update_legend_hbond(_dslist):

    # SER_20@O-SER_20@OG-HG --> SER20_O-SER20_OG-HG
    for d0 in _dslist:
        d0.legend = d0.legend.replace("_", "")
        d0.legend = d0.legend.replace("@", "_")

    for d0 in _dslist:
        if d0.legend == 'HB00000[UU]':
            d0.legend = 'total_solute_solute'


def search_hbonds_noseries(traj,
                           mask="",
                           dtype='dataset',
                           update_legend=True, *args, **kwd):
    """search hbonds for a given mask
    Parameters
    ----------
    traj : {Trajectory-like object, frame_iter object, list of traj}
    mask : str 
        Amber atom mask
    dtype : str {'list', 'pyarray', 'dataset', 'ndarray'}, default='dataset'
    *args, **kwd: optional

    Returns
    -------
    out : DataSetList | pyarray | ndarray | list | dict (depend on 'dtype')

    http://ambermd.org/doc12/Amber15.pdf (page 575)
    """

    dslist = DataSetList()
    act = adict['hbond']

    command = mask
    if "series" in command:
        raise ValueError("don't accept key `series`")
    act(command, traj, dslist=dslist, *args, **kwd)
    act.print_output()

    if update_legend:
        _update_legend_hbond(dslist)

    if dtype == 'dataframe':
        # return DataFrame.T to have better visual effect
        return dslist.to_dataframe().T
    else:
        return _get_data_from_dtype(dslist, dtype=dtype)


def search_hbonds(traj,
                  mask="",
                  dtype='dataset',
                  solventdonor=None,
                  solventacceptor=None,
                  update_legend=True, *args, **kwd):
    """search hbonds for a given mask
    Parameters
    ----------
    traj : {Trajectory-like object, frame_iter object, list of traj}
    mask : str 
        Amber atom mask
    dtype : str {'list', 'pyarray', 'dataset', 'ndarray'}, default='dataset'
    *args, **kwd: optional

    Returns
    -------
    out : DataSetList | pyarray | ndarray | list | dict (depend on 'dtype')

    Examples
    --------
    * The syntax was adapted from http://ambermd.org/doc12/Amber15.pdf (page 575)
    * The explaniation in " " is direct excerpt from this manual

    * "search for all hydrogen bonds within residues 1-22"
        dslist = search_hbonds(traj, ":1-22")

    * "search for all hydrogen bonds within residues 1-22, specifying output files"

        dslist = search_hbonds(traj, ":1-22 out nhb.dat avgout avghb.dat", dflist=dflist)
        dflist.write_all_datafile()

    * "search for all hydrogen bonds formed between donors in residue 1 and acceptors in residue 2" 

        dslist = search_hbonds(traj, "donormask :1 acceptormask :2", dtype='ndarray'))

    See Also
    --------
    http://ambermd.org/doc12/Amber15.pdf (page 575)
    """
    s_donor = "solventdonor " + str(solventdonor) if solventdonor else ""
    s_acceptor = "solventacceptor " + \
        str(solventacceptor) if solventacceptor else ""

    dslist = DataSetList()
    act = adict['hbond']

    command = " ".join(("series", mask, s_donor, s_acceptor))
    act(command, traj, dslist=dslist, *args, **kwd)
    act.print_output()

    if update_legend:
        _update_legend_hbond(dslist)
    if dtype == 'dataframe':
        # return DataFrame.T to have better visual effect
        return dslist.to_dataframe().T
    elif dtype == 'hbond_class':
        dslist_new = _get_data_from_dtype(dslist, dtype='dataset')
        return HbondAnalysisResult(dslist_new)
    else:
        return _get_data_from_dtype(dslist, dtype=dtype)


def search_nointramol_hbonds(traj,
                             mask="solventacceptor :WAT@O solventdonor :WAT",
                             dtype='dataset',
                             update_legend=True, *args, **kwd):
    """
    Search hbonds between solute and solvent, ignoring intra-hbond

    Parameters
    ----------
    traj : Trajectory-like or any iterable object that _frame_iter_mater return a Frame
    mask : str, default "solventacceptor :WAT@O solventdonor :WAT"
        cpptraj command
    dtype : str, default 'dataset'
    *args, **kwd: optional

    Examples
    --------
    >>> pyca.search_nointramol_hbonds(traj)
    >>> pyca.search_nointramol_hbonds([traj, traj2], top=traj.top)

    See Also
    --------
       search_hbonds
    """
    dslist = DataSetList()
    act = adict['hbond']
    command = "series nointramol " + mask
    act(command, traj, dslist=dslist, *args, **kwd)
    act.print_output()

    if update_legend:
        _update_legend_hbond(dslist)
    if dtype == 'hbond_class':
        return HbondAnalysisResult(dslist)
    else:
        return _get_data_from_dtype(dslist, dtype=dtype)
