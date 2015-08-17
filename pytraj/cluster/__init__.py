from __future__ import absolute_import
from pytraj._get_common_objects import _get_top, _get_data_from_dtype


def kmeans(traj=None,
           n_clusters=10,
           random_point=True,
           kseed=1,
           maxit=100,
           distance_metric='rms',
           mask='*',
           top=None,
           output_op='',
           dtype='ndarray',
           *args, **kwd):
    '''perform clustering and return cluster index for each frame

    Returns
    -------
    1D numpy array

    Examples
    --------
    >>> from pytraj.cluster import kmeans
    >>> kmeans(traj)
    >>> kmeans(traj, n_clusters=5)
    >>> kmeans(traj, n_clusters=5, mask'='@CA')
    '''

    # don't need to _get_top
    _clusters = 'kmeans clusters ' + str(n_clusters)
    _random_point = 'randompoint' if random_point else ''
    _kseed = 'kseed ' + str(kseed)
    _maxit = str(maxit)
    _distance_metric = distance_metric
    _mask = mask
    _output = output_op
    command = ' '.join((_clusters, _random_point, _kseed, _maxit,
                        _distance_metric, _mask, _output))
    return _cluster(traj, command, top=top, dtype=dtype, *args, **kwd)


def _cluster(traj=None, command="", top=None, dtype='dataset', *args, **kwd):
    """
    Parameters
    ---------
    traj : Trajectory-like | list of Trajectory-like | frame or chunk iterator
    command : cpptraj command
    top : Topology, optional
    dslist : CpptrajDatasetList, optional
    dflist : DataFileList, optional

    Notes:
    Supported algorithms: kmeans, hieragglo, and dbscan.

    Examples
    --------
        do_clustering(traj, "kmeans clusters 50 @CA")

    Returns
    -------
    CpptrajDatasetList object

    """
    from pytraj.analyses.CpptrajAnalyses import Analysis_Clustering
    from pytraj.datasets.DataSetList import DataSetList as CpptrajDatasetList
    _top = _get_top(traj, top)
    ana = Analysis_Clustering()
    # need to creat `dslist` here so that every time `do_clustering` is called,
    # we will get a fresh one (or will get segfault)
    dslist = CpptrajDatasetList()

    if traj is not None:
        dslist.add_set("coords", "__pytraj_cluster")
        #dslist[-1].top = _top
        dslist[0].top = _top
        for frame in traj:
            # dslist[-1].add_frame(frame)
            dslist[0].add_frame(frame)
        command += " crdset __pytraj_cluster"
    else:
        pass
    ana(command, _top, dslist, *args, **kwd)
    # remove frames in dslist to save memory
    dslist.remove_set(dslist['__pytraj_cluster'])
    return _get_data_from_dtype(dslist, dtype=dtype)
