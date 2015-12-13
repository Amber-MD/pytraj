from __future__ import absolute_import
from pytraj.get_common_objects import get_topology, get_data_from_dtype
from pytraj.decorators import register_pmap, register_openmp
from pytraj.get_common_objects import super_dispatch, get_fi_with_dslist
from pytraj.c_analysis import c_analysis
from pytraj.datasets.c_datasetlist import DatasetList as CpptrajDatasetList


@super_dispatch()
@register_openmp
def kmeans(traj=None,
           mask='*',
           n_clusters=10,
           random_point=True,
           kseed=1,
           maxit=100,
           metric='rms',
           top=None,
           frame_indices=None,
           options='',
           dtype='ndarray'):
    '''perform clustering and return cluster index for each frame

    Parameters
    ----------
    traj : Trajectory-like or iterable that produces Frame
    mask : str, default: * (all atoms)
    n_clusters: int, default: 10
    random_point : bool, default: True
    maxit : int, default: 100
        max iterations
    metric : str, {'rms', 'dme'}
        distance metric
    top : Topology, optional, default: None
        only need to provide this Topology if ``traj`` does not have one
    frame_indices : {None, 1D array-like}, optional
        if not None, only perform clustering for given indices. Notes that this is
        different from ``sieve`` keywords.
    options : str, optional
        extra cpptraj options controlling output, sieve, ...

    Sieve options::

        [sieve <#> [random [sieveseed <#>]]]

    Output options::

        [out <cnumvtime>] [gracecolor] [summary <summaryfile>] [info <infofile>]
        [summarysplit <splitfile>] [splitframe <comma-separated frame list>]
        [clustersvtime <filename> cvtwindow <window size>]
        [cpopvtime <file> [normpop | normframe]] [lifetime]
        [sil <silhouette file prefix>]

    Coordinate output options::

        [ clusterout <trajfileprefix> [clusterfmt <trajformat>] ]
        [ singlerepout <trajfilename> [singlerepfmt <trajformat>] ]
        [ repout <repprefix> [repfmt <repfmt>] [repframe] ]
        [ avgout <avgprefix> [avgfmt <avgfmt>] ]

    Notes
    -----
    - if the distance matrix is large (get memory Error), should add sieve number to
    ``options`` (check example)
    - install ``libcpptraj`` with ``-openmp`` flag to speed up this calculation.

    Returns
    -------
    1D numpy array of frame indices

    Examples
    --------
    >>> import pytraj as pt
    >>> from pytraj.cluster import kmeans
    >>> traj = pt.datafiles.load_tz2()
    >>> # use default options
    >>> kmeans(traj)
    array([8, 8, 6, ..., 0, 0, 0], dtype=int32)
    >>> # update n_clusters
    >>> data = kmeans(traj, n_clusters=5)
    >>> # update n_clusters with CA atoms
    >>> data = kmeans(traj, n_clusters=5, mask='@CA')
    >>> # specify distance metric
    >>> data = kmeans(traj, n_clusters=5, mask='@CA', kseed=100, metric='dme')
    >>> # add sieve number for less memory
    >>> data = kmeans(traj, n_clusters=5, mask='@CA', kseed=100, metric='rms', options='sieve 5')
    >>> # add sieve number for less memory, and specify random seed for sieve
    >>> data = kmeans(traj, n_clusters=5, mask='@CA', kseed=100, metric='rms', options='sieve 5 sieveseed 1')
    '''
    # don't need to get_topology
    _clusters = 'kmeans clusters ' + str(n_clusters)
    _random_point = 'randompoint' if random_point else ''
    _kseed = 'kseed ' + str(kseed)
    _maxit = str(maxit)
    _metric = metric
    # turn of cpptraj's cluster info
    _output = options
    command = ' '.join((_clusters, _random_point, _kseed, _maxit, _metric,
                        _output))
    return _cluster(traj, mask, frame_indices=frame_indices, top=top, dtype=dtype, options=command)


def _cluster(traj=None, mask="", frame_indices=None, dtype='dataset', top=None, options=''):
    """clustering

    Parameters
    ----------
    traj : Trajectory-like or any iterable that produces Frame
    mask : str
        atom mask
    frame_indices : {None, array-like}, optional
    dtype : str
        return data type
    top : Topology, optional
    options: str
        more cpptraj option

    Notes
    -----
    Supported algorithms: kmeans, hieragglo, and dbscan.
    """
    # Note: do not use super_dispatch here. We use get_fi_with_dslist

    ana = c_analysis.Analysis_Clustering()
    # need to creat `dslist` here so that every time `do_clustering` is called,
    # we will get a fresh one (or will get segfault)
    crdname = 'DEFAULT_NAME'
    dslist, _top, mask2 = get_fi_with_dslist(traj, mask, frame_indices, top, crdname=crdname)

    # do not output cluster info to STDOUT
    command = ' '.join((mask2, "crdset {0}".format(crdname), options, 'noinfo'))
    ana(command, dslist)

    # remove frames in dslist to save memory
    dslist.remove_set(dslist[crdname])
    return get_data_from_dtype(dslist, dtype=dtype)
