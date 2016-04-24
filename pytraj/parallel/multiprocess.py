# do not use relative import here. Treat this module as a seperated package.
from __future__ import absolute_import
import numpy as np
from functools import partial
from collections import OrderedDict
from pytraj.c_options import info as compiled_info
from pytraj import matrix
from pytraj import mean_structure, volmap
from pytraj import Frame
from pytraj import ired_vector_and_matrix, rotation_matrix
from pytraj import NH_order_parameters
from multiprocessing import cpu_count
from pytraj.tools import concat_dict
from pytraj.externals.six import string_types
from pytraj.get_common_objects import get_reference

from .base import worker_byfunc


def _pmap(func, traj, *args, **kwd):
    '''use python's multiprocessing to accelerate calculation. Limited calculations.

    Parameters
    ----------
    func : a pytraj's methods or a list of string or simply as a cpptraj' text
    traj : pytraj.TrajectoryIterator
    n_cores : int, number of cores to be used, default 2. Specify n_cores=-1 to use all available cores
    iter_options : dict, default {}
        Specify trajectory iterating option. This will be done before calling ``func``.
    frame_indices : {None, array-like}, default None, optional
        if provided, pytraj will split this frame_indices into different chunks and let
        cpptraj perform calculation for specific indices.
        frame_indices must be pickable so is can be sent to different cores.

    *args, **kwd: additional keywords

    Returns
    -------
    out : OrderedDict

    Notes
    -----
    - If you not sure about parallel's results, you should compare the output to serial run.

    - This is absolutely experimental. The syntax might be changed in future.

    Rule of thumbs: start with small number of frames (saying 10 frames), varying
    n_cores=1, 2, 3, 4 to see if the data makes sense or not.

    There are two modes in this method, use pytraj's methods (pytraj.rmsd, pytraj.radgyr,
    ...) or use cpptraj's command text syntax ('autoimage', 'rms', ...)

    If using cpptraj syntax::

        pytraj only supports limited cpptraj's Actions (not Analysis, checm Amber15 manual
        about Action and Analysis), say no  to 'matrix', 'atomicfluct', ... or any action
        that results output depending on the number of frames.


    This method only benifits you if your calculation is quite long (saying few minutes to
    few hours). For calculation that takes less than 1 minutes, you won't see the
    significant speed up (or even slower) since pytraj need to warm up and need to gather
    data when the calculation done.

    The parallel cacluation is very simple, trajectory will be split (almost equal) to
    different chunk (n_chunks = n_cores), pytraj/cpptraj perform calculation for each
    chunk in each core and then send data back to master. Note that we are using Python's
    built-in multiprocessing module, so you can use this method interactively in Ipython
    and ipython/jupyter notebook. This behavior is different from using MPI, in which you
    need to write a script, escaping ipython ession and type something like::

        mpirun -n 4 python my_script.py

    vs::

        In [1]: pt.pmap(pt.radgyr, traj, n_cores=4)
        Out[1]:
        OrderedDict([('RoG_00000',
                      array([ 18.91114428,  18.93654996,  18.84969884,  18.90449256,
                              18.8568644 ,  18.88917208,  18.9430491 ,  18.88878079,
                              18.91669565,  18.87069722]))])

    This is experimental method, you should expect its syntax, default output will be changed.

    When sending Topology to different cores, pytraj will reload Topology from
    traj.top.filename, so if you need to update Topology (in the fly), save it to disk and
    reload before using ``pytraj.pmap``

    Examples
    --------
    >>> import numpy as np
    >>> import pytraj as pt
    >>> traj = pt.load_sample_data('tz2')

    >>> # use iter_options
    >>> iter_options = {'autoimage': True, 'rmsfit': (0, '@CA')}
    >>> data = pt.pmap(pt.mean_structure, traj, iter_options=iter_options)

    >>> # cpptraj command style
    >>> data = pt.pmap(['distance :3 :7', 'vector mask :3 :12'], traj, n_cores=4)

    >>> # use reference. Need to explicitly use 'refindex', which is index of reflist
    >>> data = pt.pmap(['rms @CA refindex 0'], traj, ref=[traj[3],], n_cores=3)
    >>> data
    OrderedDict([('RMSD_00001', array([  2.68820312e-01,   3.11804885e-01,   2.58835452e-01,
             9.10475988e-08,   2.93310737e-01,   4.10197322e-01,
             3.96226694e-01,   3.66059215e-01,   3.90890362e-01,
             4.89180497e-01]))])

    >>> # use reference: if not want to use 'refindex', can use 'reference'
    >>> # the advantage is you can not specify a list of reference
    >>> data = pt.pmap(['rms @CA reference'], traj, ref=[traj[3],], n_cores=3)
    >>> data
    OrderedDict([('RMSD_00001', array([  2.68820312e-01,   3.11804885e-01,   2.58835452e-01,
             9.10475988e-08,   2.93310737e-01,   4.10197322e-01,
             3.96226694e-01,   3.66059215e-01,   3.90890362e-01,
             4.89180497e-01]))])

    >>> # use different references. Need to explicitly use 'refindex', which is index of reflist
    >>> # create a list of references
    >>> reflist = traj[3], traj[4]
    >>> # make sure to specify `refindex`
    >>> # `refindex 0` is equal to `reflist[0]`
    >>> # `refindex 1` is equal to `reflist[1]`
    >>> data = pt.pmap(['rms @CA refindex 0', 'rms !@H= refindex 1'], traj, ref=reflist, n_cores=2)
    >>> # convert to ndarray
    >>> data_arr = pt.tools.dict_to_ndarray(data)

    >>> # perform parallel calculation with given frame_indices
    >>> traj = pt.datafiles.load_tz2()
    >>> data = pt.pmap(pt.radgyr, traj, '@CA', frame_indices=range(10, 50), n_cores=4)
    >>> # serial version
    >>> data = pt.radgyr(traj, '@CA', frame_indices=range(10, 50))


    See also
    --------
    pytraj.pmap_mpi
    '''
    from multiprocessing import Pool
    from pytraj import TrajectoryIterator

    n_cores = kwd.pop('n_cores') if 'n_cores' in kwd else 2
    iter_options = kwd.pop('iter_options') if 'iter_options' in kwd else {}
    apply = kwd.pop('apply') if 'apply' in kwd else None

    if n_cores <= 0:
        # use all available cores
        n_cores = cpu_count()

    # update reference
    if 'ref' in kwd:
        kwd['ref'] = get_reference(traj, kwd['ref'])

    if isinstance(func, (list, tuple, string_types)):
        # assume using _load_batch_pmap
        from pytraj.parallel.base import _load_batch_pmap
        #check_valid_command(func)
        data = _load_batch_pmap(n_cores=n_cores,
                                traj=traj,
                                lines=func,
                                dtype='dict',
                                root=0,
                                mode='multiprocessing',
                                **kwd)
        data = concat_dict((x[1] for x in data))
        return data
    else:
        if not callable(func):
            raise ValueError('must callable argument')
        # pytraj's method
        if not hasattr(func, '_is_parallelizable'):
            raise ValueError("this method does not support parallel")
        elif not func._is_parallelizable:
            raise ValueError("this method does not support parallel")
        else:
            if hasattr(
                    func,
                    '_openmp_capability') and func._openmp_capability and 'OPENMP' in compiled_info(
            ):
                raise RuntimeError(
                    "this method supports both openmp and pmap, but your cpptraj "
                    "version was installed with openmp. Should not use both openmp and pmap at the "
                    "same time. In this case, do not use pmap since openmp is more efficient")

        if not isinstance(traj, TrajectoryIterator):
            raise ValueError('only support TrajectoryIterator')

        if 'dtype' not in kwd and func not in [
                mean_structure, matrix.dist, matrix.idea,
                ired_vector_and_matrix, rotation_matrix,
                volmap,
        ]:
            kwd['dtype'] = 'dict'

        # keyword
        if func is volmap:
            assert kwd.get('size') is not None, 'must provide "size" value'

        p = Pool(n_cores)

        pfuncs = partial(worker_byfunc,
                         n_cores=n_cores,
                         func=func,
                         traj=traj,
                         args=args,
                         kwd=kwd,
                         iter_options=iter_options,
                         apply=apply)

        data = p.map(pfuncs, [rank for rank in range(n_cores)])
        p.close()

        if func in [matrix.dist, matrix.idea, volmap]:
            mat = np.sum((val[1] * val[2] for val in data)) / traj.n_frames
            return mat
        elif func in [ired_vector_and_matrix, ]:
            # data is a list of (rank, (vectors, matrix), n_frames)
            mat = np.sum((val[1][1] * val[2] for val in data)) / traj.n_frames
            vecs = np.column_stack(val[1][0] for val in data)
            return (vecs, mat)
        elif func in [rotation_matrix, ]:
            if 'with_rmsd' in kwd.keys() and kwd['with_rmsd']:
                # data is a list of (rank, (mat, rmsd), n_frames)
                mat = np.row_stack(val[1][0] for val in data)
                rmsd_ = np.hstack(val[1][1] for val in data)
                return OrderedDict(out=(mat, rmsd_))
            else:
                mat = np.row_stack(val[1] for val in data)
                return OrderedDict(mat=mat)
        elif func == mean_structure:
            xyz = np.sum((x[2] * x[1].xyz for x in data)) / traj.n_frames
            frame = Frame(xyz.shape[0])
            frame.xyz[:] = xyz
            return frame
        else:
            return concat_dict((x[1] for x in data))


def pmap(func=None, traj=None, *args, **kwd):
    if func != NH_order_parameters:
        return _pmap(func, traj, *args, **kwd)
    else:
        if 'n_cores' in kwd.keys():
            if kwd['n_cores'] == 1:
                # use default n_cores=2 instead of 1
                kwd['n_cores'] = 2
            if kwd['n_cores'] <= 0:
                kwd['n_cores'] = cpu_count()
        else:
            # use n_cores=2 for default value
            kwd['n_cores'] = 2
        return NH_order_parameters(traj, *args, **kwd)


pmap.__doc__ = _pmap.__doc__
