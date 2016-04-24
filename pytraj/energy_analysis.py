from __future__ import absolute_import

from .shared_methods import iterframe_master
from .get_common_objects import get_data_from_dtype, super_dispatch
from .externals.six.moves import range
from .externals.six import string_types
from .decorators import register_pmap

__all__ = ['esander']


def _default_func():
    from array import array
    return array('d', [])


@register_pmap
@super_dispatch()
def esander(traj=None,
            prmtop=None,
            igb=8,
            mm_options=None,
            qm_options=None,
            mode=None,
            dtype='dict',
            frame_indices=None,
            top=None):
    """energy decomposition by calling `libsander`

    Parameters
    ----------
    traj : Trajectory-like or iterables that produce Frame
        if `traj` does not hold Topology information, `top` must be provided
    prmtop : str or Structure from ParmEd, default=None, optional
        To avoid any unexpected error, you should always provide original topology
        filename. If prmtop is None, pytraj will load Topology from traj.top.filename.

        - why do you need to load additional topology filename? Because cpptraj and sander
          use different Topology object, can not convert from one to another.
    igb : GB model, default=8 (GB-Neck2)
        If specify `mm_options`, this `igb` input will be ignored
    mm_options : InputOptions from `sander`, default=None, optional
        if `mm_options` is None, use `gas_input` with given igb.
        If `mm_options` is not None, use this
    qm_options : InputOptions from `sander` for QMMM, optional
    mode : str, default=None, optional
        if mode='minimal', get only 'bond', 'angle', 'dihedral' and 'total' energies
    top : pytraj.Topology or str, default=None, optional
        only need to specify this ``top`` if ``traj`` does not hold Topology
    dtype : str, {'dict', 'dataset', 'ndarray', 'dataframe'}, default='dict'
        return data type
    frame_indices : None or 1D array-like, default None
        if not None, only perform calculation for given frames

    Returns
    -------
    Dict of energies (to be used with DataFrame) or DatasetList

    Examples
    --------
    Examples are adapted from $AMBERHOME/test/sanderapi

    >>> import pytraj as pt
    >>> # GB energy
    >>> traj = pt.datafiles.load_ala3()
    >>> traj.n_frames
    1
    >>> data = pt.esander(traj, igb=8)
    >>> data['gb']
    array([-92.88577683])
    >>> data['bond']
    array([ 5.59350521])

    >>> # PME
    >>> import os
    >>> from pytraj.testing import amberhome
    >>> import sander
    >>> topfile = os.path.join(amberhome, "test/4096wat/prmtop")
    >>> rstfile = os.path.join(amberhome, "test/4096wat/eq1.x")
    >>> traj = pt.iterload(rstfile, topfile)
    >>> options = sander.pme_input()
    >>> options.cut = 8.0
    >>> edict = pt.esander(traj=traj, mm_options=options)
    >>> edict['vdw']
    array([ 6028.95167558])

    >>> # GB + QMMM
    >>> topfile = os.path.join(amberhome, "test/qmmm2/lysine_PM3_qmgb2/prmtop")
    >>> rstfile = os.path.join(amberhome, "test/qmmm2/lysine_PM3_qmgb2/lysine.crd")
    >>> traj = pt.iterload(rstfile, topfile)

    >>> options = sander.gas_input(8)
    >>> options.cut = 99.0
    >>> options.ifqnt = 1
    >>> qm_options = sander.qm_input()
    >>> qm_options.iqmatoms[:3] = [8, 9, 10]
    >>> qm_options.qm_theory = "PM3"
    >>> qm_options.qmcharge = 0
    >>> qm_options.qmgb = 2
    >>> qm_options.adjust_q = 0

    >>> edict = pt.esander(traj=traj, mm_options=options, qm_options=qm_options)
    >>> edict['bond']
    array([ 0.00160733])
    >>> edict['scf']
    array([-11.92177575])

    >>> # passing options to `pytraj.pmap`: need to pass string
    >>> from pytraj.testing import get_fn
    >>> fn, tn = get_fn('tz2')
    >>> traj = pt.iterload(fn, tn)
    >>> inp_str = 'mm_options = sander.pme_input()'
    >>> edict = pt.pmap(pt.esander, traj, mm_options=inp_str, n_cores=2)
    >>> edict['dihedral']
    array([ 126.39307126,  127.0460586 ,  137.26793522,  125.30521069,
            125.25110884,  137.69287326,  125.78280543,  125.14530517,
            118.41540102,  128.73535036])

    Notes
    -----
    This method does not work with `pytraj.pmap` when you specify mm_options and
    qm_options. Use `pytraj.pmap_mpi` with MPI instead.

    Work with ``pytraj.pmap``::

        pt.pmap(pt.esander, traj, igb=8, dtype='dict')

    Will NOT work with ``pytraj.pmap``::

        import sander
        inp = sander.gas_input(8)
        pt.pmap(pt.esander, traj, mm_options=inp, dtype='dict')

    Why? Because Python need to pickle each object to send to different cores and Python
    does not know how to pickle mm_options from sander.gas_input(8).

    This works with ``pytraj.pmap_mpi`` because pytraj explicitly create ``mm_options``
    in each core without pickling.
    """
    from collections import defaultdict, OrderedDict
    from pytraj.misc import get_atts
    import numpy as np

    try:
        import sander
    except ImportError:
        raise ImportError("need both `pysander` installed. Check Ambertools15")

    ddict = defaultdict(_default_func)

    if mm_options is None:
        inp = sander.gas_input(igb)
    elif igb is not None:
        inp = mm_options

    if isinstance(inp, string_types):
        # dangerous
        local_dict = {'sander': sander}
        exec(inp.lstrip(), local_dict)
        inp = local_dict['mm_options']

    if isinstance(qm_options, string_types):
        # dangerous
        local_dict = {'sander': sander}
        exec(qm_options.lstrip(), local_dict)
        qm_options = local_dict['qm_options']

    if prmtop is None:
        try:
            # try to load from file by taking top.filename
            prmtop_ = top.filename
        except AttributeError:
            raise ValueError("prmtop must be AmberParm object in ParmEd")
    else:
        # Structure, string
        prmtop_ = prmtop

    if not hasattr(prmtop_, 'coordinates') or prmtop_.coordinates is None:
        try:
            # if `traj` is Trajectory-like (not frame_iter), try to take 1st
            # coords
            coords = traj[0].xyz
        except (TypeError, AttributeError):
            # create fake list
            coords = [0. for _ in range(top.n_atoms * 3)]
    else:
        # use default coords in `AmberParm`
        coords = prmtop_.coordinates

    if top.has_box():
        box = top.box.tolist()
        has_box = True
    else:
        box = None
        has_box = False

    with sander.setup(prmtop_, coords, box, inp, qm_options):
        for frame in iterframe_master(traj):
            if has_box:
                sander.set_box(*frame.box.tolist())
            sander.set_positions(frame.xyz)
            ene, frc = sander.energy_forces()

            # potentially slow
            ene_atts = get_atts(ene)
            for att in ene_atts:
                ddict[att].append(getattr(ene, att))

    new_dict = None
    if mode == 'minimal':
        new_dict = {}
        for key in ['bond', 'angle', 'dihedral', 'tot']:
            new_dict[key] = ddict[key]
    else:
        new_dict = ddict

    for key in new_dict.keys():
        new_dict[key] = np.asarray(new_dict[key])

    if dtype == 'dict':
        return OrderedDict(new_dict)
    else:
        from pytraj.datasets.c_datasetlist import DatasetList

        dslist = DatasetList()
        size = new_dict['tot'].__len__()
        for key in new_dict.keys():
            dslist.add('double')
            dslist[-1].key = key
            dslist[-1].resize(size)
            dslist[-1].data[:] = new_dict[key]
        return get_data_from_dtype(dslist, dtype)
