from pytraj._shared_methods import iterframe_master
from pytraj._get_common_objects import _get_topology, _get_data_from_dtype
from pytraj.compat import range
from pytraj.decorators import _register_pmap

__all__ = ['get_pysander_energies']


def _default_func():
    from array import array
    return array('d', [])


@_register_pmap
def get_pysander_energies(traj=None,
                          parm=None,
                          igb=8,
                          input_options=None,
                          qmmm_options=None,
                          mode=None,
                          top=None,
                          dtype='dict',
                          verbose=False):
    # TODO: change method's name?
    """"
    Parameters
    ---------
    traj : {Traj-like object, frame, list of trajs, list of frames} from pytraj
        if `traj` does not hold Topology information, `top` must be provided
    parm : {str, Topology object from ParmEd}, default=None, optional
    igb : GB model, default=8 (GB-Neck2)
        Note: this `igb` input will be ignored if `input_options` is not None
    input_options : InputOptions object from `sander`, default=None, optional
        if `input_options` is None, use `gas_input` with `igb = 8`
        If `input_options` is not None, use this
    qmmm_options : InputOptions object from `sander` for QMMM, optional
    mode : str, default=None, optional
        if mode='minimal', get only 'bond', 'angle', 'dihedral' and 'total' energies
    top : {Topology, str}, default=None, optional
    dtype : str, {'dict', 'dataset', 'ndarray', 'dataframe'}, default='dict'
    verbose : bool, default True
        print warning message if True

    Returns
    -------
    Dict of energies (to be used with DataFrame) or DatasetList

    Examples
    --------
        # minimal input
        energy_decomposition = get_pysander_energies
        energy_decomposition(traj)

        # with option
        import sander
        inp = sander.gas_input(igb=6)
        energy_decomposition(traj, input_options=inp)

        # with list of frames, must provide Topology object
        energy_decomposition([frame0, frame1], top=my_topology_object)

        # with provided ParmEd object
        import parmed as pmd
        parm = pmd.load_file("myfile.prmtop")
        energy_decomposition(traj, parm=parm, igb=5)
    """
    from collections import defaultdict
    from pytraj.misc import get_atts
    import numpy as np

    try:
        import sander
    except ImportError:
        raise ImportError("need both `pysander` installed. Check Ambertools15")

    ddict = defaultdict(_default_func)

    _top = _get_topology(traj, top)

    if input_options is None:
        inp = sander.gas_input(igb)
    elif igb is not None:
        if verbose:
            print("inp is not None, ignore provided `igb` and use `inp`")
        inp = input_options

    if parm is None:
        try:
            # try to load from file by taking _top.filename
            if verbose:
                print("can not find `Structure` from parmed, loading %s")
            _parm = _top.filename
        except AttributeError:
            raise ValueError("parm must be AmberParm object in ParmEd")
    else:
        # Structure, string
        _parm = parm

    if not hasattr(_parm, 'coordinates') or _parm.coordinates is None:
        try:
            # if `traj` is Trajectory-like (not frame_iter), try to take 1st
            # coords
            coords = traj[0].xyz
        except (TypeError, AttributeError):
            # create fake list
            coords = [0. for _ in range(_top.n_atoms * 3)]
    else:
        # use default coords in `AmberParm`
        coords = _parm.coordinates

    if _top.has_box():
        box = _top.box.tolist()
        has_box = True
    else:
        box = None
        has_box = False

    with sander.setup(_parm, coords, box, inp, qmmm_options):
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
        return new_dict
    else:
        from pytraj.datasets.DatasetList import DatasetList

        dslist = DatasetList()
        size = new_dict['tot'].__len__()
        for key in new_dict.keys():
            dslist.add_set('double')
            dslist[-1].key = key
            dslist[-1].resize(size)
            dslist[-1].data[:] = new_dict[key]
        return _get_data_from_dtype(dslist, dtype)
