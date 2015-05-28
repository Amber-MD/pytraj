from pytraj.utils import has_
from pytraj.externals.six import string_types
from pytraj._shared_methods import _frame_iter_master
from pytraj._get_common_objects import _get_top, _get_data_from_dtype
from pytraj.compat import range

__all__ = ['get_pysander_energies']

# a copy from `pytraj.testing`
def make_random_frame(n_atoms=10000):
    import numpy as np
    from pytraj import Frame

    frame = Frame(n_atoms)
    frame.xyz[:] = np.random.randn(n_atoms, 3)
    return frame

def get_pysander_energies(traj=None, parm=None, igb=8, input_options=None, qmmm_options=None, 
                          mode=None, top=None, dtype='dict'):
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

    Returns:
    Dict of energies (to be used with DataFrame)

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
        import chemistry as chem
        parm = chem.load_file("myfile.prmtop")
        energy_decomposition(traj, parm=parm, igb=5)
    """
    from array import array as pyarray
    from collections import defaultdict
    from pytraj.misc import get_atts
    try:
        import sander
        from chemistry.amber.readparm import AmberParm
    except ImportError:
        raise ImportError("need both `pysander` and `chemistry` installed. Check Ambertools15")

    ddict = defaultdict(lambda : pyarray('d', []))
    _top = _get_top(traj, top)

    if input_options is None:
        inp = sander.gas_input(igb)
    elif igb is not None:
        print ("inp is not None, ignore provided `igb` and use `inp`")
        inp = input_options

    if isinstance(parm, string_types):
        parm = AmberParm(parm)

    if not isinstance(parm, AmberParm):
        try:
            parm = AmberParm(_top.filename)
        except:
            raise ValueError("parm must be AmberParm object in ParmEd")

    if parm.coords is None:
        try:
            parm.load_coordinates(traj[0].coords)
        except:
            # make fake coords to trick parm
            parm.load_coordinates([0. for _ in range(_top.n_atoms * 3)])

    if _top.has_box():
        box = _top.box.tolist()
        has_box = True
    else:
        box = None
        has_box = False

    with sander.setup(parm, parm.coords, box, inp, qmmm_options):
        for frame in _frame_iter_master(traj):
            if has_box:
                sander.set_box(*frame.box.tolist())
            sander.set_positions(frame.coords)
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

    if dtype == 'dict':
        return new_dict
    else:
        from pytraj import DataSetList

        dslist = DataSetList()
        size = new_dict['tot'].__len__()
        for key in new_dict.keys():
            dslist.add_set('double')
            dslist[-1].legend = key
            dslist[-1].resize(size)
            dslist[-1].data[:] = new_dict[key]
        return _get_data_from_dtype(dslist, dtype)
