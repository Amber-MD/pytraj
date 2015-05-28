from pytraj.utils import has_
from pytraj.externals.six import string_types

def get_pysander_energies(traj=None, parm=None, igb=8, input_options=None, qmmm_options=None, mode=None):
    """"
    Parameters
    ---------
    traj : Traj-like object from pytraj
    parm : {str, Topology object from ParmEd}
    igb : GB model, default=8 (GB-Neck2)
    input_options : InputOptions object from `sander`, default=None, optional
        if `input_options` is None, use `gas_input` with `igb = 8`
        If `input_options` is not None, use this
    qmmm_options : InputOptions object from `sander` for qmmm
    mode : str, default=None
        if mode='minimal', get only 'bond', 'angle', 'dihedral' and 'total' energies

    Returns:
    Dict of energies (to be used with DataFrame)

    Examples
    --------
        energy_decomposition = get_pysander_energies
        import sander
        inp = sander.gas_input(igb=6)
        energy_decomposition(traj, input_options=inp)
    """

    from pytraj.misc import get_atts
    from collections import defaultdict
    try:
        import sander
        from chemistry.amber.readparm import AmberParm
    except ImportError:
        raise ImportError("need both `pysander` and `chemistry` installed. Check Ambertools15")

    ddict = defaultdict(list, [])

    if input_options is None:
        inp = sander.gas_input(igb)
    elif igb is not None:
        print ("inp is not None, ignore provided `igb` and use `inp`")
        inp = input_options

    if isinstance(parm, string_types):
        parm = AmberParm(parm)

    if not isinstance(parm, AmberParm):
        try:
            parm = AmberParm(traj.top.filename)
        except:
            raise ValueError("parm must be AmberParm object in ParmEd")

    if parm.coords is None:
        parm.load_coordinates(traj[0].coords)

    if traj.top.has_box():
        box = traj.top.box.tolist()
        has_box = True
    else:
        box = None
        has_box = False

    with sander.setup(parm, parm.coords, box, inp, qmmm_options):
        for frame in traj:
            if has_box:
                sander.set_box(*frame.box.tolist())
            sander.set_positions(frame.coords)
            ene, frc = sander.energy_forces()

            # potentially slow
            ene_atts = get_atts(ene)
            for att in ene_atts:
                ddict[att].append(getattr(ene, att))

    if mode == 'minimal':
        new_dict = {}
        for key in ['bond', 'angle', 'dihedral', 'tot']:
            new_dict[key] = ddict[key]
        return new_dict
    else:
        return ddict
