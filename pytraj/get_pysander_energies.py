print ("method's name migh be changed")
from pytraj.utils import has_
from pytraj.externals.six import string_types

def get_pysander_energies(parm, traj, igb=8):
    """"
    Parameters
    ---------
    parm : Topology object from ParmEd
    traj : Traj-like object from pytraj
    igb : GB model, default=8 (GB-Neck2)

    Returns:
    Dict of energies (to be used with DataFrame)
    """
    import sander
    from pytraj.misc import get_atts
    from collections import defaultdict
    from chemistry.amber.readparm import AmberParm
    ddict = defaultdict(list, [])

    inp = sander.gas_input(igb)
    if isinstance(parm, string_types):
        parm = AmberParm(parm)

    if parm.coords is None:
        parm.load_coordinates(traj[0].coords)

    with sander.setup(parm, parm.coords, None, inp):
        for frame in traj:
            sander.set_positions(frame.coords)
            ene, frc = sander.energy_forces()
            ene_atts = get_atts(ene)
            for att in ene_atts:
                ddict[att].append(getattr(ene, att))
    return ddict
