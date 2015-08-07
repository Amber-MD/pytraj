import pytraj as pt

# use `iterload` to save memory
traj = pt.iterload("../tests/data/tz2.ortho.nc",
                   "../tests/data/tz2.ortho.parm7")

try:
    import pandas as pd
    # search all possible dihedral types supported by cpptraj
    dset = pt.multidihedral(traj, dtype='dataframe')
    print(dset)
except ImportError:
    pass
