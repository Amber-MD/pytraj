import pytraj as pt

try:
    import mdtraj as md
    traj = pt.load("./nvt.xtc", "nvt.gro", engine='mdtraj')
except (ImportError, RuntimeError):
    pass
