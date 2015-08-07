try:
    import pytraj as pt
    import mdtraj as md
    import MDAnalysis as MDA

    from MDAnalysisTests.datafiles import PSF, DCD

    traj0 = pt.iterload(DCD, PSF, engine='pytraj')
    traj1 = pt.iterload(DCD, PSF, engine='mdanalysis')
    traj2 = pt.load(DCD, PSF, engine='mdtraj')
    print(traj0, traj1, traj2)

except ImportError:
    print("ImportError. skip")
