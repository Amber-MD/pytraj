def test():
    import pytraj as pt
    from pytraj.utils.progress import ProgressBarTrajectory

    traj = pt.datafiles.load_tz2()

    t2 = ProgressBarTrajectory(traj)
    pt.to_pickle(t2, 'test.pk')

    t3 = pt.read_pickle('test.pk')

    for _ in t3: pass
