import pytraj as pt
from pytraj.testing import aa_eq
from pytraj.utils.progress import ProgressBarTrajectory

def test():

    traj = pt.datafiles.load_tz2()

    t2 = ProgressBarTrajectory(traj)
    pt.to_pickle(t2, 'test.pk')

    t3 = pt.read_pickle('test.pk')

    data0 = pt.rmsd(traj)
    data1 = pt.rmsd(t2)
    data2 = pt.rmsd(t3)
    aa_eq(data0, data1)
    aa_eq(data0, data2)
