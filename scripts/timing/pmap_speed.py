import pytraj as pt
from pytraj.testing import Timer
from time import time

traj = pt.iterload('md*.nc', 'prmtop', frame_slice=(0, 1000))
print(traj)

mask = '!:WAT,Na+'

t0 = time()
pt.rmsd(traj, mask=mask)
t_serial = time() - t0

for n_cores in [1, 2, 3, 4, 5, 6, 7, 8]:
    t0 = time()
    pt.pmap(pt.rmsd, traj, mask=mask, ref=traj[0], n_cores=n_cores)
    t_par = time() - t0
    efficiency = (t_serial / t_par) / n_cores
    print('n_cores = {}, time = {}, speed up = {}, efficiency = {}'.format(
        n_cores, t_par, t_serial / t_par, efficiency))
