import pytraj as pt

traj = pt.iterload('./data/tz2.ortho.nc', './data/tz2.ortho.parm7')

def worker(traj=traj):
    print(pt.rmsd(traj))

import threading

threads = []
for i in range(100):
   t = threading.Thread(target=worker)
   threads.append(t)
   t.start()
