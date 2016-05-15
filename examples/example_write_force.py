import numpy as np
import pytraj as pt
import sander

traj = pt.datafiles.load_tz2()

inp = sander.gas_input(8)

frcs = []

with sander.setup(traj.top.filename, traj[0].xyz, traj.top.box, inp):
    for frame in traj:
        sander.set_box(*frame.box.tolist())
        sander.set_positions(frame.xyz)
        ene, frc = sander.energy_forces()
        frcs.append(np.array(frc).reshape(traj.n_atoms, 3))

def get_frame_with_force(traj, forces=frcs):
    frame0 = pt.Frame()
    crdinfo = dict(has_force=True)

    frame0._allocate_force_and_velocity(traj.top, crdinfo)

    for frame, frc in zip(traj, frcs):
        frame0.xyz[:] = frame.xyz
        frame0.force[:] = frc 
        yield frame0

pt.write_traj('traj.nc',
              traj=get_frame_with_force(traj),
              top=traj.top,
              overwrite=True,
              options='force')
