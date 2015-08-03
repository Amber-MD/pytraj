import pytraj as pt

# use `iterload` to save memory
traj = pt.load("../tests/data/tz2.ortho.nc",
               "../tests/data/tz2.ortho.parm7")
print(traj)

# get new trajectory with only waters
t0 = traj[':WAT@O']

# make new Topology for Oxygen of WAT
# you can choose any water number, I randomly chose `1000@O`
wat_top = traj.top._get_new_from_mask(':1000@O')

# for each frame, make a new Trajectory
# where each water is a Frame
# get water topology first
def make_water_traj(frame, top=wat_top):
    xyz = frame.xyz.reshape(frame.shape[0], 1 , 3)
    return pt.Trajectory(xyz=xyz, top=top)

print(t0)
for f in t0:
    new_traj = make_water_traj(f)
    print(new_traj)
    #print(pt.clustering(new_traj))
