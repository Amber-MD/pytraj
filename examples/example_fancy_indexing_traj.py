from pytraj import io as mdio

traj = mdio.load("../tests/data/Tc5b.x", "../tests/data/Tc5b.top")[:]
# get x coord of 0-th atom of 0-th frame
print(traj[0, 0, 0])

# get new Trajectory with only CA atoms
print(traj['@CA'])

# get new Trajectory with only CA with several frames
print(traj[2:7, '@CA'])
print(traj[[1, 5, 7, 9], '@CA'])
