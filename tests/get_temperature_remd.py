import sys
from pytraj import io as mdio

# python ./get_temperatures_remd.py top_file md_x_file

traj = mdio.load(sys.argv[2], sys.argv[1])

print (traj.temperatures)
