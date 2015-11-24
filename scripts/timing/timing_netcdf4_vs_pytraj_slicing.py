
import pytraj as pt
from netCDF4 import Dataset
from pytraj.testing import Timer

filename = 'GAAC3.5000frames.nc'
top_name = 'GAAC.topo'
top = pt.load_topology(top_name)

@Timer()
def load_netcdf4(filename=filename):
    with Dataset(filename, 'r') as fh:
        c= fh['coordinates']
        return c[:N]

@Timer()
def load_mutable_traj(filename=filename):
    with Dataset(filename, 'r') as fh:
        c= fh['coordinates']
        traj = pt.Trajectory(xyz=c[:N], top=top.copy())
        traj.xyz[:] = traj.xyz.astype('f8')

@Timer()
def load_immutable_traj(filename=filename, top_name=top_name):
    traj = pt.iterload(filename, top_name)
    return traj[:N]

if __name__ == '__main__':
    import sys
    from collections import OrderedDict
    from matplotlib import pyplot as plt
    import pandas as pd
    import numpy as np

    data = OrderedDict()
    for N in [10, 50, 100, 200, 500, 1000, 1200]:
    #for N in [10, 50]:
        t_netcdf4 = load_netcdf4()
        timelist = np.array([load_netcdf4(), load_mutable_traj(), load_immutable_traj()])
        timelist = timelist / timelist[0]
        data['n_frames={}'.format(N)] = timelist
    df = pd.DataFrame(data)
    # transpose
    df = df.T
    df.columns = ('netcdf4', 'netcdf4-traj', 'TrajectoryIterator')
    print(df)
    #df.plot.barh()
    #plt.show()
    #plt.savefig('speed_pytraj_vs_netcdf4_slicing.png', dpi=300)

    # output: relative time (smaller is better)
    #                    netcdf4  netcdf4-traj  TrajectoryIterator
    # n_frames=10          1      2.479123           18.831928
    # n_frames=50          1      1.626648            4.657385
    # n_frames=100         1      1.584890            2.872311
    # n_frames=200         1      1.496677            1.924031
    # n_frames=500         1      1.463720            1.354975
    # n_frames=1000        1      1.446821            1.170017
    # n_frames=1200        1      1.455076            1.141454
