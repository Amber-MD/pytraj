import netCDF4 as nc
from pytraj import io
from timeit import timeit
from pytraj.testing import aa_eq

file = nc.Dataset("./md.trj")
traj = io.iterload("./md.trj", 'tc5bwat.top')
print(traj)
c = file['coordinates']

s = slice(10, 18, 2)


def test_nc():
    c[s]


def test_traj():
    traj[s]


def test_iter_nc():
    for arr in c:
        pass


def test_iter_traj():
    for frame in traj:
        pass

print("slicing")
print("netCDF4")
print(timeit(test_nc, number=5))
print("traj")
print(timeit(test_traj, number=5))  # >= 3 time faster

print("itering whole data/traj")
print("netCDF4")
print(timeit(test_iter_nc, number=5))
print("traj")
print(timeit(test_iter_traj,  number=5))  # >= 2 times faster

ncxyz = c[s]
pxyz = traj[s].xyz
aa_eq(ncxyz, pxyz)
