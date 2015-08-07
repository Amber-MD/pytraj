import pytraj as pt

# use `iterload`
# data is not yet loaded to memory
traj = pt.iterload("../tests/data/tz2.ortho.nc",
                   "../tests/data/tz2.ortho.parm7")

# iterate trajectory by chunk with autoimage On
# chunksize = 3
for chunk in traj.iterchunk(chunksize=3, autoimage=True):
    print(chunk)
