from pytraj.datasets import DataSet_integer
import numpy as np

ds = DataSet_integer() # wrapper of cpptraj' class
ds.resize(100) # allocate 100 items

# initialize some data
# must use np.int32 since we used `int` for DataSet_integer
arr0 = np.random.randint(0, 1000, 100).astype(np.int32)

# copy data to `ds`
ds.data[:] = arr0

# try to use `ds` as `numpy` array without casting
assert np.sum(ds) == np.sum(arr0)

ds2 = DataSet_integer()
ds2.resize(ds.size)
print (ds2.size, ds2.tolist())
ds2.data[:] = ds.data[:]
print (ds2.size, ds2.tolist())

ds3 = DataSet_integer()
ds3.resize(ds.size)
ds3.data[:] = arr0
print (ds3.size, ds3.tolist())

ds4 = DataSet_integer()
ds4.resize(ds.size)
ds4.data[:] = arr0[:]
print (ds4.size, ds4.tolist())
