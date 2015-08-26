import numpy as np
from memory_profiler import profile

@profile
def test():
    a = np.arange(1000000).reshape(1000, 1000)
    print(a.shape)
    a = np.concatenate(a, a[0].reshape(1, 1000))
    print(a.shape)

test()
