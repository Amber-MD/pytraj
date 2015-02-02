# test if having numpy installed
HAS_NUMPY = False
ndarray = None
try:
    import numpy as np
    HAS_NUMPY = True
    ndarray = np.ndarray
except:
    pass

