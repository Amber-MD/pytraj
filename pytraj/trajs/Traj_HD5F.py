# place holder to read/write hd5py file format
from __future__ import print_function
from pytraj.Frame import Frame
from pytraj.FrameArray import FrameArray
from pytraj.utils.check_and_assert import _import_h5py, _import_numpy

has_h5py, h5py = _import_h5py()
has_numpy, np = _import_numpy()

#assert has_h5py == True
#assert has_numpy == True

# TODO : inherit from FrameArray?
# FIXME: *** Error in `python': double free or corruption (out)
class HD5F():
    def __init__(self, filename=None, mode='r', flag='hd5f', *args, **kwd):
        print ("creating new HD5F instance")
        self.filename = filename
        self.mode = mode
        self.flag = flag # for what?

    def load_toframearray(self, filename, mode='r', top=None):
        farray = FrameArray()
        h5fh = h5py.File(filename, mode)
        farray.resize(h5fh['coordinates'].shape[0])

        for idx, arr in enumerate(h5fh['coordinates']):
            # allocate frame
            farray[idx] = Frame(arr.shape[0])
            # turn py_free_mem to False does not help avoiding memory error
            # why?
            #self[idx].py_free_mem = False
            # make sure to use np.float64 (double)
            farray[idx].set_from_crd(arr.flatten().astype(np.float64))
        return farray

    def write(self, *args, **kwd):
        raise NotImplementedError("not yet supported")
