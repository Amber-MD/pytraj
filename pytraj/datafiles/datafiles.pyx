# distutils: language = c++

from libcpp.string cimport string
from ..core.cpptraj_core cimport CpptrajState, _CpptrajState

cdef extern from "Command.h": 
    ctypedef enum RetType "Command::RetType":
        pass
    cdef cppclass _Command "Command":
        @staticmethod
        RetType ProcessInput(_CpptrajState&, const string&)
        @staticmethod
        RetType Dispatch(_CpptrajState&, const string&)


cdef class DataFile:
    def __cinit__(self, py_free_mem=True):
        self.thisptr = new _DataFile()
        self.py_free_mem = py_free_mem

    def __dealloc__(self):
        if self.py_free_mem:
            del self.thisptr

    def read_data(self, filenameIn, arglist, DatasetList datasetlist):
        return self.thisptr.ReadDataIn(filenameIn.encode(),
               ArgList(arglist).thisptr[0], datasetlist.thisptr[0])


cdef class DataFileList:
    def __cinit__(self, py_free_mem=True):
        self.thisptr = new _DataFileList()
        self.py_free_mem = py_free_mem

    def __dealloc__(self):
        if self.py_free_mem:
            del self.thisptr

    def write_all_datafiles(self):
        self.thisptr.WriteAllDF()

def load_batch(traj, txt):
    '''
    '''
    from pytraj import TrajectoryIterator
    state = CpptrajState()

    if not isinstance(traj, TrajectoryIterator):
        raise ValueError('traj must by TrajectoryIterator, '
                         'use pytraj.iterload(...)')

    txt0 = '''
    parm %s
    ''' % traj.top.filename

    for fname, frame_slice in zip(traj.filelist, traj.frame_slice_list):
        start, stop, stride = frame_slice
        if stop == -1:
            _stop  = 'last'
        elif stop < -1:
            raise RuntimeError('does not support negative stop for load_batch (except -1 (last))')
        else:
            _stop = stop
        txt0 += 'trajin {0} {1} {2} {3}'.format(fname, str(start), str(_stop), str(stride))

    lines = (txt0 + txt).lstrip().rstrip().split('\n')
    for line in lines:
        _Command.Dispatch(state.thisptr[0], line.encode())
    return state
