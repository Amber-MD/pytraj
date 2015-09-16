# distutils: language = c++
import os
import numpy as np
from ..api import Trajectory
from ..core.cpptraj_core cimport AtomMask
from ..Topology cimport Topology

from .._cyutils import get_positive_idx
from pytraj.externals.six import string_types, PY2
from .._shared_methods import my_str_method
from .._shared_methods import _xyz, _tolist
from .._shared_methods import _savetraj, _get_temperature_set
from .._shared_methods import _box
from ..utils.check_and_assert import ensure_exist
from ..utils.check_and_assert import is_array, is_range
from ..externals.six.moves import zip, range

def _split_range(int chunksize, int start, int stop):
    '''split a given range to n_chunks

    Examples
    --------
    >>> from pytraj.misc import split_range
    >>> split_range(3, 0, 10)
    [(0, 3), (3, 6), (6, 10)]
    '''
    cdef int n_chunks, i, _stop

    n_chunks = (stop - start)//chunksize

    if ((stop - start) % chunksize ) != 0:
        n_chunks += 1

    for i in range(n_chunks):
        if i < n_chunks - 1:
            _stop = start + (i + 1) * chunksize
        else:
            _stop = stop
        yield start + i * chunksize, _stop

cdef class TrajectoryCpptraj:
    def __cinit__(self):
        self.thisptr = new _TrajectoryCpptraj()
        self._top = Topology()
        self._filelist = []

    def load(self, filename=None, top=None, frame_slice=(0, -1, 1)):
        '''
        filename : a single filename or a list of filenames
        top : Topology-like object
        take_slice : add slice
        '''
        ensure_exist(filename)

        cdef Topology tmp_top
        cdef ArgList _arglist

        if top is None:
            tmp_top = <Topology> self.top
        else:
            tmp_top = <Topology> top

        # convert pytraj frame_slice to cpptraj's understandable format (str)
        # need to increase start, stop by +1 since cpptraj use +1 for cpptraj.in
        if not isinstance(frame_slice, tuple):
            if isinstance(frame_slice, list) and isinstance(frame_slice[0], tuple):
                frame_slice = frame_slice[0]
            else:
                raise ValueError("frame_slice must be a tuple or a list of tuple")
        if len(frame_slice) == 1:
            start = frame_slice[0]
            stop = -1
            stride = 1
        elif len(frame_slice) == 2:
            # no stride info
            start, stop = frame_slice
            stride = 1
        elif len(frame_slice) == 3:
            start, stop, stride = frame_slice
        else:
            raise ValueError()

        start += 1
        # don't increase stop by +1
        # slice(0, 10, None) --> python does not take last `10`
        arg = " ".join((str(start), str(stop), str(stride)))

        if isinstance(filename, string_types):
            # use absolute path so we can go to different folder
            filename = os.path.abspath(filename)
            _arglist = ArgList(arg)
            self.thisptr.AddSingleTrajin(filename.encode(), _arglist.thisptr[0], tmp_top.thisptr)
            self._filelist.append(os.path.abspath(filename))
        else:
            raise ValueError("filename must a a string")

    def __len__(self):
        return self.n_frames

    @property
    def n_frames(self):
        return self.thisptr.Size()

    @property
    def n_atoms(self):
        '''used for frame_iter'''
        return self.top.n_atoms

    def __str__(self):
        return my_str_method(self)

    def __repr__(self):
        return self.__str__()

    def __call__(self, *args, **kwd):
        return self.iterframe(*args, **kwd)

    def __iter__(self):
        '''iterately getting Frame instance
        '''
        cdef int i
        cdef int n_atoms = self.n_atoms
        cdef n_frames = self.n_frames

        # use `frame` as buffer 
        cdef Frame frame = Frame(n_atoms)

        for i in range(n_frames):
            # do not create new Frame inside this loop to reduce memory
            self.thisptr.GetFrame(i, frame.thisptr[0])
            yield frame

    property top:
        def __get__(self):
            self._top.thisptr[0] = self.thisptr.Top()
            return self._top

        def __set__(self, Topology other):
            self.thisptr.SetTopology(other.thisptr[0])

    def iterframe(self, int start=0, int stop=-1, int stride=1, mask=None):
        '''iterately get Frames with start, stop, stride 
        Parameters
        ---------
        start : int (default = 0)
        stop : int (default = max_frames)
        stride : int
        mask : str or array of interger
        '''
        cdef int i
        cdef int n_atoms = self.n_atoms
        cdef Frame frame = Frame()
        cdef AtomMask atm
        cdef int _end
        cdef int[:] int_view

        if stop == -1:
            _end = <int> self.n_frames
        else:
            _end = stop

        del frame.thisptr
        if mask is not None:
            if isinstance(mask, string_types):
                atm = self.top(mask)
            else:
                try:
                    atm = AtomMask()
                    atm.add_selected_indices(mask)
                except TypeError:
                    raise TypeError("dont know how to cast to memoryview")
            frame.thisptr = new _Frame(<int>atm.n_atoms)
        else:
            frame.thisptr = new _Frame(n_atoms)

        with self:
            i = start
            while i < _end:
                if mask is None:
                    self.thisptr.GetFrame(i, frame.thisptr[0])
                else:
                    self.thisptr.GetFrame(i, frame.thisptr[0], atm.thisptr[0])
                yield frame
                i += stride

    def iterchunk(self, int chunksize=2, int start=0, int stop=-1):
        '''iterately get Frames with start, chunk
        returning Trajectory or Frame instance depend on `chunk` value
        Parameters
        ---------
        start : int (default = 0)
        chunk : int (default = 1, return Frame instance). 
                if `chunk` > 1 : return Trajectory instance
        copy_top : bool, default=False
            if False: no Topology copy is done for new (chunk) Trajectory
        '''
        cdef int i, j, _stop
        cdef int n_frames = self.n_frames
        cdef int n_atoms = self.n_atoms
        cdef Frame frame
        cdef int _tmp_start, _tmp_stop, real_n_frames

        # check `start`
        if start < 0 or start >= n_frames:
            start = 0

        # check `stop`
        if stop <= 0 or stop >= n_frames:
            stop = <int> self.n_frames

        if chunksize <= 1:
            raise ValueError("chunk must be >= 2")

        if chunksize + start > stop:
            raise ValueError("start + chunk must be smaller than max frames")

        # only open and close file once
        with self:
            for (_tmp_start, _tmp_stop) in _split_range(chunksize, start, stop):
                # always create new Trajectory
                farray = Trajectory()
                farray.top = self.top.copy()
                real_n_frames = len(range(_tmp_start, _tmp_stop))
                farray._allocate(real_n_frames, farray.top.n_atoms)
                farray._boxes = np.empty((real_n_frames, 6), dtype='f8')

                for idx, frame in enumerate(self.iterframe(start=_tmp_start,
                    stop=_tmp_stop)):
                    farray._xyz[idx] = frame.xyz
                    farray._boxes[idx] = frame.box._get_data()
                yield farray
                    
    def __setitem__(self, idx, value):
        raise NotImplementedError("Read only Trajectory. Use Trajectory class for __setitem__")

    def __getitem__(self, idxs):
         # allocate frame for storing data
         cdef Frame frame0
         cdef Frame frame = Frame(self.top.n_atoms)
         cdef int start, stop, step
         cdef int i
         cdef int idx_1, idx_2
         cdef int[:] int_view
         cdef list tmplist
         cdef AtomMask atom_mask_obj
         cdef idxs_size
     
         if isinstance(idxs, AtomMask):
             # atm = top('@CA')
             # traj[atm]
             atom_mask_obj = <AtomMask> idxs
             _farray = Trajectory()
             _farray.top = self.top._modify_state_by_mask(atom_mask_obj)
             for i, frame in enumerate(self):
                 _frame = Frame(frame, atom_mask_obj)
                 _farray.append(_frame)
             self.tmpfarray = _farray
             # hold _farray in self.tmpfarray to avoid memory lost
             return self.tmpfarray
         elif isinstance(idxs, string_types):
             # return array with given mask
             # traj['@CA']
             mask = idxs
             try:
                 return self[self.top(mask)]
             except:
                 txt = "not supported keyword `%s`" % idxs
                 raise NotImplementedError(txt)
         elif isinstance(idxs, slice):
             start, stop, stride = idxs.indices(self.n_frames)
             self.tmpfarray = self._to_nptraj_by_indices(range(start, stop, stride))
             return self.tmpfarray
         else:
             # not is a slice
             if idxs == ():
                 return self
             elif isinstance(idxs, tuple):
                 idxs_size = len(idxs)
                 if idxs_size >= 4:
                     raise NotImplementedError("number of elements must me smaller than 4")
                 idx0 = idxs[0]
     
                 idx1 = idxs[1]
                 if isinstance(self[idx0], Frame):
                     frame = self[idx0]
                     self.tmpfarray = frame
                     if isinstance(idx1, string_types):
                         # traj[0, '@CA']
                         atm = self.top(idx1)
                         self.tmpfarray = Frame(frame, atm)
                         return self.tmpfarray
                     else:
                         frame.top = self.top
                         return self.tmpfarray[idxs[1:]]
                 elif isinstance(self[idx0], Trajectory):
                     farray = self[idx0]
                     self.tmpfarray = farray
                     if isinstance(idx1, AtomMask) or isinstance(idx1, string_types):
                         if idxs_size == 2:
                             return self.tmpfarray[idxs[1]]
                         else:
                             return self.tmpfarray[idxs[1]][idxs[2]]
                     else:
                         try:
                             return self.tmpfarray[idxs[1]]
                         except:
                             raise NotImplementedError()
             elif is_array(idxs) or isinstance(idxs, list) or is_range(idxs):
                 # traj[[2, 6, 3]]
                 # support indexing that having 'len'
                 if any(isinstance(x, bool) for x in idxs):
                     raise NotImplementedError("do not support bool indexing")
                 self.tmpfarray = self._to_nptraj_by_indices(idxs)
                 return self.tmpfarray

             else:
                 # traj[8]
                 # assuming that `idxs` is integer
                 idx_1 = <int> get_positive_idx(idxs, self.n_frames)
                 # raise index out of range
                 if idxs != 0 and idx_1 == 0:
                     raise ValueError("index is out of range")

                 with self:
                     self.thisptr.GetFrame(idx_1, frame.thisptr[0])
                 self.tmpfarray = frame
                 return self.tmpfarray

    def save(self, filename="", format='unknown', overwrite=True, *args, **kwd):
        _savetraj(self, filename, format, overwrite, *args, **kwd)

    def write(self, *args, **kwd):
        self.save(*args, **kwd)

    @property
    def unitcells(self):
        return _box(self)

    @property
    def xyz(self):
        '''return a copy of xyz coordinates (ndarray, shape=(n_frames, n_atoms, 3)
        We can not return a memoryview since Trajectory is a C++ vector of Frame object
        '''
        return _xyz(self)

    @property
    def filelist(self):
        '''return a list of input filenames. Order does matter'''
        return self._filelist

    def __enter__(self):
        return self

    def __exit__(self, *args):
        # cpptraj will take care of close/open file
        pass

    def __dealloc__(self):
        if self.thisptr:
            del self.thisptr

    def _to_nptraj_by_slice(self, int start, int stop, int stride):
        cdef int i, j
        cdef int n_atoms = self.n_atoms
        cdef Frame frame
        cdef double[:, :, ::1] xyz
        cdef int n_frames = len(range(start, stop, stride))

        traj = Trajectory()
        traj.xyz = np.zeros((n_frames, n_atoms, 3), dtype='f8')
        traj.unitcells = np.zeros((n_frames, 6), dtype='f8')
        traj.top = self.top
        xyz = traj.xyz[:]

        frame = Frame(n_atoms, xyz[0], _as_ptr=True)
        i = start
        j = 0
        while i < stop:
            # use `frame` as a pointer pointing to `xyz` memory
            # dump coords to xyz array
            frame.thisptr.SetXptr(frame.n_atoms, &xyz[j, 0, 0])
            self.thisptr.GetFrame(i, frame.thisptr[0])
            traj.unitcells[j] = frame.box._get_data()
            i += stride
            j += 1
        return traj

    def _to_nptraj_by_indices(self, indices):
        '''indices is iterable that has __len__
        '''
        cdef int i, j
        cdef int n_atoms = self.n_atoms
        cdef Frame frame
        cdef double[:, :, :] xyz
        cdef int n_frames = len(indices)

        traj = Trajectory()
        traj._allocate(n_frames, n_atoms)
        traj.unitcells = np.zeros((n_frames, 6), dtype='f8')
        traj.top = self.top
        xyz = traj.xyz

        frame = Frame(n_atoms, xyz[0], _as_ptr=True)
        for j, i in enumerate(indices):
            # use `frame` as a pointer pointing to `xyz` memory
            # dump coords to xyz array
            frame.thisptr.SetXptr(frame.n_atoms, &xyz[j, 0, 0])
            # copy coordinates of `self[i]` to j-th frame in `traj`
            self.thisptr.GetFrame(i, frame.thisptr[0])
            traj.unitcells[j] = frame.box._get_data()
        return traj

    def _iterframe_indices(self, frame_indices):
        cdef int i
        cdef Frame frame = Frame(self.n_atoms)

        for i in frame_indices:
            self.thisptr.GetFrame(i, frame.thisptr[0])
            yield frame

