# distutils: language = c++
from .._utils cimport get_positive_idx
from ..Trajectory cimport Trajectory
from ..AtomMask cimport AtomMask
from ..Topology cimport Topology

from pytraj.externals.six import string_types
from ..decorators import memoize # cache
from .._shared_methods import my_str_method
from .._shared_methods import _xyz, _tolist
from .._shared_methods import _savetraj, _get_temperature_set
from .._shared_methods import _box_to_ndarray
from ..utils.check_and_assert import _import_numpy
from ..utils.check_and_assert import is_word_in_class_name
from ..utils.check_and_assert import is_array, is_range
from .._get_common_objects import _get_top
from .Trajout import Trajout


cdef class TrajectoryCpptraj:
    def __cinit__(self):
        self.thisptr = new _TrajectoryCpptraj()
        self._top = Topology()
        self._filelist = []

    def __dealloc__(self):
        if self.thisptr:
            del self.thisptr

    def __enter__(self):
        return self

    def __exit__(self, *args):
        # cpptraj will take care of close/open file
        pass

    def load(self, filename=None, top=None, frame_slice=(0, -1, 1)):
        '''
        filename : a single filename or a list of filenames
        top : Topology-like object
        take_slice : add slice
        '''
        cdef Topology tmp_top
        cdef ArgList _arglist

        if top is None:
            tmp_top = <Topology> self.top
        else:
            tmp_top = <Topology> top

        # convert pytraj frame_slice to cpptraj's understandable format (str)
        # need to increase start, stop by +1 since cpptraj use +1 for cpptraj.in
        if not isinstance(frame_slice, tuple):
            raise ValueError("frame_slice must be a tuple")
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
            _arglist = ArgList(arg)
            self.thisptr.AddSingleTrajin(filename.encode(), _arglist.thisptr[0], tmp_top.thisptr)
            self._filelist.append(filename)
        else:
            raise ValueError("filename must a a string")

    def load_new(self, *args, **kwd):
        '''
        remove all trajectory data and load new ones
        '''
        saved_top = self.top
        del self.thisptr
        self._filelist = []
        self.thisptr = new _TrajectoryCpptraj()
        self.top = saved_top
        self.load(*args, **kwd)

    def _add_trajin(self, Trajin trajin):
        '''add memoryview for input trajin'''
        self.thisptr.AddInputTraj(trajin.baseptr_1)

    def __len__(self):
        return self.n_frames

    @property
    def n_frames(self):
        return self.size

    @property
    def size(self):
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
        return self.frame_iter(*args, **kwd)

    def __iter__(self):
        '''iterately getting Frame instance
        '''
        cdef int i
        cdef int n_atoms = self.n_atoms
        cdef n_frames = self.n_frames
        cdef Frame frame = Frame(n_atoms)

        for i in range(n_frames):
            self.thisptr.GetFrame(i, frame.thisptr[0])
            yield frame

    property top:
        def __get__(self):
            self._top.thisptr[0] = self.thisptr.Top()
            return self._top

        def __set__(self, Topology other):
            self.thisptr.SetTopology(other.thisptr[0])

    def frame_iter(self, int start=0, int stop=-1, int stride=1, mask=None):
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

    def chunk_iter(self, int chunk=2, int start=0, int stop=-1, bint copy_top=False):
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
        cdef int n_chunk, i, j, _stop
        cdef int n_frames = self.n_frames
        cdef int n_atoms = self.n_atoms
        cdef Trajectory farray
        cdef Frame frame

        # check `start`
        if start < 0 or start >= n_frames:
            start = 0

        # check `stop`
        if stop <= 0 or stop >= n_frames:
            stop = <int> self.size - 1

        if chunk <= 1:
            raise ValueError("chunk must be >= 2")

        if chunk + start > stop:
            raise ValueError("start + chunk must be smaller than max frames")

        n_chunk = int((stop- start)/chunk)
        if ((stop - start) % chunk ) != 0:
            n_chunk += 1

        # only open and close file once
        with self:
            for i in range(n_chunk):
                # always create new Trajectory
                farray = Trajectory(check_top=False)
                if copy_top:
                    farray.top = self.top.copy()
                else:
                    farray.top = self.top
                    farray.top.py_free_mem = False # let `self` do it

                if i != n_chunk - 1:
                    _stop = start + chunk*(i+1)
                else:
                    _stop = stop + 1

                for j in range(start + chunk * i,  _stop):
                    frame = Frame(n_atoms)
                    self.thisptr.GetFrame(j, frame.thisptr[0])
                    farray.append(frame, copy=False)
                yield farray

    def __setitem__(self, idx, value):
        raise NotImplementedError("Read only Trajectory. Use Trajectory class for __setitem__")

    def __getitem__(self, idxs):
         # allocate frame for storing data
         cdef Frame frame0
         cdef Frame frame = Frame(self.top.n_atoms)
         cdef Trajectory farray
         cdef int start, stop, step
         cdef int i
         cdef int idx_1, idx_2
         cdef int[:] int_view
         cdef list tmplist
         cdef AtomMask atom_mask_obj
         cdef idxs_size
     
         if isinstance(idxs, AtomMask):
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
             # mimic API of MDtraj
             if idxs == 'coordinates':
                 return self[:, :, :]
             elif idxs == 'topology':
                 return self.top
             else:
                 # return array with given mask
                 # traj[':@CA']
                 # traj[':@CA :frame']
                 # use `mask` to avoid confusion
                 mask = idxs
                 try:
                     return self[self.top(mask)]
                 except:
                     txt = "not supported keyword `%s`" % idxs
                     raise NotImplementedError(txt)

         elif isinstance(idxs, slice):
             # idxs is slice
             farray = Trajectory()
             # NOTE: MUST make a copy self.top
             # if NOT: double-free memory when using `_fast_strip_atoms`
             #farray.top = self.top
             farray.top = self.top.copy()
     
             # check comment in Trajectory class with __getitem__ method
             start, stop, step = idxs.indices(self.size)
     
             with self:
                 if start > stop and (step < 0):
                     # traj[:-1:-3]
                     is_reversed = True
                     # see comment in Trajectory (__getitem__)
                     start, stop = stop + 1, start + 1
                     step *= -1
                 else:
                     is_reversed = False
     
                 for frame in self.frame_iter(start, stop, step):
                     # add '-1' to stop 
                     # in `frame_iter`, we include `stop`
                     # but for slicing, python does not include stop
                     # always use `copy=True` since we are taking frames from 
                     # read-only Trajectory, no memview
                     farray.append(frame, copy=True)
     
                 if is_reversed:
                     # reverse vector if using negative index slice
                     # traj[:-1:-3]
                     farray.reverse()
             return farray
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
                         frame.set_top(self.top)
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
                 farray = Trajectory()
                 if hasattr(idxs, '__getitem__'):
                    if isinstance(idxs[0], bool):
                        raise NotImplementedError("don't support bool indexing")
                 # for unknown reason, this does not work, it returns a Frame (?)
                 farray.get_frames(from_traj=self, indices=idxs, update_top=True)
                 # need to use `farray` so Cython knows its type
                 self.tmpfarray = farray
                 return self.tmpfarray

             else:
                 # assuming that `idxs` is integer
                 idx_1 = <int> get_positive_idx(idxs, self.size)
                 # raise index out of range
                 if idxs != 0 and idx_1 == 0:
                     raise ValueError("index is out of range")

                 with self:
                     self.thisptr.GetFrame(idx_1, frame.thisptr[0])
                 self.tmpfarray = frame
                 return self.tmpfarray

    def _fast_slice(self, slice my_slice):
        cdef int start, stop, step
        cdef int count
        cdef int n_atoms = self.n_atoms
        cdef Trajectory farray = Trajectory(check_top=False)
        cdef _Frame* _frame_ptr

        # NOTE: MUST make a copy self.top
        # if NOT: double-free memory when using `_fast_strip_atoms`
        #farray.top = self.top
        farray.top = self.top.copy()

        start, stop, step  = my_slice.indices(self.size)
        count = start
        with self:
            while count < stop:
                _frame_ptr = new _Frame(n_atoms)
                self.thisptr.GetFrame(count, _frame_ptr[0])
                farray.frame_v.push_back(_frame_ptr)
                count += step
        return farray

    def save(self, filename="", format='unknown', overwrite=True, *args, **kwd):
        _savetraj(self, filename, format, overwrite, *args, **kwd)

    def write(self, *args, **kwd):
        self.save(*args, **kwd)

    def box_to_ndarray(self):
        return _box_to_ndarray(self)

    def to_mutable_traj(self):
        '''same as self[:] but more explicit'''
        return self[:]

    @property
    def xyz(self):
        '''return a copy of xyz coordinates (ndarray, shape=(n_frames, n_atoms, 3)
        We can not return a memoryview since Trajectory is a C++ vector of Frame object
        '''
        # xyz = traj.xyz[:]
        # xyz += 1.
        # print (xyz[0, 0])
        # print (traj.xyz[0, 0])
        return _xyz(self)

    @property
    def filelist(self):
        return self._filelist
