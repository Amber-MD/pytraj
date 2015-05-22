# distutils: language = c++
from __future__ import absolute_import
import os
cimport cython
from cython.parallel cimport prange
from cpython.array cimport array as pyarray
from .._utils cimport get_positive_idx
from ..Trajectory cimport Trajectory
from ..AtomMask cimport AtomMask

from ..decorators import memoize # cache
from ..externals.six import string_types
from .._shared_methods import my_str_method
from .._shared_methods import _xyz, _tolist
from .._shared_methods import _savetraj, _get_temperature_set
from .._shared_methods import _box_to_ndarray
from ..utils.check_and_assert import _import_numpy
from ..utils.check_and_assert import is_word_in_class_name
from ..utils.check_and_assert import is_array, is_range
from pytraj.externals.six.moves import range
from .Trajout import Trajout


cdef class Trajin (TrajectoryFile):

    def __cinit__(self):
        self.baseptr_1 = <_Trajin*> self.baseptr0
        self.debug = False
        pass

    def __dealloc__(self):
        pass

    def __enter__(self):
        self._begin_traj()
        return self

    def __exit__(self, arg1, arg2, arg3):
        self._end_traj()

    def __iter__(self):
        """call `with Trajin_instace` before using this iteration"""
        cdef Frame frame = Frame(self.top.n_atoms)
        cdef int i

        with self:
           for i in range(self.baseptr_1.TotalFrames()):
               # don't use Python method to avoid overhead
               #self._get_next_frame(frame)
               self.baseptr_1.GetNextFrame(frame.thisptr[0])
               yield frame

    def __call__(self, *args, **kwd):
        return self.frame_iter(*args, **kwd)

    def frame_iter(self, int start=0, int stop=-1, int stride=1, mask=None):
        """iterately get Frames with start, stop, stride 
        Parameters
        ---------
        start : int (default = 0)
        stop : int (default = max_frames - 1)
        stride : int
        mask : str or array of interger
        """
        cdef int i
        cdef int n_atoms = self.n_atoms
        cdef Frame frame
        cdef AtomMask atm
        cdef int _end
        cdef int[:] int_view

        if stop == -1:
            _end = <int> self.n_frames
        else:
            _end = stop + 1

        if mask is not None:
            frame2 = Frame() # just make a pointer
            if isinstance(mask, string_types):
                atm = self.top(mask)
            else:
                try:
                    atm = AtomMask()
                    atm.add_selected_indices(mask)
                except TypeError:
                    raise TypeError("dont know how to cast to memoryview")
            frame2.thisptr = new _Frame(<int>atm.n_atoms)
        else:
            frame = Frame(n_atoms)

        # only open/close file once
        with self:
            i = start
            while i < _end:
                self.baseptr_1.ReadTrajFrame(i, frame.thisptr[0])
                if mask is not None:
                    frame2.thisptr.SetCoordinates(frame.thisptr[0], atm.thisptr[0])
                    yield frame2
                else:
                    yield frame
                i += stride

    def chunk_iter(self, int chunk=2, int start=0, int stop=-1, bint copy_top=False):
        """iterately get Frames with start, chunk
        returning Trajectory or Frame instance depend on `chunk` value
        Parameters
        ---------
        start : int (default = 0)
        chunk : int (default = 1, return Frame instance). 
                if `chunk` > 1 : return Trajectory instance
        copy_top : bool, default=False
            if False: no Topology copy is done for new (chunk) Trajectory
        """
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
                    self.baseptr_1.ReadTrajFrame(j, frame.thisptr[0])
                    farray.append(frame, copy=False)
                yield farray

    def __str__(self):
        return my_str_method(self)

    def __repr__(self):
        return self.__str__()

    def __len__(self):
        return self.size

    @property
    def n_frames(self):
        """self.n_frames == self.size
        Use this to match MDTraj API
        """
        return self.size

    @property
    def n_atoms(self):
        return self.top.n_atoms

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
        elif not isinstance(idxs, slice):
            if isinstance(idxs, tuple):
                idx_0 = idxs[0]
    
                all_are_slice_instances = True
                for tmp in idxs:
                    if not isinstance(tmp, slice): all_are_slice_instances = False
    
                has_numpy, _np = _import_numpy()
                # got Segmentation fault if using "is_instance3 and not has_numpy"
                # TODO : Why?
                #if is_instance3 and not has_numpy:
                if all_are_slice_instances:
                    # traj[:, :, :]
                    # traj[1:2, :, :]
                    tmplist = []
                    for frame in self[idxs[0]]:
                        tmplist.append(frame[idxs[1:]])
                    if has_numpy:
                        return _np.asarray(tmplist)
                    else:
                        return tmplist
                    #raise NotImplementedError("not yet supported if all indcies are slices")
                if isinstance(self[idx_0], Frame):
                    frame = self[idx_0]
                    self.tmpfarray = frame
                    return self.tmpfarray[idxs[1:]]
                elif isinstance(self[idx_0], Trajectory):
                    farray = self[idx_0]
                    self.tmpfarray = farray
                    return self.tmpfarray[idxs[1:]]
            elif is_array(idxs) or isinstance(idxs, list) or is_range(idxs):
                farray = Trajectory()
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
                    self._read_traj_frame(idx_1, frame)
                self.tmpfarray = frame
                return self.tmpfarray
        else:
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
    
                for frame in self.frame_iter(start, stop-1, step):
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
            # I am not really sure about below comment (happend before but 
            # not sure if this is the reason.)
            # use tmpfarray to hold farray for nested indexing
            # if not, Python will free memory for sub-Trajectory 
            #self.tmpfarray = farray
            #return self.tmpfarray

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
                self.baseptr_1.ReadTrajFrame(count, _frame_ptr[0])
                farray.frame_v.push_back(_frame_ptr)
                count += step
        return farray

    def _fast_slice_openmp(self, slice my_slice):
        """
        don't try to use this method. will get segfault
        """
        cdef int start, stop, step
        cdef int count, i
        cdef int n_atoms = self.n_atoms
        cdef Trajectory farray = Trajectory(check_top=False)
        cdef _Frame* _frame_ptr
        cdef _Frame* _frame_ptr2
        cdef int new_size

        farray.top = self.top

        start, stop, step  = my_slice.indices(self.size)
        new_size = len(range(start, stop, step)) # better way?

        farray.frame_v.reserve(new_size)

        with self:
            # segfault?
            # YES: need to compile parallel version?
            for i in prange(start, stop, step, nogil=True):
                _frame_ptr = new _Frame(n_atoms)
                self.baseptr_1.ReadTrajFrame(i, _frame_ptr[0])
                _frame_ptr2 = farray.frame_v[i]
                _frame_ptr2  = &_frame_ptr[0]
        return farray


    def __setitem__(self, idx, value):
        raise NotImplementedError("Read only Trajectory. Use Trajectory class for __setitem__")

    def is_empty(self):
        return self.max_frames == 0

    def _check_allocated(self):
        # decorator?
        if self.is_empty():
            raise MemoryError("not allocate pointer yet or have empty traj")

    @classmethod
    def _check_frame_args(cls, ArgList argIn, int maxFrames):
        cdef int startArg, stopArg, offsetArg 
        startArg = stopArg = offsetArg = 0
        _Trajin.CheckFrameArgs(argIn.thisptr[0], maxFrames, startArg, stopArg, offsetArg)
        return startArg, stopArg, offsetArg

    def _get_next_frame(self, Frame frame):
        #cdef Frame frame = Frame()
        self.baseptr_1.GetNextFrame(frame.thisptr[0])
        #return frame

    property max_frames:
        def __get__(self):
            if self.baseptr_1:
                return self.baseptr_1.TotalFrames()
            else:
                return 0

        def __set__(self, int value):
            self.baseptr_1.SetTotalFrames(value)

    @property
    @memoize
    def size(self):
        # alias of max_frames
        return self.max_frames

    def load(self, tnameIn, Topology tparmIn, ArgList argIn):
        """
        Load trajectory from file.

        Parameters:
        filename :: string (trajectory file's name)
        ArgList instance
        Topology instance
        """
        tnameIn = tnameIn.encode("UTF-8")
        if os.path.isfile(tnameIn):
            return self.baseptr_1.SetupTrajRead(tnameIn, argIn.thisptr[0], tparmIn.thisptr)
        else:
            raise ValueError("File does not exist")

    def _begin_traj(self, bint showProgress=False):
        return self.baseptr_1.BeginTraj(showProgress)

    def _end_traj(self):
        self.baseptr_1.EndTraj()

    def close(self):
        self._end_traj()

    def _read_traj_frame(self, int currentFrame, Frame frameIn):
        # TODO : add checking frame.n_atoms == self.top.n_atoms?
        return self.baseptr_1.ReadTrajFrame(currentFrame, frameIn.thisptr[0])

    def save(self, filename="", fmt='unknown', overwrite=True):
        _savetraj(self, filename, fmt, overwrite)

    def write(self, *args, **kwd):
        self.save(*args, **kwd)

    def get_subframes(self, mask=None, indices=None):
        raise NotImplementedError()

    @property
    def temperatures(self):
        """return a Python array of temperatures"""
        cdef pyarray tarr = pyarray('d', [])

        for frame in self:
            tarr.append(frame.temperature)
        return tarr

    @property
    def temperature_set(self):
        return _get_temperature_set(self)

    def fit_to(self, ref=None):
        txt = """
        This is immutatble class. You can not use with fit_to
        You Trajectory class or you can iterate to get Frame (mutable)

        >>> farray = Trajectory()
        >>> for frame in traj:
        >>>     frame.fit_to(ref)
        >>>     farray.append(frame)
        >>>     # or do anything interesting with `frame`
        """
        __doc__ = txt
        raise NotImplementedError(txt)

    @property
    @memoize
    def shape(self):
        return (self.size, self[0].n_atoms, 3)

    @property
    #@memoize
    def xyz(self):
        """return a copy of xyz coordinates (ndarray, shape=(n_frames, n_atoms, 3)
        We can not return a memoryview since Trajectory is a C++ vector of Frame object
        """
        # NOTE: turn off `memoize`
        # xyz = traj.xyz[:]
        # xyz += 1.
        # print (xyz[0, 0])
        # print (traj.xyz[0, 0])
        return _xyz(self)

    @memoize
    def tolist(self):
        return _tolist(self)

    @property
    @memoize
    def coordinfo(self):
        cdef CoordinateInfo cinfo = CoordinateInfo()
        cinfo.thisptr[0] = self.baseptr_1.TrajCoordInfo()
        return cinfo

    @memoize
    def box_to_ndarray(self):
        return _box_to_ndarray(self)
