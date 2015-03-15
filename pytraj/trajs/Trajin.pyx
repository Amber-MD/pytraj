# distutils: language = c++
from __future__ import absolute_import
import os
cimport cython
from cpython.array cimport array as pyarray
from pytraj._utils cimport get_positive_idx
from pytraj.FrameArray cimport FrameArray
from pytraj.AtomMask cimport AtomMask

from pytraj.utils.check_and_assert import _import_numpy
from .Trajout import Trajout
from pytraj._save_traj import _save
from pytraj.externals.six import string_types


cdef class Trajin (TrajectoryFile):

    def __cinit__(self):
        self.baseptr_1 = <_Trajin*> self.baseptr0
        self.debug = False
        pass

    def __dealloc__(self):
        pass

    def __enter__(self):
        self.begin_traj()
        return self

    def __exit__(self, arg1, arg2, arg3):
        self.end_traj()

    def __iter__(self):
        """call `with Trajin_instace` before using this iteration"""
        cdef Frame frame = Frame(self.top.n_atoms)
        cdef int i

        with self:
           for i in range(self.baseptr_1.TotalFrames()):
               # don't use Python method to avoid overhead
               #self.get_next_frame(frame)
               self.baseptr_1.GetNextFrame(frame.thisptr[0])
               yield frame

    def __call__(self, *args, **kwd):
        return self.frame_iter(*args, **kwd)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.profile(True)
    @cython.infer_types(True)
    def frame_iter(self, int start=0, int stop=-1, int stride=1, mask=None):
    #def frame_iter(self, int start=0, int stop=-1, int stride=1):
        """iterately get Frames with start, stop, stride 
        Parameters
        ---------
        start : int (default = 0)
        chunk : int (default = 1)
        stop : int (default = max_frames - 1)
        mask : str (default=None)
            take coords with given mask
        """
        cdef int i
        cdef Frame frame = Frame(self.n_atoms)
        cdef Frame frame2
        cdef AtomMask atm
        cdef int _end

        if stop == -1:
            _end = <int> self.n_frames 
        else:
            _end = stop + 1

        i = start
        with self:
            while i < _end:
                self.baseptr_1.ReadTrajFrame(i, frame.thisptr[0])
                if mask is not None:
                    atm = self.top(mask)
                    frame2 = Frame(atm.n_atoms)
                    frame2.thisptr.SetCoordinates(frame.thisptr[0], atm.thisptr[0])
                    yield frame2
                else:
                    yield frame
                i += stride

    def chunk_iter(self, int chunk=2, int start=0, int stop=-1):
        """iterately get Frames with start, chunk
        returning FrameArray or Frame instance depend on `chunk` value
        Parameters
        ---------
        start : int (default = 0)
        chunk : int (default = 1, return Frame instance). 
                if `chunk` > 1 : return FrameArray instance
        """
        cdef int newstart
        cdef int n_chunk, i 

        # check `start`
        if start < 0 or start >= self.n_frames:
            start = 0

        # check `stop`
        if stop <= 0 or stop >= self.n_frames:
            stop = <int> self.size - 1

        if chunk <= 1:
            raise ValueError("chunk must be >= 2")

        if chunk + start > stop:
            raise ValueError("start + chunk must be smaller than max frames")

        n_chunk = int((stop- start)/chunk)
        if ((stop - start) % chunk ) != 0:
            n_chunk += 1

        for i in range(n_chunk):
            if i != n_chunk - 1:
                yield self[start + chunk*i : start + chunk*(i+1)]
            else:
                # use `stop + 1` since Python ignore last index
                yield self[start + chunk*i : stop+1]

    def __str__(self):
        name = self.__class__.__name__
        n_atoms = 0 if self.top.is_empty() else self.top.n_atoms
        tmps = """%s instance with %s frames, %s atoms/frame
                  ID = %s
               """ % (
                name, self.size, n_atoms,
                hex(id(self))
                )
        return tmps

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

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def __getitem__(self, idxs):
        # allocate frame for storing data
        cdef Frame frame0
        cdef Frame frame = Frame(self.top.n_atoms)
        cdef FrameArray farray
        cdef int start, stop, step
        cdef int i
        cdef int idx_1, idx_2
        cdef list tmplist

        if isinstance(idxs, string_types):
            # mimic API of MDtraj
            if idxs == 'coordinates':
                return self[:, :, :]
            elif idxs == 'topology':
                return self.top
            else:
                # return array with given mask
                # traj[':@CA']
                # traj[':@CA :frame']
                try:
                    # use `mask` to avoid confusion
                    mask = idxs
                    if ':frame' not in mask:
                        # return numpy array
                        has_numpy, np = _import_numpy()
                        N = self.top(mask).n_atoms
                        if has_numpy:
                            arr0 = np.empty(N*self.size*3).reshape(self.size, N, 3)
                            for i, frame in enumerate(self):
                                arr0[i] = frame[self.top(mask)]
                            return arr0
                        else:
                            # create 3D list with shape of (n_frames, n_atoms, 3)
                            _coord_list = []
                            for frame in self:
                                _coord_list.append(frame[self.top(mask)])
                            return _coord_list
                    else:
                        _farray = FrameArray()
                        _farray.top = self.top.modify_state_by_mask(self.top(mask))
                        for i, frame in enumerate(self):
                            _frame = frame.get_subframe(mask, self.top)
                            _farray.append(_frame)
                        self.tmpfarray = _farray
                        # hold _farray in self.tmpfarray to avoid memory lost
                        return self.tmpfarray
                except:
                    txt = "not supported keyword `%s`" % idxs
                    raise NotImplementedError(txt)

        if not isinstance(idxs, slice):
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
                elif isinstance(self[idx_0], FrameArray):
                    farray = self[idx_0]
                    self.tmpfarray = farray
                    return self.tmpfarray[idxs[1:]]
            else:
                # assuming that `idxs` is integer
                idx_1 = get_positive_idx(idxs, self.size)
                # raise index out of range
                if idxs != 0 and idx_1 == 0:
                    raise ValueError("index is out of range")

                with self:
                    self.read_traj_frame(idx_1, frame)
                self.tmpfarray = frame
                return self.tmpfarray
        else:
            if self.debug:
                print idxs
            farray = FrameArray()
            # should we copy self.top or use memview?
            farray.top = self.top.copy()

            # check comment in FrameArray class with __getitem__ method
            start, stop, step = idxs.indices(self.size)
            if self.debug:
                print (start, stop, step)

            with self:
                if start > stop and (step < 0):
                    # traj[:-1:-3]
                    is_reversed = True
                    # see comment in FrameArray (__getitem__)
                    start, stop = stop + 1, start + 1
                    step *= -1
                else:
                    is_reversed = False

                for i in range(start, stop, step):
                    self.read_traj_frame(i, frame)
                    farray.append(frame)

                if is_reversed:
                    # reverse vector if using negative index slice
                    # traj[:-1:-3]
                    farray.reverse()

            # use tmpfarray to hold farray for nested indexing
            # if not, Python will free memory for sub-FrameArray 
            self.tmpfarray = farray
            return self.tmpfarray

    def __setitem__(self, idx, value):
        raise NotImplementedError("Read only Trajectory. Use FrameArray class for __setitem__")

    def is_empty(self):
        return self.max_frames == 0

    def check_allocated(self):
        # decorator?
        if self.is_empty():
            raise MemoryError("not allocate pointer yet or have empty traj")

    @classmethod
    def check_frame_args(cls, ArgList argIn, int maxFrames):
        cdef int startArg, stopArg, offsetArg 
        startArg = stopArg = offsetArg = 0
        _Trajin.CheckFrameArgs(argIn.thisptr[0], maxFrames, startArg, stopArg, offsetArg)
        return startArg, stopArg, offsetArg

    def get_next_frame(self, Frame frame):
        #cdef Frame frame = Frame()
        self.baseptr_1.GetNextFrame(frame.thisptr[0])
        #return frame

    def prepare_for_read(self,bint b):
        self.check_allocated()
        self.baseptr_1.PrepareForRead(b)

    property max_frames:
        def __get__(self):
            if self.baseptr_1:
                return self.baseptr_1.TotalFrames()
            else:
                return 0

        def __set__(self, int value):
            self.baseptr_1.SetTotalFrames(value)

    @property
    def size(self):
        # alias of max_frames
        return self.max_frames

    @property
    def total_read_frames(self):
        self.check_allocated()
        return self.baseptr_1.TotalReadFrames()

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

    def begin_traj(self, bint showProgress=False):
        return self.baseptr_1.BeginTraj(showProgress)

    def end_traj(self):
        self.baseptr_1.EndTraj()

    def read_traj_frame(self, int currentFrame, Frame frameIn):
        # TODO : add checking frame.n_atoms == self.top.n_atoms?
        return self.baseptr_1.ReadTrajFrame(currentFrame, frameIn.thisptr[0])

    def save(self, filename="", fmt='unknown', overwrite=False):
        _save(self, filename, fmt, overwrite)

    def write(self, *args, **kwd):
        self.save(*args, **kwd)

    def get_subframes(self, mask, indices=None):
        cdef FrameArray farray = FrameArray()
        raise NotImplementedError("not yet")

    @property
    def temperatures(self):
        """return a Python array of temperatures"""
        cdef pyarray tarr = pyarray('d', [])

        for frame in self:
            tarr.append(frame.temperature)
        return tarr

    def fit_to(self, ref=None):
        txt = """
        This is immutatble class. You can not use with fit_to
        You FrameArray class or you can iterate to get Frame (mutable)

        >>> farray = FrameArray()
        >>> for frame in traj:
        >>>     frame.fit_to(ref)
        >>>     farray.append(frame)
        >>>     # or do anything interesting with `frame`
        """
        __doc__ = txt
        raise NotImplementedError(txt)

    @property
    def shape(self):
        return (self.n_frames, self[0].n_atoms, 3)
