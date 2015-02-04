#print print  distutils: language = c++
from cpython.array cimport array as pyarray
cimport cython
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as incr
from pytraj.Topology cimport Topology
from pytraj._utils cimport get_positive_idx

# python level
from pytraj.externals.six import string_types
from pytraj.TrajReadOnly import TrajReadOnly
from pytraj.utils.check_and_assert import _import_numpy
from pytraj.utils.check_and_assert import file_exist
from pytraj.trajs.Trajout import Trajout
from pytraj._save_traj import _save

# we don't allow sub-class in Python level since we will mess up with memory
@cython.final
cdef class FrameArray (object):
    def __cinit__(self, filename='', top=None, indices=None, 
                  bint warning=False, n_frames=None, flag=None):
        if isinstance(top, string_types):
            self.top = Topology(top)
        elif isinstance(top, Topology):
            self.top = top.copy()
        else:
            # create empty topology
            self.top = Topology()

        if n_frames is not None:
            # reserve n_frames
            self.resize(n_frames)

        self.oldtop = None

        self.warning = warning

        # since we are using memoryview for slicing this class istance, we just need to 
        # let `parent` free memory
        # this variable is intended to let FrameArray control 
        # freeing memory for Frame instance but it's too complicated
        #self.is_mem_parent = True
        if filename != "" and flag != 'hd5f':
            # TODO : check if file exist
            if not file_exist(filename):
                raise ValueError("There is not file having this name")
            if self.top.is_empty():
                raise ValueError("Need to have non-empty Topology")
            self.load(filename=filename, indices=indices)

    def copy(self):
        "Return a copy of FrameArray"
        cdef FrameArray other = FrameArray()
        cdef _Frame _frame
        for _frame in self.frame_v:
            other.frame_v.push_back(_frame)
        # copy self.top too
        other.top = self.top.copy()
        return other

    def __dealloc__(FrameArray self):
        """should we free memory for Frame instances here?
        (we set frame.py_free_mem = False in __getitem__)
        """
        #print "Test FrameArray exiting"
        pass
        #cdef Frame frame
        #if self.is_mem_parent:
        #    for frame in self:
        #        # we don't __dealloc__ here.
        #        # just turn py_free_mem on to let Frame class frees memory
        #        # work?
        #        # NO : Error in `python': double free or corruption (out)`
        #        # --> don't need this method. We still have the commented code here to 
        #        # remind not need to add in future.
        #        #frame.py_free_mem = True
        #        del frame.thisptr

    def __del__(self):
        cdef Frame frame
        for frame in self:
            del frame.thisptr

    def load(self, filename='', Topology top=None, indices=None):
        # TODO : add more test cases
        # should we add hdf5 format here?
        cdef Trajin_Single ts
        cdef int idx

        if isinstance(filename, string_types):
            # we don't use UTF-8 here since ts.load(filename) does this job
            #filename = filename.encode("UTF-8")
            ts = Trajin_Single()
            if top is not None:
                ts.top = top.copy()
                ts.load(filename)
                # update top for self too
                if not self.top.is_empty():
                    print "updating FrameArray topology"
                self.top = top.copy()
            else:
                # use self.top
                ts.top = self.top.copy()
                # this does not load whole traj into disk, just "prepare" to load
                ts.load(filename)

            if indices is None:
                # load all frames
                self.join(ts[:])
            elif isinstance(indices, slice):
                self.join(ts[indices])
            else:
                # indices is tuple, list, ...
                # we loop all traj frames and extract frame-ith in indices 
                # TODO : check negative indexing?
                # increase size of vector
                old_size = self.size
                indices_len = len(indices)

                # make more room for vector
                self.frame_v.resize(old_size + indices_len, _Frame())

                for idx in range(ts.size):
                    # TODO : make it faster
                    # really need?
                    if idx in indices:
                        # only work with `list` or `tuple`
                        idx_idx = indices.index(idx)
                        self[old_size + idx_idx] = ts[idx]

        elif isinstance(filename, (list, tuple)):
            # list of files
            for fh in filename:
                # recursive
                self.load(fh, top, indices)
        else:
            raise ValueError("can not load file/files")

    @property
    def shape(self):
        return (self.n_frames, self[0].n_atoms, 3)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def __getitem__(self, idxs):
        # TODO : same as Trajin class
        # should combine or inherit or ?
        #"""Return a reference of FrameArray[idx]
        #To get a copy
        #>>>frame = FrameArray_instance[10].copy()

        # TODO : why not using existing slice of list?

        cdef Frame frame = Frame(self.top.n_atoms)
        cdef FrameArray farray
        cdef int start, stop, step
        cdef int i
        cdef int idx_1, idx_2
        #cdef list tmplist

        # test memoryview for traj[:, :, :]
        cdef double[:, :, :] arr3d

        frame.py_free_mem = False

        if self.warning:
            print "return a Frame or sub-Framearray view of this instance"
            print "Use with care. For safetype, use `copy` method"

        if len(self) == 0:
            raise ValueError("Your FrameArray is empty, how can I index it?")

        if isinstance(idxs, string_types):
            # mimic API of MDtraj
            if idxs == 'coordinates':
                return self[:, :, :]
            elif idxs == 'topology':
                return self.top
            else:
                # return array with given mask
                # traj[':@CA']
                try:
                    # use `mask` to avoid confusion
                    mask = idxs
                    has_numpy, np = _import_numpy()
                    if not has_numpy:
                        raise NotImplementedError("must have numpy")
                    N = self.top(mask).n_selected
                    arr0 = np.empty(N*self.size*3).reshape(self.size, N, 3)
                    for i, frame in enumerate(self):
                        arr0[i] = frame[self.top(mask)]
                    return arr0
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
                # TODO : make memoryview for traj[:, :, :]
                if all_are_slice_instances:
                    # return 3D array or list of 2D arrays?
                    # traj[:, :, :]
                    # traj[1:2, :, :]
                    tmplist = []
                    for frame in self[idxs[0]]:
                        tmplist.append(frame[idxs[1:]])
                    if has_numpy:
                        # test memoryview, does not work yet.
                        # don't delete those line to remind we DID work on this
                        #arr3d = _np.empty(shape=_np.asarray(tmplist).shape)
                        #for i, frame in enumerate(self[idxs[0]]):
                        #    for j, f0 in enumerate(frame[idxs[1]]):
                        #        arr3d[i][j] = f0[:]
                        #return arr3d

                        return _np.asarray(tmplist)
                    else:
                        return tmplist
                if isinstance(self[idx_0], Frame):
                    frame = self[idx_0]
                    return frame[idxs[1:]]
                elif isinstance(self[idx_0], FrameArray):
                    farray = self[idx_0]
                    return farray[idxs[1:]]
                #return frame[idxs[1:]]
            else:
                idx_1 = get_positive_idx(idxs, self.size)
                # raise index out of range
                if idxs != 0 and idx_1 == 0:
                    # need to check if array has only 1 element. 
                    # arr[0] is  arr[-1]
                    if idxs != -1:
                        raise ValueError("index is out of range")
                # get memoryview
                #frame.thisptr = &(self.frame_v[idx_1])
                # change [] to `at` for bound-checking
                frame.thisptr = &(self.frame_v.at(idx_1))
                return frame
        else:
            if self.warning:
                print "FrameArray slice"
            # creat a subset array of `FrameArray`
            farray = FrameArray()
            # farray.is_mem_parent = False

            # should we make a copy of self.top or get memview?
            farray.top = self.top.copy()
            # create positive indexing for start, stop if they are None
            start, stop, step  = idxs.indices(self.size)
            
            # mimic negative step in python list
            # debug
            #print "before updating (start, stop, step) = (%s, %s, %s)" % (start, stop, step)
            if start > stop and (step < 0):
                # since reading TRAJ is not random access for large file, we read from
                # begining to the end and append Frame to FrameArray
                # we will reverse later after getting all needed frames
                # traj[:-1:-3]
                is_reversed = True
                # swap start and stop but adding +1 (Python does not take last index)
                # a = range(10) # a[5:1:-1] = [5, 4, 3, 2]
                # a[2:5:1] = [2, 3, 4, 5]
                start, stop = stop + 1, start + 1
                step *= -1
            else:
                is_reversed = False

            # debug
            #print "after updating (start, stop, step) = (%s, %s, %s)" % (start, stop, step)
      
            for i in range(start, stop, step):
                # turn `copy` to `False` to have memoryview
                # turn `copy` to `True` to make a copy
                farray.append(self[i], copy=True)
            if is_reversed:
                # reverse vector if using negative index slice
                # traj[:-1:-3]
                farray.reverse()

            # hold farray by self.tmpfarray object
            # so self[:][0][0] is still legit
            self.tmpfarray = farray
            #if self.tmpfarray.size == 1:
            #    return self.tmpfarray[0]
            return self.tmpfarray

    def __setitem__(self, idx, other):
        # TODO : add slice
        # make a copy
        # to make thing simple, we don't use fancy slicing here
        cdef Frame frame
        if len(self) == 0:
            raise ValueError("Your FrameArray is empty, how can I index it?")
        if isinstance(idx, (long, int)) and isinstance(other, Frame):
            frame = <Frame> other
            self.frame_v[idx] = frame.thisptr[0]
        # TODO : check this
        #self[idx].py_free_mem = False
        else:
            # example: self[0, 0, 0] = 100.
            self[idx[0]][idx[1:]] = other
            #txt = "not yet implemented. Try using framearray[idx1][idx2, idx3] = value"
            #raise NotImplementedError(txt)
        
    def __delitem__(FrameArray self, int idx):
        self.erase(idx)

    def __str__(FrameArray self):
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

    def reverse(self):
        # should we just create a fake operator?
        cpp_reverse(self.frame_v.begin(), self.frame_v.end())

    def erase(self, idxs):
        cdef int idx
        # dealloc frame pointer too?
        if isinstance(idxs, (int, long)):
            idx = idxs
            self.frame_v.erase(self.frame_v.begin() + idx)
        else:
            # assume : list, slice, iteratable object
            for idx in idxs:
                self.erase(idx)
        
    @property
    def size(self):
        return self.frame_v.size()

    def is_empty(self):
        return self.size == 0

    @property
    def n_frames(self):
        """same as self.size"""
        return self.size

    def __len__(self):
        return self.size

    def __iter__(self):
        """return a reference of Frame instance
        >>> for frame in FrameArray_instance:
        >>>     pass
                
        """
        cdef vector[_Frame].iterator it  = self.frame_v.begin()
        cdef Frame frame 

        while it != self.frame_v.end():
            frame = Frame()
            # use memoryview, don't let python free memory of this instance
            frame.py_free_mem = False
            frame.thisptr = &(deref(it))
            yield frame
            incr(it)

    def __add__(self, FrameArray other):
        self += other
        return self

    def __iadd__(self, FrameArray other):
        """
        append `other`'s frames to `self`
        frame0 += other
        """
        cdef _Frame _frame
        if self.top.n_atoms != other.top.n_atoms:
            raise ValueError("n_atoms of two arrays do not match")

        for _frame in other.frame_v:
            # make a copy
            self.frame_v.push_back(_Frame(_frame))
        return self

    def append(self, Frame framein, copy=True):
        cdef Frame frame = Frame() 
        if copy:
            frame = framein.copy()
            self.frame_v.push_back(frame.thisptr[0])
        else:
            # create memory view
            # need to set `py_free_mem` to False
            framein.py_free_mem = False
            self.frame_v.push_back(framein.thisptr[0])

    def join(self, *args):
        """join multiple FrameArray instances"""
        for arg in args:
            self._join(arg)

    cdef void _join(self, FrameArray other, mask=None):
        # TODO : do we need this method when we have `get_frames`
        if self.top.n_atoms != other.top.n_atoms:
            raise ValueError("n_atoms of two arrays do not match")
        if mask:
            raise NotImplementedError("not yet")
        self.frame_v.reserve(self.frame_v.size() + other.frame_v.size())
        self.frame_v.insert(self.frame_v.end(), 
                            other.frame_v.begin(), other.frame_v.end())

    def resize(self, int n_frames):
        self.frame_v.resize(n_frames)

    @property
    def temperatures(self):
        """return a Python array of temperatures"""
        cdef pyarray tarr = pyarray('d', [])

        for frame in self:
            tarr.append(frame.temperature)
        return tarr

    def get_frames(self, from_traj=None, indices=None, update_top=False, copy=True):
        # TODO : fater loading?
        """get frames from Trajin instance
        def get_frames(from_traj=None, indices=None, update_top=False, copy=True)
        Parameters:
        ----------
        from_traj : TrajReadOnly or FrameArray, default=None
        indices : default=None
        update_top : bool, default=False
        copy : bool, default=True

        Note:
        ----
        Have not support indices yet. Get max_frames from trajetory
        """
        
        cdef int i
        cdef int start, stop, step
        cdef Frame frame

        if from_traj is not None:
            ts = from_traj
            # append new frames to `self`
            if update_top:
                self.top = ts.top.copy()

            if not update_top:
                if self.top.n_atoms != ts.top.n_atoms:
                    raise ValueError("FrameArray.top.n_atoms should be equal to Trajin_Single.top.n_atoms or set update_top=True")

            if isinstance(ts, Trajin_Single) or isinstance(ts, TrajReadOnly):
                if indices is not None:
                    # slow method
                    # TODO : use `for idx in leng(indices)`?
                    if isinstance(indices, slice):
                        # use slice for saving memory
                        start, stop, step = indices.start, indices.stop, indices.step
                        for i in range(start, stop, step):
                            self.append(ts[i], copy=copy)
                    else:
                        # regular list, tuple, array,...
                        for i in indices:
                            print "debug FrameArray.get_frames"
                            self.append(ts[i], copy=copy)
                else:    
                    # get whole traj
                    frame = Frame()
                    frame.set_frame_v(ts.top, ts.has_vel(), ts.n_repdims)
                    ts.begin_traj()
                    for i in range(ts.max_frames):
                        ts.get_next_frame(frame)
                        self.append(frame, copy=copy)
                    ts.end_traj()

            #elif isinstance(ts, FrameArray2) or isinstance(ts, FrameArray):
            elif isinstance(ts, FrameArray):
                # TODO : rename FrameArray2
                # use try and except?
                if indices is None:
                    for i in range(ts.size):
                        # TODO : make indices as an array?
                        self.append(ts[i], copy=copy)
                else:
                    for i in indices:
                        # TODO : make indices as an array?
                        self.append(ts[i], copy=copy)

        else:
            # if from_traj is None, return new FrameArray
            newfarray = FrameArray()
            if update_top:
                newfarray.top = self.top.copy()
            for i in indices:
                newfarray.append(self[i], copy=copy)
            return newfarray

    def strip_atoms(self, mask=None, update_top=True, bint has_box=False):
        """if you use memory for numpy, you need to update after resizing Frame
        >>> arr0 = np.asarray(frame.buffer)
        >>> frame.strip_atoms(top,"!@CA")
        >>> # update view
        >>> arr0 = np.asarray(frame.buffer)
        """

        cdef vector[_Frame].iterator it
        cdef Frame frame = Frame()
        cdef Topology tmptop = Topology()

        if mask == None: 
            raise ValueError("Must provide mask to strip")
        mask = mask.encode("UTF-8")

        # do not dealloc since we use memoryview for _Frame
        frame.py_free_mem = False
        it = self.frame_v.begin()
        while it != self.frame_v.end():
            frame.thisptr = &(deref(it))
            # we need to update topology since _strip_atoms will modify it
            tmptop = self.top.copy()
            frame._strip_atoms(tmptop, mask, update_top, has_box)
            incr(it)
        if update_top:
            self.top = tmptop.copy()

    # taking from Trajin_Single
    @classmethod
    def write_options(cls):
        TrajReadOnly.write_options()

    @classmethod
    def read_options(cls):
        TrajReadOnly.read_options()

    def save(self, filename="", fmt='unknown', overwrite=False):
        _save(self, filename, fmt, overwrite)
        #if fmt == 'unknown':
        #    # convert to "UNKNOWN_TRAJ"
        #    fmt = fmt.upper() + "_TRAJ"
        #else:
        #    fmt = fmt.upper()

        #with Trajout(filename=filename, top=self.top, fmt=fmt, 
        #             overwrite=overwrite, more_args=None) as trajout:
        #    for idx, frame in enumerate(self):
        #        trajout.writeframe(idx, frame, self.top)
