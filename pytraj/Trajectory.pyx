# distutils: language = c++
from __future__ import absolute_import, division
cimport cython
from cpython.array cimport array as pyarray
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as incr
from cython.parallel cimport prange, parallel
from libc.string cimport memcpy
from .Topology cimport Topology
from .AtomMask cimport AtomMask
from .Frame cimport Frame
from .trajs.Trajin cimport Trajin
from .math.Matrix_3x3 cimport Matrix_3x3
from .cpp_algorithm cimport iter_swap

# python level
from ._cyutils import get_positive_idx, _int_array1d_like_to_memview
from ._set_silent import set_error_silent
from .trajs.Trajin_Single import Trajin_Single
from .externals.six import string_types, PY2
from .externals._load_pseudo_parm import load_pseudo_parm
from .TrajectoryIterator import TrajectoryIterator
from .utils.check_and_assert import (_import_numpy, is_int, is_frame_iter,
                                     file_exist, is_mdtraj, is_pytraj_trajectory,
                                     is_word_in_class_name,
                                     is_array, is_range)
from .trajs.Trajout import Trajout
from ._get_common_objects import _get_top, _get_data_from_dtype
from ._shared_methods import (_savetraj, _get_temperature_set,
                              _xyz, _tolist, _split_and_write_traj)
from ._shared_methods import my_str_method, _box_to_ndarray
from ._xyz import XYZ

from . import common_actions as pyca
from .hbonds import search_hbonds

__all__ = ['Trajectory']


cdef class Trajectory (object):
    def __cinit__(self, filename=None, top=None, indices=None, 
            bint warning=False, n_frames=None, check_top=True):
        """
        Parameters
        ----------
        filename: str or Trajectory-like or array-like
            str : filename
            Trajectory-like: pytraj's Trajectory, TrajectoryIterator, mdtraj's Trajectory,
                DataSet_Coords_CRD, DataSet_Coords_TRJ
        top : str or Topology, default=None
        indices : array-like, frames to take, default=None
        warning : bool, default=False
            for debuging
        n_frames : int, optional, default=None
            preallocate n_frames
        check_top : bool, optional, default=True, don't check Topology

        Examples
        --------
            traj = Trajectory()
            traj = Trajectory("md.x", "prmtop")
            traj = Trajectory("md.x", t2.top)
            traj = Trajectory(xyz, t2.top) # create new Trajectory with given `xyz` array
            traj = Trajectory(n_frames=100) # preallocate 100 frames
            traj = Trajectory(check_top=False) # don't check any Topology to save time
        """
        
        cdef Frame frame

        if check_top:
            try:
                self.top = _get_top(filename, top)
            except:
                if is_mdtraj(filename):
                    self.top = load_pseudo_parm(filename.top)
                else:
                    raise ValueError()
            if self.top is None:
                self.top = Topology()
        else:
            self.top = Topology()

        if n_frames is not None:
            # reserve n_frames
            self.resize(n_frames)

        self.oldtop = None
        self.warning = warning

        # since we are using memoryview for slicing this class istance, we just need to 
        # let `parent` free memory
        # this variable is intended to let Trajectory control 
        # freeing memory for Frame instance but it's too complicated
        #self.is_mem_parent = True
        if filename is not None:
            self.load(filename, self.top, indices)

    def copy(self):
        "Return a copy of Trajectory"
        cdef Trajectory other = Trajectory()
        cdef Frame frame

        other.top = self.top.copy()

        for frame in self:
            other.append(frame, copy=True)
        return other

    def __dealloc__(self):
        """should we free memory for Frame instances here?
        (we set frame.py_free_mem = False in __getitem__)
        """
        #print "Test Trajectory exiting"
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

    def __array__(self):
        raise NotImplementedError("pytraj.Trajectory does not have buffer interface")

    def __del__(self):
        """deallocate all frames"""
        cdef Frame frame
        for frame in self:
            del frame.thisptr

    def __call__(self, *args, **kwd):
        """return frame_iter"""
        return self.frame_iter(*args, **kwd)

    def load(self, filename='', Topology top=None, indices=None):
        # TODO : add more test cases
        # should we add hdf5 format here?
        #cdef Trajin_Single ts
        cdef int idx
        cdef Frame frame
        cdef Trajin trajin

        if top is not None:
            if self.top.is_empty():
                self.top = top.copy()
            else:
                pass
            # don't update top if not self.top.is_empty()
        else:
            if self.top.is_empty():
                # if both top and self.top are empty, need to raise ValueError
                try:
                    tmpobj = filename
                    if hasattr(tmpobj, 'top'):
                        self.top = tmpobj.top.copy()
                    elif hasattr(tmpobj[0], 'top'):
                        self.top = tmpobj[0].top.copy()
                except:
                    raise ValueError("need to have non-empty Topology")

        # always use self.top
        if isinstance(filename, string_types):
            # load from single filename
            # we don't use UTF-8 here since ts.load(filename) does this job
            #filename = filename.encode("UTF-8")
            ts = Trajin_Single()
            ts.top = self.top.copy()
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
                for idx in indices:
                    self.append(ts[idx], copy=True) # copy=True because we load from disk
        elif isinstance(filename, Frame):
            self.append(filename)
        elif isinstance(filename, (list, tuple)):
            # load from a list/tuple of filenames
            # or a list/tuple of numbers
            _f0 = filename[0]
            if isinstance(_f0, string_types) or hasattr(_f0, 'n_frames'):
                # need to check `string_types` since we need to load list of numbers too.
                # list of filenames
                list_of_files_or_trajs = filename
                for fh in list_of_files_or_trajs:
                    if self.warning:
                        print ("Loading from list/tuple. Ignore `indices`")
                    # recursive
                    self.load(fh, self.top, indices)
            else:
                # load xyz
                try:
                    _xyz = filename
                    self.append_xyz(_xyz)
                except:
                    raise ValueError("must be a list/tuple of either filenames/Traj/numbers")
        elif hasattr(filename, 'n_frames') and not is_mdtraj(filename):
            # load from Traj-like object
            # make temp traj to remind about traj-like
            traj = filename
            if indices is None:
                for frame in traj:
                    self.append(frame)
            else:
                for idx, frame in enumerate(traj):
                    # slow method.
                    if idx in indices:
                        self.append(frame)
        elif is_frame_iter(filename):
            # load from frame_iter
            _frame_iter = filename
            for frame in _frame_iter:
                self.append(frame)
        elif is_mdtraj(filename):
            _traj = filename
            # add "10 *" since mdtraj use 'nm' while pytraj use 'Angstrom'
            #self.append_ndarray(10 * _traj.xyz)
            # keep original coorsd, don't cast
            self.append_ndarray(_traj.xyz)
        elif is_word_in_class_name(filename, 'DataSetList'):
            # load DataSetList
            # iterate all datasets and get anything having frame_iter
            dslist = filename
            for _d0 in dslist:
                if hasattr(_d0, 'frame_iter'):
                    _d0.top = self.top.copy()
                    # don't let _d0 free memory since we use Topology 'view'
                    for frame in _d0.frame_iter():
                        self.append(frame)
        else:
            try:
                # load from array
                _xyz = filename
                self.append_xyz(_xyz)
            except:
                raise ValueError("filename must be str, traj-like or numpy array")

    @cython.infer_types(True)
    @cython.cdivision(True)
    def append_xyz(self, xyz_in):
        cdef int n_atoms = self.top.n_atoms
        cdef int natom3 = n_atoms * 3
        cdef int n_frames, i 
        """Try loading xyz data with 
        shape=(n_frames, n_atoms, 3) or (n_frames, n_atoms*3) or 1D array

        If using numpy array with shape (n_frames, n_atoms, 3),
        try "append_ndarray" method (much faster)
        """

        if n_atoms == 0:
            raise ValueError("n_atoms = 0: need to set Topology or use `append_ndarray`'")

        has_np, np = _import_numpy()
        if has_np:
            xyz = np.asarray(xyz_in)
            if len(xyz.shape) == 1:
                n_frames = int(xyz.shape[0]/natom3)
                _xyz = xyz.reshape(n_frames, natom3) 
            elif len(xyz.shape) in [2, 3]:
                _xyz = xyz
            else:
                raise NotImplementedError("only support array/list/tuples with ndim=1,2,3")
            for arr0 in _xyz:
                frame = Frame(n_atoms)
                # flatten either 1D or 2D array
                frame.set_from_crd(arr0.flatten())
                self.append(frame)
        else:
            if isinstance(xyz_in, (list, tuple)):
                xyz_len = len(xyz_in)
                if xyz_len % (natom3) != 0:
                    raise ValueError("Len of list must be n_frames*n_atoms*3")
                else:
                    n_frames = int(xyz_len / natom3)
                    for i in range(n_frames):
                        frame = Frame(n_atoms)
                        frame.set_from_crd(xyz_in[natom3 * i : natom3 * (i + 1)])
                        self.append(frame)
            elif hasattr(xyz_in, 'memview'):
                    frame = Frame(n_atoms)
                    for i in range(xyz_in.shape[0]):
                        frame.append_xyz(xyz_in[i]) 
                        self.append(frame)
            else:
                raise NotImplementedError("must have numpy or list/tuple must be 1D")

    def append_ndarray(self, xyz):
        """load ndarray with shape=(n_frames, n_atoms, 3)"""
        cdef Frame frame
        cdef int i
        cdef double[:, :] myview
        cdef int n_frames = xyz.shape[0]
        cdef int n_atoms = xyz.shape[1]
        cdef int oldsize = self.frame_v.size()
        cdef int newsize = oldsize + n_frames
        import numpy as np

        # need to use double precision
        if xyz.dtype != np.float64:
            _xyz = xyz.astype(np.float64)
        else:
            _xyz = xyz
        self.frame_v.resize(newsize)

        for i in range(n_frames):
            # make memoryview for ndarray
            myview = _xyz[i]
            # since we use `vector[Frame*]`, we need to allocate Frame's size
            self[i + oldsize] = Frame(n_atoms)
            frame = self[i + oldsize]
            # copy coords
            frame._fast_copy_from_xyz(myview[:])

    @property
    def shape(self):
        return (self.n_frames, self[0].n_atoms, 3)

    @property
    def xyz(self):
        """return a copy of xyz coordinates (wrapper of ndarray, shape=(n_frames, n_atoms, 3)
        We can not return a memoryview since Trajectory is a C++ vector of Frame object

        Notes
        -----
            read-only
        """
        cdef bint has_numpy
        cdef int i
        cdef int n_frames = self.n_frames
        cdef int n_atoms = self.n_atoms
        cdef Frame frame

        has_numpy, np = _import_numpy()
        myview = np.empty((n_frames, n_atoms, 3), dtype='f8')

        if self.n_atoms == 0:
            raise NotImplementedError("need to have non-empty Topology")
        if has_numpy:
            for i, frame in enumerate(self):
                myview[i] = frame.buffer2d
            return XYZ(myview)
        else:
            raise NotImplementedError("must have numpy")

    @property
    def _xyz(self):
        """return a copy of xyz coordinates (wrapper of ndarray, shape=(n_frames, n_atoms, 3)
        We can not return a memoryview since Trajectory is a C++ vector of Frame object

        Notes
        -----
            read-only
        """
        import numpy as np
        cdef:
            int i
            int n_frames = self.n_frames
            int n_atoms = self.n_atoms
            int count
            _Frame _frame
            double *ptr_src
            double *ptr_dest

        cdef double[:, :] memview = cython.view.array(shape=(n_frames, n_atoms * 3),
                                             itemsize=sizeof(double),
                                             format='d')

        for i in range(n_frames):
            _frame = self.frame_v[i][0]
            ptr_src = _frame.xAddress()
            ptr_dest = &memview[i, 0]
            count = n_atoms * 3 * sizeof(double)
            memcpy(<void*> ptr_dest, <void*> ptr_src, count)

        return np.asarray(memview).reshape(n_frames, n_atoms, 3)

    @property
    def coordinates(self):
        """return 3D numpy.ndarray, same as `TrajectoryIterator.xyz`
        """
        return self.xyz

    def update_coordinates(self, double[:, :, :] xyz):
        '''update coords from 3D xyz array, dtype=f8'''
        cdef int idx, n_frames
        cdef double* ptr_src
        cdef double* ptr_dest
        cdef size_t count

        n_frames = xyz.shape[0]
        n_atoms = xyz.shape[1]
        count = sizeof(double) * n_atoms * 3

        for idx in range(n_frames):
            ptr_dest = self.frame_v[idx].xAddress()
            ptr_src = &(xyz[idx, 0, 0])
            memcpy(<void*> ptr_dest, <void*> ptr_src, count)

    def update_xyz(self, double[:, :, :] xyz):
        self.update_coordinates(xyz)

    def tolist(self):
        return _tolist(self)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def __getitem__(self, idxs):
        # TODO : same as Trajin class
        # should combine or inherit or ?
        #"""Return a reference of Trajectory[idx]
        #To get a copy
        #>>>frame = Trajectory_instance[10].copy()

        # TODO : why not using existing slice of list?

        cdef Frame frame = Frame(self.top.n_atoms) # need to allocate here?
        cdef Frame _frame # used for AtomMask selection. will allocate mem later
        cdef Trajectory farray
        cdef int start, stop, step, count
        cdef int i, j
        cdef int idx_1, idx_2
        cdef int[:] int_view
        cdef AtomMask atom_mask_obj
        cdef pyarray list_arr
        cdef int idxs_size

        # test memoryview for traj[:, :, :]
        cdef double[:, :, :] arr3d

        frame.py_free_mem = False

        if self.warning:
            print "return a Frame or sub-Trajectory view of this instance"
            print "Use with care. For safetype, use `copy` method"

        if len(self) == 0:
            raise ValueError("Your Trajectory is empty, how can I index it?")

        elif isinstance(idxs, AtomMask):
            atom_mask_obj = <AtomMask> idxs
            _farray = Trajectory(check_top=False) # just create naked Trajectory
            set_error_silent(True) # turn off cpptraj' verbose
            _farray.top = self.top._modify_state_by_mask(atom_mask_obj)
            set_error_silent(False)
            for i, frame in enumerate(self):
                _frame = Frame(frame, atom_mask_obj) # 1st copy
                _frame.py_free_mem = False #
                _farray.append(_frame, copy=False) # 2nd copy if using `copy=True`
            return _farray

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
                    atom_mask_obj = self.top(mask)
                    return self[atom_mask_obj]
                except:
                    txt = "not supported keyword `%s` or there's proble with your topology" % idxs
                    raise NotImplementedError(txt)

        elif isinstance(idxs, slice):
            # is slice
            # creat a subset array of `Trajectory`
            #farray = Trajectory()
            # farray.is_mem_parent = False

            # should we make a copy of self.top or get memview?
            #farray.top = self.top
            #farray.top.py_free_mem = False # let `master` Trajectory do freeing mem
            # create positive indexing for start, stop if they are None
            start, stop, step  = idxs.indices(self.size)
            
            # mimic negative step in python list
            # debug
            #print "before updating (start, stop, step) = (%s, %s, %s)" % (start, stop, step)
            if start > stop and (step < 0):
                # since reading TRAJ is not random access for large file, we read from
                # begining to the end and append Frame to Trajectory
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
      
            farray = self._fast_slice((start, stop, step))
            #i = start
            #while i < stop:
            #    # turn `copy` to `False` to have memoryview
            #    # turn `copy` to `True` to make a copy
            #    farray.append(self[i], copy=False)
            #    i += step
            if is_reversed:
                # reverse vector if using negative index slice
                # traj[:-1:-3]
                farray.reverse()

            # hold farray by self.tmpfarray object
            # so self[:][0][0] is still legit (but do we really need this with much extra memory?)
            #self.tmpfarray = farray
            #if self.tmpfarray.size == 1:
            #    return self.tmpfarray[0]
            #return self.tmpfarray
            return farray

        else:
            # not slice
            if idxs == ():
                # empty tuple
                return self
            if isinstance(idxs, tuple):
                idxs_size = len(idxs)
                idx1 = idxs[1]
                if idxs_size > 3:
                    raise NotImplementedError("support indexing up to 3 elements")
                idx0 = idxs[0]

                if isinstance(self[idx0], Frame):
                    frame = self[idx0]
                    frame.py_free_mem = False
                    if isinstance(idx1, string_types):
                        # traj[0, '@CA']
                        frame.set_top(self.top)
                    # TODO: need to check memory
                    if idxs_size == 2:
                        return frame[idxs[1]]
                    else:
                        return frame[idxs[1]][idxs[2:]]

                elif isinstance(self[idx0], Trajectory):
                    farray = self[idx0]
                    # place holder to avoid memory free
                    # atm = traj.top("@CA")
                    # traj[0, atm]
                    #if isinstance(idx1, AtomMask) or isinstance(idx1, string_types):
                    #    if idxs_size == 2:
                    #        return farray[idxs[1]]
                    #    else:
                    #        return farray[idxs[1]][idxs[2]]
                    #else:
                    #    try:
                    #        return farray[idxs[1:]]
                    #    except:
                    #        raise NotImplementedError("")
                    if idxs_size == 2:
                        return farray[idxs[1]]
                    else:
                        return farray[idxs[1]][idxs[2]]

            elif is_array(idxs) or isinstance(idxs, list) or is_range(idxs):
                _farray = Trajectory(check_top=False)
                _farray.top = self.top # just make a view, don't need to copy Topology
                # check if there are bool index
                if isinstance(idxs, list):
                    if isinstance(idxs[0], bool):
                        import numpy as np
                        idxs = np.array(idxs)
                if hasattr(idxs, 'dtype') and idxs.dtype.name == 'bool':
                    for i in range(idxs.shape[0]):
                        if idxs[i] == True:
                            frame.thisptr = self.frame_v[i] # point to i-th item
                            frame.py_free_mem = False # don't free mem
                            _farray.frame_v.push_back(frame.thisptr) # just copy pointer
                else:
                    for i in idxs:
                        frame.thisptr = self.frame_v[i] # point to i-th item
                        frame.py_free_mem = False # don't free mem
                        _farray.frame_v.push_back(frame.thisptr) # just copy pointer
                return _farray
            else:
                idx_1 = get_positive_idx(idxs, self.size)
                # raise index out of range
                if idxs != 0 and idx_1 == 0:
                    # need to check if array has only 1 element. 
                    # arr[0] is  arr[-1]
                    if idxs != -1:
                        raise ValueError("index is out of range")
                #print ("get memoryview")
                #frame.thisptr = &(self.frame_v[idx_1])
                frame.py_free_mem = False
                frame.thisptr = self.frame_v[idx_1]
                return frame

    def _fast_slice(self, my_slice):
        """only positive indexing

        Examples
        --------
            traj._fast_slice(slice(0, 8, 2))
            traj._fast_slice((0, 8, 2))

        """
        cdef int start, stop, step
        cdef int count
        cdef Trajectory myview = Trajectory(check_top=False)
        cdef _Frame* _frame_ptr

        myview.top = self.top

        if isinstance(my_slice, slice):
            start, stop, step  = my_slice.indices(self.size)
        else:
            try:
                start, stop, step  = my_slice
            except:
                raise ValueError("don't know how to unpack start, stop and step")

        count = start
        while count < stop:
            _frame_ptr = self.frame_v[count]
            myview.frame_v.push_back(_frame_ptr)
            count += step

        return myview

    @cython.boundscheck(False)
    @cython.wraparound(False)
    @cython.initializedcheck(False)
    def __setitem__(self, idx, other):
        # TODO : add slice
        # make a copy
        # to make thing simple, we don't use fancy slicing here
        cdef Frame frame = Frame() # create _Frame pointer
        frame.py_free_mem = False
        cdef AtomMask atm
        cdef double[:, :, :] view3d
        cdef double* ptr
        cdef int[:] int_view
        cdef int i, j, k
        cdef int n_frames = self.n_frames
        cdef Trajectory other_traj

        if len(self) == 0:
            raise ValueError("Your Trajectory is empty, how can I index it?")

        if other is None:
            raise ValueError("why bothering assign None?")
        if is_int(idx):
            if isinstance(other, Frame):
                frame = <Frame> other.copy()
                frame.py_free_mem = False
                if frame.n_atoms != self.n_atoms:
                    raise ValueError("don't have the same n_atoms")
                self.frame_v[idx] = frame.thisptr
            else:
                # xyz
                try:
                    self[<int> idx]._fast_copy_from_xyz(other)
                except:
                    msg = "`other` must be a Frame or an array xzy with shape=(natoms, 3), dtype=float64"
                    raise ValueError(msg)
        elif idx == '*':
            # update all atoms, use fast version
            self.update_xyz(other) # xyz
        elif isinstance(idx, AtomMask) or isinstance(idx, string_types):
            if isinstance(idx, AtomMask):
                atm = <AtomMask> idx
            else:
                atm = self.top(idx)
            if isinstance(other, Trajectory):
                # TODO: not use numpy?
                other_traj = other
                for i in range(n_frames):
                    self[i][atm] = other_traj[i].xyz
            else:
                view3d = other
                try:
                    int_view = atm.indices.astype('i4')
                except ValueError:
                    int_view = atm.indices
                # loop all frames
                for i in range(view3d.shape[0]):
                    # don't use pointer: frame.thisptr = self.frame_v[i]
                    # (got segfault)
                    #frame = self[i]
                    frame.thisptr = self.frame_v[i]
                    # loop all selected atoms
                    for j in range(view3d.shape[1]):
                        # take atom index
                        k = int_view[j]
                        # update coords for each atoms
                        # take pointer position
                        ptr = frame.thisptr.xAddress() + 3 * k
                        # assignment
                        ptr[0] = view3d[i, j, 0]
                        ptr[1] = view3d[i, j, 1]
                        ptr[2] = view3d[i, j, 2]
        else:
            # example: self[0, 0, 0] = 100.
            self[idx[0]][idx[1:]] = other
            #txt = "not yet implemented. Try using framearray[idx1][idx2, idx3] = value"
            #raise NotImplementedError(txt)
        
    def __delitem__(self, int idx):
        self.erase(idx)

    def __str__(self):
        return my_str_method(self)

    def __repr__(self):
        return self.__str__()
    
    def __enter__(self):
        return self

    def __exit__(self, *args):
        # we don't do anythin here. Just create the same API for TrajectoryIterator
        pass

    def frame_iter(self, start=0, stop=None, stride=1, mask=None, autoimage=False, rmsfit=None):
        # TODO: combined with TrajectoryIterator
        from pytraj.core.frameiter import FrameIter

        if mask is None:
            _top = self.top
        else:
            _top = self.top._get_new_from_mask(mask)

        if rmsfit is not None:
            if len(rmsfit) != 2:
                raise ValueError("rmsfit must be a tuple of two elements (frame, mask)")
            if is_int(rmsfit[0]):
                index = rmsfit[0]
                rmsfit = tuple([self[index], rmsfit[1]])

        # check how many frames will be calculated
        if stop is None or stop >= self.n_frames:
            stop = self.n_frames
        elif stop < 0:
            stop = get_positive_idx(stop, self.n_frames)

        # make sure `range` return iterator
        n_frames = len(range(start, stop, stride))

        frame_iter_super = self._frame_iter(start, stop, stride)

        return FrameIter(frame_iter_super,
                         original_top=self.top,
                         new_top=_top,
                         start=start,
                         stop=stop,
                         stride=stride,
                         mask=mask,
                         autoimage=autoimage,
                         rmsfit=rmsfit,
                         n_frames=n_frames)

    def _frame_iter(self, int start=0, int stop=-1, int stride=1, mask=None):
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
        cdef Frame frame = Frame()
        frame.py_free_mem = False # dont let Python free mem
        cdef AtomMask atm
        cdef int _end
        cdef int[:] int_view

        if stop == -1 or stop >= self.n_frames:
            _end = <int> self.n_frames
        else:
            _end = stop

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
            #frame = Frame(n_atoms)
            # don't need to allocate frame here
            pass

        i = start
        while i < _end:
            frame.thisptr = self.frame_v[i]
            if mask is not None:
                frame2.thisptr.SetCoordinates(frame.thisptr[0], atm.thisptr[0])
                yield frame2
            else:
                yield frame
            i += stride

    def reverse(self):
        # should we just create a fake operator?
        cpp_reverse(self.frame_v.begin(), self.frame_v.end())

    def swap(self, arr0, arr1):
        """swap one or more pairs of frames"""
        cdef int i, j
        cdef int[:] i_view
        cdef int[:] j_view
        cdef int k

        if is_int(arr0) and is_int(arr1):
            i = <int> arr0
            j = <int> arr1
            iter_swap(self.frame_v.begin() + i, self.frame_v.begin() + j)
        else:
            if hasattr(arr0, 'itemsize') and arr0.itemsize != 4:
                raise ValueError("must be int32")
            elif hasattr(arr1, 'itemsize') and arr1.itemsize != 4:
                raise ValueError("must be int32")
            elif isinstance(arr0, (list, tuple)) and isinstance(arr1, (list, tuple)):
                # convert to memview
                i_view = _int_array1d_like_to_memview(arr0)
                j_view = _int_array1d_like_to_memview(arr1)
                self._swap_from_array(i_view, j_view)
            else:
                try:
                    # use memview
                    i_view = arr0
                    j_view = arr1
                    self._swap_from_array(i_view, j_view)
                except:
                    raise NotImplementedError()

    def _swap_from_array(self, cython.integral[:] i_view, cython.integral[:] j_view):
        cdef int i, j
        cdef int k

        assert i_view.shape[0] == j_view.shape[0]
        for k in range(i_view.shape[0]):
            i = i_view[k]
            j = j_view[k]
            iter_swap(self.frame_v.begin() + i, self.frame_v.begin() + j)

    def erase(self, idxs):
        cdef int idx
        # dealloc frame pointer too?
        if is_int(idxs):
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

    @property
    def n_atoms(self):
        return self.top.n_atoms

    def __len__(self):
        return self.size

    def __iter__(self):
        """return a reference of Frame instance
        >>> for frame in Trajectory_instance:
        >>>     pass
                
        """
        cdef vector[_Frame*].iterator it  = self.frame_v.begin()
        cdef Frame frame 

        # use memoryview, don't let python free memory of this instance
        while it != self.frame_v.end():
            frame = Frame()
            frame.py_free_mem = False
            frame.thisptr = deref(it)
            yield frame
            incr(it)

    def insert(self, Frame frame, int idx, copy=True):
        """insert a Frame at position idx with/without copying

        Parameters
        ----------
        frame : Frame object
        idx : int, position
        copy : bool, default=True
            make a copy of Frame or not

        Examples
        --------
            traj.insert(frame, 3, copy=False)
        """
        cdef Frame tmp_frame
        cdef vector[_Frame*].iterator it = self.frame_v.begin()

        if copy:
            tmp_frame = frame.copy()
            tmp_frame.py_free_mem = False # keep lifetime longer
        else:
            tmp_frame = frame

        self.frame_v.insert(it + idx, tmp_frame.thisptr)
        
    def append(self, Frame framein, copy=True):
        """append new Frame

        Parameters
        ---------
        framein : Frame object
        copy : bool, default=True
            if 'True', make a copy of Frame. If 'False', create a view
        """
        cdef Frame frame
        # Note: always use `copy=True`
        # use `copy = False` if you want to get memoryview for `self`
        # need to set `py_free_mem = False`
        if copy:
            frame = Frame(framein)
            frame.py_free_mem = False
            self.frame_v.push_back(frame.thisptr)
        else:
            framein.py_free_mem = False
            self.frame_v.push_back(framein.thisptr)

    def join(self, traj, copy=True):
        """traj.join(traj2) with/without copy
        """
        cdef Trajectory other, farray
        cdef Frame frame

        if traj is self:
            raise ValueError("why do you join your self?")
        if is_pytraj_trajectory(traj) or is_frame_iter(traj):
            if self.top.n_atoms != traj.top.n_atoms:
                raise ValueError("n_atoms of two arrays do not match")
            for frame in traj:
                self.append(frame, copy=copy)
        elif isinstance(traj, (list, tuple)):
            # assume a list or tuple of Trajectory
            for farray in traj:
                self.join(farray, copy=copy)
        else:
            raise ValueError("traj must a Trajectory-like for frame iter")

    def resize(self, int n_frames):
        self.frame_v.resize(n_frames)

    @property
    def temperatures(self):
        """return a Python array of temperatures
        """
        cdef pyarray tarr = pyarray('d', [])

        for frame in self:
            tarr.append(frame.temperature)
        return tarr

    @property
    def temperature_set(self):
        return _get_temperature_set(self)

    def get_frames(self, from_traj=None, indices=None, update_top=False, copy=False):
        """get frames from Trajin instance
        def get_frames(from_traj=None, indices=None, update_top=False, copy=True)
        Parameters:
        ----------
        from_traj : TrajectoryIterator or Trajectory, default=None
            if `from_traj` is None, return a new Trajectory (view or copy)
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

        msg = """Trajectory.top.n_atoms should be equal to Trajin_Single.top.n_atoms 
               or set update_top=True"""

        if from_traj is not None:
            ts = from_traj
            # append new frames to `self`
            if update_top:
                self.top = ts.top.copy()

            if not update_top:
                if self.top.n_atoms != ts.top.n_atoms:
                    raise ValueError(msg)

            if isinstance(ts, Trajin_Single) or isinstance(ts, TrajectoryIterator):
                # alway make a copy
                if indices is not None:
                    # slow method
                    # TODO : use `for idx in leng(indices)`?
                    if isinstance(indices, slice):
                        # use slice for saving memory
                        start, stop, step = indices.start, indices.stop, indices.step
                        for i in range(start, stop, step):
                            self.append(ts[i], copy=True)
                    else:
                        # regular list, tuple, array,...
                        for i in indices:
                            #print "debug Trajectory.get_frames"
                            self.append(ts[i], copy=True)
                else:    
                    # get whole traj
                    frame = Frame()
                    #frame.set_frame_v(ts.top, ts.has_vel(), ts.n_repdims)
                    frame.set_frame_v(ts.top)
                    with ts:
                        for i in range(ts.n_frames):
                            self.append(ts[i])

            elif isinstance(ts, Trajectory):
                # can return a copy or no-copy based on `copy` value
                # use try and except?
                if indices is None:
                    for i in range(ts.size):
                        # TODO : make indices as an array?
                        # create `view`
                        self.append(ts[i], copy=copy)
                else:
                    for i in indices:
                        # TODO : make indices as an array?
                        self.append(ts[i], copy=copy)

        else:
            # if from_traj is None, return new Trajectory
            newfarray = Trajectory()
            if update_top:
                newfarray.top = self.top.copy()
            for i in indices:
                newfarray.append(self[i], copy=copy)
            return newfarray

    def strip_atoms(self, mask=None, update_top=True, bint has_box=False):
        """
        Notes
        -----
        if you use memory for numpy, you need to update after resizing Frame
            arr0 = np.asarray(frame.buffer)
            frame.strip_atoms(top,"!@CA")
            # update view
            arr0 = np.asarray(frame.buffer)
        """
        cdef AtomMask atm = self.top(mask)
        if atm.n_atoms == 0:
            raise ValueError("number of stripped atoms must be > 1")
        # read note about `_strip_atoms`
        atm.invert_mask()

        cdef vector[_Frame*].iterator it
        cdef Frame frame = Frame()
        cdef Topology tmptop = Topology()

        if mask == None: 
            raise ValueError("Must provide mask to strip")
        #mask = mask.encode("UTF-8")

        # do not dealloc since we use memoryview for _Frame
        frame.py_free_mem = False
        it = self.frame_v.begin()
        while it != self.frame_v.end():
            frame.thisptr = deref(it)
            # we need to update topology since _strip_atoms will modify it
            tmptop = self.top.copy()
            frame._strip_atoms(tmptop, atm, update_top, has_box)
            incr(it)
        if update_top:
            self.top = tmptop.copy()

    def _fast_strip_atoms(self, mask=None, bint update_top=True, bint has_box=False):
        """
        Paramters
        ---------
        mask : str
        update_top : bool, default=True
            'True' : automatically update Topology
        has_box : bool, default=False (does not work with `True` yet)
        Notes
        -----
        * Known bug: 
        * if you use memory for numpy, you need to update after resizing Frame
        >>> arr0 = np.asarray(frame.buffer)
        >>> frame.strip_atoms(top,"!@CA")
        >>> # update view
        >>> arr0 = np.asarray(frame.buffer)
        """

        cdef Frame frame = Frame()
        cdef Topology newtop = Topology() # make a pointer
        cdef _Frame _tmpframe
        cdef int i 
        cdef int n_frames = self.frame_v.size()
        cdef AtomMask atm

        if mask is None: 
            raise ValueError("Must provide mask to strip")

        atm = self.top(mask)
        if atm.n_atoms == 0:
            raise ValueError("number of stripped atoms must be > 1")
        # use `invert_mask` to construct new frames
        atm.invert_mask() # read note in `Frame._strip_atoms`
        # don't try using bare pointer for newtop.thisptr
        # kept getting segfault
        newtop = self.top.copy()
        newtop.strip_atoms(mask)
        # do not dealloc since we use memoryview for _Frame
        frame.py_free_mem = False

        for i in range(n_frames):
            #newtop.thisptr = self.top.thisptr.modifyStateByMask(atm.thisptr[0])
            # point to i-th _Frame
            frame.thisptr = self.frame_v[i]
            # make a copy and modify it. (Why?)
            # why do we need this?
            # update new _Topology
            # need to copy all other informations
            # allocate
            _tmpframe.SetupFrameV(newtop.thisptr.Atoms(), newtop.thisptr.ParmCoordInfo())
            _tmpframe.SetFrame(frame.thisptr[0], atm.thisptr[0])
            # make a copy: coords, vel, mass...
            # if only care about `coords`, use `_fast_copy_from_frame`
            frame.thisptr[0] = _tmpframe
        if update_top:
            # C++ assignment
            self.top.thisptr[0] = newtop.thisptr[0]
            print (self.top)

    def save(self, filename="", fmt='unknown', overwrite=True, *args, **kwd):
        _savetraj(self, filename, fmt, overwrite, *args, **kwd)

    def write(self, *args, **kwd):
        """same as `save` method"""
        self.save(*args, **kwd)

    def split_and_write_traj(self, *args, **kwd):
        _split_and_write_traj(self, *args, **kwd)

    def split(self, n_chunks):
        chunksize = self.n_frames // n_chunks
        return (fa for fa in self.chunk_iter(chunksize=chunksize))

    def set_frame_mass(self):
        """update mass for each Frame from self.top"""
        cdef Frame frame
        for frame in self:
            frame.set_frame_mass(self.top)

    def rmsfit(self, ref=None, mask="*", mode='pytraj'):
        """do the fitting to reference Frame by rotation and translation
        Parameters
        ----------
        ref : {Frame object, int, str}, default=None 
            Reference
        mask : str or AtomMask object, default='*' (fit all atoms)
        mode : 'cpptraj' (faster but can not use AtomMask)| 'pytraj'

        Examples
        --------
            traj.rmsfit(0) # fit to 1st frame
            traj.rmsfit('last', '@CA') # fit to last frame using @CA atoms
        """
        # not yet dealed with `mass` and box
        from pytraj.actions.CpptrajActions import Action_Rmsd
        cdef Frame frame
        cdef AtomMask atm
        cdef Frame ref_frame
        cdef int i


        if isinstance(ref, Frame):
            ref_frame = <Frame> ref
        elif is_int(ref):
            i = <int> ref
            ref_frame = self[i]
        elif isinstance(ref, string_types):
            if ref.lower() == 'first':
                i = 0
            if ref.lower() == 'last':
                i = -1
            ref_frame = self[i]
        else:
            # first
            ref_frame = self[0]

        if mode == 'pytraj':
            if isinstance(mask, string_types):
                atm = self.top(mask)
            elif isinstance(mask, AtomMask):
                atm = <AtomMask> mask
            else:
                raise ValueError("mask must be string or AtomMask object")

            for frame in self:
                _, mat, v1, v2 = frame.rmsd(ref_frame, atm, get_mvv=True)
                frame.trans_rot_trans(v1, mat, v2)

        elif mode == 'cpptraj':
            # switch to fast speed
            # we still use mode 'pytraj' so we can use AtomMask (for what?)
            act = Action_Rmsd()
            act(mask, [ref_frame, self], top=self.top)
        else:
            raise ValueError("mode = pytraj | cpptraj")

    # start copy and paste from "__action_in_traj.py"
    def calc_distance(self, mask="", *args, **kwd):
        return pyca.calc_distance(self, mask, *args, **kwd)

    def calc_distrmsd(self, mask="", *args, **kwd):
        return pyca.calc_distrmsd(self, mask, *args, **kwd)

    def calc_rmsd(self, ref=None, mask='', mass=False, fit=True, *args, **kwd):
        """
        Examples
        --------
        >>> from pytraj import io
        >>> traj = io.load_sample_data('tz2')[:]
        >>> traj.calc_rmsd() 
        >>> traj.calc_rmsd(0)
        >>> traj.calc_rmsd(-1)
        >>> traj.calc_rmsd(-1, '@CA', True, True, dtype='dataset')
        >>> traj.calc_rmsd(-1, '@CA', True, True, dtype='pyarray')
        """
        if is_int(ref):
            # index
            _ref = self[ref]
        else:
            _ref = ref
        return pyca.calc_rmsd(mask=mask, traj=self, ref=_ref,
                              mass=mass, fit=fit, *args, **kwd)

    def rmsd(self, *args, **kwd):
            return self.calc_rmsd(*args, **kwd)

    def calc_bfactors(self, mask="", *args, **kwd):
        return pyca.calc_bfactors(self, mask, *args, **kwd)

    def calc_radgyr(self, mask="", *args, **kwd):
        return pyca.calc_radgyr(self, mask, *args, **kwd)

    def calc_angle(self, mask="", *args, **kwd):
        return pyca.calc_angle(self, mask, *args, **kwd)

    def calc_matrix(self, mask="", *args, **kwd):
        return pyca.calc_matrix(self, mask, *args, **kwd)

    def calc_dssp(self, mask="", *args, **kwd):
        return pyca.calc_dssp(self, mask, *args, **kwd)

    def calc_dihedral(self, mask="", *args, **kwd):
        return pyca.calc_dihedral(self, mask, *args, **kwd)

    def calc_multidihedral(self, mask="", *args, **kwd):
        return pyca.calc_multidihedral(self, mask, *args, **kwd)

    def calc_molsurf(self, mask="", *args, **kwd):
        return pyca.calc_molsurf(self, mask, *args, **kwd)

    def calc_center_of_mass(self, mask="", *args, **kwd):
        return pyca.calc_center_of_mass(self, mask, *args, **kwd)

    def calc_COM(self, mask="", *args, **kwd):
        return pyca.calc_center_of_mass(self, mask, *args, **kwd)

    def calc_center_of_geometry(self, mask="", *args, **kwd):
        return pyca.calc_center_of_geometry(self, mask, *args, **kwd)

    def calc_COG(self, mask="", *args, **kwd):
        return pyca.calc_center_of_geometry(self, mask, *args, **kwd)

    def calc_vector(self, mask="", dtype='dataset', *args, **kwd):
        return pyca.calc_vector(self, mask, dtype=dtype, *args, **kwd)

    def search_hbonds(self, mask="*", *args, **kwd):
        return pyca.search_hbonds(self, mask, *args, **kwd)

    def get_average_frame(self, mask="", *args, **kwd):
        return pyca.get_average_frame(self, mask, *args, **kwd)

    def calc_watershell(self, mask="", *args, **kwd):
        return pyca.calc_watershell(self, mask, *args, **kwd)

    def autoimage(self, command=""):
        # NOTE: I tried to used cpptraj's Action_AutoImage directly but
        # there is no gain in speed. don't try.
        pyca.do_autoimage(self, command)

    def rotate(self, mask="", matrix=None):
        cdef Frame frame
        cdef Matrix_3x3 mat
        cdef AtomMask atm

        if matrix is None:
            pyca.do_rotation(self, mask)
        else:
            try:
                mat = Matrix_3x3(mask) # cheap to copy
                atm = self.top(mask)
                for frame in self:
                    frame.rotate_with_matrix(mat, atm)
            except:
                raise ValueError("require string or Matrix-like object")

    def translate(self, mask=""):
        pyca.do_translation(self, mask)

    def center(self, mask=""):
        """
        Examples
        --------
            traj.center() # all atoms, center to box center (x/2, y/2, z/2)
            traj.center('@CA origin') # center at origin, use @CA
            traj.center('mass') # center to box center, use mass weighted.
            traj.center(':1 mass') # residue 1, use mass weighted.

        Notes
        -----
            if using 'mass', make sure to `set_frame_mass` before calling `center`

        See Also
        --------
            Amber15 manual (http://ambermd.org/doc12/Amber15.pdf, page 546)

        """
        from pytraj.actions.CpptrajActions import Action_Center
        act = Action_Center()
        act(mask, self)

    def set_nobox(self):
        cdef Frame frame

        for frame in self:
            frame.set_nobox()

    def _allocate(self, int n_frames, int n_atoms):
        """pre-allocate (n_atoms, n_atoms, 3)
        """
        cdef Frame frame

        self.frame_v.resize(n_frames)
        for i in range(n_frames):
            self.frame_v[i] = new _Frame(n_atoms)

    def box_to_ndarray(self):
        return _box_to_ndarray(self)

    # math
    def __tmpidiv__(self, value):
        cdef Frame frame
        cdef int i
        cdef int size = self.size

        if isinstance(value, Trajectory) or is_frame_iter(value):
            for i, frame in enumerate(value):
                # frame /= other_frame
                self[i] /= frame
        else:
            # numpy
            for frame in self:
                try:
                    frame.xyz.__idiv__(value)
                except:
                    frame.xyz.__itruediv__(value)

    def __idiv__(self, value):
        self.__tmpidiv__(value)
        return self

    def __itruediv__(self, value):
        self.__tmpidiv__(value)
        return self

    def __iadd__(self, value):
        cdef Frame frame
        cdef Trajectory tmp_traj
        cdef int i
        cdef int size = self.size
        cdef int n_atoms = self.n_atoms
        cdef Frame other_frame

        if isinstance(value, Frame):
            other_frame = value
            for i in range(size):
                # frame += other_frame
                self.frame_v[i][0] += other_frame.thisptr[0]
        elif is_pytraj_trajectory(value) or is_frame_iter(value):
            # nogain with OPENMP
            for i,frame in enumerate(value):
                # frame += other_frame
                self.frame_v[i][0] += frame.thisptr[0]
        else:
            # numpy
            for frame in self:
                frame.xyz.__iadd__(value)
        return self

    def __isub__(self, value):
        cdef Frame frame
        cdef Trajectory tmp_traj
        cdef int i
        cdef int size = self.size
        cdef Frame other_frame

        if is_pytraj_trajectory(value) or is_frame_iter(value):
            # nogain with OPENMP
            for i, frame in enumerate(value):
                # frame -= other_frame
                self.frame_v[i][0] -= frame.thisptr[0]
        elif isinstance(value, Frame):
            other_frame = value
            for i in range(size):
                # frame -= other_frame
                self.frame_v[i][0] -= other_frame.thisptr[0]
        else:
            # numpy
            for frame in self:
                frame.xyz.__isub__(value)
        return self

    def __imul__(self, value):
        cdef Frame frame
        cdef Trajectory tmp_traj
        cdef int i
        cdef int size = self.size
        cdef Frame other_frame

        if isinstance(value, Frame):
            other_frame = value
            for i in range(size):
                # frame *= other_frame
                self.frame_v[i][0] *= other_frame.thisptr[0]
        elif is_pytraj_trajectory(value) or is_frame_iter(value):
            # nogain with OPENMP
            for i, frame in enumerate(value):
                # frame *= other_frame
                self.frame_v[i][0] *= frame.thisptr[0]
        else:
            # numpy
            for frame in self:
                frame.xyz.__imul__(value)
        return self

    def apply(self, func=None, args=None, indices_or_mask=None):
        """apply `func` to traj's coords"""
        cdef Frame frame

        if isinstance(indices_or_mask, string_types):
            indices = self.top(indices_or_mask).indices
        else:
            indices = indices_or_mask

        # TODO: make this shorter
        if args is None:
            for frame in self:
                if indices is not None:
                    frame[indices] = func(frame.xyz[indices])
                else:
                    frame.xyz[:] = func(frame.xyz)
        else:
            for frame in self:
                if indices is not None:
                    frame[indices] = func(frame.xyz[indices], args)
                else:
                    frame.xyz[:] = func(frame.xyz, args)

    def __contains__(self, Frame other):
        """check if frame is in self"""
        cdef Frame frame

        for frame in self:
            if other.is_(frame):
                return True
        return False

    def average(self, mask=None):
        cdef AtomMask atm
        cdef Frame frame
        cdef _Frame* _frame
        cdef _Frame* frame_ptr
        cdef bint has_mask
        cdef int n_atoms

        if mask is not None or (isinstance(mask, string_types) and not mask.is_empty()):
            atm = self.top(mask)
            frame = Frame(atm.n_atoms)
            has_mask = True
        else:
            frame = Frame(self.n_atoms)
            has_mask = False
        frame.zero_coords()

        for frame_ptr in self.frame_v:
            if has_mask: 
                _frame = new _Frame(frame_ptr[0], atm.thisptr[0])
            else:
                _frame = new _Frame(frame_ptr[0])
            frame.thisptr[0] += deref(_frame)
        frame /= self.n_frames
        return frame

    def chunk_iter(self, int chunksize=2, int start=0, int stop=-1):
        """iterately get Frames with start, chunk
        returning Trajectory

        Parameters
        ---------
        start : int (default = 0)
        chunk : int (default = 1, return Frame instance). 
                if `chunk` > 1 : return Trajectory instance
        copy_top : bool, default=False
            if False: no Topology copy is done for new (chunk) Trajectory
        """
        # use `chunk` so we don't need to change `chunk` to `chunksize`
        cdef int chunk = chunksize
        cdef int n_chunk, i, _stop
        cdef int n_frames = self.n_frames
        cdef int n_atoms = self.n_atoms
        cdef Trajectory farray
        cdef int size = self.n_frames
        cdef vector[_Frame*].iterator it = self.frame_v.begin()

        # check `start`
        if start < 0 or start >= n_frames:
            start = 0

        # check `stop`
        if stop <= 0 or stop >= n_frames:
            stop = size - 1

        if chunk <= 1:
            raise ValueError("chunk must be >= 2")

        if chunk + start > stop:
            raise ValueError("start + chunk must be smaller than max frames")

        n_chunk = int((stop- start)/chunk)
        if ((stop - start) % chunk ) != 0:
            n_chunk += 1

        for i in range(n_chunk):
            # always create new Trajectory
            farray = Trajectory(check_top=False)
            farray.top = self.top

            if i != n_chunk - 1:
                _stop = start + chunk*(i+1)
            else:
                _stop = stop + 1
            farray.frame_v.resize(_stop - i * chunk)
            farray.frame_v.assign(it + i * chunk, it + _stop)
            yield farray

    @classmethod
    def from_iterable(cls, iteratable, top=None, copy=True):
        """return a new Trajectory from `iteratable` object
        """
        if top is None:
            try:
                top = iteratable.top
            except AttributeError:
                raise AttributeError("require a Topology")

        traj = cls(top=top)
        for f in iteratable:
            traj.append(f, copy=copy)
        return traj

    @property
    def _estimated_MB(self):
        """esimated MB of data will be loaded to memory
        """
        return self.n_frames * self.n_atoms * 3 * 8 / 1E6
