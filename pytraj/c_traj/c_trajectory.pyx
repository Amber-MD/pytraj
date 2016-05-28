# distutils: language = c++
import os
import numpy as np
from ..trajectory import Trajectory
from ..core.c_core cimport AtomMask
from ..core.box cimport Box
from ..topology cimport Topology

from ..cyutils import get_positive_idx
from pytraj.externals.six import string_types
from ..shared_methods import (my_str_method, _xyz, _box)
from ..utils.check_and_assert import ensure_exist
from ..utils.check_and_assert import is_array, is_range

# do not use compat for range here. Let Cython handle
#from ..externals.six.moves import range


def _split_range(int chunksize, int start, int stop):
    '''split a given range to n_chunks

    Examples
    --------
    >>> _split_range(3, 0, 10)
    [(0, 3), (3, 6), (6, 10)]
    '''
    cdef int n_chunks, i, _stop

    n_chunks = (stop - start)//chunksize

    if ((stop - start) % chunksize) != 0:
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
        # we use TopPtr here, self._top is a binding
        self._top._own_memory = False
        self._filelist = []
        self._own_memory = True
        self._being_transformed = False
        self._being_superposed = False
        self._transform_commands = []
        self._max_count_to_reset = 50000 # for reset superpose memory

    def _load(self, filename=None, top=None, frame_slice=(0, -1, 1)):
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
            step = 1
        elif len(frame_slice) == 2:
            # no step info
            start, stop = frame_slice
            step = 1
        elif len(frame_slice) == 3:
            start, stop, step = frame_slice
        else:
            raise ValueError()

        start += 1
        # don't increase stop by +1
        # slice(0, 10, None) --> python does not take last `10`
        arg = " ".join((str(start), str(stop), str(step)))

        if isinstance(filename, string_types):
            # use absolute path so we can go to different folder
            filename = os.path.abspath(filename)
            _arglist = ArgList(arg)
            self.thisptr.AddSingleTrajin(
                filename.encode(),
                _arglist.thisptr[0],
                tmp_top.thisptr)
            self._filelist.append(os.path.abspath(filename))
        else:
            raise ValueError("filename must a a string")

        self._initialize_actionlist()

    def _initialize_actionlist(self):
        assert self.top.n_atoms > 0, 'must set Topology'

        if self._being_superposed:
            # if self._being_superposed, its _cdslist and _actionlist were already maded
            # free memory
            self._cdslist.thisptr.Clear()
            self._cdslist.thisptr.ClearTop()
            self._cdslist.thisptr.ClearRef()
            self._actionlist.thisptr.Clear()
            del self._cdslist.thisptr
            self._cdslist.thisptr = new _CpptrajDatasetList()
            del self._actionlist.thisptr
            self._actionlist.thisptr = new _ActionList()

        self._actionlist = ActionList(top=self.top)
        self._cdslist = CpptrajDatasetList()

    def __len__(self):
        return self.n_frames

    @property
    def n_frames(self):
        """total snapshots
        """
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
        cdef int n_frames = self.n_frames

        # use `frame` as buffer
        cdef Frame frame = Frame()

        # del frame.thisptr # do not do this.
        frame.thisptr[0] = self.thisptr.AllocateFrame()

        for i in range(n_frames):
            # do not create new Frame inside this loop to reduce memory
            self.thisptr.GetFrame(i, frame.thisptr[0])
            if self._being_transformed:
                self._do_transformation(frame)
            yield frame

    property top:
        """Topology
        """
        def __get__(self):
            #self._top.thisptr[0] = self.thisptr.Top()
            self._top.thisptr = self.thisptr.TopPtr()
            return self._top

        def __set__(self, Topology other):
            # self.thisptr.SetTopology(other.thisptr[0])
            self.thisptr.CoordsSetup(other.thisptr[0], self.thisptr.CoordsInfo())

    def iterframe(self, int start=0, int stop=-1, int step=1, mask=None):
        '''iterately get Frames with start, stop, step
        Parameters
        ---------
        start : int (default = 0)
        stop : int (default = max_frames)
        step : int
        mask : str or array of interger
        '''
        cdef int i
        cdef int n_atoms = self.n_atoms
        cdef Frame frame = Frame()
        cdef AtomMask atm
        cdef int _end
        cdef int[:] int_view
        cdef unsigned int max_frame = self.n_frames

        if stop == -1:
            _end = <int> max_frame
        else:
            _end = stop

        if mask is not None:
            #    del frame.thisptr
            if isinstance(mask, string_types):
                atm = self.top(mask)
            else:
                try:
                    atm = AtomMask()
                    atm.add_selected_indices(mask)
                except TypeError:
                    raise TypeError("dont know how to cast to memoryview")
        #    frame.thisptr = new _Frame(<int>atm.n_atoms)
        # else:
        #    frame.thisptr[0] = self.thisptr.AllocateFrame()

        frame.thisptr[0] = self.thisptr.AllocateFrame()

        with self:
            i = start
            while i < _end:
                if mask is None:
                    self.thisptr.GetFrame(i, frame.thisptr[0])
                else:
                    self.thisptr.GetFrame(i, frame.thisptr[0], atm.thisptr[0])
                if self._being_transformed:
                    self._do_transformation(frame)
                yield frame
                i += step

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
                # do not make Topology copy here to save time.
                # we don't use farray.top assignment to avoid extra copy
                farray._top = self.top
                real_n_frames = len(range(_tmp_start, _tmp_stop))
                farray._allocate(real_n_frames, farray.top.n_atoms)
                farray._boxes = np.empty((real_n_frames, 6), dtype='f8')

                for idx, frame in enumerate(self.iterframe(start=_tmp_start,
                                                           stop=_tmp_stop)):
                    if self._being_transformed:
                        self._do_transformation(frame)
                    farray._xyz[idx] = frame.xyz
                    farray._boxes[idx] = frame.box._get_data()
                yield farray

    def __setitem__(self, idx, value):
        raise NotImplementedError(
            "Read only Trajectory. Use Trajectory class for __setitem__")

    def __getitem__(self, idxs):
        # allocate frame for storing data
        cdef Frame frame0
        cdef Frame frame = Frame()
        cdef int start, stop, step
        cdef int i
        cdef int idx_1, idx_2
        cdef int[:] int_view
        cdef list tmplist
        cdef AtomMask atom_mask_obj
        cdef idxs_size

        frame.thisptr[0] = self.thisptr.AllocateFrame()

        if isinstance(idxs, AtomMask):
            # atm = top('@CA')
            # traj[atm]
            atom_mask_obj = <AtomMask> idxs
            _farray = Trajectory()
            _farray.top = self.top._modify_state_by_mask(atom_mask_obj)
            for i, frame in enumerate(self):
                if self._being_transformed:
                    self._do_transformation(frame)
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
            start, stop, step = idxs.indices(self.n_frames)
            self.tmpfarray = self._load_traj_by_indices(range(start, stop, step))
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
                    if self._being_transformed:
                        self._do_transformation(frame)
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
                self.tmpfarray = self._load_traj_by_indices(idxs)
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
                    if self._being_transformed:
                        self._do_transformation(frame)
                self.tmpfarray = frame
                return self.tmpfarray

    @property
    def unitcells(self):
        '''return unitcells (Box) with shape=(n_frames, 6)
        '''
        return _box(self)

    @property
    def xyz(self):
        '''return a copy of xyz coordinates (ndarray, shape=(n_frames, n_atoms, 3))

        Notes
        -----
            - It will be very expensive to call this attribute since pytraj will load all coordinates
        to memory. This attribute is good for small immutable trajectory.

            -
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
        if self.thisptr and self._own_memory:
            del self.thisptr

    def _load_traj_by_indices(self, indices):
        '''indices is iterable that has __len__
        '''
        cdef int i, j
        cdef int n_atoms = self.n_atoms
        cdef Frame frame
        cdef double[:, :, :] xyz
        cdef int n_frames = len(indices)

        traj = Trajectory()
        traj.top = self.top
        if n_frames == 0:
            # return empty traj
            return traj
        traj._allocate(n_frames, n_atoms)
        traj.unitcells = np.zeros((n_frames, 6), dtype='f8')
        traj.top = self.top
        xyz = traj.xyz

        if not self.thisptr.CoordsInfo().HasVel():
            # faster
            frame = Frame(n_atoms, xyz[0], _as_ptr=True)
            for j, i in enumerate(indices):
                # use `frame` as a pointer pointing to `xyz` memory
                # dump coords to xyz array
                frame.thisptr.SetXptr(frame.n_atoms, &xyz[j, 0, 0])
                # copy coordinates of `self[i]` to j-th frame in `traj`
                self.thisptr.GetFrame(i, frame.thisptr[0])
                if self._being_transformed:
                    self._do_transformation(frame)
                traj.unitcells[j] = frame.box._get_data()
            return traj

        else:
            # slower
            frame = Frame()
            frame.thisptr[0] = self.thisptr.AllocateFrame()
            for j, i in enumerate(indices):
                # use `frame` as a pointer pointing to `xyz` memory
                # dump coords to xyz array
                # copy coordinates of `self[i]` to j-th frame in `traj`
                self.thisptr.GetFrame(i, frame.thisptr[0])
                if self._being_transformed:
                    self._do_transformation(frame)
                traj.xyz[j] = frame.xyz
                traj.unitcells[j] = frame.box._get_data()
            return traj

    def _iterframe_indices(self, frame_indices):
        cdef int i
        cdef Frame frame = Frame()
        cdef unsigned int max_frame = self.n_frames

        frame.thisptr[0] = self.thisptr.AllocateFrame()

        for i in frame_indices:
            assert 0 <= i < max_frame, 'frame index must be between 0 and max_frame - 1'
            self.thisptr.GetFrame(i, frame.thisptr[0])
            if self._being_transformed:
                self._do_transformation(frame)
            yield frame

    def translate(self, command):
        return self._add_transformation('translate', command)

    def scale(self, command):
        return self._add_transformation('scale', command)

    def center(self, command=''):
        return self._add_transformation('center', command)

    def rotate(self, command):
        return self._add_transformation('rotate', command)

    def autoimage(self, command=''):
        return self._add_transformation('autoimage', command)

    def principal(self, command):
        return self._add_transformation('principal', command)

    def superpose(self, mask='', Frame ref=None, ref_mask="", mass=False):
        command = dict(mask=mask, ref=ref, ref_mask=ref_mask, mass=mass)
        return self._add_transformation('superpose', command)

    def _align(self, mask='', Frame ref=None, ref_mask="", bint mass=False):
        """register to superpose to reference frame when iterating. 

        Notes
        -----
        This method is different from ``superpose`` in pytraj.Trajectory.
        It does not change the coordinates of TrajectoryCpptraj/TrajectoryIterator itself but 
        the copy of Frame.
        """
        cdef Frame frame = ref

        refset = self._cdslist.add('reference')
        refset.top = self.top if ref.top is None else ref.top
        refset.add_frame(frame)
        refset.name = 'myref' + str(len(self._cdslist))

        mass_ = 'mass' if mass else ''

        command = '{mask} {refmask} ref {refname} {mass}'.format(refname=refset.name,
                                                            mask=mask,
                                                            refmask=ref_mask,
                                                            mass=mass_)
        self._actionlist.add('align', command, dslist=self._cdslist)

        self._being_transformed = True
        self._being_superposed = True

        # trick cpptraj to re-setup
        # why? if not traj_on_disk().autoimage().superpose() won't perform superpose
        self._actionlist.is_setup = False
        return self

    def _add_transformation(self, name, command):
        '''

        Parameters
        ----------
        name : str, cpptraj action name
        command : str, what do you want.
        '''
        if name not in ['superpose']:
            self._actionlist.add(name, command)
            self._being_transformed = True
        else:
            self._align(**command)
        self._transform_commands.append((name, command))
        return self

    def _remove_transformations(self):
        self._actionlist = ActionList()
        self._cdslist = CpptrajDatasetList()
        self._being_transformed = False
        self._being_superposed = False
        self._transform_commands = []

    def _reset_transformation(self):
        old_commands = self._transform_commands[:]
        self._transform_commands = []

        # self._cdslist.thisptr.Clear()
        self._initialize_actionlist()

        for name, command in old_commands:
             self._add_transformation(name, command)

    def _do_transformation(self, Frame frame):
        '''wrapper for self._actionlist.compute

        To avoid memory leak, we need to reset Dataset that holds RMSD value
        '''
        self._actionlist.compute(frame)
 
    @property
    def metadata(self):
        '''return a dict of general information

        Examples
        --------
        >>> traj.metadata
        {'box_type': 'ortho',
         'has_box': True,
         'has_force': False,
         'has_replcica_dims': False,
         'has_temperature': False,
         'has_time': True,
         'has_velocity': False,
         'n_atoms': 5293,
         'n_frames': 10}
        '''
        return self._crdinfo

    property _crdinfo:
        def __get__(self):
            cdef _CoordinateInfo cinfo

            cinfo = self.thisptr.CoordsInfo()
            return {'has_velocity': cinfo.HasVel(),
                    'has_temperature': cinfo.HasTemp(),
                    'has_time': cinfo.HasTime(),
                    'has_force': cinfo.HasForce(),
                    'has_box': cinfo.HasBox(),
                    'has_replcica_dims': cinfo.HasReplicaDims(),
                    'n_frames': self.n_frames,
                    'n_atoms': self.n_atoms,
                    'box_type': self.top.box.type}

        def __set__(self, dict crdinfo):
            '''value is a dict
            '''
            cdef CoordinateInfo cinfo
            cdef Box box

            crdinfo2 = dict((k, v) for k, v in crdinfo.item())

            if 'box' not in crdinfo2:
                crdinfo2['box'] = self.top.box

            cinfo = CoordinateInfo(crdinfo2)

            self.thisptr.CoordsSetup(self.thisptr.Top(), cinfo.thisptr[0])
