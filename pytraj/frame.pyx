# distutils: language = c++

#from __future__ import absolute_import, division
# turn off `division` for automatically casting
# NOTE: need to import numpy locally in some methods to avoid error + segfault (hell
# knows, cython 0.23.1)
#
from __future__ import absolute_import
import math
cimport cython
from libc.math cimport sqrt
from cython cimport view
from cpython cimport array as cparray  # for extend python array
from cpython.array cimport array as pyarray
from cython.parallel import prange
from cython.operator cimport dereference as deref
from libcpp.vector cimport vector
from libc.string cimport memcpy

import numpy as np
from pytraj.utils.check_and_assert import is_int
from pytraj.core.c_core import ArgList
from pytraj.c_traj.c_trajout import TrajectoryWriter
from pytraj.externals.six import string_types

DEF RADDEG = 57.29577951308232

cdef extern from "TorsionRoutines.h" nogil:
    double cpptorsion "Torsion" (const double *, const double *, const double *, const double *)
    double cppangle "CalcAngle" (const double*, const double*, const double*)

cdef extern from "DistRoutines.h" nogil:
    double DIST2_NoImage(double*, double*)

__all__ = ['Frame']

cdef class Frame (object):
    """A snapshot of trajectory, hodling coordinates, unicell (box), vector (optional), force (optional), ...
    This class is for internal use.

    Parameters
    ----------
    n_atoms : int, default=0
        create a new Frame with n_atoms
    frame : a Frame, default=None
        make a copy from `frame`
    atommask : AtomMask, default=None
        make a copy from `frame` with atommask

    Examples
    --------
    >>> import pytraj as pt
    >>> # create an empty Frame
    >>> frame = pt.Frame()
    >>> frame.n_atoms
    0

    >>> # create a Frame with 304 atoms
    >>> frame = pt.Frame(304)
    >>> frame.n_atoms
    304

    >>> # create an empty Frame then append coordinate
    >>> traj = pt.datafiles.load_tz2()
    >>> f0 = traj[0]
    >>> frame = pt.Frame()
    >>> frame.append_xyz(f0.xyz)
    <Frame with 223 atoms>
    >>> frame.n_atoms == f0.n_atoms
    True

    >>> # copy from other Frame
    >>> frame2 = pt.Frame(frame)

    >>> # create a Frame as a pointer (does not own any memory), used for fast iterating
    >>> frame = pt.Frame(f0.n_atoms, f0.xyz, _as_ptr=True)
    >>> frame.n_atoms == f0.n_atoms
    True
    """

    def __cinit__(self, *args, _as_ptr=False):
        """Constructor for Frame instance

        Examples
        --------
        >>> from pytraj import Frame
        >>> # created empty Frame instance
        >>> frame0 = Frame()
        >>> # create Frame instance with 304 atoms
        >>> frame1 = Frame(304)
        >>> # create a copy of frame1
        >>> frame2 = Frame(frame1)
        >>> # create a copy of frame with atommask
        >>> frame3 = Frame(frame1, atm_instance)
        """
        import numpy as np
        cdef Frame frame
        cdef AtomMask atmask
        cdef vector[_Atom] vt
        cdef Topology top
        cdef Atom at
        cdef list atlist
        cdef int natom, natom3
        cdef int i
        cdef double[:, ::1] view

        self._own_memory = True
        self._as_view = False

        if not args:
            self.thisptr = new _Frame()
        else:
            if not _as_ptr:
                if len(args) == 2:
                    frame, atmask = args
                    self.thisptr = new _Frame(frame.thisptr[0], atmask.thisptr[0])
                elif len(args) == 1:
                    # copy Frame
                    if isinstance(args[0], Frame):
                        frame = args[0]
                        self.thisptr = new _Frame(frame.thisptr[0])
                    # creat a new Frame instance with natom
                    elif isinstance(args[0], int):
                        natom = <int> args[0]
                        self.thisptr = new _Frame(natom)
                    else:
                        # Create Frame from list of atom mask
                        atlist = args[0]
                        for at in atlist:
                            vt.push_back(at.thisptr[0])
                        self.thisptr = new _Frame(vt)
                else:
                    raise ValueError()
            else:
                # create Frame as a view.
                natom = args[0]
                view = args[1]
                self._as_view = True
                self.thisptr = new _Frame(natom, &view[0, 0])

    def __dealloc__(self):
        if self._own_memory and self.thisptr:
            del self.thisptr

    def __del__(self):
        del self.thisptr

    def copy(self):
        """return a deep copy of Frame"""
        cdef Frame frame = Frame(self)
        return frame

    def append_xyz(self, double[:, :] xyz):
        """append 3D array and return itself
        """
        cdef int i
        cdef int N = xyz.shape[0]

        for i in range(N):
            self.thisptr.AddXYZ( &xyz[i, 0])
        return self

    cdef void _append_xyz_2d(self, double[:, :] xyz):
        # for internal use
        # TODO: add assert
        cdef int i
        cdef int N = xyz.shape[0]

        for i in range(N):
            self.thisptr.AddXYZ( &xyz[i, 0])

    cdef void _append_xyz_1d(self, double[:] xyz):
        # TODO: add assert
        # for internal use
        cdef int i
        cdef int N = <int> xyz.shape[0] / 3

        for i in range(N):
            self.thisptr.AddXYZ( &xyz[i*3])

    def swap_atoms(self, cython.integral[:, :] int_view):
        """swap coordinates for an array of atom pairs

        Parameters
        ----------
        int_view: 2D-int array-like, shape=(2, n_atoms)
        """
        cdef int i

        for i in range(int_view.shape[1]):
            self.thisptr.SwapAtoms(int_view[0, i], int_view[1, i])

    def __str__(self):
        tmp = "<%s with %s atoms>" % (
            self.__class__.__name__,
            self.n_atoms,
            )
        return tmp

    def __repr__(self):
        return self.__str__()

    def is_(self, Frame other):
        '''check if Frame is itself or not
        '''
        return self.thisptr == other.thisptr

    property shape:
        def __get__(self):
            '''return shape of its coordinates
            '''
            return self.xyz.shape

    property n_frames:
        def __get__(self):
            '''always return 1. Used for tricking action.
            '''
            return 1

    def __getitem__(self, idx):
        """slicing

        Examples
        --------
        >>> from pytraj import io
        >>> traj = io.load_sample_data('tz2')
        >>> f0 = traj[0]
        >>> data = f0[0]
        >>> data = f0[0, 0]
        >>> data = f0[:,  0]
        >>> data = f0.top = traj.top
        >>> data = f0['@CA']
        """
        cdef AtomMask atm
        cdef cython.view.array cy_arr
        cdef int new_size
        cdef int i, j
        cdef int[:] int_view

        if isinstance(idx, pyarray):
            return self.xyz[idx]
        elif isinstance(idx, AtomMask):
            # return a sub-array copy with indices got from
            # idx.selected_indices()
            # TODO : add doc
            if idx.n_atoms == 0:
                raise ValueError("emtpy mask")
            return self.xyz[idx.indices]
        elif isinstance(idx, tuple) and isinstance(idx[0], AtomMask):
            # (AtomMask, )
            if len(idx) == 1:
                return self[idx[0]]
            elif len(idx) == 2:
                return self[idx[0]][idx[1]]
            elif len(idx) == 3:
                return self[idx[0]][(idx[1], idx[2])]
            else:
                raise NotImplementedError()
        elif isinstance(idx, dict):
            # Example: frame[dict(top=top, mask='@CA')]
            # return a sub-array copy with indices got from
            # idx as a `dict` instance
            atm = AtomMask(idx['mask'])
            idx['top'].set_integer_mask(atm)
            return self[atm.indices]
        elif isinstance(idx, string_types):
            # Example: frame['@CA']
            if self.top is not None and not self.top.is_empty():
                return self[<AtomMask> self.top(idx)]
            else:
                raise ValueError('must have non-empty topology. Use self.set_top'
                                 ' or use self[AtomMask]')

        elif isinstance(idx, tuple) and isinstance(idx[0], string_types):
            # (AtomMask, )
            if len(idx) == 1:
                return self[idx[0]]
            elif len(idx) == 2:
                return self[idx[0]][idx[1]]
            elif len(idx) == 3:
                return self[idx[0]][(idx[1], idx[2])]
            else:
                raise NotImplementedError()
        else:
            return self.xyz[idx]

    def __setitem__(self, idx, value):
        if isinstance(idx, AtomMask):
            self.xyz[idx.indices] = value
        elif isinstance(value, string_types):
            # assume this is atom mask
            if self.top is None:
                raise ValueError("must set Topology for frame")
            else:
                self[self.top(idx)] = value
        else:
            self.xyz[idx] = value

    def __iter__(self):
        cdef int i
        for i in range(self.n_atoms):
            yield np.asarray(self._buffer2d[i])

    def __array__(self):
        """
        arr0 = np.asarray(frame)

        (== (arr0 = np.asarray(frame[:])))
        """
        import numpy as np
        cdef double* ptr = self.thisptr.xAddress()
        return np.asarray(<double[:self.thisptr.Natom(), :3]> ptr, dtype='f8')

    def _fast_copy_from_frame(self, Frame other):
        """only copy coords"""
        # no boundchecking
        # a bit faster than: self.thisptr[0] = other.thisptr[0]
        # (since we copy only coords)
        cdef double *ptr_src
        cdef double *ptr_dest
        cdef int count

        ptr_src = other.thisptr.xAddress()
        ptr_dest = self.thisptr.xAddress()
        count = self.thisptr.Natom() * 3 * sizeof(double)
        memcpy(<void*> ptr_dest, <void*> ptr_src, count)

    def _fast_copy_from_xyz(self, double[:, :] xyz, indices=None):
        """only copy coords

        Parameters
        ----------
        xyz : 2D array-like, dtype='double', must have buffer interface
        indices : 1D array-like, dtype='i4', must have buffer interface
            default=None
        """

        cdef double *ptr_src
        cdef double *ptr_dest
        cdef int count
        cdef int i, j
        cdef int[:] int_view
        # no boundchecking

        if indices is None:
            # copy all
            ptr_src = &xyz[0, 0]
            ptr_dest = self.thisptr.xAddress()
            count = self.thisptr.Natom() * 3 * sizeof(double)
            memcpy(<void*> ptr_dest, <void*> ptr_src, count)
        else:
            count = 3 * sizeof(double)
            try:
                int_view = indices  # create `view`
            except:
                int_view = indices.astype('i4')
            # NOTE: try `prange` with different `schedule` but no gain
            for i in range(int_view.shape[0]):
                j = int_view[i]
                ptr_dest = self.thisptr.xAddress() + j * 3
                ptr_src = &xyz[i, 0]
                # copy coords of each Atom
                memcpy(<void*> ptr_dest, <void*> ptr_src, count)

    def __len__(self):
        return self.size

    property _buffer1d:
        def __get__(self):
            """return memory view for Frame coordinates
            TODO : rename?
            """
            # debug
            # print "from calling _buffer1d: _own_memory = ", self._own_memory
            # end debug
            def _buffer(N):
                cdef double* ptr = self.thisptr.xAddress()
                cdef view.array my_arr
                my_arr = <double[:N]> ptr
                return my_arr
            return _buffer(self.size)

    property _buffer2d:
        def __get__(self):
            """return memory view for Frame coordinates but reshape
            (just like self._buffer3 = self.buffer.reshape())
            TODO : rename?
            """
            # debug
            # print "from calling buffer: _own_memory = ", self._own_memory
            # end debug
            def _buffer(N):
                cdef double* ptr = self.thisptr.xAddress()
                cdef view.array my_arr
                my_arr = <double[:N, :3]> ptr
                return my_arr
            return _buffer(self.n_atoms)

    property xyz:
        def __get__(self):
            """return numpy array as a view of Frame xyz coords"""
            return np.asarray(self._buffer2d)

        def __set__(self, value):
            if not hasattr(value, 'shape') and value.shape != self.shape:
                raise ValueError("shape mismatch")
            self.xyz[:] = value

    property velocity:
        def __get__(self):
            """velocity, default None
            """
            cdef int n_atoms = self.n_atoms
            if self.has_velocity():
                return np.asarray(<double[:n_atoms, :3]> self.thisptr.vAddress())
            else:
                return None

    property force:
        """force, default None
        """
        def __get__(self):
            cdef int n_atoms = self.n_atoms
            if self.has_force():
                return np.asarray(<double[:n_atoms, :3]> self.thisptr.fAddress())
            else:
                return None

    property n_atoms:
        '''n_atoms
        '''
        def __get__(self):
            '''
            '''
            return self.thisptr.Natom()

    property size:
        def __get__(self):
            """don't change the name of this method"""
            return self.thisptr.size()

    property temperature:
        '''temperature
        '''

        def __get__(self):
            return self.thisptr.Temperature()

        def __set__(self, double tin):
            self.thisptr.SetTemperature(tin)

    property time:
        '''time'''
        def __get__(self):
            return self.thisptr.Time()

        def __set__(self, double timein):
            self.thisptr.SetTime(timein)

    def _update_atoms(self, indices, xyz):
        """TODO: add doc"""
        # xyz : 1D array
        if len(indices) != len(xyz)/3:
            raise ValueError("TODO: add doc")
        self._cy_update_atoms(pyarray('i', indices), pyarray('d', xyz), len(indices))

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef _cy_update_atoms(self, int[:] indices, double[:] xyz, int N):
        """Update atom coordinates with indices and given xyz array
        Parameters:
        ----------
        indices : int[:]
        xyz : double[:]
        N : int (len of indices)
        """
        cdef double* ptr
        cdef int i, idx

        for i in range(N):
            idx = indices[i]
            ptr = self.thisptr.xAddress() + 3 * idx
            ptr[0] = xyz[3*i + 0]
            ptr[1] = xyz[3*i + 1]
            ptr[2] = xyz[3*i + 2]

    def atom(self, int atomnum):
        """return xyz coordinates of idx-th atom

        Examples
        --------
        >>> import pytraj as pt
        >>> traj = pt.datafiles.load_tz2()
        >>> frame = traj[0]
        >>> frame.atom(10)
        array('d', [0.27399998903274536, 11.727999687194824, 8.701000213623047])
        """
        # return XYZ for atomnum
        # cpptraj: return double*
        # use python array to store double*
        # can we change it?

        cdef int i
        cdef pyarray arr = pyarray('d', [])

        if atomnum >= self.n_atoms or atomnum < 0:
            raise ValueError("Index is out of range")

        # TODO: check if this is not empty Frame
        arr.append(self.thisptr.XYZ(atomnum)[0])
        arr.append(self.thisptr.XYZ(atomnum)[1])
        arr.append(self.thisptr.XYZ(atomnum)[2])
        return arr

    def set_nobox(self):
        '''set nobox

        Examples
        --------
        >>> import pytraj as pt
        >>> traj = pt.datafiles.load_tz2_ortho()
        >>> frame = traj[0]
        >>> frame.box
        <Box: ortho, (x, y, z, alpha, beta, gamma) = (35.262779662258, 41.845547679864616, 36.16862952899312, 90.0, 90.0, 90.0)>
        >>> frame.set_nobox()
        >>> frame.box
        <Box: nobox, (x, y, z, alpha, beta, gamma) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)>
        '''
        self._boxview[:] = pyarray('d', [0. for _ in range(6)])

    def has_box(self):
        return self.box.has_box()

    def has_velocity(self):
        return self.thisptr.HasVelocity()

    def has_force(self):
        return self.thisptr.HasForce()

    property coordinates:
        '''return a copy of Frame's coordinates
        '''
        def __get__(self):
            return np.array(self.xyz)

    property mass:
        """return mass array"""
        def __get__(self):
            cdef double[:] arr = np.empty(self.n_atoms, dtype='f8')
            cdef int i

            for i in range(self.thisptr.Natom()):
                arr[i] = self.thisptr.Mass(i)
            return np.array(arr)

    property box:
        '''unitcell

        Example
        -------
        >>> import pytraj as pt
        >>> frame = pt.Frame(304)
        >>> frame.box
        >>> traj = pt.datafiles.load_tz2_ortho()
        >>> frame = traj[0]
        >>> frame.box
        '''
        def __get__(self):
            '''return unitcell'''
            cdef Box box = Box()
            box.thisptr.SetBox(self.thisptr.bAddress())
            return box

        def __set__(self, other):
            """
            other : {Box, array-like}
            """
            _box = Box(other)
            self._boxview[:] = _box[:]

    property _boxview:
        def __get__(self):
            """return a memoryview of box array"""
            cdef double* ptr = self.thisptr.bAddress()
            cdef view.array my_arr
            my_arr = <double[:6]> ptr
            return my_arr

    def set_mass(self, Topology top):
        '''set mass for Frame, requires a Topology
        '''
        self.thisptr.SetMass(top.thisptr.Atoms())

    def _set_mass_from_array(self, double[:] arr):
        '''mostly for pickling
        '''
        cdef Atom atom
        cdef vector[_Atom] va
        cdef unsigned int i

        for i in range(arr.shape[0]):
            # just try to create dummy atom (name, type, charge, mass)
            atom = Atom(name='X', type='X', charge=0.0, mass=arr[i])
            va.push_back(atom.thisptr[0])
        self.thisptr.SetMass(va)

    def __iadd__(Frame self, value):
        cdef Frame other
        # +=
        # either of two methods are correct
        #self.thisptr[0] = self.thisptr[0].addequal(other.thisptr[0])
        if isinstance(value, Frame):
            other = <Frame> value
            self.thisptr[0] += other.thisptr[0]
        else:
            self.xyz[:] += value
        return self

    def __add__(self, value):
        cdef Frame other
        cdef Frame frame

        frame = Frame(self.n_atoms)
        if isinstance(value, Frame):
            other = value
            frame.xyz = self.xyz + other.xyz
        else:
            frame.xyz = self.xyz + value
        return frame

    def __sub__(Frame self, value):
        cdef Frame other
        cdef Frame frame = Frame()

        if isinstance(value, Frame):
            other = <Frame> value
            frame.thisptr[0] = self.thisptr[0] - other.thisptr[0]
        else:
            frame = Frame(self)
            frame.xyz -= value
        return frame

    # NOTE: nogain with openmp for iadd, isub, ...
    # (and not faster than numpy, even with 500K atoms)
    def __isub__(Frame self, value):
        cdef Frame other
        # -=
        # either of two methods are correct
        #self.thisptr[0] = self.thisptr[0].subequal(other.thisptr[0])
        if isinstance(value, Frame):
            other = value
            self.thisptr[0] -= other.thisptr[0]
        else:
            self.xyz[:] -= value
        return self

    # NOTE: nogain with openmp for iadd, isub, ...
    # (and not faster than numpy, even with 500K atoms)
    def __imul__(Frame self, value):
        cdef Frame other
        # *=
        # either of two methods are correct
        #self.thisptr[0] = self.thisptr[0].mulequal(other.thisptr[0])
        if isinstance(value, Frame):
            other = value
            self.thisptr[0] *= other.thisptr[0]
        else:
            self.xyz[:] *= value
        return self

    def __mul__(Frame self, value):
        cdef Frame frame = Frame()
        cdef Frame other

        if isinstance(value, Frame):
            other = value
            frame.thisptr[0] = self.thisptr[0] * other.thisptr[0]
        else:
            frame = Frame(self)
            frame.xyz[:] *= value
        return frame

    def __tmp_idiv__(self, value):
        cdef Frame other
        cdef int i
        cdef int natom3 = self.n_atoms * 3
        cdef double* self_ptr
        cdef double* other_ptr

        if isinstance(value, Frame):
            other = value
            self_ptr = self.thisptr.xAddress()
            other_ptr = other.thisptr.xAddress()
            for i in range(natom3):
                self_ptr[i] /= other_ptr[i]
        else:
            self.xyz /= value

    def __idiv__(self, value):
        self.__tmp_idiv__(value)
        return self

    def __itruediv__(self, value):
        self.__tmp_idiv__(value)
        return self

    def __div__(self, value):
        cdef Frame other, frame

        frame = Frame(self.n_atoms)
        if isinstance(value, Frame):
            other = value
            frame.xyz = self.xyz / other.xyz
        else:
            frame.xyz = self.xyz / value
        return frame

    def __truediv__(self, value):
        cdef Frame other, frame

        frame = Frame(self.n_atoms)
        if isinstance(value, Frame):
            other = value
            frame.xyz = self.xyz / other.xyz
        else:
            frame.xyz = self.xyz / value
        return frame

    def center_of_mass(self, AtomMask atmask):
        """return Vec3 storing center of mass"""
        cdef Vec3 v3 = Vec3()
        v3.thisptr[0] = self.thisptr.VCenterOfMass(atmask.thisptr[0])
        return v3

    def center_of_geometry(self, AtomMask atmask):
        """return Vec3 storing center of geometry"""
        cdef Vec3 v3 = Vec3()
        v3.thisptr[0] = self.thisptr.VGeometricCenter(atmask.thisptr[0])
        return v3

    def _trans_rot_trans(self, Vec3 vec3, Matrix_3x3 m3, Vec3 vec3_2):
        # TODO : add doc, make test case
        self.thisptr.Trans_Rot_Trans(vec3.thisptr[0], m3.thisptr[0], vec3_2.thisptr[0])

    def rmsd(self, Frame frame, AtomMask atommask=None,
             mask=None, top=None,
             bint mass=False, get_mvv=False):
        """Calculate rmsd betwen two frames

        Parameters
        ----------
        frame : Frame instance
        mass : bool, default = False
        get_mvv : bool
            if True: return rmsd, Matrix_3x3, Vec3, Vec3
            if False: return rmsd
        """

        cdef Matrix_3x3 m3
        cdef Vec3 v1, v2
        if top is not None and mask is not None and atommask is None:
            atm = AtomMask(mask)
            top.set_integer_mask(atm)
            new_self = Frame(self, atm)
            new_ref = Frame(frame, atm)
        if top is None and mask is None and atommask is not None:
            new_self = Frame(self, atommask)
            new_ref = Frame(frame, atommask)
        if top is None and mask is None and atommask is None:
            # we need to make a copy since cpptraj update coords of frame after rmsd calc
            # all atoms
            new_self = Frame(self)
            new_ref = Frame(frame)

        if not get_mvv:
            return new_self.thisptr.RMSD(new_ref.thisptr[0], mass)
        else:
            m3 = Matrix_3x3()
            v1, v2 = Vec3(), Vec3()
            rmsd_ = new_self.thisptr.RMSD(new_ref.thisptr[0], m3.thisptr[0],
                                          v1.thisptr[0], v2.thisptr[0], mass)
            return rmsd_, m3, v1, v2

    def rmsd_nofit(self, Frame frame, bint mass=False):
        """Calculate rmsd betwen two frames without fitting

        Parameters
        ----------
        frame : Frame instance
        mass : bool, default = False
        """
        return self.thisptr.RMSD_NoFit(frame.thisptr[0], mass)

    def rmsfit(self, ref=None, AtomMask atm=None):
        """do the fitting to reference Frame by rotation and translation
        """
        # not yet dealed with `mass` and box
        cdef Matrix_3x3 mat
        cdef Vec3 v1

        _, mat, v1, v2 = self.rmsd(ref, atm, get_mvv=True)
        self._trans_rot_trans(v1, mat, v2)

    def _set_axis_of_rotation(self, int atom1, int atom2):
        cdef Vec3 vec = Vec3()
        vec.thisptr[0] = self.thisptr.SetAxisOfRotation(atom1, atom2)
        return vec

    def _calc_inertia(self, AtomMask atm, Matrix_3x3 Inertia):
        cdef Vec3 vec = Vec3()
        vec.thisptr[0] = self.thisptr.CalculateInertia(atm.thisptr[0], Inertia.thisptr[0])
        return vec

    def strip(self, AtomMask atm):
        """strip

        Parameters
        ----------
        atm: AtomMask

        Returns
        -------
        self
        """
        atm.invert_mask()
        cdef Frame frame = Frame(atm.n_atoms)

        frame.thisptr.SetFrame(self.thisptr[0], atm.thisptr[0])
        # deallocate old coordinates
        del self.thisptr
        # point to the new one
        self.thisptr = frame.thisptr
        # do not let ``frame`` deallocate, let ``self`` do it
        frame._own_memory = False
        return self

    property top:
        def __get__(self):
            return self._top

        def __set__(self, value):
            self._top = value

    def _dihedral(self, cython.integral[:, :] int_arr):
        """return python array of dih angle for four atoms with indices idx1-4
        Parameters
        ----------
        int_arr : 2D array of int with shape=(n_groups, 4)

        Returns
        -------
        1D python array
        """
        cdef int id0, idx1, idx2, idx3
        cdef int n_arr = int_arr.shape[0]
        cdef int i
        cdef pyarray arr0 = cparray.clone(pyarray('d', []), n_arr, zero=False)
        cdef double[:] arr0_view = arr0

        for i in range(n_arr):
            idx0 = int_arr[i, 0]
            idx1 = int_arr[i, 1]
            idx2 = int_arr[i, 2]
            idx3 = int_arr[i, 3]
            arr0_view[i] = RADDEG * cpptorsion(
                self.thisptr.XYZ(idx0),
                self.thisptr.XYZ(idx1),
                self.thisptr.XYZ(idx2),
                self.thisptr.XYZ(idx3))
        return arr0

    def _angle(self, cython.integral[:, :] int_arr):
        """return python array of angles for three atoms with indices idx1-3
        Parameters
        ----------
        int_arr : 2D array of int with shape=(n_groups, 3)

        Returns
        -------
        1D python array
        """
        cdef int idx0, idx1, idx2
        cdef int n_arr = int_arr.shape[0]
        cdef int i
        cdef pyarray arr0 = cparray.clone(pyarray('d', []), n_arr, zero=False)
        cdef double[:] arr0_view = arr0

        for i in range(n_arr):
            idx0 = int_arr[i, 0]
            idx1 = int_arr[i, 1]
            idx2 = int_arr[i, 2]
            arr0_view[i] = RADDEG * cppangle(self.thisptr.XYZ(idx0),
                                             self.thisptr.XYZ(idx1), self.thisptr.XYZ(idx2))
        return arr0

    def _distance(self, arr, parallel=False):
        # TODO: use `cdef _calc_distance`
        # need to make nested function to use default `parallel` value (=False)
        # if not, will get error: Special method __defaults__
        # has wrong number of arguments
        return self._calc_distance(arr, parallel)

    @cython.boundscheck(False)
    def _calc_distance(self, cython.integral[:, :] int_arr, bint parallel):
        """return python array of distance for two atoms with indices idx0, idx1
        Parameters
        ----------
        int_arr : 2D array of int with shape=(n_groups, 2)

        Returns
        -------
        1D python array
        """
        cdef int idx0, idx1
        cdef int n_arr = int_arr.shape[0]
        cdef int i
        cdef pyarray arr0 = cparray.clone(pyarray('d', []), n_arr, zero=False)
        cdef double[:] arr0_view = arr0

        if parallel:
            for i in prange(n_arr, nogil=True):
                idx0 = int_arr[i, 0]
                idx1 = int_arr[i, 1]
                arr0_view[i] = sqrt(
                    DIST2_NoImage(
                        self.thisptr.XYZ(idx0),
                        self.thisptr.XYZ(idx1)))
        else:
            for i in range(n_arr):
                # just duplicate code (ugly)
                idx0 = int_arr[i, 0]
                idx1 = int_arr[i, 1]
                arr0_view[i] = sqrt(
                    DIST2_NoImage(
                        self.thisptr.XYZ(idx0),
                        self.thisptr.XYZ(idx1)))
        return arr0

    def to_ndarray(self):
        """return a ndarray as a view of Frame's coordinates"""
        import numpy as np
        return np.asarray(self._buffer2d)

    def _to_dataframe(self, top=None):
        import pandas as pd
        import numpy as np
        if pd:
            if top is None:
                labels = list('xyzm')
                arr = np.vstack((self.xyz.T, self.mass)).T
            else:
                labels = ['resid', 'resname', 'atomname', 'mass',
                          'x', 'y', 'z']
                mass_arr = np.array(self.mass, dtype='f4')
                resid_arr = np.empty(mass_arr.__len__(), dtype='i')
                resname_arr = np.empty(mass_arr.__len__(), dtype='U4')
                atomname_arr= np.empty(mass_arr.__len__(), 'U4')

                for idx, atom in enumerate(top.atoms):
                    # TODO: make faster?
                    resid_arr[idx] = atom.resid
                    resname_arr[idx] = top._get_residue(atom.resid).name
                    atomname_arr[idx] = atom.name

                arr = np.vstack((resid_arr, resname_arr, atomname_arr,
                                 mass_arr, self.xyz.T)).T
            return pd.DataFrame(arr, columns=labels)
        else:
            raise ValueError("must have pandas")

    def __setstate__(self, state):
        # when pickle, python return an empty frame
        cdef int n_atoms
        del self.thisptr
        coords = state['coordinates']
        # allocate xyz and mass too
        n_atoms = coords.shape[0]
        self.thisptr = new _Frame(n_atoms)
        self.xyz[:] = coords
        self._set_mass_from_array(state['mass'])

    def __getstate__(self):
        # need to make a copy of xyz to avoid memory free
        # (need for parallel)
        # TODO: velocity?
        return {'coordinates': np.array(self.xyz, dtype='f8'),
                'mass': self.mass}

    def _allocate_force_and_velocity(self, Topology top, crdinfo):
        '''allocate_memory for velocity, force arrays

        Parameters
        ----------
        top : Topology
        crdinfo : dict

        Examples
        --------
        >>> import pytraj as pt
        >>> top = pt.tools.make_fake_topology(n_atoms=100)
        >>> # make empty frame
        >>> frame = pt.Frame()
        >>> frame.has_force()
        False
        >>> frame.has_velocity()
        False
        >>> frame.n_atoms
        0

        # allocate
        >>> frame._allocate_force_and_velocity(top, crdinfo={'has_force': True, 'has_velocity': True})
        >>> frame.has_force()
        True
        >>> frame.has_velocity()
        True
        >>> frame.n_atoms
        100
        '''
        cdef CoordinateInfo crdinfo_

        crdinfo2 = dict((k, v) for k, v in crdinfo.items())
        if 'box' not in crdinfo2:
            crdinfo2['box'] = top.box

        crdinfo_ = CoordinateInfo(crdinfo2)
        self.thisptr.SetupFrameV(top.thisptr.Atoms(), crdinfo_.thisptr[0])
