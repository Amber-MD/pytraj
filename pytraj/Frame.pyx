# distutils: language = c++

#from __future__ import absolute_import, division
# turn off `division` for automatically casting
from __future__ import absolute_import
cimport cython
from libc.math cimport sqrt
from cython cimport view
from cpython cimport array as cparray # for extend python array
from cpython.array cimport array as pyarray
from cython.parallel import prange
from cython.operator cimport dereference as deref
from libcpp.vector cimport vector
from libc.string cimport memcpy
from cpython.buffer cimport Py_buffer
import math
from pytraj.decorators import for_testing, iter_warning
from pytraj.decorators import name_will_be_changed
from pytraj.utils.check_and_assert import _import_numpy
from pytraj.utils.check_and_assert import is_int
from pytraj.ArgList import ArgList
from pytraj.trajs.Trajout import Trajout
from pytraj.externals.six import string_types
from pytraj.exceptions import *

DEF RADDEG       =   57.29577951308232

cdef extern from "TorsionRoutines.h" nogil:
    double cpptorsion "Torsion" (const double *, const double *, const double *, const double *)
    double cppangle "CalcAngle" (const double*, const double*, const double*)

cdef extern from "DistRoutines.h" nogil:
    double DIST2_NoImage(double*, double*)


# TODO : reogarnize memory view, there are too many ways to assess
# need to finalize

__all__ = ['Frame']

def check_instance(inst, clsname):
    if not isinstance(inst, clsname):
        raise ValueError("Must be instance of %s") % clsname.__name__

cdef class Frame (object):
    """
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
    >>> pt.Frame.from_ndarray(xyz)
    >>> frame = pt.Frame(304)
    >>> frame.append_xyz(xyz)
    >>> frame2 = pt.Frame(frame)
    """
    def __cinit__(self, *args):
        # Should I include topology in Frame?
        # May by not: memory
        # Include topology in Trajectory instance? 
        """Constructor for Frame instance
        >>> from pytraj.Frame import Frame
        >>> # created empty Frame instance
        >>> frame0 = Frame()
        >>> # create Frame instance with 304 atoms
        >>> frame1 = Frame(304)
        >>> # create a copy of frame1
        >>> frame2 = Frame(frame1)
        >>> # create a copy of frame with atommask
        >>> frame3 = Frame(frame1, atm_instance)
        """
        cdef Frame frame
        cdef AtomMask atmask
        cdef vector[_Atom] vt
        cdef Topology top
        cdef Atom at
        cdef list atlist
        cdef int natom, natom3
        cdef int i

        self.py_free_mem = True

        if not args:
            self.thisptr = new _Frame()
        else:
            if len(args) == 2:
                frame, atmask = args
                self.thisptr = new _Frame(frame.thisptr[0], atmask.thisptr[0])
            elif len(args) == 1:
                _, _np = _import_numpy()
                # copy Frame
                if isinstance(args[0], Frame):
                    frame = args[0]
                    self.thisptr = new _Frame(frame.thisptr[0])
                # creat a new Frame instance with natom
                elif isinstance(args[0], int):
                    natom = <int> args[0]
                    self.thisptr = new _Frame(natom)
                elif isinstance(args[0], (list, tuple, pyarray, _np.ndarray)):
                    # TODO : specify all things similar to list or array
                    natom3 = <int> len(args[0])
                    self.thisptr = new _Frame(natom3/3)
                    for i in range(natom3):
                        self.set_from_crd(pyarray('d', args[0]))
                else:
                    # Create Frame from list of atom mask
                    atlist = args[0]
                    for at in atlist:
                        vt.push_back(at.thisptr[0])
                    self.thisptr = new _Frame(vt)
            else:
                raise ValueError()

    def __dealloc__(self):
        if self.py_free_mem and self.thisptr:
            del self.thisptr

    def __del__(self):
        del self.thisptr

    def same_coords_as(self, Frame other):
        """check if two frames having the same coords"""
        return (self.coords == other.coords)

    def copy(self):
        """return a copy"""
        cdef Frame frame = Frame(self)
        return frame

    def is_empty(self):
        return self.size == 0

    def set_from_crd(self, CRDtype farray, *args):
        """"""
        cdef int numCrd, numBoxCrd
        cdef bint hasVel
        cdef AtomMask mask

        if len(args) == 3:
            numCrd, numBoxCrd, hasVel = args
            self.thisptr.SetFromCRD(farray, numCrd, numBoxCrd, hasVel)
        elif len(args) == 4:
            mask, numCrd, numBoxCrd, hasVel = args
            self.thisptr.SetFromCRD(farray, mask.thisptr[0], numCrd, numBoxCrd, hasVel)
        else:
            # TODO : shape checking
            numCrd = len(farray)
            self.thisptr.SetFromCRD(farray, numCrd, 0, False)

    def info(self, char* msg=''):
        self.thisptr.Info(msg)

    def clear_atoms(self):
        self.thisptr.ClearAtoms()

    def append_xyz(self, double[:, :] xyz):
        """append 3D array and return itself
        """
        cdef int i
        cdef int N = xyz.shape[0]

        for i in range(N):
            self.thisptr.AddXYZ(&xyz[i, 0])
        return self

    cdef void _append_xyz_2d(self, double[:, :] xyz):
        # for internal use
        # TODO: add assert
        cdef int i
        cdef int N = xyz.shape[0]

        for i in range(N):
            self.thisptr.AddXYZ(&xyz[i, 0])

    cdef void _append_xyz_1d(self, double[:] xyz):
        # TODO: add assert
        # for internal use
        cdef int i
        cdef int N = <int> xyz.shape[0] / 3

        for i in range(N):
            self.thisptr.AddXYZ(&xyz[i*3])

    def append_vec3(self, Vec3 vec):
        self.thisptr.AddVec3(vec.thisptr[0])

    def swap_atoms(self, int atom1, int atom2):
        """
        Parameters
        ----------
        atom1 : int
        atom2 : int

        """
        self.thisptr.SwapAtoms(atom1, atom2)

    def swap_atom_array(self, cython.integral[:, :] int_view):
        """
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
        return self.thisptr == other.thisptr

    @property
    def shape(self):
        return self.buffer2d[:].shape

    @property
    def n_frames(self):
        return 1

    def __getitem__(self, idx):
        """
        Examples
        --------
        >>> from pytraj import io
        >>> traj = io.load_sample_data('tz2')
        >>> f0 = traj[0]
        >>> f0[0]
        >>> f0[0, 0]
        >>> f0[:,  0]
        >>> f0.top = traj.top
        >>> f0['@CA']
        >>> atm = traj.top.select("@CB")
        >>> f0[atm]
        >>> f0[atm, 0]
        """
        cdef AtomMask atm
        cdef cython.view.array cy_arr
        cdef int new_size
        cdef int i, j
        cdef int[:] int_view

        has_numpy, np = _import_numpy()
        if isinstance(idx, pyarray):
            if not has_numpy:
                # create memoryview
                int_view = idx
                new_size = int_view.shape[0]
                cy_arr = cython.view.array(shape=(new_size, 3), 
                         itemsize=sizeof(double), format='d')
                for i in range(new_size):
                    # get index for `self`
                    j = int_view[i]
                    cy_arr[i] = self[j]
                return cy_arr
            else:
                return self.xyz[idx]
        elif isinstance(idx, AtomMask):
            # return a sub-array copy with indices got from 
            # idx.selected_indices()
            # TODO : add doc
            if idx.n_atoms == 0:
                raise ValueError("emtpy mask")
            if not has_numpy:
                return self[idx.indices]
            else:
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
            if has_numpy:
                return self.xyz[idx]
            else:
                return self.buffer2d[idx]

    def __setitem__(self, idx, value):
        has_np, np = _import_numpy()

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

        has_numpy, np = _import_numpy()
        for i in range(self.n_atoms):
            if has_numpy:
                yield np.asarray(self.buffer2d[i])
            else:
                yield self.buffer2d[i]

    def __array__(self):
        """
        arr0 = np.asarray(frame)

        (== (arr0 = np.asarray(frame[:])))
        """
        _, np = _import_numpy()
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
                int_view = indices # create `view`
            except:
                int_view = indices.astype('i4')
            # NOTE: try `prange` with different `schedule` but no gain
            for i in range(int_view.shape[0]):
                j = int_view[i]
                ptr_dest = self.thisptr.xAddress() + j * 3
                ptr_src = &xyz[i, 0]
                # copy coords of each Atom
                memcpy(<void*> ptr_dest, <void*> ptr_src, count)


    def frame_iter(self):
        """
        """
        yield self

    def __len__(self):
        return self.size

    @property
    def buffer1d(self):
        """return memory view for Frame coordinates
        TODO : rename?
        """
        # debug
        #print "from calling buffer1d: py_free_mem = ", self.py_free_mem
        # end debug
        #"name_will_be_changed, this is for development"
        def _buffer(N):
            cdef double* ptr = self.thisptr.xAddress()
            cdef view.array my_arr
            my_arr = <double[:N]> ptr
            return my_arr
        return _buffer(self.size)

    @property
    def buffer2d(self):
        """return memory view for Frame coordinates but reshape
        (just like self._buffer3 = self.buffer.reshape())
        TODO : rename?
        """
        # debug
        #print "from calling buffer: py_free_mem = ", self.py_free_mem
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
            has_np, np = _import_numpy()
            if has_np:
                return np.asarray(self.buffer2d)
            else:
                raise NotImplementedError("need numpy. Use `buffer2d` instead")
        def __set__(self, value):
            if not hasattr(value, 'shape') and value.shape != self.shape:
                raise ValueError("shape mismatch")
            self.xyz[:] = value
        
    def is_empty(self):
        return self.thisptr.empty()

    def has_vel(self):
        return self.thisptr.HasVelocity()

    @property
    def n_atoms(self):
       return self.thisptr.Natom()

    @property
    def size(self):
        """don't change the name of this method"""
        return self.thisptr.size()

    @property
    def n_repdims(self):
        return self.thisptr.NrepDims()

    property temperature:
        def __get__(self):
            return self.thisptr.Temperature()
        def __set__(self, double tin):
            self.thisptr.SetTemperature(tin)

    property time:
        def __get__(self):
            return self.thisptr.Time()
        def __set__(self, double timein):
            self.thisptr.SetTime(timein)

    def update_atom(self, int idx, double[:] xyz):
        cdef double* ptr = self.thisptr.xAddress() + 3 * idx
        if len(xyz) > 3:
            raise IndexError("")
        ptr[0] = xyz[0]
        ptr[1] = xyz[1]
        ptr[2] = xyz[2]

    def update_atoms(self, indices, xyz):
        """TODO: add doc"""
        # xyz : 1D array
        if len(indices) != len(xyz)/3:
            raise ValueError("TODO: add doc")
        self._update_atoms(pyarray('i', indices), pyarray('d', xyz), len(indices))

    @cython.boundscheck(False)
    @cython.wraparound(False)
    cdef _update_atoms(self, int[:] indices, double[:] xyz, int N):
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

    def atoms(self, int atomnum):
        """return xyz coordinates of idx-th atom"""
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

    @property
    def coordinates(self):
        return self.xyz

    @property
    def coords(self):
        """
        return 1D-coords (copy) of Frame
        Notes
        -----
        same as `Frame.coords`. We use `coordinates` to be the same as in
        `parmed`
        """
        cdef pyarray arr = cparray.clone(pyarray('d', []), 
                            self.n_atoms*3, zero=False)
        cdef int i
        cdef double* ptr = self.thisptr.xAddress()
        cdef int natom3 = 3 * self.thisptr.Natom()

        for i in range(natom3):
            arr[i] = ptr[i]
        return arr

    def v_xyz(self, int atnum):
        """return a copy of velocity"""
        cdef int i
        cdef pyarray arr = pyarray('d', [])

        arr.append(self.thisptr.VXYZ(atnum)[0])
        arr.append(self.thisptr.VXYZ(atnum)[1])
        arr.append(self.thisptr.VXYZ(atnum)[2])
        return arr

    @property
    def mass(self):
        """return mass array"""
        cdef pyarray arr = pyarray('d', [])
        cdef int i

        for i in range(self.thisptr.Natom()):
            arr.append(self.thisptr.Mass(i))
        return arr

    def set_nobox(self):
        self.boxview[:] = pyarray('d', [0. for _ in range(6)])

    def box_crd(self):
        cdef Box box = Box()
        box.thisptr[0] = self.thisptr.BoxCrd()
        return box.tolist()

    def v_address(self):
        return <double[:self.thisptr.Natom(), :3]> self.thisptr.vAddress()

    property box:
        def __get__(self):
            cdef Box box = Box()
            box.thisptr.SetBox(self.thisptr.bAddress())
            return box
        def __set__(self, other):
            """
            other : {Box, array-like}
            """
            _box = Box(other)
            self.boxview[:] = _box[:]

    def has_box(self):
        return self.box.has_box()

    @property
    def boxview(self):
        """return a memoryview of box array"""
        cdef double* ptr = self.thisptr.bAddress()
        cdef view.array my_arr
        my_arr = <double[:6]> ptr
        return my_arr

    def _t_address(self):
        # cpptraj: return double*
        raise NotImplementedError()

    def _i_address(self):
        # cpptraj: return int*
        raise NotImplementedError()

    def set_box_angles(self, double[:] ain):
        self.thisptr.SetBoxAngles(&ain[0])

    def set_frame(self, *args):
        cdef int atomnum 
        cdef AtomMask atm
        cdef Frame frame

        if isinstance(args[0], Frame) and isinstance(args[1], AtomMask):
            frame = args[0]
            atm = args[1]
            self.thisptr.SetFrame(frame.thisptr[0], atm.thisptr[0])
        elif is_int(args):
            atomnum = <int> args[0]
            self.thisptr.SetupFrame(atomnum)

    def set_frame_mass(self, Topology top):
        return self.thisptr.SetupFrameM(top.thisptr.Atoms())

    def set_frame_x_m(self, vector[double] Xin, vector[double] massIn):
        return self.thisptr.SetupFrameXM(Xin, massIn)

    def set_frame_v(self, Topology top):
        """TODO: add doc
        """
        return self.thisptr.SetupFrameV(top.thisptr.Atoms(), top.thisptr.ParmCoordInfo())

    def set_frame_from_mask(self, mask, atomlist_or_top):
        cdef Atom atom
        cdef vector[_Atom] v 
        cdef AtomMask atm
        cdef Topology top

        if isinstance(mask, string_types):
            # if providing mask string
            # create AtomMask instance
            atm = AtomMask(mask)
            # assume atht atomlist_or_top is Topology instance
            atomlist_or_top.set_integer_mask(atm)
        elif isinstance(mask, AtomMask):
            atm = mask 
        else:
            raise NotImplementedError("must be `str` or AtomMask")

        if isinstance(atomlist_or_top, (list, tuple)):
            # list of atoms
            for atom in atomlist_or_top:
                v.push_back(atom.thisptr[0])
            return self.thisptr.SetupFrameFromMask(atm.thisptr[0], v)
        elif isinstance(atomlist_or_top, Topology):
            top = atomlist_or_top
            return self.thisptr.SetupFrameFromMask(atm.thisptr[0], top.thisptr.Atoms())
        else:
            raise NotImplementedError("must be list of Atom objects or Topology")

    def set_coords(self, Frame frame, *args):
        """set coords for Frame with given mask
        Paramters:
        ---------
        frame : Frame object
        atmask : AtomMask object (optional)
        """
        cdef AtomMask atmask 
        if not args:
            self.thisptr.SetCoordinates(frame.thisptr[0])
        else:
            atmask = args[0]
            self.thisptr.SetCoordinates(frame.thisptr[0], atmask.thisptr[0])

    def set_coords_by_map(self, Frame frame, vector[int] mapIn):
        self.thisptr.SetCoordinatesByMap(frame.thisptr[0], mapIn)

    def zero_coords(self):
        self.thisptr.ZeroCoords()

    # NOTE: nogain with openmp for iadd, isub, ...
    # (and not faster than numpy, even with 500K atoms) 
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

    def divide(self, double divisor, *args):
        # TODO : check
        cdef Frame frame
        if not args:
            self.thisptr.Divide(divisor)
        else:
            frame = args[0]
            return self.thisptr.Divide(frame.thisptr[0], divisor)

    def add_by_mask(self, Frame frame, AtomMask atmask):
        """Increment atoms in `self` by selected atoms from `frame`
        """
        self.thisptr.AddByMask(frame.thisptr[0], atmask.thisptr[0])

    def check_coords_invalid(self):
        return self.thisptr.CheckCoordsInvalid()

    def center_of_mass(self, AtomMask atmask):
        """return Vec3"""
        cdef Vec3 v3 = Vec3()
        v3.thisptr[0] = self.thisptr.VCenterOfMass(atmask.thisptr[0])
        return v3

    def center_of_geometry(self, AtomMask atmask):
        """return Vec3"""
        cdef Vec3 v3 = Vec3()
        v3.thisptr[0] = self.thisptr.VGeometricCenter(atmask.thisptr[0])
        return v3

    def translate(self, *args):
        # TODO : doc
        """
        Paramters:
        ---------
               (Vec3, first_atom_idx, last_atom_idx)
            or (Vec3, atom_idx) 
            or (Vec3)
        TODO: add doc
        """
        cdef firstAtom, lastAtom, atom
        cdef Vec3 vec3

        if not args:
            raise ValueError()

        if isinstance(args[0], Vec3):
            vec3 = args[0]
        else:
            # try to convert to Vec3. no warranty :D
            vec3 = Vec3(args[0])

        if len(args) == 3:
            vec3, firstAtom, lastAtom = args
            self.thisptr.Translate(vec3.thisptr[0], firstAtom, lastAtom)
        elif len(args) == 2:
            vec3, atom = args
            check_instance(vec3, Vec3)
            self.thisptr.Translate(vec3.thisptr[0], atom)
        elif len(args) == 1:
            self.thisptr.Translate(vec3.thisptr[0])
        else:
            raise ValueError()

    def neg_translate(self, Vec3 vec):
        self.thisptr.NegTranslate(vec.thisptr[0])

    def rotate_with_matrix(self, mat, *args):
        """
        Parameters
        ----------
        mat : Matrix-like, shape=(3,3)
            3x3 matrix (pytraj or numpy)
        """
        cdef AtomMask atm
        cdef Matrix_3x3 _mat

        has_numpy, np = _import_numpy()
        if not has_numpy:
            assert isinstance(mat, Matrix_3x3)
        if isinstance(mat, np.matrix) or isinstance(mat, np.ndarray):
            _mat = Matrix_3x3(mat)
        else:
            # assume Matrix_3x3
            _mat = mat
        if args:
            atm = <AtomMask> args[0]
            self.thisptr.Rotate(_mat.thisptr[0], atm.thisptr[0])
        else:
            self.thisptr.Rotate(_mat.thisptr[0])

    def rotate(self, int x=0, int y=0, int z=0, atommask=None):
        """rotate(Matrix_3x3 m3, *args)
        Paramters:
        m3 : Matrix_3x3
        *args : optional
            atmmaks : AtomMaks instance
         or (mask (str), Topology instance)

        """
        cdef Matrix_3x3 m3 = Matrix_3x3()
        cdef AtomMask atmask
        cdef string mask
        cdef Topology top

        m3.thisptr.CalcRotationMatrix(math.radians(x), 
                                      math.radians(y), 
                                      math.radians(z))

        if atommask is None:
            self.thisptr.Rotate(m3.thisptr[0])
        else:
            atmask = <AtomMask> atommask
            self.thisptr.Rotate(m3.thisptr[0], atmask.thisptr[0])

    def trans_rot_trans(self, Vec3 vec3, Matrix_3x3 m3, Vec3 vec3_2):
        # TODO : add doc, make test case
        self.thisptr.Trans_Rot_Trans(vec3.thisptr[0], m3.thisptr[0], vec3_2.thisptr[0])

    def scale(self, AtomMask atm, double sx, double sy, double sz):
        # TODO : add doc, make test case
        self.thisptr.Scale(atm.thisptr[0], sx, sy, sz)

    def center_on_origin(self,bint useMassIn):
        # TODO : add doc, make test case
        cdef Vec3 v = Vec3()
        v.thisptr[0] = self.thisptr.CenterOnOrigin(useMassIn)
        return v

    def rmsd(self, Frame frame, AtomMask atommask=None, 
             mask=None, top=None,
             bint use_mass=False, get_mvv=False):
        # TODO : use_mass does not work properly
        """Calculate rmsd betwen two frames
        rmsd(Frame frame, bint use_mass=False, get_mvv=False):
        Parameters:
        ----------
        frame : Frame instance
        use_mass : bool, default = False
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
            return new_self.thisptr.RMSD(new_ref.thisptr[0], use_mass)
        else:
            m3 = Matrix_3x3()
            v1, v2 = Vec3(), Vec3()
            rmsd_ = new_self.thisptr.RMSD(new_ref.thisptr[0], m3.thisptr[0], 
                                          v1.thisptr[0], v2.thisptr[0], use_mass)
            return rmsd_, m3, v1, v2

    def rmsd_centered_ref(self, Frame ref, bint use_mass=False, *args):
        """Calculate rmsd betwen two frames
        Parameters:
        ----------
        frame : Frame instance
        use_mass : bool, default = False
        *args :  optional, 2 args (Matrix_3x3 instance, Vec3 instance)
        """
        cdef Matrix_3x3 mat 
        cdef Vec3 v

        if not args:
            return self.thisptr.RMSD_CenteredRef(ref.thisptr[0], use_mass)
        else:
            mat, v = args
            assert isinstance(mat, Matrix_3x3) == True
            assert isinstance(v, Vec3) == True
            return self.thisptr.RMSD_CenteredRef(ref.thisptr[0], mat.thisptr[0], v.thisptr[0], use_mass)

    def rmsd_nofit(self, Frame frame, bint use_mass=False):
        """Calculate rmsd betwen two frames without fitting
        Parameters:
        ----------
        frame : Frame instance
        use_mass : bool, default = False
        """
        return self.thisptr.RMSD_NoFit(frame.thisptr[0], use_mass)

    def dist_rmsd(self, Frame frame, atommask=None):
        """Calculate dist_rmsd betwen two frames
        Parameters:
        ----------
        frame : Frame instance
        """
        cdef Frame f1, f2

        if atommask is None:
            return self.thisptr.DISTRMSD(frame.thisptr[0])
        else:
            f1 = Frame(self, <AtomMask>atommask)
            f2 = Frame(frame, <AtomMask>atommask)
            return f1.thisptr.DISTRMSD(f2.thisptr[0])

    def rmsfit(self, ref=None, AtomMask atm=None):
        """do the fitting to reference Frame by rotation and translation
        TODO : add assert test
        """
        # not yet dealed with `mass` and box
        cdef Matrix_3x3 mat
        cdef Vec3 v1

        _, mat, v1, v2 = self.rmsd(ref, atm, get_mvv=True)
        self.trans_rot_trans(v1, mat, v2)

    def set_axis_of_rotation(self, int atom1, int atom2):
        cdef Vec3 vec = Vec3()
        vec.thisptr[0] = self.thisptr.SetAxisOfRotation(atom1, atom2)
        return vec

    def calc_inertia(self, AtomMask atm, Matrix_3x3 Inertia):
        cdef Vec3 vec = Vec3()
        vec.thisptr[0] = self.thisptr.CalculateInertia(atm.thisptr[0], Inertia.thisptr[0])
        return vec

    def calc_temperature(self, AtomMask mask, int deg_of_freedom):
        return self.thisptr.CalcTemperature(mask.thisptr[0], deg_of_freedom)

    cdef void _strip_atoms(Frame self, Topology top, AtomMask atm, bint update_top, bint has_box):
        """this method is too slow vs cpptraj
        if you use memory for numpy, you need to update after resizing Frame
        >>> arr0 = np.asarray(frame.buffer)
        >>> frame.strip_atoms(top,"!@CA")
        >>> # update view
        >>> arr0 = np.asarray(frame.buffer)
        """
        # NOTE:`atm` here is the KEPT-atommaks (for performance)
        # we will do `atm.invert_mask` laster in `strip_atoms` method
        # Important: don't try to use raw pointers. Make sure to create Python's objects
        # (Frame, Topology) to control cpptraj' objects' lifetime.

        cdef Topology newtop = Topology()
        newtop.py_free_mem = False
        cdef Frame tmpframe = Frame() 

        del tmpframe.thisptr

        tmpframe.thisptr = new _Frame(self.thisptr[0])
        
        newtop.thisptr = top.thisptr.modifyStateByMask(atm.thisptr[0])
        tmpframe.thisptr.SetupFrameV(newtop.thisptr.Atoms(), newtop.thisptr.ParmCoordInfo())
        tmpframe.thisptr.SetFrame(self.thisptr[0], atm.thisptr[0])

        # make a copy: coords, vel, mass...
        # if only care about `coords`, use `_fast_copy_from_frame`
        self.thisptr[0] = tmpframe.thisptr[0]
        if not has_box:
            self.set_nobox()
        if update_top:
            top.thisptr[0] = newtop.thisptr[0]

    def strip_atoms(Frame self, mask=None, Topology top=Topology(), 
                    bint update_top=False, bint has_box=False, bint copy=False):
        """strip_atoms(string mask, Topology top=Topology(), 
                       bint update_top=False, bint has_box=False)

        Return:
          None : if copy=False
          new striped Frame instance if copy=True

        Parameters:
        ----------
        mask : str, mask, non-default
        top : Topology, default=Topology()
        update_top : bint, default=False
        has_box : bint, default=False
        copy : bint, default=False
        """
        cdef AtomMask atm = top(mask)
        atm.invert_mask()

        if mask is None or top.is_empty():
            raise ValueError("need non-empty mask and non-empty Topology")
        cdef Frame frame

        mask = mask.encode("UTF-8")

        if not copy:
            self._strip_atoms(top, atm, update_top, has_box)
        else:
            frame = Frame(self)
            frame._strip_atoms(top, atm, update_top, has_box)
            return frame

    def get_subframe(self, mask=None, top=None):
        cdef AtomMask atm

        if isinstance(mask, string_types):
            assert top is not None
            atm = AtomMask(mask)
            top.set_integer_mask(atm)
        elif isinstance(mask, AtomMask):
            atm = <AtomMask> mask
            assert top is None
        else:
            raise NotImplementedError('mask mut string  or AtomMask object')
        return Frame(self, atm)

    property top:
        def __get__(self):
            return self._top
        def __set__(self, value):
            self._top = value

    def save(self, filename="", top=None, format='unknown', 
             overwrite=False, *args, **kwd):
        if format== 'unknown':
            # convert to "UNKNOWN_TRAJ"
            format= format.upper() + "_TRAJ"
        with Trajout(filename=filename, top=top, format=format, 
                     overwrite=overwrite, *args, **kwd) as trajout:
            trajout.write(0, self, top)

    def calc_dihedral(self, cython.integral[:, :] int_arr):
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
            arr0_view[i] = RADDEG * cpptorsion(self.thisptr.XYZ(idx0), self.thisptr.XYZ(idx1),
                          self.thisptr.XYZ(idx2), self.thisptr.XYZ(idx3))
        return arr0

    def calc_angle(self, cython.integral[:, :] int_arr):
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
            arr0_view[i] = RADDEG * cppangle(self.thisptr.XYZ(idx0), self.thisptr.XYZ(idx1),
                        self.thisptr.XYZ(idx2))
        return arr0

    def calc_distance(self, arr, parallel=False):
        # TODO: use `cdef _calc_distance`
        # need to make nested function to use default `parallel` value (=False)
        # if not, will get error: Special method __defaults__ 
        # has wrong number of arguments
        return self._calc_distance(arr, parallel)

    def _calc_distance(self, cython.integral [:, :] int_arr, bint parallel):
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
                arr0_view[i] = sqrt(DIST2_NoImage(self.thisptr.XYZ(idx0), self.thisptr.XYZ(idx1)))
        else:
            for i in range(n_arr):
                # just duplicate code (ugly)
                idx0 = int_arr[i, 0]
                idx1 = int_arr[i, 1]
                arr0_view[i] = sqrt(DIST2_NoImage(self.thisptr.XYZ(idx0), self.thisptr.XYZ(idx1)))
        return arr0

    def tolist(self):
        """return a list of 2D coords"""
        return [list(x) for x in self]

    def to_ndarray(self):
        """return a ndarray as a view of self.buffer2d"""
        has_np, np = _import_numpy()
        if has_np:
            return np.asarray(self.buffer2d)
        else:
            raise PytrajNumpyError()

    @classmethod
    def from_ndarray(cls, xyz):
        """create new Frame from a numpy.ndarray
        """
        return Frame().append_xyz(xyz)

    def to_dataframe(self, top=None):
        from pytraj.utils import _import_pandas
        _, pd = _import_pandas()
        if pd:
            _, np = _import_numpy()
            if top is None:
                labels = list('xyzm')
                arr = np.vstack((self.xyz.T, self.mass)).T
            else:
                labels = ['resnum', 'resname', 'atomname', 'mass',
                          'x', 'y', 'z']
                mass_arr = np.array(self.mass, dtype='f4')
                resnum_arr = np.empty(mass_arr.__len__(), dtype='i')
                resname_arr = np.empty(mass_arr.__len__(), dtype='U4')
                atomname_arr= np.empty(mass_arr.__len__(), 'U4')

                for idx, atom in enumerate(top.atoms):
                    # TODO: make faster?
                    resnum_arr[idx] = atom.resnum
                    resname_arr[idx] = top._get_residue(atom.resnum).name
                    atomname_arr[idx] = atom.name

                arr = np.vstack((resnum_arr, resname_arr, atomname_arr, 
                                 mass_arr, self.xyz.T)).T
            return pd.DataFrame(arr, columns=labels)
        else:
            raise ValueError("must have pandas")

    def as_3darray(self):
        return self.xyz.reshape((1, self.n_atoms, 3))

    def __setstate__(self, state):
        self.append_xyz(state['coordinates'])

    def __getstate__(self):
        return {'coordinates' : self.xyz}
