# distutils: language = c++

cimport cython
from cython cimport view
from cpython.array cimport array as pyarray
from cython.operator cimport dereference as deref
from libcpp.vector cimport vector
from cpython.buffer cimport Py_buffer
from pytraj._utils cimport _get_buffer1D

from pytraj.decorators import for_testing, iter_warning
from pytraj.decorators import name_will_be_changed
from pytraj.utils.check_and_assert import _import_numpy
from pytraj.ArgList import ArgList
from pytraj.trajs.Trajout import Trajout
from pytraj.externals.six import string_types

# TODO : reogarnize memory view, there are too many ways to assess
# need to finalize

def check_instance(inst, clsname):
    if not isinstance(inst, clsname):
        raise ValueError("Must be instance of %s") % clsname.__name__

cdef class Frame (object):
    """Original cpptraj doc (Frame.h) (written by Daniel R. Roe)
    (pytraj doc will be updated) 
    Class: Frame
        Hold coordinates, perform various operations/transformations on them.
        Intended to hold coordinates e.g. from a trajectory or reference frame,
        along with box coordinates (used in imaging calculations), mass information,
        and optionally velocity information. Frame can be set up coords only (all 
        masses set to 1.0), coords and masses, or coords/masses/velocities. Mass is 
        stored since several functions (like COM, RMSD, Inertia etc) have the option
        to factor in the mass of the atoms involved, and this avoids having to pass
        a mass pointer in, which takes the burden of keeping track of mass away from 
        actions etc. Mass is stored when the frame is initially created, and is 
        modified if necessary by SetFrame (which is the case when e.g. calculating
        per-residue RMSD).
        
        - Implementation Details:
        
        In addition to the constructors, there are two classes of routine that
        can be used to set up Frames. The SetupX routines do any memory allocation,
        and assign masses, and the SetX routines assign coordinates/velocities. The
        SetX routines will dynamically adjust the size of the frame up to maxnatom,
        but no reallocation will occur so the frame should be set up for the largest
        possible # of atoms it will hold. This avoids expensive reallocations.
        The representation of coordinates (X) and velocities (V) are double*
        instead of STL vectors so as to easily interface with the FileIO routines
        which tend to be much faster than iostream ops. 

        pytraj doc
        ============
        Should Frame hold topology info? (may be NOT, it's expesive)
        TODO : should we really need all methods here?
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

    def same_coords_as(self, Frame other):
        """check if two frames having the same coords"""
        return (self.coords == other.coords)

    def copy(self):
        """return a copy"""
        cdef Frame frame = Frame()
        frame.thisptr[0] = self.thisptr[0]
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

    def add_xyz(self, double[:] xyz):
        self.thisptr.AddXYZ(&xyz[0])

    def add_vec3(self, Vec3 vec):
        self.thisptr.AddVec3(vec.thisptr[0])

    def swap_atoms(self, int atom1, int atom2):
        self.thisptr.SwapAtoms(atom1, atom2)

    def __str__(self):
        tmp = "%s instance with %s atoms. ID = %s" % (
                self.__class__.__name__,
                self.n_atoms,
                hex(id(self))
                )
        return tmp

    @property
    def shape(self):
        return self.buffer3d[:].shape

    def __getitem__(self, idx):
        # always return memoryview 
        has_numpy, np = _import_numpy()
        if isinstance(idx, AtomMask):
            # return a sub-array copy with indices got from 
            # idx.selected_indices()
            # TODO : add doc
            if not has_numpy:
                raise NotImplementedError("supported if having numpy installed")
            arr0 = np.asarray(self.buffer3d[:])
            if isinstance(idx, AtomMask):
                if idx.n_selected == 0:
                    raise ValueError("emtpy mask")
                return arr0[np.array(idx.selected_indices())]
            else:
                # TODO : double-check if we can use np.ndarray or pyarray here?
                # isinstance(idx, pyarray)
                return arr0[idx]
        elif isinstance(idx, dict):
            # Example: frame[dict(top=top, mask='@CA')]
            # return a sub-array copy with indices got from 
            # idx as a `dict` instance
            # TODO : add doc
            if not has_numpy:
                raise NotImplementedError("supported if having numpy installed")
            atm = AtomMask(idx['mask'])
            idx['top'].set_integer_mask(atm)
            arr0 = np.asarray(self.buffer3d[:])
            return arr0[np.array(atm.selected_indices())]
        elif isinstance(idx, string_types):
            # Example: frame['@CA']
            if self.top is not None and not self.top.is_empty():
                return self[self.top(idx)]
            else:
                raise ValueError("must have non-empty topology")
        else:
            if has_numpy:
                return np.asarray(self.buffer3d[idx])
            else:
                return self.buffer3d[idx]

    def __setitem__(self, idx, value):
        # TODO : should we use buffer. Kind of dangerous
        # TODO : add examples hereo
        #if not isinstance(idx, (list, tuple)):
        #    self.buffer[idx] = value
        #else:
        #    if isinstance(value, (list, tuple)):
        #        value = pyarray('d', value)
        #    self.buffer3d[idx] = value
        if isinstance(value, (tuple, list)):
            value = pyarray('d', value)
            self.buffer3d[idx] = value
        elif isinstance(value, AtomMask):
            raise NotImplementedError("not yet support AtomMask indexing")
        else:
            self.buffer3d[idx] = value

    def __iter__(self):
        cdef int i

        has_numpy, np = _import_numpy()
        for i in range(self.n_atoms):
            if has_numpy:
                yield np.asarray(self.buffer3d[i])
            else:
                yield self.buffer3d[i]

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
    def buffer3d(self):
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

    @property 
    def temperature(self):
        return self.thisptr.Temperature()

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
    def coords(self):
        """
        return a copy of frame coords (python array)
        """
        cdef pyarray arr = pyarray('d', [])
        cdef int i

        for i in range(3 * self.thisptr.Natom()):
            arr.append(deref(self.thisptr.CRD(i)))
        return arr

    def v_xyz(self, int atnum):
        """return a copy of velocity"""
        cdef int i
        cdef pyarray arr = pyarray('d', [])

        arr.append(self.thisptr.VXYZ(atnum)[0])
        arr.append(self.thisptr.VXYZ(atnum)[1])
        arr.append(self.thisptr.VXYZ(atnum)[2])
        return arr

    def mass(self,int atnum):
        """return `mass` of atnum-th atom"""
        return self.thisptr.Mass(atnum)

    def set_nobox(self):
        self.boxview[:] = pyarray('d', [0. for _ in range(6)])

    def box_crd(self):
        cdef Box box = Box()
        box.thisptr[0] = self.thisptr.BoxCrd()
        return box.tolist()

    def v_address(self):
        # cpptraj: return double*
        raise NotImplementedError()

    def get_box(self):
        cdef Box box = Box()
        box.thisptr.SetBox(self.thisptr.bAddress())
        return box

    def has_box(self):
        box = self.get_box()
        return box.has_box()

    @property
    def boxview(self):
        """return a memoryview of box array"""
        cdef double* ptr = self.thisptr.bAddress()
        cdef view.array my_arr
        my_arr = <double[:6]> ptr
        return my_arr

    def t_address(self):
        # cpptraj: return double*
        raise NotImplementedError()

    def i_address(self):
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
        elif isinstance(args, (int, long)):
            atomnum = <int> args[0]
            self.thisptr.SetupFrame(atomnum)

    def set_frame_x_m(self, vector[double] Xin, vector[double] massIn):
        return self.thisptr.SetupFrameXM(Xin, massIn)

    # BROKEN WITH NEW CPPTRAJ
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

    def set_coord(self, Frame frame, *args):
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

    def __iadd__(Frame self, Frame other):
        # += 
        # either of two methods are correct
        #self.thisptr[0] = self.thisptr[0].addequal(other.thisptr[0])
        self.thisptr[0] += other.thisptr[0]
        return self

    def __sub__(Frame self, Frame other):
        cdef Frame frame = Frame()
        frame.thisptr[0] = self.thisptr[0] - other.thisptr[0]
        return frame

    def __isub__(Frame self, Frame other):
        # -= 
        # either of two methods are correct
        #self.thisptr[0] = self.thisptr[0].subequal(other.thisptr[0])
        self.thisptr[0] -= other.thisptr[0]
        return self

    def __imul__(Frame self, Frame other):
        # *=
        # either of two methods are correct
        #self.thisptr[0] = self.thisptr[0].mulequal(other.thisptr[0])
        self.thisptr[0] *= other.thisptr[0]
        return self

    def __mul__(Frame self, Frame other):
        cdef Frame frame = Frame()
        frame.thisptr[0] = self.thisptr[0] * other.thisptr[0]
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
        self.thisptr.AddByMask(frame.thisptr[0], atmask.thisptr[0])

    def check_coords_invalid(self):
        return self.thisptr.CheckCoordsInvalid()

    def VCenterOfMass(self, AtomMask atmask):
        # return Vec3 instance
        cdef Vec3 v3 = Vec3()
        v3.thisptr[0] = self.thisptr.VCenterOfMass(atmask.thisptr[0])
        return v3

    def VGeometricCenter(self, AtomMask atmask):
        # return Vec3 instance
        cdef Vec3 v3 = Vec3()
        v3.thisptr[0] = self.thisptr.VGeometricCenter(atmask.thisptr[0])
        return v3

    def translate(self, *args):
        """
        Args:
               (Vec3, first_atom_idx, last_atom_idx)
            or (Vec3, atom_idx) 
            or (Vec3)
        TODO: add doc
        """
        cdef firstAtom, lastAtom, atom
        cdef Vec3 vec3

        if len(args) == 3:
            vec3, firstAtom, lastAtom = args
            self.thisptr.Translate(vec3.thisptr[0], firstAtom, lastAtom)
        elif len(args) == 2:
            vec3, atom = args
            check_instance(vec3, Vec3)
            self.thisptr.Translate(vec3.thisptr[0], atom)
        elif len(args) == 1:
            vec3 = args[0]
            self.thisptr.Translate(vec3.thisptr[0])
        else:
            raise ValueError()

    def neg_translate(self, Vec3 vec):
        self.thisptr.NegTranslate(vec.thisptr[0])

    def rotate(self, Matrix_3x3 m3, *args):
        """rotate(Matrix_3x3 m3, *args)
        Paramters:
        m3 : Matrix_3x3
        *args : optional
            atmmaks : AtomMaks instance
         or (mask (str), Topology instance)

        """
        cdef AtomMask atmask
        cdef string mask
        cdef Topology top

        if not args:
            self.thisptr.Rotate(m3.thisptr[0])
        elif len(args) == 1:
            atmask = args[0]
            self.thisptr.Rotate(m3.thisptr[0], atmask.thisptr[0])
        elif len(args) == 2:
            mask, top = args
            atmask = AtomMask(mask)
            top.set_integer_mask(atmask)
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
            # all atoms
            new_self = self
            new_ref = frame

        if not get_mvv:
            return new_self.thisptr.RMSD(new_ref.thisptr[0], use_mass)
        else:
            m3 = Matrix_3x3()
            v1, v2 = Vec3(), Vec3()
            rmsd_ = new_self.thisptr.RMSD(new_ref.thisptr[0], m3.thisptr[0], v1.thisptr[0], v2.thisptr[0], use_mass)
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
            return self.thisptr.RMSD_CenteredRef(ref.thisptr[0], mat.thisptr[0], v.thisptr[0], use_mass)

    def rmsd_nofit(self, Frame frame, bint use_mass=False):
        """Calculate rmsd betwen two frames without fitting
        Parameters:
        ----------
        frame : Frame instance
        use_mass : bool, default = False
        """
        return self.thisptr.RMSD_NoFit(frame.thisptr[0], use_mass)

    def dist_rmsd(self, Frame frame):
        """Calculate dist_rmsd betwen two frames
        Parameters:
        ----------
        frame : Frame instance
        """
        return self.thisptr.DISTRMSD(frame.thisptr[0])

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

    # use Action_Strip here?
    # TODO : BROKEN
    cdef void _strip_atoms(Frame self, Topology top, string mask, bint update_top, bint has_box):
        """this method is too slow vs cpptraj
        if you use memory for numpy, you need to update after resizing Frame
        >>> arr0 = np.asarray(frame.buffer)
        >>> frame.strip_atoms(top,"!@CA")
        >>> # update view
        >>> arr0 = np.asarray(frame.buffer)
        """

        cdef Topology newtop = Topology()
        newtop.py_free_mem = False
        cdef AtomMask atm = AtomMask()
        cdef Frame tmpframe = Frame() 

        del tmpframe.thisptr

        tmpframe.thisptr = new _Frame(self.thisptr[0])
        
        atm.thisptr.SetMaskString(mask)
        atm.thisptr.InvertMask()
        top.thisptr.SetupIntegerMask(atm.thisptr[0])
        newtop.thisptr = top.thisptr.modifyStateByMask(atm.thisptr[0])
        if not has_box:
            newtop.thisptr.SetParmBox(_Box())
        tmpframe.thisptr.SetupFrameV(newtop.thisptr.Atoms(), newtop.thisptr.ParmCoordInfo())
        tmpframe.thisptr.SetFrame(self.thisptr[0], atm.thisptr[0])
        self.thisptr[0] = tmpframe.thisptr[0]
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
        if mask is None or top.is_empty():
            raise ValueError("need non-empty mask and non-empty Topology")
        cdef Frame frame

        mask = mask.encode("UTF-8")

        if not copy:
            self._strip_atoms(top, mask, update_top, has_box)
        else:
            frame = Frame(self)
            frame._strip_atoms(top, mask, update_top, has_box)
            return frame

    def __getbuffer__(self, Py_buffer* viewout, int flags):
        # NotImplementedYet
        # idea: use np.asarray(frame) rather using np.asarray(frame.buffer)
        pass

    def get_subframe(self, mask="", Topology top=Topology(), AtomMask atommask=AtomMask()):
        cdef AtomMask atm = AtomMask(mask)

        if top.is_empty() and atommask.n_selected == 0:
            raise ValueError("topology and atommask can not be both empty")
        elif not top.is_empty() and atommask.n_selected != 0:
            raise ValueError("topology and atommask can not both exist")
        else:
            if not top.is_empty():
                top.set_integer_mask(atm)
            else:
                atm = atommask
        return Frame(self, atm)

    def set_top(self, value):
        self.top = value

    def save(self, filename="", top=None, fmt='unknown', overwrite=False):
        if fmt == 'unknown':
            # convert to "UNKNOWN_TRAJ"
            fmt = fmt.upper() + "_TRAJ"
        with Trajout(filename=filename, top=top, fmt=fmt, 
                     overwrite=overwrite, more_args=None) as trajout:
            trajout.writeframe(0, self, top)
