# distutils: language = c++
from libcpp.vector cimport vector
from libcpp.string cimport string
from pytraj.Atom cimport _Atom, Atom
from pytraj.AtomMask cimport _AtomMask, AtomMask
from pytraj.Box cimport _Box, Box, BoxType
from pytraj.Topology cimport _Topology, Topology
from pytraj.Vec3 cimport _Vec3, Vec3
from pytraj.Matrix_3x3 cimport _Matrix_3x3, Matrix_3x3
from pytraj.CoordinateInfo cimport _CoordinateInfo, CoordinateInfo

ctypedef vector[float] CRDtype
ctypedef vector[double] Darray

cdef extern from "Frame.h": 
    # Frame.h
    ctypedef enum CenterMode "Frame::CenterMode":
        ORIGIN "Frame::ORIGIN"
        BOXCTR "Frame::BOXCTR"
        POINT "Frame::POINT"
    cdef cppclass _Frame "Frame":
        _Frame() 
        #~_Frame() 
        _Frame(int)
        _Frame(const vector[_Atom]&)
        _Frame(const _Frame&, const _AtomMask&)
        _Frame(const _Frame&)
        #_Frame& operator =(_Frame)
        void SetFromCRD(const CRDtype&, int, int, bint)
        void SetFromCRD(const CRDtype&, const _AtomMask&, int, int, bint)
        CRDtype ConvertToCRD(int, bint) const 
        void printAtomCoord(int) const 
        void Info(const char *) const 
        void ClearAtoms() 
        void AddXYZ(const double *)
        void AddVec3(const _Vec3&)
        void SwapAtoms(int, int)
        double& index_opr "operator[]"(int idx)
        const double& operator[](int idx) const 
        bint empty() const 
        bint HasVelocity() const 
        int Natom() const 
        int size() const 
        int NrepDims() const 
        double Temperature() const 
        const double * XYZ(int atnum) const 
        const double * CRD(int idx) const 
        const double * VXYZ(int atnum) const 
        double Mass(int atnum) const 
        const _Box& BoxCrd() const 
        inline double * xAddress() 
        inline double * vAddress() 
        inline double * bAddress() 
        inline double * tAddress() 
        inline int * iAddress() 
        #inline const double * xAddress() const 
        #inline const double * vAddress() const 
        #inline const double * bAddress() const 
        #inline const double * tAddress() const 
        inline const int * iAddress() const 
        inline void SetBoxAngles(const double *)
        int SetupFrame(int)
        int SetupFrameM(const vector[_Atom]&)
        int SetupFrameXM(const vector[double]&, const vector[double]&)
        int SetupFrameV(const vector[_Atom]&, const _CoordinateInfo&)
        int SetupFrameFromMask(const _AtomMask&, const vector[_Atom]&)
        void SetCoordinates(const _Frame&, const _AtomMask&)
        void SetCoordinates(const _Frame&)
        void SetFrame(const _Frame&, const _AtomMask&)
        void SetCoordinatesByMap(const _Frame&, const vector[int]&)
        void StripUnmapped_Atoms(const _Frame&, const vector[int]&)
        void ModifyByMap(const _Frame&, const vector[int]&)
        void ZeroCoords() 
        _Frame& addequal "operator +=" (const _Frame&)
        #_Frame& operator +=(const _Frame&)
        _Frame& subequal "operator -=" (const _Frame&)
        #_Frame& operator -= (const _Frame&)
        _Frame& mulequal "operator *=" (const _Frame&)
        #_Frame& operator *= (const _Frame&)
        const _Frame operator *(const _Frame&) const 
        const _Frame operator -(const _Frame&) const 
        int Divide(const _Frame&, double)
        void Divide(double)
        int AddByMask(const _Frame&, const _AtomMask&)
        inline bint CheckCoordsInvalid() const 
        inline _Vec3 VCenterOfMass(const _AtomMask&) const 
        inline _Vec3 VGeometricCenter(const _AtomMask&) const 
        inline _Vec3 VCenterOfMass(int, int) const 
        inline _Vec3 VGeometricCenter(int, int) const 
        inline void Translate(const _Vec3&, int, int)
        inline void Translate(const _Vec3&, int)
        inline void Translate(const _Vec3&)
        inline void NegTranslate(const _Vec3&)
        inline void Rotate(const _Matrix_3x3&)
        inline void Rotate(const _Matrix_3x3&, const _AtomMask&)
        inline void Trans_Rot_Trans(const _Vec3&, const _Matrix_3x3&, const _Vec3&)
        void Scale(const _AtomMask&, double, double, double)
        # Not in cpptraj anymore
        void Center(const _AtomMask&, CenterMode, const _Vec3&, bint)

        _Vec3 CenterOnOrigin(bint)
        double RMSD(_Frame&, bint)
        double RMSD(_Frame&, _Matrix_3x3&, _Vec3&, _Vec3&, bint)
        double RMSD_CenteredRef(const _Frame&, bint)
        double RMSD_CenteredRef(const _Frame&, _Matrix_3x3&, _Vec3&, bint)
        double RMSD_NoFit(const _Frame&, bint) const 
        double RMSD_FitToRef(const _Frame&, const _Vec3&)
        double DISTRMSD(const _Frame&) const 
        _Vec3 SetAxisOfRotation(int, int)
        _Vec3 CalculateInertia(const _AtomMask&, _Matrix_3x3&) const 
        double CalcTemperature(const _AtomMask&, int) const 


cdef class Frame:
    cdef _Frame* thisptr
    cdef public bint py_free_mem
    cdef void _strip_atoms(Frame self, Topology top, string m, bint update_top, bint has_box)
    cdef _update_atoms(self, int[:], double[:], int)
    # create and object as alias to Topology instance
    cdef object top

cdef inline int get_positive_idx(idx, size):
    # TODO : do we need this method?
    # we can we memoryview to get slicing too
    """Used for negative indexing"""
    if idx < 0:
        idx = size + idx
        if idx < 0:
            raise ValueError("index is out of range")
    if idx >= size:
        raise ValueError("index is out of range")
    return idx
