# distutils: language = c++
from __future__ import absolute_import
from libcpp.map cimport map
from libcpp.vector cimport vector
from libcpp.string cimport string 
from ..math.cpp_math cimport _Grid, _Vec3, Vec3, _Matrix_3x3, Matrix_3x3, _Matrix
from .base cimport _DataSet, DataSet, DataType
from ..Frame cimport _Frame, Frame
from ..Topology cimport _Topology, Topology
from ..core.cpptraj_core cimport _ArgList, ArgList, _AtomMask, AtomMask


cdef extern from "DataSet_1D.h": 
    cdef cppclass _DataSet_1D "DataSet_1D" (_DataSet):
        _DataSet_1D() 
        _DataSet_1D(_DataSet)
        # virtual methods
        #virtual ~_DataSet_1D() 
        int Allocate1D(size_t)
        double Dval(size_t) const
        double Xcrd(size_t) const
        # end virtual methods
        inline bint IsTorsionArray() const 
        double Avg() const 
        double Avg(double& sd) const 
        double Min() const 
        double Max() const 
        int CrossCorr(const _DataSet_1D&, _DataSet_1D&, int, bint, bint) const 
        double CorrCoeff(const _DataSet_1D&) const 


cdef class DataSet_1D (DataSet):
    # baseptr0 is from DataSet
    cdef _DataSet_1D* baseptr_1
# distutils: language = c++
cdef extern from "DataSet_double.h": 
    cdef cppclass _DatasetDouble "DataSet_double" (_DataSet_1D):
        _DatasetDouble() 
        @staticmethod
        _DataSet * Alloc() 
        double& operator[](size_t idx)
        double& index_opr "operator[]"(size_t idx)
        const vector[double]& Data() const 
        void assign_opr "operator =" (const vector[double]& rhs)
        void AddElement(double d)
        void Resize(size_t sizeIn)
        size_t Size()
        int Sync() 
        void Info() const 
        int Allocate1D(size_t)
        void Add(size_t, const void *)
        double Dval(size_t idx) const 
        double Xcrd(size_t idx) const 
        void Append(const _DatasetDouble&)
        void SetNOE(double b, double bh, double r)
        double NOE_bound() const 
        double NOE_boundH() const 
        double NOE_rexp() const 
        void ShiftTorsions(double, double)


cdef class DatasetDouble (DataSet_1D):
    cdef _DatasetDouble* thisptr
    cdef bint py_free_mem 

cdef extern from "DataSet_float.h": 
    cdef cppclass _DatasetFloat "DataSet_float" (_DataSet_1D):
        _DatasetFloat() 
        @staticmethod
        _DataSet * Alloc() 
        float& operator[](size_t idx)
        float& index_opr "operator[]"(size_t idx)
        int Size()
        void Resize(size_t)

cdef class DatasetFloat (DataSet_1D):
    cdef _DatasetFloat* thisptr
    cdef bint py_free_mem 

cdef extern from "DataSet_integer.h": 
    cdef cppclass _DatasetInteger "DataSet_integer" (_DataSet_1D):
        _DatasetInteger() 
        @staticmethod
        _DataSet * Alloc() 
        int& operator[](size_t idx)
        int& index_opr "operator[]"(size_t idx)
        void AddElement(int i)
        int Size()
        void Resize(size_t)
        void Add( size_t, const void* )

cdef class DatasetInteger (DataSet_1D):
    cdef _DatasetInteger* thisptr
    cdef bint py_free_mem 

cdef extern from "DataSet_string.h": 
    cdef cppclass _DatasetString "DataSet_string" (_DataSet_1D):
        _DatasetString()
        _DataSet * Alloc() 
        string& index_opr "operator[]"(size_t idx)
        void AddElement(const string& s)
        void Resize(size_t sizeIn)
        int Size()

cdef class DatasetString(DataSet_1D):
    cdef _DatasetString* thisptr
    cdef bint py_free_mem

cdef extern from "DataSet_Vector.h": 
    cdef cppclass _DatasetVector "DataSet_Vector" (_DataSet_1D):
        _DatasetVector() 
        _DataSet * Alloc() 
        void SetIred() 
        bint IsIred() const 
        void reset() 
        void Resize(size_t s)
        void Resize(size_t s, const _Vec3& v)
        bint Empty() const 
        #const _Vec3& operator[](int i) const 
        _Vec3& index_opr "operator[]" (int i)
        const _Vec3& OXYZ(int i) const 
        void ReserveVecs(size_t n)
        void AddVxyz(const _Vec3& v)
        void AddVxyz(const _Vec3& v, const _Vec3& c)
        #const_iterator begin() const 
        #const_iterator end() const 
        const _Vec3& Back() const 
        int CalcSphericalHarmonics(int)
        #const _ComplexArray& SphericalHarmonics(int) const 
        double SphericalHarmonicsNorm(int)


cdef class DatasetVector (DataSet_1D):
    cdef _DatasetVector* thisptr
    cdef bint py_free_mem


cdef extern from "DataSet_2D.h": 
    # DataSet_2D.h
    ctypedef enum MatrixType "DataSet_2D::MatrixType":
        pass
    ctypedef enum MatrixKind "DataSet_2D::MatrixKind":
        pass
    cdef cppclass _DataSet_2D "DataSet_2D" (_DataSet):
        _DataSet_2D() 
        _DataSet_2D(DataType tIn, int wIn, int pIn)
        # virtual methods
        int Allocate2D(size_t, size_t) 
        int AllocateHalf(size_t) 
        int AllocateTriangle(size_t) 
        double GetElement(size_t, size_t) const  
        size_t Nrows() const  
        size_t Ncols() const  
        double * MatrixArray() const  
        MatrixKind Kind "MatrixKind"() const  
        # end virtual methods

        void Add(size_t, const void *)
        const char * MatrixTypeString(MatrixType m)
        const char * MatrixOutputString(MatrixType m)

cdef class DataSet_2D (DataSet):
    cdef _DataSet_2D* baseptr_1


#ctypedef Matrix[double].iterator iterator
ctypedef vector[double] Darray

cdef extern from "DataSet_MatrixDbl.h": 
    cdef cppclass _DatasetMatrixDouble "DataSet_MatrixDbl" (_DataSet_2D):
        _DatasetMatrixDouble() 
        double& index_opr "operator[]"(size_t idx)
        @staticmethod
        _DataSet * Alloc() 
        size_t Size() const 
        int Sync() 
        void Info() const 
        int Allocate2D(size_t x, size_t y)
        int AllocateHalf(size_t x)
        int AllocateTriangle(size_t x)
        #void Write2D(_CpptrajFile&, int, int) const 
        double GetElement(size_t x, size_t y) const 
        size_t Nrows() const 
        size_t Ncols() const 
        #double * MatrixArray() const # not implemented
        MatrixKind Kind() const 
        # make alias to avoid naming conflict with DataSet (DataType)
        MatrixType matType "Type"() const 
        unsigned int Nsnapshots() const 
        void IncrementSnapshots() 
        double& Element(size_t x, size_t y)
        int AddElement(double d)
        void SetElement(size_t x, size_t y, double d)
        #iterator begin() 
        #iterator end() 
        const Darray& Vect() const 
        Darray& V1() 
        void AllocateVector(size_t vsize)
        #Darray.iterator v1begin() 
        #Darray.iterator v1end() 
        void SetTypeAndKind(MatrixType tIn, MatrixKind kIn)
        void StoreMass(const Darray& mIn)
        const Darray& Mass() const 


cdef class DatasetMatrixDouble (DataSet_2D):
    cdef _DatasetMatrixDouble* thisptr
    cdef bint py_free_mem


cdef extern from "DataSet_MatrixFlt.h": 
    cdef cppclass _DatasetMatrixFloat  "DataSet_MatrixFlt" (_DataSet_2D):
        _DataSet_MatrixFlt() 
        float& index_opr "operator[]" (size_t idx)
        @staticmethod
        _DataSet * Alloc() 


cdef class DatasetMatrixFloat(DataSet_2D):
    cdef _DatasetMatrixFloat * thisptr
    cdef bint py_free_mem


cdef extern from "DataSet_3D.h": 
    cdef cppclass _DataSet_3D "DataSet_3D" (_DataSet):
        _DataSet_3D() 
        _DataSet_3D(DataType tIn, int wIn, int pIn)
        void Add(size_t, const void *)
        int Allocate_N_O_D(size_t, size_t, size_t, const _Vec3&, const _Vec3&)
        int Allocate_N_C_D(size_t, size_t, size_t, const _Vec3&, const _Vec3&)
        int Allocate_X_C_D(const _Vec3&, const _Vec3&, const _Vec3&)
        inline bint CalcBins(double, double, double, int&, int&, int&) const 
        inline double DX() const 
        inline double DY() const 
        inline double DZ() const 
        inline double OX() const 
        inline double OY() const 
        inline double OZ() const 
        inline double MX() const 
        inline double MY() const 
        inline double MZ() const 
        inline _Vec3 BinCorner(int, int, int)
        inline _Vec3 BinCenter(int, int, int)

cdef class DataSet_3D (DataSet):
    cdef _DataSet_3D* baseptr_1

cdef extern from "DataSet_GridFlt.h": 
    cdef cppclass _DatasetGridFloat "DataSet_GridFlt" (_DataSet_3D):
        _DatasetGridFloat()
        float& index_opr "operator[]"(size_t idx)
        _DataSet * Alloc() 
        const _Grid[float]& InternalGrid() const 
        size_t Size() const 
        int Sync() 
        void Info() const 
        int Allocate3D(size_t x, size_t y, size_t z)
        double GetElement(int x, int y, int z) const 
        void SetElement(int x, int y, int z, float v)
        double operator[](size_t idx) const 
        size_t NX() const 
        size_t NY() const 
        size_t NZ() const 
        #iterator begin() 
        #iterator end() 
        inline long int Increment(const _Vec3&, float)
        inline long int Increment(const double *, float)
        inline long int Increment(int, int, int, float)
        float GridVal(int x, int y, int z) const 
        long int CalcIndex(int i, int j, int k) const 

cdef class DatasetGridFloat (DataSet_3D):
    cdef _DatasetGridFloat* thisptr
    cdef public bint py_free_mem

cdef extern from "DataSet_Modes.h": 
    cdef cppclass _DataSet_Modes "DataSet_Modes" (_DataSet):
        _DataSet_Modes() 
        _DataSet * Alloc() 
        size_t Size() const 
        int Sync() 
        void Info() const 
        void Add(size_t, const void *)
        const Darray& AvgCrd() const 
        const Darray& Mass() const 
        int NavgCrd() const 
        double * Avg_FramePtr() 
        const double * Avg_FramePtr() const 
        void AllocateAvgCoords(int n)
        void SetAvgCoords(const _DataSet_2D&)
        int SetModes(bint, int, int, const double *, const double *)
        int CalcEigen(const _DataSet_2D&, int)
        void PrintModes() 
        int EigvalToFreq(double)
        int MassWtEigvect(Darray&)
        int ReduceVectors() 
        double Eigenvalue(int i) const 
        const double * Eigenvectors() const 
        const double * Eigenvector(int i) const 
        int Nmodes() const 
        int VectorSize() const 
        #MatrixType Type() const 
        bint IsReduced() const 

cdef class DataSet_Modes (DataSet):
    cdef _DataSet_Modes* thisptr

ctypedef map[double, int] TmapType
cdef extern from "DataSet_RemLog.h": 
    cdef cppclass _DataSet_RemLog "DataSet_RemLog":
        _DataSet_RemLog() 
        _DataSet * Alloc() 
        void AllocateReplicas(int)
        void AddRepFrame(int rep, const _ReplicaFrame& frm)
        const _ReplicaFrame& RepFrame(int exch, int rep) const 
        int NumExchange() const 
        bint ValidEnsemble() const 
        void TrimLastExchange() 
        size_t Size() const 
        int Sync() 
        void Info() const 
        void Add(size_t, const void *)


    cdef cppclass _ReplicaFrame "DataSet_RemLog::ReplicaFrame":
        _Replica_Frame() 
        int SetTremdFrame(const char *, const TmapType&)
        int SetHremdFrame(const char *, const vector[int]&)
        int ReplicaIdx() const 
        int PartnerIdx() const 
        int CoordsIdx() const 
        bint Success() const 
        double Temp0() const 
        double PE_X1() const 
        double PE_X2() const 


cdef class DataSet_RemLog:
    cdef _DataSet_RemLog* thisptr

cdef class ReplicaFrame:
    cdef _ReplicaFrame* thisptr



cdef extern from "DataSet_Mat3x3.h": 
    ctypedef vector[_Matrix_3x3].iterator mat_iterator
    cdef cppclass _DatasetMatrix3x3 "DataSet_Mat3x3" (_DataSet_1D):
        _DatasetMatrix3x3()
        @staticmethod
        _DataSet * Alloc() 
        bint Empty()
        void AddMat3x3(_Matrix_3x3)
        mat_iterator begin()
        mat_iterator end()
        _Matrix_3x3& operator[](int i)


cdef class DatasetMatrix3x3(DataSet_1D):
    cdef _DatasetMatrix3x3* thisptr
    cdef bint py_free_mem 


cdef extern from "DataSet_Mesh.h": 
    cdef cppclass _DatasetMesh "DataSet_Mesh" (_DataSet_1D):
        _DatasetMesh()
        _DatasetMesh(int, double, double)
        _DataSet * Alloc() 
        size_t Size() const 
        int Sync() 
        void Info() const 
        int Allocate1D(size_t)
        void Add(size_t, const void *)
        double Dval(size_t idx) const 
        double Xcrd(size_t idx) const 
        inline void AddXY(double, double)
        double X(int i) const 
        double Y(int i) const 
        void CalculateMeshX(int, double, double)
        int SetMeshXY(const _DataSet_1D&)
        double Integrate_Trapezoid(_DatasetMesh&) const 
        double Integrate_Trapezoid() const 
        int SetSplinedMeshY(const vector[double]&, const vector[double]&)
        int SetSplinedMesh(const _DataSet_1D&)
        int LinearRegression(double&, double&, double&, bint) const 

cdef class DatasetMesh(DataSet_1D):
    cdef _DatasetMesh* thisptr
    cdef public bint py_free_mem


cdef extern from "DataSet_Coords.h": 
    cdef cppclass _DataSet_Coords "DataSet_Coords" (_DataSet):
        _DataSet_Coords() 
        _DataSet_Coords(DataType)
        #virtual ~_DataSet_Coords() 
        _Frame AllocateFrame() const 
        
        # virtual methods
        void AddFrame(const _Frame&) 
        void SetCRD(int, const _Frame&) 
        void GetFrame(int, _Frame&) 
        void GetFrame(int, _Frame&, const _AtomMask&) 
        # end virtual methods

        void SetTopology(const _Topology&)
        inline const _Topology& Top() const 


cdef class DataSet_Coords (DataSet):
    # DataSet has baseptr0
    cdef _DataSet_Coords* baseptr_1
    cdef Topology _top
    cdef bint py_free_mem

    # use tmpfarray object to hold Frame or Trajectory 
    # (if we want to use dset[0][0] correctly)
    cdef object tmpfarray
# distutils: language = c++

cdef extern from "DataSet_Coords_CRD.h": 
    cdef cppclass _DataSet_Coords_CRD "DataSet_Coords_CRD" (_DataSet_Coords):
        _DataSet_Coords_CRD() 
        @staticmethod
        _DataSet * Alloc() 
        size_t Size() const 
        int Sync() 
        void Info() const 
        int Allocate1D(size_t)
        void Add(size_t, const void *)
        double Dval(size_t)const 
        double Xcrd(size_t idx)const 
        inline void AddFrame(const _Frame& fIn)
        inline void GetFrame(int idx, _Frame & fIn)
        inline void GetFrame(int idx, _Frame & fIn, const _AtomMask& mIn)
        inline void SetCRD(int idx, const _Frame& fIn)


cdef class DataSet_Coords_CRD (DataSet_Coords):
    cdef _DataSet_Coords_CRD* thisptr
cdef extern from "DataSet_Coords_REF.h": 
    cdef cppclass _DataSet_Coords_REF "DataSet_Coords_REF" (_DataSet_Coords):
        _DataSet_Coords_REF() 

        # turn off those methods since they are in parent class
        @staticmethod
        _DataSet * Alloc() 
        size_t Size() const 
        #int Sync() 
        #void Info() const 
        #void Add(size_t, const void *)
        #int AllocateCoords(size_t)
        #inline void AddFrame(const _Frame& fIn)
        #inline void GetFrame(int idx, _Frame& fIn)
        #inline void GetFrame(int idx, _Frame& fIn, const _AtomMask& mIn)
        #inline void SetCRD(int idx, const _Frame& fIn)

        int LoadRef(const string&, const _Topology&, int)
        int SetupRef_Frame(const string&, const string&, const _Topology&, _ArgList&, int)
        int SetupRef_Frame(_DataSet_Coords *, const string&, int, int)
        int StripRef(const string&)
        int StripRef(const _AtomMask&)
        const _Frame& RefFrame() const 
        int RefIndex() const 

cdef class DataSet_Coords_REF (DataSet_Coords):
    cdef _DataSet_Coords_REF* thisptr
