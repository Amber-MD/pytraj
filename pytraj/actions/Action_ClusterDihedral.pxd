# distutils: language = c++
#from libcpp.vector cimport vector
#from libcpp.string cimport string
from pytraj.actions.Action cimport _Action, Action, FunctPtr, _DispatchObject


cdef extern from "Action_ClusterDihedral.h": 
#    cdef cppclass _Action_ClusterDihedral::DCnode "Action_ClusterDihedral::DCnode" (_Action):
#        DCnode() 
#        DCnode(vector[int]& binIn, int frameIn)
#        DCnode(const DCnode& rhs)
#        #DCnode& operator =(const DCnode& rhs)
#        void Increment() 
#        void Add_Frame(int fIn)
#        bint operator[( const DCnode& rhs) const 
#        bint operator](const DCnode& rhs) const 
#        #bint operator = =(const DCnode& rhs) const 
#        bint BinMatch(vector[int]& binIn)
#        long int Count() 
#        bin_it binbegin() 
#        bin_it binend() 
#        frame_it framebegin() 
#        frame_it frameend() 
#        int Num_Frames() 


    cdef cppclass _Action_ClusterDihedral "Action_ClusterDihedral" (_Action):
        _Action_ClusterDihedral() 
        _DispatchObject * Alloc() 
        void Help() 


#    cdef cppclass _Action_ClusterDihedral::DCmask "Action_ClusterDihedral::DCmask" (_Action):
#        DCmask() 
#        DCmask(int a1, int a2, int a3, int a4, int bins, double min)
#        int A1() 
#        int A2() 
#        int A3() 
#        int A4() 
#        int Bins() 
#        double Step() 
#        double Min() 


cdef class Action_ClusterDihedral (Action):
    cdef _Action_ClusterDihedral* thisptr

