# distutils: language = c++

from libcpp.string cimport string
from Trajin cimport *

cdef cppclass Trajin2(_Trajin):
    int SetupTrajRead(const string& test, _ArgList& argin, _Topology *top):
        return 0
    int ReadTrajFrame(int idx, _Frame& frame):
        return 0
    int BeginTraj(bint b):
        return 0
    void EndTraj():
        print "None"
    void PrintInfo(int id):
        print "None"
    bint HasVelocity():
        return False
    int NreplicaDimension():
        return 0

#cdef Trajin2* ptr = new Trajin2()
cdef Trajin2 ptr 

ptr.BeginTraj(False)
ptr.EndTraj()
ptr.PrintInfo(12)
print ptr.HasVelocity()
print ptr.NreplicaDimension()
