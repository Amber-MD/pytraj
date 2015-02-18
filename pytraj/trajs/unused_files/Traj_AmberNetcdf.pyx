# distutils: language = c++
include "config.pxi"

IF BINTRAJ:
    cdef class Traj_AmberNetcdf:
        def __cinit__(self):
            self.thisptr = new _Traj_AmberNetcdf()
    
        def __dealloc__(self):
            del self.thisptr
    
        def ID_TrajFormat(self, CpptrajFile fileIn):
            return self.thisptr.ID_TrajFormat(fileIn.thisptr[0])
    
        def setupTrajin(self, string filename, Topology trajParm):
            return self.thisptr.setupTrajin(filename, trajParm.thisptr)
    
        def setupTrajout(self, string filename, Topology trajParm, int NframesToWrite, bint append):
            return self.thisptr.setupTrajout(filename, trajParm.thisptr, NframesToWrite, append)
    
        def openTrajin(self):
            return self.thisptr.openTrajin()
    
        def closeTraj(self):
            self.thisptr.closeTraj()
    
        def readFrame(self,int myset, Frame frameIn):
            return self.thisptr.readFrame(myset, frameIn.thisptr[0])
    
        def readVelocity(self,int myset, Frame frameIn):
            return self.thisptr.readVelocity(myset, frameIn.thisptr[0])
    
        def writeFrame(self,int myset, Frame frameOut):
            return self.thisptr.writeFrame(myset, frameOut.thisptr[0])
    
        def Info(self):
            self.thisptr.Info()
    
        def processWriteArgs(self, ArgList arglist):
            return self.thisptr.processWriteArgs(arglist.thisptr[0])
    
        def processReadArgs(self,ArgList arglist):
            return self.thisptr.processReadArgs(arglist.thisptr[0])
    
        # not yet interested yet
        #def  int createReservoir(self,bint, double, int):
    
        #def writeReservoir(self,int, Frame, double, int):
    
