# distutils: language = c++
include "config.pxi"

cdef class NetcdfFile:
    def __cinit__(self):
        self.thisptr = new _NetcdfFile()

    def __dealloc__(self):
        del self.thisptr

    #def GetNetcdfConventions(self, *args):
    #    cdef char* fname
    #    if not args:
    #        return self.thisptr.GetNetcdfConventions()
    #    elif len(args) == 1:
    #        fname = args[0]
    #        return self.thisptr.GetNetcdfConventions(fname)
    #    else:
    #        raise NotImplemented()
    IF BINTRAJ:
        def GetNetcdfConventions(self, char* fname):
            return self.thisptr.GetNetcdfConventions(fname)

        def NetcdfDebug(self):
            self.thisptr.NetcdfDebug()

        def GetAttrText(self, char* attribute):
            return self.thisptr.GetAttrText(attribute)

        def NC_openRead(self, string name):
            return self.thisptr.NC_openRead(name)

        def NC_openWrite(self, string name):
            return self.thisptr.NC_openWrite(name)

        #def int NC_createReservoir(self,bint, double, int, int, int):

        #def int NC_create(self, string, NCTYPE, int, bint, bint, bint, bint, bint, bint, ReplicaDimArray, string):

        def NC_close(self):
            self.thisptr.NC_close()

        def SetupFrame(self):
            return self.thisptr.SetupFrame()

        def SetupCoordsVelo(self,bint b):
            return self.thisptr.SetupCoordsVelo(b)

        def SetupTime(self):
            return self.thisptr.SetupTime()

        #def int SetupBox(self,double *, NCTYPE):

        def SetupTemperature(self):
            return self.thisptr.SetupTemperature()

        #def int SetupMultiD(self,ReplicaDimArray):

        #def void FloatToDouble(self,double *, float *):

        #def void DoubleToFloat(self,float *, double *):
            pass
        
        @property
        def Ncid(self):
            return self.thisptr.Ncid()

        property Ncatom:
            def __get__(self):
                return self.thisptr.Ncatom()
            def __set__(self, natomIn):
                self.thisptr.SetNcatom(natomIn)

        @property
        def Ncatom3(self):
            return self.thisptr.Ncatom3()

        @property
        def Ncframe(self):
            return self.thisptr.Ncframe()

        def HasVelocities(self):
            return self.thisptr.HasVelocities()

        def HasCoords(self):
            return self.thisptr.HasCoords()

