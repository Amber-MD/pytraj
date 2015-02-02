# distutils: language = c++
from pycpptraj.decorators import makesureABC


cdef class TrajectoryIO (BaseIOtype):
    def __cinit__(self):
        # _TrajectoryIO class has virtual method
        #  --> create instance of derived class and cast to _TrajectoryIO
        #self.baseptr_1 = <_TrajectoryIO*> new _TrajectoryFile()
        pass

    def __dealloc__(self):
        pass
        # __dealloc__ here?
        #if self.baseptr_1 != NULL:
        #    del self.baseptr_1

    #def virtual bint ID_TrajFormat(self,CpptrajFile):

    @makesureABC("TrajectoryIO")
    def total_frames(self, string filename, Topology top):
        return self.baseptr_1.setupTrajin(filename, top.thisptr)

    #def virtual int setupTrajout(self, string, Topology *, int, bint):

    #def virtual int openTrajin(self) :

    #def virtual int readFrame(self,int, Frame):

    #def virtual int readVelocity(self,int, Frame):

    #def virtual int writeFrame(self,int, Frame):

    #def virtual void closeTraj(self) :

    #def virtual void Info(self) :

    #def virtual int processWriteArgs(self,ArgList):

    #def virtual int processReadArgs(self,ArgList):

    @makesureABC("TrajectoryIO")
    def has_box(self):
        return self.baseptr_1.HasBox()

    @makesureABC("TrajectoryIO")
    def traj_box(self):
        cdef Box box = Box()
        box.thisptr[0] = self.baseptr_1.TrajBox()
        return box

    @makesureABC("TrajectoryIO")
    def has_velocity(self):
        return self.baseptr_1.HasV()

    @makesureABC("TrajectoryIO")
    def has_temperature(self):
        return self.baseptr_1.HasT()

    #@property
    #def title(self):
    #    return self.baseptr_1.Title()

    #title.setter
    #def title(self, string tIn):
    #    self.baseptr_1.SetTitle(tIn)

    @makesureABC("TrajectoryIO")
    def replica_dim(self):
        cdef ReplicaDimArray repdim = ReplicaDimArray()
        repdim.thisptr[0] = self.baseptr_1.ReplicaDimensions()
        return repdim

    @makesureABC("TrajectoryIO")
    def set_debug(self,int dIn):
        self.baseptr_1.SetDebug(dIn)

    @makesureABC("TrajectoryIO")
    def set_box(self, Box bIn):
        self.baseptr_1.SetBox(bIn.thisptr[0])

    @makesureABC("TrajectoryIO")
    def set_velocity(self,bint vIn):
        self.baseptr_1.SetVelocity(vIn)

    @makesureABC("TrajectoryIO")
    def set_temperature(self,bint tIn):
        self.baseptr_1.SetTemperature(tIn)

    @makesureABC("TrajectoryIO")
    def set_replica_dims(self, ReplicaDimArray rIn):
        self.baseptr_1.SetReplicaDims(rIn.thisptr[0])
