# distutils: language = c++


cdef class CoordinateInfo:
    def __cinit__(self):
        self.thisptr = new _CoordinateInfo()

    def __dealloc__(self):
        del self.thisptr

    #def CoordinateInfo(self, Box b, bint v, bint t, bint m):
    #def CoordinateInfo(self, ReplicaDimArray r, Box b, bint v, bint t, bint m, bint f):

    def has_box(self):
        return self.thisptr.HasBox()

    def traj_box(self):
        cdef Box box = Box()
        box.thisptr[0] = self.thisptr.TrajBox()
        return box

    def has_vel(self):
        return self.thisptr.HasVel()

    def has_temp(self):
        return self.thisptr.HasTemp()

    def has_time(self):
        return self.thisptr.HasTime()

    def has_force(self):
        return self.thisptr.HasForce()

    def has_replica_dims(self):
        return self.thisptr.HasReplicaDims()

    #def ReplicaDimArray ReplicaDimensions(self):

    #def void SetTime(self,bint m):

    #def void SetTemperature(self,bint t):

    #def void SetVelocity(self,bint v):

    #def void SetBox(self, Box b):

