# distutils: language = c++

cdef class CoordinateInfo:
    def __cinit__(self, dict crdinfo={}):
        cdef Box box
        cdef bint has_velocity, has_force, has_time
        cdef bint has_temperature

        self.thisptr = new _CoordinateInfo()

        box = crdinfo.get('box', Box())
        has_velocity = crdinfo.get('has_velocity', False)
        has_time = crdinfo.get('has_time', False)
        has_temperature = crdinfo.get('has_temperature', False)
        has_force = crdinfo.get('has_force', False)

        self.thisptr.SetBox(box.thisptr[0])
        self.thisptr.SetVelocity(has_velocity)
        self.thisptr.SetForce(has_force)
        self.thisptr.SetTime(has_time)
        self.thisptr.SetTemperature(has_temperature)

    def __dealloc__(self):
        if self.thisptr:
            del self.thisptr
