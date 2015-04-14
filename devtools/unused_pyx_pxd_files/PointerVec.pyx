# distutils: language = c++
from cython.operator cimport dereference as deref

cdef class TopVec:
    def __init__(self):
        pass

    def append(self, Topology top):
        self.ptrvec.push_back(top.thisptr)

    def __getitem__(self, int idx):
        cdef Topology top = Topology()
        top.thisptr[0] = deref(self.ptrvec[idx])
        return top

    def __iter__(self):
        cdef vector[_Topology*].iterator it
        cdef Topology top 

        it = self.ptrvec.begin()
        while it != self.ptrvec.end():
            top = Topology()
            # make a copy or a pointer?
            top.thisptr[0] = deref(deref(it))
            yield top

    @property
    def size(self):
        return self.ptrvec.size()


# Do we really need thi class while having FrameArray?
cdef class FrameVec:
    def __init__(self):
        pass

    def append(self, Frame frame):
        self.ptrvec.push_back(frame.thisptr)

    def __getitem__(self, int idx):
        cdef Frame frame = Frame()
        frame.thisptr[0] = deref(self.ptrvec[idx])
        return frame

    def __iter__(self):
        cdef vector[_Frame*].iterator it
        cdef Frame frame 

        it = self.ptrvec.begin()
        while it != self.ptrvec.end():
            frame = Frame()
            # make a frame copy or use pointer?
            frame.thisptr[0] = deref(deref(it))
            yield frame

    @property
    def size(self):
        return self.ptrvec.size()
