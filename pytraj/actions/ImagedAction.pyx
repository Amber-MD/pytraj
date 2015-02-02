# distutils: language = c++


cdef class ImagedAction:
    def __cinit__(self):
        self.thisptr = new _ImagedAction()

    def __dealloc__(self):
        del self.thisptr

    #def ImagedAction(self):

    #def void InitImaging(self,bint imageIn):

    #def void SetupImaging(self,Box::BoxType parmboxtype):

    #def bint ImagingEnabled(self):

    #def bint UseImage(self):

    #def ImagingType ImageType(self):

