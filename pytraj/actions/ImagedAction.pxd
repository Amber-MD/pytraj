# distutil: language = c++

from pytraj.DistRoutines cimport *
from pytraj.Box cimport _Box, Box, BoxType

cdef extern from "ImagedAction.h":
    cdef cppclass _ImagedAction "ImagedAction":
        _ImagedAction()
        void InitImaging(bint)
        void SetupImaging(BoxType)
        bint ImagingEnabled(self)
        bint UseImage()
        ImageType()

cdef class ImagedAction:
    cdef _ImagedAction* thisptr
