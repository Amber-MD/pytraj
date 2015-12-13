# distutils: language = c++
from __future__ import absolute_import
import numpy as np
from cython cimport view
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as incr

from pytraj.c_dict import BoxTypeDict, get_key
from pytraj.externals.six import string_types


cdef class Box(object):
    def __cinit__(self, *args):
        cdef double[:] boxIn
        cdef Box rhs
        cdef Matrix_3x3 mat

        if not args:
            self.thisptr = new _Box()
        elif len(args) == 1:
            if isinstance(args[0], Box):
                rhs = args[0]
                self.thisptr = new _Box(rhs.thisptr[0])
            elif isinstance(args[0], Matrix_3x3):
                mat = args[0]
                self.thisptr = new _Box(mat.thisptr[0])
            elif isinstance(args[0], string_types):
                # set box based on type
                self.thisptr = new _Box()
                self.type = args[0]
            else:
                try:
                    # if args[0] has buffer interface
                    boxIn = args[0]
                except:
                    boxIn = np.array([item for item in args[0]], dtype='f8')
                self.thisptr = new _Box( &boxIn[0])
        else:
            raise ValueError("")

    def copy(self):
        cdef Box other = Box()
        other.thisptr = new _Box(self.thisptr[0])
        return other

    def __str__(self):
        boxlisttxt = ", ".join([str(tmp) for tmp in self.tolist()])
        boxlisttxt = "(" + boxlisttxt + ")"
        txt = "<Box: %s, (x, y, z, alpha, beta, gamma) = %s>" % (self.type, boxlisttxt)
        return txt

    def __repr__(self):
        return self.__str__()

    property data:
        def __get__(self):
            """memoryview for box array"""
            return <double[:6]> self.thisptr.boxPtr()

    def _get_data(self):
        """memoryview for box array"""
        return <double[:6]> self.thisptr.boxPtr()

    def __dealloc__(self):
        del self.thisptr

    def __getitem__(self, idx):
        """add fancy indexing?"""
        return self._get_data()[idx]

    def __setitem__(self, idx, value):
        self.data[idx] = value

    def __iter__(self):
        for x in self._get_data():
            yield x

    @classmethod
    def all_box_types(cls):
        return [x.lower() for x in BoxTypeDict.keys()]

    def _update_box_type(self):
        """trick to let cpptraj correctly set boxtype
        Example
        -------
        >>> from pytraj.core import Box
        >>> box = Box()
        >>> box.alpha, box.beta, box.gamma = 90., 90., 90.
        >>> print (box.type)
        nobox
        >>> box._update_box_type()
        >>> print (box.type)
        ortho
        """
        cdef double[:] data = self.data

        self.thisptr.SetBox( &data[0])

    def set_beta_lengths(self, double beta, double xin, double yin, double zin):
        self.thisptr.SetBetaLengths(beta, xin, yin, zin)

    def set_nobox(self):
        self.thisptr.SetNoBox()

    def to_recip(self):
        cdef Matrix_3x3 ucell = Matrix_3x3()
        cdef Matrix_3x3 recip = Matrix_3x3()
        self.thisptr.ToRecip(ucell.thisptr[0], recip.thisptr[0])
        return np.array(ucell), np.array(recip)

    property type:
        def __get__(self):
            if self.tolist()[3:] == [60., 90., 60]:
                # cpptraj does not handle this correctly
                # so we introduce our own version
                return 'rhombic'
            else:
                return get_key(self.thisptr.Type(), BoxTypeDict).lower()

        def __set__(self, value):
            all_box_types = self.all_box_types()
            value = value.lower()
            if value == 'ortho':
                self.alpha, self.beta, self.gamma = 90., 90., 90.
            elif value == 'truncoct':
                # use cpptraj' method
                self.thisptr.SetTruncOct()
            elif value == 'rhombic':
                # check cpptraj' code to know why
                self.alpha, self.beta, self.gamma = 0., 60., 0.
            elif value == 'nobox':
                self.set_nobox()
            else:
                msg = "supported boxtype is ortho | truncoct | rhombic | nobox\n"
                msg2 = """use box.alpha, box.beta, box.gamma to explicitly assign values
                          and use `_update_box_type() method`"""
                raise ValueError(msg + msg2)
            # need to update all info so cpptraj will `SetBoxType` (private method)
            # sounds dummy to set your box to yourself to do this trick :D
            # should update cpptraj code
            self._update_box_type()

    property x:
        def __get__(self):
            return self.thisptr.BoxX()

        def __set__(self, double value):
            self.thisptr.SetX(value)

    property y:
        def __get__(self):
            return self.thisptr.BoxY()

        def __set__(self, double value):
            self.thisptr.SetY(value)

    property z:
        def __get__(self):
            return self.thisptr.BoxZ()

        def __set__(self, double value):
            self.thisptr.SetZ(value)

    property alpha:
        def __get__(self):
            return self.thisptr.Alpha()

        def __set__(self, double value):
            self.thisptr.SetAlpha(value)

    property beta:
        def __get__(self):
            return self.thisptr.Beta()

        def __set__(self, double value):
            self.thisptr.SetBeta(value)

    property gamma:
        def __get__(self):
            return self.thisptr.Gamma()

        def __set__(self, double value):
            self.thisptr.SetGamma(value)

    def has_box(self):
        return self.thisptr.HasBox()

    @property
    def center(self):
        cdef Vec3 vec = Vec3()
        vec.thisptr[0] = self.thisptr.Center()
        return np.array(vec.values)

    @property
    def lengths(self):
        cdef Vec3 vec = Vec3()
        vec.thisptr[0] = self.thisptr.Lengths()
        return vec

    def tolist(self):
        return list([x for x in self.data])

    @property
    def values(self):
        """return a copy of Box's data: [x, y, z, alpha, beta, gamma]"""
        return np.array(self._get_data()[:])
