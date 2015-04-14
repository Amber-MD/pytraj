# distutils: language = c++


cdef class ComplexArray:
    def __cinit__(self):
        self.thisptr = new _ComplexArray()

    def __dealloc__(self):
        del self.thisptr

    #def ComplexArray(self):

    #def ComplexArray(self,int):

    #def ComplexArray(self, ComplexArray):

    #def void Allocate(self,int):

    #def void Assign(self, ComplexArray):

    #def void PadWithZero(self,int):

    #def void Normalize(self,double):

    #def void SquareModulus(self):

    #def void ComplexConjTimes(self, ComplexArray):

    #def double * CAptr(self):

    #def int size(self):

    #def double operator[](self,int idx):

    #def  double operator[](self,int idx):

    #def  iterator begin(self):

    #def  iterator end(self):

