# distutils: language = c++


cdef class Grid:
    def __cinit__(self):
        self.thisptr = new _Grid[float]()

    def __dealloc__(self):
        del self.thisptr

    #def Grid(self):

    #def Grid(self, Grid):

    #def T operator[](self,size_t idx):

    #def  T operator[](self,size_t idx):

    #def size_t size(self):

    #def int resize(self,size_t, size_t, size_t):

    #def size_t NX(self):

    #def size_t NY(self):

    #def size_t NZ(self):

    #def long int incrementBy(self,int, int, int, T):

    #def void setGrid(self,int, int, int, T):

    #def  T element(self,int, int, int):

    #def long int CalcIndex(self,int x, int y, int z):

    #def iterator begin(self):

    #def iterator end(self):

