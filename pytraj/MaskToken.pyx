# distutils: language = c++


cdef class MaskToken:
    def __cinit__(self):
        self.thisptr = new _MaskToken()

    def __dealloc__(self):
        del self.thisptr

    def TypeName(self):
        return self.thisptr.TypeName()

    def Print(self):
        self.thisptr.Print()

    def SetToken(self,MaskTokenType tktype, string s):
        return self.thisptr.SetToken(tktype, s)

    def SetDistance(self,string s):
        return self.thisptr.SetDistance(s)

    def SetOperator(self,MaskTokenType tktype):
        self.thisptr.SetOperator(tktype)

    def Type(self):
        return self.thisptr.Type()

    def Res1(self):
        return self.thisptr.Res1()

    def Res2(self):
        return self.thisptr.Res2()

    def Name(self):
        cdef NameType nt = NameType()
        del nt.thisptr
        nt.thisptr[0] = self.thisptr.Name()
        return nt 
    
    def OnStack(self): 
        return self.thisptr.OnStack()
    
    def Within(self): 
        return self.thisptr.Within()
            
    def ByAtom(self):
        return self.thisptr.ByAtom()

    def Distance(self):
        return self.thisptr.Distance()

    def SetOnStack(self):
        self.thisptr.SetOnStack()

