# distutils: language = c++

#do I need to wrap this ReadLine class? Python is very good for text processing
cdef extern from "ReadLine.h":
    char* duplicate_string(const char*)
    char* command_generator(const char*, int)
    char** cpptraj_completion(const char*, int, int)
    cdef cppclass _ReadLine "ReadLine":
        _ReadLine()
        int GetInput()
        bint YesNoPrompt(char *)
        const char* c_str()
        bint empty()
    char* duplicate_string(const char*)
    char* command_generator(const char*, int)
    char** cpptraj_completion(const char*, int, int)

cdef class ReadLine:
    cdef _ReadLine* thisptr
    def __cinit__(self):
        self.thisptr = new _ReadLine()

    def GetInput(self):
        self.thisptr.GetInput()

    def YesNoPrompt(self, mystring):
        return self.thisptr.YesNoPrompt(mystring)

    def c_str(self):
        return self.thisptr.c_str()

    def empty(self):
        return self.thisptr.empty()
