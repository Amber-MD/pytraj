# distutil: language = c++

cdef import from "ArrayIterator.h":
    cdef cppclass ArrayIterator[T]:
        ArrayIterator()
        ArrayIterator(const ArrayIterator&)
        ArrayIterator(T*)
        
        #Relations
        bint operator==(const ArrayIterator&)
        bint operator!=(const ArrayIterator&)

        #Increment
        ArrayIterator& operator++()
        ArrayIterator& operator++(int)

        #Value
        T& operator*()
        #not yet supported
        #T* operator->()
        #ArrayIterator& operator+=(int)
        ArrayIterator operator+(int)

#cdef class ArrayIterator:
#    cdef _ArrayIterator[T]* thisptr
