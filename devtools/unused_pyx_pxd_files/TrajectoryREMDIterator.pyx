# disutils: language = c++

# a trick to change to appropriate name for REMD
# why not using "TrajectoryIterator"? I don't know how.
cdef class TrajectoryREMDIterator(Trajin_Single):
    pass
