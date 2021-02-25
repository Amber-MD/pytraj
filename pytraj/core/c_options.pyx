# distutils: language = c++
from libcpp.string cimport string
from libcpp.vector cimport vector

cdef extern from "Version.h":
    cdef string CPPTRAJ_VERSION_STRING
    cdef string CPPTRAJ_INTERNAL_VERSION

cdef extern from "Cpptraj.h":
    cdef cppclass Cpptraj:
        @staticmethod
        string Defines()

__cpptraj_version__ = CPPTRAJ_VERSION_STRING.decode()
__cpptraj_internal_version__ = CPPTRAJ_INTERNAL_VERSION.decode()


def info():
    cdef string s

    s = Cpptraj.Defines()
    return s.decode()

# These are functions in cpptraj specially wrapped for pytraj
cdef extern from "ExternalFxn.h":
    void EXT_SetDefaultRng(int)

## set_default_rng
#
# This calls Cpptraj::Random_Number::SetDefaultRng() to change default
# random number generator. Primarily used for tests to set the RNG
# back to the old Marsaglia generator (rtype 0).
def set_default_rng(rtype):
    EXT_SetDefaultRng(rtype) 

cdef extern from "CpptrajStdio.h":
    void SuppressErrorMsg(bint)
    void SetWorldSilent(bint)


def set_error_silent(turnoff=True):
    SuppressErrorMsg(turnoff)


def set_world_silent(turnoff=True):
    SetWorldSilent(turnoff)

def set_cpptraj_verbose(cm=True):
    if cm:
        set_world_silent(False)
    else:
        set_world_silent(True)
