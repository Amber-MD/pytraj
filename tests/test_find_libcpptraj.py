import unittest
from pytraj.misc import find_libcpptraj, find_library
from ctypes import cdll
import unittest

p = find_libcpptraj()

#print(p)
plib = cdll.LoadLibrary(p[-1])
#print(plib)
#print(find_libcpptraj(unique=True))
#print(find_library('netcdf', unique=True))
