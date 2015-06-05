from pytraj.misc import find_libcpptraj
from ctypes import cdll
import unittest

p = find_libcpptraj()

print (p)
plib = cdll.LoadLibrary(p[-1])
print (plib)
