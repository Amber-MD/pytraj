# distutils: language = c++
# do we really need this routine?

from libcpp.string cimport string
from libcpp.vector cimport vector

cdef extern from "StringRoutines.h": 
    string tildeExpansion(const string&)
    ctypedef vector[string] StrArray
    StrArray ExpandToFilenames(const string&)
    bint fileExists(const string&)

    string NumberFilename(const string, int)
    int DigitWidth(long int)
    int FloatWidth(double)

    int convertToInteger(const string&)
    double convertToDouble(const string&)
    void RemoveTrailingWhitespace(string&)
    string integerToString(int)
    string integerToString(int, int)
    string doubleToString(double)

    bint validInteger(const string&)
    bint validDouble(const string&)

    string SetDoubleFormatString(int, int, int)
    string SetStringFormatString(int, bint)
    string SetIntegerFormatString(int)
