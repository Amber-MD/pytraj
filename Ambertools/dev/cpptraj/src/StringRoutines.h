#ifndef INC_STRINGROUTINES_H
#define INC_STRINGROUTINES_H
#include <string>
#include <vector>
/*! \file StringRoutines.h
    \brief Collection of useful string routines.

    Any commonly used routines that make use of the string class. This
    includes converting to/from numbers, string modification, setting
    printf-type format strings, etc.
 */
std::string tildeExpansion(std::string const&);
typedef std::vector<std::string> StrArray;
StrArray ExpandToFilenames(std::string const&);
bool fileExists(std::string const&);

std::string NumberFilename(std::string const &, int);
int DigitWidth(long int);
int FloatWidth(double);

int convertToInteger(std::string const &);
double convertToDouble(std::string const &);
void RemoveTrailingWhitespace(std::string &);
std::string NoTrailingWhitespace(std::string const&);
std::string integerToString(int);
std::string integerToString(int,int);
std::string doubleToString(double);
/// Brief check that the passed in string begins with a digit or '-'
bool validInteger(std::string const&);
/// Brief check that the passed in string begins with a digit, '-', or '.'
bool validDouble(std::string const&);

std::string SetDoubleFormatString(int, int, int);
std::string SetStringFormatString(int, bool);
std::string SetIntegerFormatString(int);
/// \return the current date/time with format 'mm/dd/yy  hh:mm:ss'
std::string TimeString();
// NOTE: Not really string routines, but here for convenience.
long long AvailableMemory();
double AvailableMemory_MB();
# ifdef __APPLE__
long long TotalGlobalMemory();
# endif
#endif
