#include <cstdio>    // fprintf, fopen, fclose, sprintf
#include <cmath>     // log10
#include <ctime>     // for TimeString()
#include <cerrno>
#include <cstring>   // strerror
#include <sstream>   // istringstream, ostringstream
#include <locale>    // isspace
#include <stdexcept> // BadConversion
#include <vector>
#ifndef __PGI
#  include <glob.h>  // For tilde expansion
#endif
#include "StringRoutines.h"
#include "CpptrajStdio.h"
#ifdef _MSC_VER
# include <windows.h>
#else
# include <unistd.h>
#endif
#ifdef __APPLE__
# include <sys/sysctl.h>
# include <mach/mach_host.h>
#endif

// tildeExpansion()
/** Use glob.h to perform tilde expansion on a filename, returning the
  * expanded filename. If the file does not exist or globbing fails return an
  * empty string. Do not print an error message if the file does not exist
  * so that this routine and fileExists() can be used to silently check files.
  */
std::string tildeExpansion(std::string const& filenameIn) {
  if (filenameIn.empty()) {
    mprinterr("Error: tildeExpansion: null filename specified.\n");
    return std::string("");
  }
# ifdef __PGI
  // NOTE: It seems some PGI compilers do not function correctly when glob.h
  //       is included and large file flags are set. Just disable globbing
  //       for PGI and return a copy of filenameIn.
  // Since not globbing, check if file exists before returning filename.
  FILE *infile = fopen(filenameIn.c_str(), "rb");
  if (infile == 0) return std::string("");
  fclose(infile);
  return filenameIn;
# else
  glob_t globbuf;
  globbuf.gl_offs = 1;
  std::string returnFilename;
  int err = glob(filenameIn.c_str(), GLOB_TILDE, NULL, &globbuf);
  if ( err == GLOB_NOMATCH )
    //mprinterr("Error: '%s' does not exist.\n", filenameIn); // Make silent
    return returnFilename;
  else if ( err != 0 )
    mprinterr("Internal Error: glob() failed for '%s' (%i)\n", filenameIn.c_str(), err);
  else {
    returnFilename.assign( globbuf.gl_pathv[0] );
    globfree(&globbuf);
  }
  return returnFilename;
# endif
}

// ExpandToFilenames()
StrArray ExpandToFilenames(std::string const& fnameArg) {
  StrArray fnames;
  if (fnameArg.empty()) return fnames;
# ifdef __PGI
  // NOTE: It seems some PGI compilers do not function correctly when glob.h
  //       is included and large file flags are set. Just disable globbing
  //       for PGI and return a copy of filenameIn.
  // Check for any wildcards in fnameArg
  if ( fnameArg.find_first_of("*?[]") != std::string::npos )
    fprintf(stdout,"Warning: Currently wildcards in filenames not supported with PGI compilers.\n");
  fnames.push_back( fnameArg );
# else
  glob_t globbuf;
  int err = glob(fnameArg.c_str(), GLOB_TILDE, NULL, &globbuf );
  //printf("DEBUG: %s matches %zu files.\n", fnameArg.c_str(), (size_t)globbuf.gl_pathc);
  if ( err == 0 ) { 
    for (unsigned int i = 0; i < (size_t)globbuf.gl_pathc; i++)
      fnames.push_back( globbuf.gl_pathv[i] );
  } else if (err == GLOB_NOMATCH )
    fprintf(stderr,"Error: %s matches no files.\n", fnameArg.c_str());
  else
    fprintf(stderr,"Error: occurred trying to find %s\n", fnameArg.c_str());
  if ( globbuf.gl_pathc > 0 ) globfree(&globbuf);
# endif
  return fnames;
}

// fileExists()
/** Return true if file can be opened "r".  */
bool fileExists(std::string const& filenameIn) {
  // Perform tilde expansion
  std::string fname = tildeExpansion(filenameIn);
  if (fname.empty()) return false;
  FILE *infile = fopen(fname.c_str(), "rb");
  if (infile==0) {
    mprinterr("Error: File '%s': %s\n", fname.c_str(), strerror( errno ));
    return false;
  }
  fclose(infile);
  return true;
}

// NumberFilename()
/** Given a filename and a number, append number to filename, i.e.
  * filename.number.
  */
std::string NumberFilename(std::string const &fname, int number) {
  std::ostringstream oss;
  oss << fname << "." << number;
  return oss.str();
}

// DigitWidth()
/** \return the number of characters necessary to express the given digit. */
int DigitWidth(long int numberIn) {
  double numf;
  int minusSign = 0;

  if (numberIn == 0L) return 1;
  if (numberIn < 0L) {
    numf = (double)(-numberIn);
    minusSign = 1;
  } else
    numf = (double) numberIn;

  numf = log10( numf );
  ++numf;
  // The cast back to long int implicitly rounds down
  int numi = (int)numf;
  return (minusSign + numi);
}

// FloatWidth()
/** \return the number of characters necessary to express given float. */
int FloatWidth(double floatIn) {
  double float_exponent = fabs( log10( floatIn ) );
  ++float_exponent;
  return (int)float_exponent; // Cast to int implicitly rounds down
}

// ---------- STRING CONVERSION ROUTINES --------------------------------------- 
/*! \class: BadConversion
    \brief Runtime exception for catching bad conversions from the convertToX routines.
  */
class BadConversion : public std::runtime_error {
public:
  BadConversion(std::string const &s)
    : std::runtime_error(s)
    { }
};

// convertToInteger()
/// Convert the input string to an integer.
int convertToInteger(std::string const &s) {
  std::istringstream iss(s);
  long int i;
  if (!(iss >> i))
    throw BadConversion("convertToInteger(\"" + s + "\")");
  return (int)i;
}

// convertToDouble()
/// Convert the input string to a double.
double convertToDouble(std::string const &s) {
  std::istringstream iss(s);
  double d;
  if (!(iss >> d))
    throw BadConversion("convertToDouble(\"" + s + "\")");
  return d;
}

// RemoveTrailingWhitespace()
/// Remove any trailing whitespace from string.
void RemoveTrailingWhitespace(std::string &line) {
  std::locale loc;

  std::string::iterator p = line.end();
  --p;
  for (; p != line.begin(); p--)
    if (!isspace( *p, loc) && *p!='\n' && *p!='\r') break;
  size_t lastSpace = (size_t)(p - line.begin()) + 1;
  //mprintf("lastSpace = %zu\n",lastSpace);
  if (lastSpace==1)
    line.clear();
  else
    line.resize( lastSpace );
}

std::string NoTrailingWhitespace(std::string const& line) {
  std::string duplicate(line);
  RemoveTrailingWhitespace(duplicate);
  return duplicate;
}

// integerToString()
std::string integerToString(int i) {
  std::ostringstream oss;
  oss << i;
  return oss.str();
}

// integerToString()
std::string integerToString(int i, int width) {
  std::ostringstream oss;
  oss.fill('0');
  oss.width( width );
  oss << std::right << i;
  return oss.str();
}

// doubleToString()
std::string doubleToString(double d) {
  std::ostringstream oss;
  oss << d;
  return oss.str();
}

// validInteger()
bool validInteger(std::string const &argument) {
  if (argument.empty()) return false;
  std::locale loc;
  std::string::const_iterator c;
  if (argument[0]=='-') {
    c = argument.begin()+1;
    if (c == argument.end()) return false;
  } else
    c = argument.begin();
  for (; c != argument.end(); ++c)
    if (!isdigit(*c,loc)) return false;
  return true;
}

// validDouble()
bool validDouble(std::string const &argument) {
  if (argument.empty()) return false;
  std::locale loc;
  std::string::const_iterator c;
  bool hasDecPt = (argument[0]=='.');
  if (argument[0]=='-' || hasDecPt) {
    c = argument.begin()+1;
    if (c == argument.end()) return false;
  } else
    c = argument.begin();
  if (!isdigit(*c,loc)) return false;
  for (; c != argument.end(); ++c)
  {
    if (*c == 'e' || *c == 'E')
      return validInteger( argument.substr(c - argument.begin() + 1) );
    bool isDecPt = (*c == '.');
    if (!isdigit(*c,loc) && (!isDecPt || (hasDecPt && isDecPt))) return false;
    if (isDecPt) hasDecPt = true;
  }
  return true;
}

// ---------- STRING FORMAT ROUTINES -------------------------------------------
// SetDoubleFormatString()
/** Set up a printf-style format string for float/double of given width, 
  * precision, and alignment, e.g. '%8.3lf'.
  */
std::string SetDoubleFormatString(int width, int precision, int type)
{
  std::string format;
  std::string width_arg;
  std::string prec_arg;
  std::string type_arg; // Will be f, lf, or E.

  // Type: 1 = float, 2 = scientific (E), otherwise double
  switch (type) {
    case 1:  type_arg = "f"; break;
    case 2:  type_arg = "E"; break;
    default: type_arg = "lf"; break;
  }
  // Set width and/or precision if applicable.
  if (width > 0)
    width_arg = integerToString( width );
  if (precision > -1)
    prec_arg = "." + integerToString( precision );
  // Set format string.
  format.append( "%" + width_arg + prec_arg + type_arg );
  return format; 
}

// SetStringFormatString()
/** Set up a printf-style format string for string (char*) of given
  * width and alignment, e.g. '%20s'.
  */
std::string SetStringFormatString(int width, bool leftAlign)
{
  std::string format;
  std::string width_arg;
  // If not left-aligned, need leading space.
  if (!leftAlign) 
    format.assign("%");
  else
    format.assign("%-");
  // Set width if applicable
  if (width > 0)
    width_arg = integerToString( width );
  // Set format string.
  format.append( width_arg + "s" );
  return format;
}

// SetIntegerFormatString()
/** Set up a printf-style format string for integer of given width
  * and alignment, e.g. '%8i'.
  */
std::string SetIntegerFormatString(int width)
{
  std::string format;
  std::string width_arg;
  // Set width if applicable
  if (width > 0)
    width_arg = integerToString( width );
  // Set format string.
  format.append( "%" + width_arg + "i" );
  return format;
}

// -----------------------------------------------------------------------------
std::string TimeString() {
  char buffer[20];
  time_t rawtime;

  time( &rawtime );
  struct tm* timeinfo = localtime( &rawtime );
  // Use snprintf for safety here. Speed not a factor for this routine.
  snprintf(buffer, 20, "%02i/%02i/%02i  %02i:%02i:%02i",
           timeinfo->tm_mon+1,timeinfo->tm_mday,timeinfo->tm_year%100,
           timeinfo->tm_hour,timeinfo->tm_min,timeinfo->tm_sec);
  return std::string( buffer );
}
// -----------------------------------------------------------------------------
long long AvailableMemory() {
#ifdef _MSC_VER
  MEMORYSTATUS status;
  GlobalMemoryStatus(&status);
  if (status.dwLength != sizeof(status))
    return -1;
  return (long long)status.dwAvailPhys;
#elif defined(__APPLE__)
  mach_msg_type_number_t count = HOST_VM_INFO_COUNT;
  vm_statistics_data_t vmstat;
  if (KERN_SUCCESS == host_statistics(mach_host_self(), HOST_VM_INFO,
                                      (host_info_t)&vmstat, &count)) {
    double total = vmstat.wire_count + vmstat.active_count +
                   vmstat.inactive_count + vmstat.free_count;
    double free = vmstat.free_count / total; // fraction
    return (long long)(TotalGlobalMemory() * free);
  }
  return -1;
#elif defined(_SC_AVPHYS_PAGES) && defined(_SC_PAGE_SIZE)
  long pages = sysconf(_SC_AVPHYS_PAGES);
  long page_size = sysconf(_SC_PAGE_SIZE);
  if (pages < 0L || page_size < 0L) return -1;
  return (long long)(pages * page_size);
#else
  return -1;
#endif
}

double AvailableMemory_MB() { 
  double avail_in_bytes = AvailableMemory();
  if (avail_in_bytes < 0.0) 
    return -1.0;
  else
    return (double)AvailableMemory() / (1024 * 1024);
}

#ifdef __APPLE__
long long TotalGlobalMemory() {
    int mib[] = {CTL_HW, HW_MEMSIZE};
    int64_t size = 0;
    size_t len = sizeof(size);
    if (sysctl(mib, 2, &size, &len, NULL, 0) == 0)
        return (long long) size;
    return 0ll;
}
#endif
