#include <cstdio> // sscanf
#include "SDFfile.h"
#include "StringRoutines.h" // RemoveTrailingWhitespace

// CONSTRUCTOR
SDFfile::SDFfile() : debug_(0), Natoms_(0), Nbonds_(0) {}


bool SDFfile::ID_SDF(CpptrajFile& fileIn) {
  // NOTE: ASSUMES FILE IS ALREADY SETUP!
  if (fileIn.OpenFile()) return false;
  // Search for V2000 somewhere in line 4
  const char* ptr = 0;
  for (int i = 0; i < 4; i++)
    if ( (ptr = fileIn.NextLine()) == 0 ) {
      fileIn.CloseFile();
      return false;
    }
  fileIn.CloseFile();
  std::string line( ptr ); // Line 4, Connection table
  if ( line.find( "V2000" ) != std::string::npos ) return true;
  return false;
}

bool SDFfile::ReadHeader() {
  if (!IsOpen()) return true;
  // Read title
  title_ = GetLine();
  RemoveTrailingWhitespace( title_ );
  // Read past info and comment
  if (NextLine() == 0) return true;
  if (NextLine() == 0) return true;
  // Read Connection table
  const char* ptr = NextLine();
  if (ptr == 0) return true;
  if (sscanf(ptr, "%i %i", &Natoms_, &Nbonds_) != 2) return true;
  return false;
}

int SDFfile::SDF_XYZ(double* XYZ) {
  const char* ptr = NextLine();
  if (ptr == 0) return 1;
  if (sscanf(ptr, "%lf %lf %lf %s", XYZ, XYZ+1, XYZ+2, Name_) != 4) return 1;
  return 0;
}

Atom SDFfile::SDF_Atom() {
  return Atom( NameType(Name_), ' ', Name_ );
}

int SDFfile::SDF_Bond(int& at1, int& at2) {
  const char* ptr = NextLine();
  if (ptr == 0) return 1;
  if (sscanf(ptr, "%i %i", &at1, &at2) != 2) return 1;
  return 0;
}
