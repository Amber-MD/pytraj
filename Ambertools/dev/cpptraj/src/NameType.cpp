#include "NameType.h"

NameType::NameType() :
  NameSize_(6)
{
  c_array_[0]=' ';
  c_array_[1]=' ';
  c_array_[2]=' ';
  c_array_[3]=' ';
  c_array_[4]=' ';
  c_array_[5]='\0';
}

NameType::NameType(const NameType &rhs) :
 NameSize_(6)
{
  c_array_[0] = rhs.c_array_[0];
  c_array_[1] = rhs.c_array_[1];
  c_array_[2] = rhs.c_array_[2];
  c_array_[3] = rhs.c_array_[3];
  c_array_[4] = rhs.c_array_[4];
  // 5 is always null 
  c_array_[5]='\0';
}

NameType::NameType(const char *rhs) :
  NameSize_(6)
{
  const char *ptr = rhs;
  for (unsigned int j = 0; j < NameSize_-1; j++) {
    c_array_[j] = *ptr;
    if (*ptr=='\0') break;
    ++ptr;
  }
  FormatName();
}

NameType::NameType(std::string const& str) :
  NameSize_(6)
{
  unsigned int ns1 = NameSize_ - 1;
  unsigned int strend = (unsigned int)str.size();
  if (strend > ns1)
    strend = ns1;
  for (unsigned int j = 0; j < strend; j++) 
    c_array_[j] = str[j];
  c_array_[strend] = '\0';
  FormatName();
}
 
NameType &NameType::operator=(const NameType &rhs) {
  if (&rhs==this) return *this;
  c_array_[0] = rhs.c_array_[0];
  c_array_[1] = rhs.c_array_[1];
  c_array_[2] = rhs.c_array_[2];
  c_array_[3] = rhs.c_array_[3];
  c_array_[4] = rhs.c_array_[4];
  // 5 is always null 
  return *this;
}

// NameType::ToBuffer()
/** For interfacing with old C stuff. Only set 1st 4 chars. */
void NameType::ToBuffer(char *buffer) const {
  buffer[0] = c_array_[0];
  buffer[1] = c_array_[1];
  buffer[2] = c_array_[2];
  buffer[3] = c_array_[3];
  buffer[4] = '\0';
}

bool NameType::Match(NameType const& maskName) const { 
  int c = 0;
  for (int m = 0; m < 5; m++) {
    if (maskName.c_array_[m] == '\0' && c_array_[c] == ' ')
      // At end of mask and whitespace in name: OK
      break;
    if (maskName.c_array_[m] == '\\') { 
      // Backslash: match literal next char in mask
      ++m;
      if (maskName.c_array_[m] != c_array_[c])
        return false;
    } else if (maskName.c_array_[m] == '*') { 
      // Mask wildcard: instant match
      return true;
    } else if (maskName.c_array_[m] != '?' && 
               maskName.c_array_[m] != c_array_[c]) { 
      // Not mask single wildcard and mismatch
      return false;
    }
    //mprintf("(%c,%c)",maskName.c_array_[m],c_array_[c]);
    ++c;
  }
  return true;
}

bool NameType::operator==(const NameType &rhs) const {
  if (c_array_[0] != rhs.c_array_[0]) return false;
  if (c_array_[1] != rhs.c_array_[1]) return false;
  if (c_array_[2] != rhs.c_array_[2]) return false;
  if (c_array_[3] != rhs.c_array_[3]) return false;
  return true;
}

bool NameType::operator==(const char *rhs) const {
  NameType tmp(rhs);
  return (*this == tmp);
}

bool NameType::operator!=(const NameType &rhs) const {
  if (c_array_[0] != rhs.c_array_[0]) return true;
  if (c_array_[1] != rhs.c_array_[1]) return true;
  if (c_array_[2] != rhs.c_array_[2]) return true;
  if (c_array_[3] != rhs.c_array_[3]) return true;
  return false;
}

bool NameType::operator!=(const char *rhs) const {
  NameType tmp(rhs);
  return (*this != tmp);
}

char NameType::operator[](int idx) const {
  if (idx < 0 || idx >= (int)NameSize_) return '\0';
  return c_array_[idx];
}

std::string NameType::Truncated() const {
  std::string outname( c_array_ );
  // Remove trailing spaces
  if (outname[3] == ' ') outname.resize(3);
  if (outname[2] == ' ') outname.resize(2);
  if (outname[1] == ' ') outname.resize(1);
  return outname;
}

/** Replace asterisks with a single quote */
void NameType::ReplaceAsterisk() {
  if (c_array_[0]=='*') c_array_[0]='\'';
  if (c_array_[1]=='*') c_array_[1]='\'';
  if (c_array_[2]=='*') c_array_[2]='\'';
  if (c_array_[3]=='*') c_array_[3]='\'';
}

// NameType::FormatName()
/** For consistency with Amber names, replace any null in the first 4 chars
  * with spaces. Remove any leading whitespace.
  */
void NameType::FormatName() 
{
  // Ensure at least 4 chars long.
  if (c_array_[0]=='\0') { // 0 chars
    c_array_[0]=' ';
    c_array_[1]=' ';
    c_array_[2]=' ';
    c_array_[3]=' ';
    c_array_[4]='\0';
  } else if (c_array_[1]=='\0') { // 1 char
    c_array_[1]=' ';
    c_array_[2]=' ';
    c_array_[3]=' ';
    c_array_[4]='\0';
  } else if (c_array_[2]=='\0') { // 2 chars
    c_array_[2]=' ';
    c_array_[3]=' ';
    c_array_[4]='\0';
  } else if (c_array_[3]=='\0') { // 3 chars
    c_array_[3]=' ';
    c_array_[4]='\0';
  }
  // Remove leading whitespace.
  if        (c_array_[0]==' ') { // Some leading whitespace
    if (c_array_[1]!=' ') {        // [_XXX]
      c_array_[0]=c_array_[1];
      c_array_[1]=c_array_[2];
      c_array_[2]=c_array_[3];
      c_array_[3]=' ';
    } else if (c_array_[2]!=' ') { // [__XX]
      c_array_[0]=c_array_[2];
      c_array_[1]=c_array_[3];
      c_array_[2]=' ';
      c_array_[3]=' ';
    } else if (c_array_[3]!=' ') { // [___X]
      c_array_[0]=c_array_[3];
      c_array_[1]=' ';
      c_array_[2]=' ';
      c_array_[3]=' ';
    }
  }
}
