#include <cstring> //strtok, strchr
#include "ArgList.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" //validDouble, validInteger, convertToX

const std::string ArgList::emptystring = "";

ArgList::ArgList(const char* input) {
  if (input != 0) SetList( std::string(input), " ");
}

ArgList::ArgList(std::string const& input) {
  SetList( input, " " );
}

ArgList::ArgList(std::string const& input, const char* sep) {
  SetList( input, sep );
}

// COPY CONSTRUCTOR
ArgList::ArgList(const ArgList &rhs) :
  argline_(rhs.argline_),
  arglist_(rhs.arglist_),
  marked_(rhs.marked_)
{}

// ArgList::operator=()
ArgList &ArgList::operator=(const ArgList &rhs) {
  if (&rhs==this) return *this;
  // Allocate and copy. Assignment ops should automatically deallocate.
  argline_ = rhs.argline_;
  arglist_ = rhs.arglist_;
  marked_ = rhs.marked_;
  return *this;
}

// ArgList::operator[]
std::string const& ArgList::operator[](int idx) const {
  if (idx < 0 || idx >= (int)arglist_.size()) {
    mprinterr("Internal Error: Position %i out of range for Argument List.\n",idx);
    return emptystring;
  }
  return arglist_[idx];
}

// ArgList::ClearList()
void ArgList::ClearList() {
  argline_.clear();
  arglist_.clear();
  marked_.clear();
}

static inline char stringBack(std::string const& str) {
  if (str.empty()) return 0;
  return str[ str.size() - 1];
}

// ArgList::SetList()
/** Separate input by the characters in separator and store as separate args.
  * This overwrites any existing args and completely resets the list.
  * \param inputString null-terminated string to be converted to arguments
  * \param separator string containing characters used to separate arguments
  * \return 0 if arglist successfully set up, 1 if not.
  */
int ArgList::SetList(std::string const& inputString, const char *separator) {
  if (inputString.empty() || separator==0) return 1;
  // Copy inputString to temp since it is destroyed by tokenize.
  size_t inputStringSize = inputString.size();
  if (inputStringSize < 1) return 1;
  char* tempString = new char[ inputStringSize+1 ];
  inputString.copy( tempString, inputStringSize, 0 );
  tempString[ inputStringSize ] = '\0'; // copy() does not append null 
  // Remove newline char from tempString if present
  if ( tempString[ inputStringSize - 1 ] == '\n' )
    tempString[ inputStringSize - 1 ] = '\0';
  // Free existing arglist
  arglist_.clear();
  marked_.clear();
  // Store inputString
  argline_.assign(inputString);
  // Begin tokenization
  char* pch = strtok(tempString, separator);
  if (pch != 0) {
    while (pch != 0) {
      // If the argument is not quoted add it to the list
      if (pch[0] != '\"' && pch[0] != '\'') 
        arglist_.push_back( std::string(pch) );
      else {
        char quotechar = pch[0];
        // If the argument begins with a quote, place this and all subsequent
        // arguments this argument until closing quote reached.
        std::string argument(pch);
        // Erase beginning quote.
        argument.erase(0, 1);
        // Make sure this argument doesnt just end with a quote.
        if (stringBack(argument) != quotechar) {
          while (pch != 0) {
            argument.append(" ");
            pch = strtok(0, separator);
            // If pch is 0 at this point there was no closing quote.
            if (pch == 0) {
              mprintf("Warning: argument list closing quote [%c] missing or misplaced\n", 
                      quotechar);
              break;
            }
            argument.append(pch);
            // Check if this argument has the closing quote.
            if (stringBack(argument) == quotechar) break;
          }
        }
        // Remove final quote from the argument. If there was no closing quote
        // this removes the final space.
        if (!argument.empty()) argument.resize( argument.size() - 1 );
        if (!argument.empty()) arglist_.push_back(argument);
      }
      pch = strtok(0,separator); // Next argument
    } // END while loop
    // Set up marked array
    marked_.resize( arglist_.size(), false );
  }
  delete[] tempString;
  return 0;
}

// ArgList::RemainingArgs()
ArgList ArgList::RemainingArgs() {
  ArgList remain;
  for (unsigned int arg = 0; arg < arglist_.size(); ++arg) {
    if ( !marked_[arg] ) {
      remain.arglist_.push_back( arglist_[arg] );
      if (!remain.argline_.empty()) remain.argline_.append(" ");
      remain.argline_.append( arglist_[arg] );
      marked_[arg] = true;
    }
  }
  remain.marked_.resize( remain.arglist_.size(), false );
  return remain;
}

// ArgList::AddArg()
/** \param input string of space-delimited args to add to argument list.
  */
void ArgList::AddArg(std::string const& input) {
  ArgList inputArgs( input );
  for (int i = 0; i < inputArgs.Nargs(); i++) {
    arglist_.push_back( inputArgs[i] );
    argline_.append(" ");
    argline_.append( inputArgs[i] );
    marked_.push_back( false );
  }
}

// ArgList::MarkArg()
void ArgList::MarkArg(int arg) {
  if (arg < 0 || arg >= (int) marked_.size()) return;
  marked_[arg]=true;
}

// ArgList::CheckForMoreArgs()
/** Check if all arguments have been processed. If not print a warning along
  * with all unprocessed arguments.
  */
bool ArgList::CheckForMoreArgs() const {
  std::string notmarked;
  for (unsigned int arg=0; arg < arglist_.size(); arg++) {
    if (!marked_[arg]) 
      notmarked.append(arglist_[arg] + " ");
  }
  if (!notmarked.empty()) { 
    mprinterr("Error: [%s] Not all arguments handled: [ %s]\n",
              arglist_[0].c_str(), notmarked.c_str());
    return true;
  }
  return false;
}

// ArgList::PrintList()
void ArgList::PrintList() const {
  for (unsigned int arg = 0; arg < arglist_.size(); arg++) 
    mprintf("  %u: %s\n",arg+1,arglist_[arg].c_str());
}

// ArgList::PrintDebug()
void ArgList::PrintDebug() const {
  mprintf("ArgLine: %s\n",argline_.c_str());
  for (unsigned int arg = 0; arg < arglist_.size(); arg++)
    mprintf("\tArg %u: %s (%i)\n",arg+1,arglist_[arg].c_str(),(int)marked_[arg]);
}

// ArgList::RemoveFirstArg()
void ArgList::RemoveFirstArg() {
  if (arglist_.empty()) return;
  arglist_.erase( arglist_.begin() );
  marked_.erase( marked_.begin() );
}

// ArgList::Command()
/* \return pointer to the first argument
 */
const char *ArgList::Command() const {
  if (arglist_.empty()) return 0;
  return arglist_[0].c_str();
}

// ArgList::CommandIs()
/** \param key Key to check first argument against
  * \return true if first argument matches key
  */
bool ArgList::CommandIs(const char *key) const {
  if (arglist_.empty()) return false;
  if (arglist_[0].compare( key )==0) return true;
  return false;
}

// ArgList::GetStringNext()
std::string const& ArgList::GetStringNext() {
  for (unsigned int arg = 0; arg < arglist_.size(); ++arg)
    if (!marked_[arg]) {
      marked_[arg]=true;
      return arglist_[arg];
    }
  return emptystring;
}

// ArgList::GetMaskNext()
/** Return next unmarked Mask. A mask MUST include one of the following: 
  *   ':' residue
  *   '@' atom
  *   '*' everything
  * \return the next unmarked atom mask expression
  */
std::string const& ArgList::GetMaskNext() {
  for (unsigned int arg = 0; arg < arglist_.size(); ++arg) {
    if (!marked_[arg]) {
      size_t found = arglist_[arg].find_first_of(":@*");
      if (found != std::string::npos) {
        marked_[arg] = true;
        return arglist_[arg];
      }
    }
  }
  return emptystring;
}

// ArgList::getNextTag()
/** Return the next unmarked tag. A tag is defined as a character string
  * bounded by brackets, e.g. [tag].
  */
std::string const& ArgList::getNextTag() {
  for (unsigned int arg = 0; arg < arglist_.size(); arg++) {
    if (!marked_[arg]) {
      std::string::reverse_iterator lastchar  = arglist_[arg].rbegin();
      std::string::iterator         firstchar = arglist_[arg].begin();
      if (*firstchar=='[' && *lastchar==']') {
        marked_[arg]=true;
        return arglist_[arg];
      }
    }
  }
  return emptystring;
}

// ArgList::getNextInteger()
/** \param def Value to return if no integer args found
  * \return Next unmarked integer argument or def
  */
int ArgList::getNextInteger(int def) {
  for (unsigned int arg=0; arg < arglist_.size(); arg++)
    if (!marked_[arg]) {
      // Check that first char is indeed an integer or '-', if not then continue
      if (validInteger(arglist_[arg])) {
        int ival = convertToInteger(arglist_[arg]);
        marked_[arg]=true;
        return ival;
      }
    }
  return def;
}

// ArgList::getNextDouble()
/** \param def Value to return if no double args found
  * \return Next unmarked double argument or def
  */
double ArgList::getNextDouble(double def) {
  for (unsigned int arg=0; arg < arglist_.size(); arg++)
    if (!marked_[arg]) {
      // Check that first char is indeed a digit, '.', or '-', if not then continue
      if (validDouble(arglist_[arg])) {
        double dval = convertToDouble(arglist_[arg]);
        marked_[arg]=true;
        return dval;
      }
    }
  return def;
}

// ArgList::GetStringKey()
/** Search the argument list for key, return the argument following key
  * as a string if found, otherwise return 0.
  * \param key string to search for
  */
std::string const& ArgList::GetStringKey(const char *key) {
  int nargs = (int)arglist_.size() - 1;
  for (int arg=0; arg < nargs; arg++)
    if (!marked_[arg]) {
      if (arglist_[arg].compare(key)==0) {
        marked_[arg]=true;
        arg++;
        marked_[arg]=true;
        return arglist_[arg];
      }
    }
  return emptystring;
}

// ArgList::getKeyInt()
/** Search the argument list for key, return the argument following key
  * as an integer if found, otherwise return def.
  * \param key string to search for
  * \param def Value to return if key not found.
  */
int ArgList::getKeyInt(const char *key, int def) {
  int nargs = (int)arglist_.size() - 1;
  for (int arg=0; arg < nargs; arg++)
    if (!marked_[arg]) {
      if (arglist_[arg].compare(key)==0) {
        if (validInteger(arglist_[arg+1])) {
          marked_[arg]=true;
          arg++;
          int ival = convertToInteger(arglist_[arg]);
          marked_[arg]=true;
          return ival;
        }
      }
    }
  return def;
}

// ArgList::getKeyDouble()
/** Search the argument list for key, return the argument following key
  * as a double if found, otherwise return def.
  * \param key string to search for
  * \param def Value to return if key not found.
  */
double ArgList::getKeyDouble(const char *key, double def) {
  int nargs = (int)arglist_.size() - 1;
  for (int arg=0; arg < nargs; arg++)
    if (!marked_[arg]) {
      if (arglist_[arg].compare(key)==0) {
        if (validDouble(arglist_[arg+1])) {
          marked_[arg]=true;
          arg++;
          double dval = convertToDouble(arglist_[arg]);
          marked_[arg]=true;
          return dval;
        }
      }
    }
  return def;
}

// ArgList::hasKey()
/** Search the argument list for key, mark and return true if found.
  * \param key string to search for
  * \return true if key is found, false if not.
  */
bool ArgList::hasKey(const char *key) {
  for (unsigned int arg = 0; arg < arglist_.size(); arg++) 
    if (!marked_[arg]) {
      if (arglist_[arg].compare(key)==0) {
        marked_[arg]=true;
        return true;
      }
    }
  return false;
}

// ArgList::Contains()
/** \param key string to search for
  * \return true if key is found, false if not.
  */
// NOTE: Should this be ignoring previously marked strings?
bool ArgList::Contains(const char *key) const {
  for (unsigned int arg = 0; arg < arglist_.size(); arg++) 
    if (!marked_[arg]) {
      if (arglist_[arg].compare(key)==0) {
        return true;
      }
    }
  return false;
}
