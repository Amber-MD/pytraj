// Range
#include <algorithm> // find
#include "Range.h"
#include "ArgList.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Range::Range() { }

// CONSTRUCTOR
/// Takes argument string as input
Range::Range( std::string const& argIn ) {
  if (!argIn.empty())
    SetRange( argIn );
}

Range::Range( std::string const& argIn, int offsetIn) {
  if (!argIn.empty()) {
    SetRange( argIn );
    ShiftBy( offsetIn );
  }
}

// COPY CONSTRUCTOR
Range::Range(const Range &rhs) :
  rangeArg_(rhs.rangeArg_),
  rangeList_(rhs.rangeList_)
{}

// ASSIGNMENT OPERATOR
Range &Range::operator=(const Range &rhs) {
  // Check for self assignment
  if ( this == &rhs ) return *this;
  rangeArg_ = rhs.rangeArg_;
  rangeList_ = rhs.rangeList_;
  return *this;
}

// Range::SetRange()
/** Given an argument containing numbers separated by "," (concatentation), and 
  * "-" (number range), construct an ordered list of numbers corresponding to 
  * the argument. Remove any duplicate numbers.
  * \return 0 on success, 1 on error.
  */
int Range::SetRange(std::string const& ArgIn) {
  std::string arg;
  int R[2], upper;
  ArgList DashList;

  //mprintf("DEBUG: SetRange(%s)\n",ArgIn);

  if (ArgIn.empty()) return 1;
  rangeList_.clear();

  // Set rangeArg
  rangeArg_.assign(ArgIn);
  // Check if ArgIn is a mask expression
  size_t maskcharpos = rangeArg_.find_first_of(":@*");
  if (maskcharpos != std::string::npos) {
    mprinterr("Error: Using a mask expression for range (%s)\n",ArgIn.c_str());
    mprinterr("Error: Ranges should only contain digits, dashes, and commas (e.g. 3-5,8-10)\n");
    return 1;
  } 
  // Split range by comma
  ArgList CommaList(rangeArg_, ",");
  int err = 0;
  while ( err == 0 ) {
    arg = CommaList.GetStringNext();
    if (arg.empty()) break; // Exit the while loop.
    // Then split by dash
    DashList.SetList(arg, "-");
    R[0] = DashList.getNextInteger(-1);
    R[1] = DashList.getNextInteger(-1);
    if (R[0]==-1) {
      mprinterr("Error: Range::SetRange(%s): Range is -1 for %s\n",ArgIn.c_str(), 
                DashList.ArgLine());
      err=1;
      break;
    }
    upper = R[1];
    if (upper==-1) upper=R[0];
    ++upper; // Want up to and including the upper argument
    if ( this->SetRange(R[0], upper) )
      mprintf("Warning: Converting %s to range [%i-%i] is not valid.\n",
              ArgIn.c_str(), R[0], R[1]);
  }

  // Dont return an empty list
  if ( err>0 || rangeList_.empty() ) 
    return 1;
  
  // Sort frames using default comparison
  rangeList_.sort();
  //for (it=rangeList_.begin(); it!=rangeList_.end(); it++)
  //  fprintf(stdout,"RangeList= %i\n",*it); 
  // Remove duplicates
  err=-1;
  //fprintf(stdout,"Size of RangeList is %lu\n",rangeList_.size());
  std::list<int>::iterator it = rangeList_.begin();
  while (!rangeList_.empty()) {
    //fprintf(stdout,"     List= %i  Last= %i",*it,err);
    upper=*it;
    if (*it == err) {
      //fprintf(stdout,"REMOVING."); 
      // Erasing effectively increments the iterator
      it = rangeList_.erase(it);
    } else {
      ++it;
    }
    err=upper;
    //fprintf(stdout,"\n");
    // If we are past the last element exit now
    if (it == rangeList_.end()) break;
  }

  return 0;
}

// Range::SetRange()
/** Given a start and end number, set up a range from start to (not 
  * including) end.  
  */
int Range::SetRange(int start, int end) {
  //Check that end is greater than start so that the range list is
  if (end <= start) {
    mprintf("Error: Range::SetRange: end (%i) <= start (%i)\n",end,start);
    return 1;
  }
  for (int range=start; range < end; range++)
    rangeList_.push_back(range);

  return 0;
}

// Range::ShiftBy()
/** Shift all numbers in range by specified value. */
void Range::ShiftBy(int val) {
  for (std::list<int>::iterator rangeNum = rangeList_.begin(); 
                                rangeNum != rangeList_.end(); ++rangeNum)
    *rangeNum += val;
}

// Range::RemoveFromRange()
/** Remove all instances of num from the range. */
void Range::RemoveFromRange(int num) {
  std::list<int>::iterator it=rangeList_.begin();
  while (it!=rangeList_.end()) {
    if (*it == num) 
      it = rangeList_.erase(it);
    else
      it++;
  }
}

// Range::PrintRange()
/** Print all numbers in the range to a line. Increment by offset. */
void Range::PrintRange(const char* header, int offset) const {
  if (header!=0)
    mprintf("%s",header);
  for (std::list<int>::const_iterator it=rangeList_.begin(); it!=rangeList_.end(); it++)
    mprintf(" %i",(*it)+offset);
  //mprintf("\n");
}

bool Range::InRange(int idx) const {
  return ( std::find( rangeList_.begin(), rangeList_.end(), idx ) != rangeList_.end() );
}
