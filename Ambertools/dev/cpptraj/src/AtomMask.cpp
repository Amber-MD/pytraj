#include <algorithm> // sort, unique
#include <locale> // isspace
#include <stack>
#include "AtomMask.h"
#include "CpptrajStdio.h"
#include "ArgList.h"

// CONSTRUCTOR
AtomMask::AtomMask() :
  debug_(0),
  maskChar_('T'),
  Natom_(0),
  nselected_(0)
{}

// CONSTRUCTOR
AtomMask::AtomMask(std::string const& maskstring) :
  debug_(0),
  maskChar_('T'),
  Natom_(0),
  nselected_(0)
{
  SetMaskString(maskstring);
}

// CONSTRUCTOR
AtomMask::AtomMask(int beginAtom, int endAtom) :
  debug_(0),
  maskChar_('T'),
  Natom_(0),
  nselected_(0)
{
  AddAtomRange(beginAtom, endAtom);
}

AtomMask::AtomMask(int atomNum) : debug_(0), maskChar_('T'), Natom_(1),
  nselected_(1), Selected_(1, atomNum)
{}

// COPY CONSTRUCTOR
AtomMask::AtomMask(const AtomMask &rhs) : 
  debug_(rhs.debug_),
  CharMask_(rhs.CharMask_),
  maskChar_(rhs.maskChar_),
  maskString_(rhs.maskString_),
  Natom_(rhs.Natom_),
  nselected_(rhs.nselected_),
  Selected_(rhs.Selected_),
  maskTokens_(rhs.maskTokens_)
{}

// AtomMask::operator=()
/// AtomMask assignment
AtomMask& AtomMask::operator=(const AtomMask &rhs) {
  // Check for self assignment
  if ( this == &rhs ) return *this;
  // Deallocate
  // Allocate and copy
  debug_ = rhs.debug_;
  nselected_ = rhs.nselected_;
  Natom_ = rhs.Natom_;
  maskChar_ = rhs.maskChar_;
  maskString_ = rhs.maskString_;
  Selected_ = rhs.Selected_;
  CharMask_ = rhs.CharMask_;
  maskTokens_ = rhs.maskTokens_;
  // Return *this
  return *this;
}

static bool IsOperator(char op) {
  if (op=='!') return true;
  if (op=='&') return true;
  if (op=='|') return true;
  if (op=='<') return true;
  if (op=='>') return true;
  return false;
}

static bool IsOperand(char op) {
  std::locale loc;
  if (op=='*')  return true;
  if (op=='/')  return true;
  if (op=='\\')  return true;
  if (op=='%')  return true;
  if (op=='-')  return true;
  if (op=='?')  return true;
  if (op==',')  return true;
  if (op=='\'') return true;
  if (op=='.')  return true;
  if (op=='=')  return true;
  if (op=='+')  return true;
  if (isalnum(op, loc)) return true;
  return false;
}

static int OperatorPriority(char op) {
  if (op == '>') return(6);
  if (op == '<') return(6);
  if (op == '!') return(5);
  if (op == '&') return(4);
  if (op == '|') return(3);
  if (op == '(') return(2);
  if (op == '_') return(1);

  mprinterr("OperatorPriority(): unknown operator ==%c== on stack when processing atom mask",op);
  return(0);
}

// AtomMask::Tokenize()
/** STEP1: preprocess the input string:
  * - remove spaces 
  * - isolate 'operands' into brackets [...]
  * - split expressions of the type :1-10@CA,CB into two parts;
  *   the two parts are joined with '&' operator and (for the sake
  *   of preserving precedence of other operators) enclosed into (..)
  *   :1-10@CA,CB    is split into  (:1-10 & @CA,CB)
  * - do basic error checking
  *
  * STEP2: convert to RPN
  * - operands (part enclosed in [..]) are copied directly to 'postfix'
  * - left parentheses are always pushed onto the stack
  * - when a right parenthesis is encountered the symbol at the top
  *   of the stack is popped off the stack and copied to 'postfix'
  *   until the symbol at the top of the stack is a left parenthesis.
  *   When that occurs both parentheses are discarded
  * - if the symbol scanned from 'infix' has a higher precedence 
  *   then the symbol at the top of the stack, the symbol being 
  *   scanned is pushed onto the stack
  * -  if the precedence of the symbol being scanned is lower than 
  *   or equal to the precedence of the symbol at the top of the 
  *   stack, the stack is popped to 'postfix' until the condition
  *   holds
  * - when the terminating symbol '_' is reached on the input scan
  *   the stack is popped to 'postfix' until the terminating symbol
  *   is also reached on the stack. Then the algorithm terminates.
  * - if the top of the stack is '(' and the terminating symbol '_'
  *   is scanned, or ')' is scanned when '_' is at the top of the
  *   stack, the parentheses of the atom expression were unbalanced
  *   and an unrecoverable error has occurred.
  *
  * \authors Daniel R. Roe, based on the tokenize() routine from PTRAJ
  *          by Viktor Hornak (2003).
  */
int AtomMask::Tokenize() {
  std::string infix;
  std::string buffer;
  std::locale loc;
  std::string postfix;
  std::stack<char> Stack;

  // 0 means new operand or operand was just completed, and terminated with ']', 
  // 1 means operand with ":" read,
  // 2 means operand with "@" read
  // 3 means '<' or '>' read, waiting for numbers.
  int flag = 0;

  for (std::string::iterator p = maskString_.begin(); p != maskString_.end(); p++)
  {
    // Skip spaces and newlines
    if ( isspace(*p, loc) ) continue;
    
    if ( IsOperator(*p) || *p == '(' || *p == ')' ) {
      if (flag > 0) {
        buffer += "])";
        flag = 0;
        infix += buffer;
      }

      infix += *p;

      if ( *p == '>' || *p == '<' ) {
        buffer.assign("([");
        buffer += *p;
        ++p;
        buffer += *p;
        flag = 3;
        if ( *p != ':' && *p != '@' ) {
          --p;
          mprinterr("Error: Tokenize: Wrong syntax for distance mask [%c]\n",*p);
          return 1;
        }
      }
    } else if ( IsOperand(*p) || isalnum(*p, loc) ) {
      if (flag==0) {
        buffer.assign("([");
        flag = 1;
        if ( *p != '*') {
          mprinterr("Error: Tokenize: Wrong syntax [%c]\n",*p);
          return 1;
        }
      }
      if (*p == '=') { // The AMBER9 definition of wildcard '=' is equivalent to '*'.
        if (flag > 0)
          *p = '*';
        else {
          mprinterr("Error: Tokenize: '=' not in name list syntax\n");
          return 1;
        }
      }
      buffer += *p;
    } else if ( *p == ':' ) {
      if (flag == 0) {
        buffer.assign("([:");
        flag = 1;
      } else {
        buffer += "])|([:";
        flag = 1;
      }
    } else if ( *p == '@' ) {
      if (flag == 0) {
        buffer.assign("([@");
        flag = 2;
      } else if (flag == 1) {
        buffer += "]&[@";
        flag = 2;
      } else if (flag == 2) {
        buffer += "])|([@";
        flag = 2;
      }
    } else {
      mprinterr("Error: Tokenize: Unknown symbol (%c) expression when parsing atom mask [%s]\n",
                *p, maskString_.c_str());
      return 1;
    }
  } // END for loop over maskString
  // Terminate buffer if necessary
  if (flag > 0) {
    buffer += "])";
    infix += buffer;
  }

  // NOTE: Check for malformed tokens?

  // Add terminal symbol '_', needed for RPN conversion
  infix += "_";

  if (debug_ > 0)
    mprintf("DEBUG: NEW_INFIX ==%s==\n",infix.c_str());

  // -----------------------------------
  // Convert to RPN
  //postfix.clear(); 
  // push terminal symbol '_' to stack
  Stack.push('_');

  // 1 when start with "[", 0 when finished.
  flag = 0;
  char pp = ' ';
  for (std::string::const_iterator p = infix.begin(); p != infix.end(); p++) {
    if (*p == '[') {
      postfix += *p;
      flag = 1;
    } else if (*p == ']') {
      postfix += *p;
      flag = 0;
    } else if (flag == 1) {
      postfix += *p;
    } else if (*p == '(') {
      Stack.push(*p);
    } else if (*p == ')') {
      while ( (pp = Stack.top()) != '(') {
        Stack.pop();
        if (pp == '_') {
          mprinterr("Error: Mask::ToRPN: unbalanced parentheses in expression\n");
          return 1;
        }
        postfix += pp;
      }
      Stack.pop(); // Discard '('
      // At this point both parentheses are discarded
    } else if (*p == '_') {
      while ( (pp = Stack.top()) != '_') {
        Stack.pop();
        if (pp == '(') {
          mprinterr("Error: Mask::ToRPN: unbalanced parentheses in expression\n");
          return 1;
        }
        postfix += pp;
      }
      Stack.pop(); // Discard '_'
    } else if ( IsOperator(*p) ) {
      int P1 = OperatorPriority( *p );
      int P2 = OperatorPriority( Stack.top() );
      if ( P1==0 || P2==0 ) return 1; // 0 indicates error in op
      if (P1 > P2) {
        Stack.push( *p );
      } else {
        while ( P1 <= P2 ) {
          pp = Stack.top();
          Stack.pop();
          postfix += pp;
          P1 = OperatorPriority( *p );
          P2 = OperatorPriority( Stack.top() );
          if ( P1==0 || P2==0 ) return 1; // 0 indicates error in op
        }
        Stack.push( *p );
      }
    } else {
      mprinterr("Error: ToRPN: Unknown symbol in atom mask (%c)\n", *p);
      return 1;
    } 
  } // END for loop over infix
  if (debug_ > 0)
  mprintf("DEBUG: NEW_POSTFIX ==%s==\n",postfix.c_str());

  // Convert to MaskTokens in same order. The postfix expression is composed
  // of operands enclosed within brackets, and single character operators.
  // The exception is the distance operator, which is also an operand. An
  // operand can have multiple entries separated by commas (e.g. [:1,2-7,5] 
  // has 3). Once an operand is complete the OnStack bit of the token is set
  // to indicate the mask should go on the stack for processing by operators.
  // Operators store the result of their operation on the mask on top of
  // the stack so they dont need to be pushed.
  std::string tokenString;
  MaskToken token;
  maskTokens_.clear();
  for (std::string::const_iterator p = postfix.begin(); p != postfix.end(); p++) 
  {  // Operand begins here
    if (*p == '[')
      buffer.clear();
    // Operand is completed
    else if (*p == ']') {
      //mprintf("PROCESSING MASK OPERAND [%s]\n",buffer.c_str());
      if (buffer[0]=='<' || buffer[0]=='>') {
        // Distance criterion operand/operator
        // Since operator, resulting mask doesnt need to go on stack.
        token.SetDistance( buffer );
        maskTokens_.push_back( token );
      } else if (buffer[0]=='*') {
        // Select all operand - result goes on stack
        token.SetOperator( MaskToken::SelectAll );
        maskTokens_.push_back( token );
        maskTokens_.back().SetOnStack();
      } else {
        // Basic Operand list. After last entry in list processed result goes
        // on stack.
        // Determine type from first char. Default to Num; MaskToken::SetToken
        // will convert to Name if appropriate.
        MaskToken::MaskTokenType tokenType = MaskToken::OP_NONE;
        if (buffer[0]==':') // Residue
          tokenType = MaskToken::ResNum; 
        else if (buffer[0]=='@') { // Atom
          tokenType = MaskToken::AtomNum; 
          if      (buffer[1]=='%') tokenType = MaskToken::AtomType;
          else if (buffer[1]=='/') tokenType = MaskToken::AtomElement;
        }
        if (tokenType==MaskToken::OP_NONE) {
          mprinterr("Error: Unrecognized token type.\n");
          maskTokens_.clear();
          return 1;
        }
        // Create new string without type character(s)
        if (tokenType==MaskToken::ResNum || tokenType==MaskToken::AtomNum)
          tokenString.assign( buffer.begin()+1, buffer.end() );
        else
          tokenString.assign( buffer.begin()+2, buffer.end() );
        if (tokenString.empty()) {
          mprinterr("Error: empty token for '%c'\n",buffer[0]);
          return 1;
        }
        // DEBUG
        //mprintf("DEBUG: buffer=[%s]  tokenString=[%s]\n",buffer.c_str(),tokenString.c_str());
        // Split operand by comma
        ArgList commaList(tokenString, ",");
        //commaList.PrintList();
        // Assign each comma-separated arg to a new token
        for (ArgList::const_iterator arg = commaList.begin(); arg != commaList.end(); ++arg) {
          if (token.SetToken( tokenType, *arg ))
            return 1;
          maskTokens_.push_back( token );
        }
        // Indicate that after last token is processed the resulting mask should 
        // go on the stack.
        maskTokens_.back().SetOnStack();
      }
    // operand is a part inside [...]
    } else if ( IsOperand( *p ) || *p == ':' || *p == '@' || *p == '<' || *p == '>' ) {
      buffer += *p;
    // Operators
    } else if (*p == '!' ) {
      token.SetOperator( MaskToken::OP_NEG );
      maskTokens_.push_back( token );
    } else if (*p == '&' ) {
      token.SetOperator( MaskToken::OP_AND );
      maskTokens_.push_back( token );
    } else if (*p == '|' ) {
      token.SetOperator( MaskToken::OP_OR );
      maskTokens_.push_back( token );
    // Distance operator; No longer used, operand is the operator
    } else if (*p == '<' || *p == '>') {
      continue;
    } else {
      mprinterr("Error: Unknown symbol while evaluating mask (%c)\n",*p);
      maskTokens_.clear();
      return 1;
    }
  } // END loop over postfix
  // Test that operators will work correctly.
  if (!maskTokens_.empty()) {
    std::stack<char> tempStack;
    bool validMask = true;
    std::vector<MaskToken>::const_iterator T = maskTokens_.begin();
    for ( ; T != maskTokens_.end(); ++T)
    {
      if (T->Type() == MaskToken::OP_AND ||
          T->Type() == MaskToken::OP_OR) // Requires 2 operands
      {
        if (tempStack.empty()) { validMask = false; break; }
        tempStack.pop();
        if (tempStack.empty()) { validMask = false; break; }
      } else if (T->Type() == MaskToken::OP_NEG  ||
                 T->Type() == MaskToken::OP_DIST) // Requires 1 operand
      {
        if (tempStack.empty()) { validMask = false; break; }
      }
      if ( T->OnStack() )
        tempStack.push('T');
    } 
    if ( !validMask ) {
      mprinterr("Error: Misplaced operator %s.\n", T->TypeName());
      maskTokens_.clear();
      return 1;
    } 
  }
  if (debug_ > 0)
    for (std::vector<MaskToken>::const_iterator T = maskTokens_.begin(); 
                                                T != maskTokens_.end(); T++)
      T->Print();

  return 0;
}

// AtomMask::ResetMask()
void AtomMask::ResetMask() {
  nselected_ = 0;
  Natom_ = 0;
  maskChar_ = 'T';
  maskString_.clear();
  Selected_.clear();
  CharMask_.clear();
  maskTokens_.clear();
}

// AtomMask::ClearSelected()
void AtomMask::ClearSelected() {
  nselected_ = 0;
  Selected_.clear();
  if (!CharMask_.empty()) {
    if (maskChar_ == 'T') CharMask_.assign(Natom_, 'F');
    else                  CharMask_.assign(Natom_, 'T');
  }
}

// AtomMask::InvertMask()
/** Currently for integer masks only. Reverse the selected mask char for 
  * next selection. By default atoms in the Selected_ array are those marked 
  * by 'T'. After a call to InvertMask the atoms in the Selected_ array will
  * be those marked by 'F'. Subsqeuent calls flip the character back and 
  * forth.
  */
void AtomMask::InvertMask() {
  if (maskChar_ == 'T')
    maskChar_ = 'F';
  else
    maskChar_ = 'T';
  if (!Selected_.empty()) {
    // Invert the integer mask.
    std::vector<int> invert;
    invert.reserve( Natom_ - nselected_ );
    const_iterator selected_atom = Selected_.begin();
    for (int idx = 0; idx < Natom_; idx++) {
      if (idx == *selected_atom) // Atom was selected, ignore.
        ++selected_atom;
      else                       // Atom was not selected, add.
        invert.push_back( idx );
    }
    Selected_ = invert;
    nselected_ = (int)Selected_.size();
  }
  if (!CharMask_.empty()) {
    // Invert the character mask.
    for (std::vector<char>::iterator atchar = CharMask_.begin();
                                     atchar != CharMask_.end(); ++atchar)
      if ( *atchar == 'T' )
        *atchar = 'F';
      else
        *atchar = 'T';
    nselected_ = Natom_ - nselected_;
  }
}

// AtomMask::NumAtomsInCommon()
/** Given an atom mask, determine how many selected atoms this mask
  * has in common with it.
  */
int AtomMask::NumAtomsInCommon(AtomMask const& maskIn) {
  std::vector<int> intersect;
  std::vector<int>::iterator intersect_end;

  if (Selected_.empty() || maskIn.Selected_.empty()) return 0;
  // Max size of the intersection is the min size of either array
  intersect.resize( Selected_.size() );
  // Create copies of arrays so they can be sorted
  std::vector<int> selected_1 = Selected_;
  std::vector<int> selected_2 = maskIn.Selected_;
  // Sort the arrays
  std::sort(selected_1.begin(), selected_1.end());
  std::sort(selected_2.begin(), selected_2.end());
  // Set intersect to the intersection of selected_1 and selectd_2
  intersect_end = std::set_intersection(selected_1.begin(), selected_1.end(),
                                   selected_2.begin(), selected_2.end(),
                                   intersect.begin());
  // DEBUG:
  //mprintf("DBG:\tIntersection of [%s] and [%s] is:",maskString_.c_str(),maskIn.maskString_.c_str());
  //for ( std::vector<int>::iterator atom = intersect.begin();
  //                                 atom != intersect_end;
  //                                 atom++)
  //  mprintf(" %i",*atom);
  //mprintf("\n");
  return int(intersect_end - intersect.begin());
}

// AtomMask::AddAtom()
/** Attempt to enforce some sorting by looking for the atom in the mask;
  * as soon as an atom # is found larger than atomIn, insert it at the
  * previous spot.
  */
// 32 33 34 48 49 50
void AtomMask::AddAtom(int atomIn) {
  // Ensure atom is not already in mask
  for (std::vector<int>::iterator atom = Selected_.begin(); atom != Selected_.end(); atom++) {
    if ( *atom == atomIn) return;
    if ( *atom > atomIn) {
      // Insert at the current position, which is the first atom # > atomIn
      Selected_.insert(atom, atomIn);
      nselected_ = (int) Selected_.size();
      return;
    }
  }

  // Add atom to mask
  Selected_.push_back(atomIn);
  nselected_ = (int) Selected_.size();
}

// AtomMask::AddAtoms()
/** Given an array, add the atom numbers in array to the Selected_ array.
  * The resulting array is sorted and any duplicates are removed.
  */
void AtomMask::AddAtoms(std::vector<int> const& atomsIn) {
  std::vector<int>::const_iterator atom;
  // Make room for atomsIn in Selected_
  //Selected_.reserve( Selected_.size() + atomsIn.size() );
  // Put every atom in atomsIn in Selected_ array
  for (atom = atomsIn.begin(); atom != atomsIn.end(); atom++) 
    Selected_.push_back( *atom ); 
  // Sort Selected_
  std::sort( Selected_.begin(), Selected_.end() );
  // Remove duplicates
  atom = unique( Selected_.begin(), Selected_.end() );
  Selected_.resize( atom - Selected_.begin() );
  nselected_ = (int) Selected_.size();
}

// AtomMask::AddAtomRange()
/** Add atoms in range from minAtom up to but not including maxAtom to 
  * Selected_ array. The resulting array is sorted and duplicates are removed.
  */
void AtomMask::AddAtomRange(int minAtom, int maxAtom) {
  //mprintf("DEBUG:\t\tAdding atoms %i to %i\n",minAtom,maxAtom);
  if (minAtom >= maxAtom) return;
  for (int atom = minAtom; atom < maxAtom; atom++)
    Selected_.push_back( atom );
  // Sort Selected_
  std::sort( Selected_.begin(), Selected_.end() );
  // Remove duplicates
  std::vector<int>::iterator atomit = unique( Selected_.begin(), Selected_.end() );
  Selected_.resize( atomit - Selected_.begin() );
  //mprintf("\t\t[");
  //for (std::vector<int>::iterator da = Selected_.begin(); da != Selected_.end(); da++)
  //  mprintf(" %i",*da);
  //mprintf("]\n");
  nselected_ = (int) Selected_.size();
}

// AtomMask::AddMaskAtPosition()
/** Add the atoms in Mask to this mask starting at the specified positon.
  * Currently only used by Closest for modifying the original stripped
  * mask with the calculated closest waters.
  */
void AtomMask::AddMaskAtPosition(AtomMask const& maskIn, int idx) {
  // NOTE: NO BOUNDS CHECK!
  for (const_iterator atom = maskIn.begin(); atom != maskIn.end(); atom++) 
    Selected_[idx++] = *atom;
}

// AtomMask::PrintMaskAtoms()
void AtomMask::PrintMaskAtoms(const char *header) const {
  mprintf("%s=",header);
  if (this->None()) 
    mprintf("No atoms selected.");
  else if (!Selected_.empty()) {
    for (std::vector<int>::const_iterator atom = Selected_.begin(); 
                                          atom != Selected_.end(); ++atom)
      mprintf(" %i",*atom + 1);
  } else if (!CharMask_.empty()) {
    int atomnum = 0;
    for (std::vector<char>::const_iterator mc = CharMask_.begin(); 
                                           mc != CharMask_.end(); mc++) {
      if (*mc == maskChar_) mprintf(" %i",atomnum);
      ++atomnum;
    }
  } else 
    mprintf("Warning: Mask [%s] has not been set up yet.\n",maskString_.c_str());
  mprintf("\n");
}

// AtomMask::SetMaskString()
/** Take the given mask expression and preprocess it for subsequent use
  * with the mask parser. Convert to infix, then postfix notation.
  */
int AtomMask::SetMaskString(const char* maskString_In) {
  if (maskString_In!=0) 
    maskString_.assign( maskString_In );
  else
    maskString_.assign( "*" );

  if (debug_ > 0) mprintf("expression: ==%s==\n", maskString_.c_str());

  // Convert mask expression to maskTokens 
  if (Tokenize()) return 1;

  return 0;
}

// AtomMask::SetMaskString()
/** If the input string is not empty, set the AtomMask expression to 
  * input string.
  * \return 0 if string was set and successfully tokenized.
  * \return 1 if tokenization failed.
  */
int AtomMask::SetMaskString( std::string const& maskStringIn ) {
  if (!maskStringIn.empty())
    return ( SetMaskString( maskStringIn.c_str() ) );
  else
    return ( SetMaskString( 0 ) );
}

// AtomMask::SetupIntMask()
/** Set up an atom mask containing selected atom numbers given a char
  * array of size natom with T for selected atoms and F for unselected
  * atoms. The actual mask parser is called from Topology. 
  * maskChar_ is used to determine whether atoms denoted by 'T' or 'F' will
  * be selected (the latter is the case e.g. with stripped atoms). 
  */
void AtomMask::SetupIntMask(const char *charmask,int natom,int debugIn) {
  debug_ = debugIn;
  // Wipe out previous char mask if allocated
  CharMask_.clear();

  // Set up an integer list of the selected atoms. 
  // NOTE: For large selections this will use 4x the memory of the char atom
  //       mask. Could check to see which will be bigger.
  nselected_=0;
  Natom_ = natom;
  Selected_.clear();
  for (int atom=0; atom < natom; atom++) {
    if (charmask[atom]==maskChar_) {
      Selected_.push_back( atom );
      ++nselected_;
    }
  }

  if (debug_>0) {
    if (maskChar_=='F')
      mprintf("          Inverse of Mask %s corresponds to %i atoms.\n",
              maskString_.c_str(), natom - nselected_);
    else
      mprintf("          Mask %s corresponds to %i atoms.\n",maskString_.c_str(),
              nselected_);
  }
}

// AtomMask::SetupCharMask()
/** Given an input char mask of size natom, set up a corresponding char mask.
  * Useful for cases where we need to know both atoms in and out of mask.
  */
void AtomMask::SetupCharMask(const char *charmask, int natom, int debugIn) {
  debug_ = debugIn;
  // Wipe out previous Selected_ mask if allocated
  Selected_.clear();

  // Allocate atom mask - free mask if already allocated
  CharMask_.clear();

  nselected_=0;
  Natom_ = natom;
  CharMask_.reserve( natom );
  for (int i = 0; i < natom; i++) {
    CharMask_.push_back( charmask[i] );
    // Determine number of selected atoms
    if (charmask[i] == maskChar_) ++nselected_;
  }
}

// AtomMask::AtomInCharMask()
bool AtomMask::AtomInCharMask(int atom) const {
  if (CharMask_.empty()) return false;
  if (atom < 0) return false;
  if (atom >= (int)CharMask_.size()) return false;
  if (CharMask_[atom]==maskChar_) return true;
  return false;
}

bool AtomMask::AtomsInCharMask(int startatom, int endatom) const {
  if (CharMask_.empty()) return false;
  if (startatom > endatom) return false;
  if (startatom < 0) return false;
  if (endatom > (int)CharMask_.size()) return false;
  for (int idx = startatom; idx < endatom; ++idx)
    if (CharMask_[idx] == maskChar_) return true;
  return false;
}

// AtomMask::ConvertToCharMask()
int AtomMask::ConvertToCharMask() {
  if (Selected_.empty()) return 0;
  // If Natom_ is 0 (which can happen if mask was set up with AddAtomX
  // functions) this cant work.
  if (Natom_==0) {
    mprinterr("Error: ConvertToCharMask(): Natom for integer mask is 0.\n");
    return 1;
  }
  CharMask_.assign( Natom_, 'F' );
  for (std::vector<int>::iterator maskatom = Selected_.begin();
                                  maskatom != Selected_.end(); ++maskatom)
    CharMask_[*maskatom]='T';
  Selected_.clear();
  return 0;
}

// AtomMask::ConvertToIntMask()
int AtomMask::ConvertToIntMask() {
  if (CharMask_.empty()) return 0;
  Selected_.reserve( nselected_ );
  for (int atom = 0; atom < Natom_; atom++) {
    if (CharMask_[atom]=='T')
      Selected_.push_back( atom );
  }
  CharMask_.clear();
  return 0;
}

void AtomMask::MaskInfo() const {
  mprintf("\tMask [%s] corresponds to %i atoms.\n", maskString_.c_str(), nselected_);
}

void AtomMask::BriefMaskInfo() const {
  mprintf(" [%s](%i)", maskString_.c_str(), nselected_);
}
