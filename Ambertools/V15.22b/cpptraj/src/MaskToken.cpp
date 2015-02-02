#include <locale>
#include "MaskToken.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // convertToInteger, convertToDouble

MaskToken::MaskToken() :
  type_(OP_NONE),
  res1_(-1),
  res2_(-1),
  name_(""),
  onStack_(false),
  d_within_(false),
  d_atom_(false),
  distance_(0.0)
{ }

const char* MaskToken::MaskTypeString[] = {
  "OP_NONE", "ResNum", "ResName", "AtomNum", "AtomName", "AtomType", "AtomElement", "SelectAll",
  "OP_AND", "OP_OR", "OP_NEG", "OP_DIST"
};

void MaskToken::Print() const {
  mprintf("TOKEN: [%s]",MaskTypeString[type_]);
  switch (type_) {
    case ResName:
    case AtomName: mprintf(" Name=[%s]",*name_); break;
    case ResNum:
    case AtomNum: mprintf(" First=%i  Second=%i",res1_,res2_); break;
    case OP_DIST: 
      mprintf(" within=%i  d_atom=%i  distance^2=%lf",
              (int)d_within_, (int)d_atom_, distance_);
      break;
    default: mprintf(" ");
  }
  mprintf(" OnStack=%i\n",(int)onStack_);
/*
  mprintf("TOKEN: [%s] Res1=%i  Res2=%i  Name=[%s]  OnStack=%i\n",
          MaskTypeString[type_], res1_, res2_, *name_, (int)onStack_);
  mprintf("            within=%i  d_atom=%i  distance^2=%lf\n",
          (int)d_within_, (int)d_atom_, distance_);*/
} 

const char *MaskToken::TypeName() const {
  return MaskTypeString[type_];
}

void MaskToken::MakeNameType() {
  if (type_ == ResNum)
    type_ = ResName;
  else if (type_ == AtomNum)
    type_ = AtomName;
}

// Basic : or @ operand
int MaskToken::SetToken( MaskTokenType typeIn, std::string const& tokenString ) {
  std::locale loc;
  if (tokenString.empty()) return 1;
  // Set initial token type
  type_ = typeIn;
  onStack_ = false;
  // Does this token argument have an asterisk? If its at the first position
  // make this an ALL token, otherwise make this a name token.
  size_t asteriskPosition = tokenString.find_first_of("*");
  if ( asteriskPosition != std::string::npos ) {
    if (asteriskPosition == 0) {
      type_ = SelectAll;
      return 0;
    } else {
      MakeNameType();
    }
  }
  // Check that all chars are digits or - for number range 
  if (type_ == ResNum || type_ == AtomNum) {
    for (std::string::const_iterator p = tokenString.begin(); p != tokenString.end(); ++p) {
      if (*p != '-' && isalpha(*p, loc)) {
        //mprintf("DEBUG: making name type because of %c\n",*p);
        MakeNameType();
        break;
      } 
    }
  }
  if (type_ == ResNum || type_ == AtomNum) {
    // Does this token argument have a dash? Only valid for number ranges.
    size_t dashPosition = tokenString.find_first_of("-");
    if (dashPosition != std::string::npos) {
      // Get first and second args. If first arg is blank negative number specified.
      std::string arg1(tokenString.begin(), tokenString.begin()+dashPosition);
      if (arg1.empty()) {
        mprinterr("Error: Mask expressions cannot contain negative numbers (%s)\n",
                  tokenString.c_str());
        return 1;
      }
      std::string arg2(tokenString.begin()+dashPosition+1, tokenString.end());
      if (arg2.empty()) {
        mprinterr("Error: Incomplete number range given (%s).\n", tokenString.c_str());
        return 1;
      }
      res1_ = convertToInteger( arg1 );
      res2_ = convertToInteger( arg2 );
    } else {
      // Get the number arg
      res1_ = convertToInteger( tokenString );
      res2_ = res1_;
    }
    // Ensure that res1 and res2 are valid
    if (res2_ < res1_) {
      mprinterr("Error: Mask range, second num (%i) less than first (%i).\n",res2_,res1_);
      return 1;
    }
    // It is expected that number args will start from 1
    if (res1_ < 1 || res2_ < 1) {
      mprinterr("Error: One or both numbers of mask arg (%s) < 1 (%i, %i)\n",
                tokenString.c_str(), res1_,res2_);
      return 1;
    }
  } else {
    // This is a string arg.
    // Use AssignNoFormat so that * not converte to '
    //name_.AssignNoFormat( tokenString.c_str() ); // TODO: Convert directly from string
    name_ = tokenString;
  }
  return 0;
}

// [<|>][@|:]<dist>
int MaskToken::SetDistance(std::string &distop) {
  if (distop.empty()) return 1;
  type_ = OP_DIST;
  onStack_ = false;
  // Min size is 3 chars
  if (distop.size() < 3) {
    mprinterr("Error: Malformed distance operator [%s]\n",distop.c_str());
    return 1;
  }
  // 1st char indicates within (<) or without (>)
  if (distop[0]=='<')
    d_within_ = true;
  else if (distop[0]=='>')
    d_within_ = false;
  else {
    mprinterr("Error: Malformed distance operator: expected '<' or '>' (%c)\n",distop[0]);
    return 1;
  }
  // 2nd char indidcates atoms (@) or residues (:)
  if (distop[1]=='@')
    d_atom_ = true;
  else if (distop[1]==':')
    d_atom_ = false;
  else {
    mprinterr("Error: Malformed distance operator: expected ':' or '@' (%c)\n",distop[1]);
    return 1;
  }
  // 3rd char onwards is the distance argument
  std::string distarg(distop.begin()+2, distop.end());
  distance_ = convertToDouble(distarg);
  // Pre-square the distance
  distance_ *= distance_;
  return 0;
}

void MaskToken::SetOperator(MaskTokenType typeIn) {
  type_ = typeIn;
  onStack_ = false;
}

void MaskToken::SetOnStack() {
  onStack_ = true;
}

