#include "MapAtom.h"
// TODO: Just set the char to whatever the atomic weight is?
const char MapAtom::AtomicElementChar[Atom::NUMELEMENTS] = { 0,
  'H',  'B',  'C',  'N',  'O', 'F',  
  'P',  'S',  'X',  'Y',  'f', 'c',
  'I',  'M',  'U',  'L',  'K', 'R',  
  'E',  'Z',  'n',    0,    0,   0,
    0,    0,    0,    0,    0,   0,
    0,    0,    0,    0,    0,   0,
    0,    0,    0,    0,    0,   0,
    0,    0,    0,    0,    0,   0,
    0,    0,    0,    0,    0,   0,
    0,    0,    0,    0,    0,   0,
    0,    0,    0,    0,    0,   0,
    0,    0,    0,    0,    0,   0,
    0,    0,
    0
};

/// CONSTRUCTOR
MapAtom::MapAtom() :
  isChiral_(false),
  boundToChiral_(false), 
  isMapped_(false), 
  complete_(false), 
  Nduplicated_(0),
  name_(0) 
{}

// COPY CONSTRUCTOR
MapAtom::MapAtom(const MapAtom& rhs) : Atom(rhs),
   isChiral_(rhs.isChiral_),
   boundToChiral_(rhs.boundToChiral_), 
   isMapped_(rhs.isMapped_), 
   complete_(rhs.complete_),
   atomID_(rhs.atomID_), 
   unique_(rhs.unique_), 
   Nduplicated_(rhs.Nduplicated_),
   name_(rhs.name_) 
{}

// COPY CONSTRUCTOR
/// Copy base atom to this MapAtom
MapAtom::MapAtom(const Atom& rhs) : Atom(rhs),
  isChiral_(false),
  boundToChiral_(false),
  isMapped_(false),
  complete_(false),
  Nduplicated_(0),
  name_(AtomicElementChar[Element()])
{}

// Assignment
MapAtom& MapAtom::operator=(const MapAtom& rhs) {
  if (&rhs == this) return *this;
  Atom::operator=(rhs);
  isChiral_ = rhs.isChiral_;
  boundToChiral_ = rhs.boundToChiral_;
  isMapped_ = rhs.isMapped_;
  complete_ = rhs.complete_;
  atomID_   = rhs.atomID_;
  unique_   = rhs.unique_;
  Nduplicated_ = rhs.Nduplicated_;
  name_ = rhs.name_;
  return *this;
}
