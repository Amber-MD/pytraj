#include <cctype> // isalpha
#include <algorithm> // sort
#include "Atom.h"
#include "CpptrajStdio.h"

const int Atom::AtomicElementNum[NUMELEMENTS] = { 0,
 1,  5,  6,  7,  8,  9,  
 15, 16, 17, 35, 26, 20, 
 53, 12, 29, 3,  19, 37, 
 55, 30, 11, 13, 18, 33,
 47, 79, 85, 4,  56, 83,
 24, 27, 48, 87, 31, 32,
  2, 72, 80, 49, 77, 36,
 25, 42, 10, 28, 41, 76,
 46, 78, 82, 84, 44, 45,
 75, 86, 88, 14, 21, 34,
 38, 50, 51, 22, 43, 52,
 73, 81, 23, 74, 54, 40,
 39, 71,
 0
};

/// Atom names corresponding to AtomicElementType.
// 2 chars + null.
const char* Atom::AtomicElementName[NUMELEMENTS] = { "??",
  "H",  "B",  "C",  "N",  "O",   "F",  
  "P",  "S",  "CL", "BR", "FE", "CA",
  "I",  "MG", "CU", "LI", "K",  "RB", 
  "CS", "ZN", "NA", "AL", "AR", "AS",
  "AG", "AU", "AT", "BE", "BA", "BI",
  "CR", "CO", "CD", "FR", "GA", "GE",
  "HE", "HF", "HG", "IN", "IR", "KR",
  "MN", "MO", "NE", "NI", "NB", "OS",
  "PD", "PT", "PB", "PO", "RU", "RH",
  "RE", "RN", "RA", "SI", "SC", "SE",
  "SR", "SN", "SB", "TI", "TC", "TE",
  "TA", "TL", "V",  "W",  "XE", "ZR",
  "Y",  "LU",
  "XP"
};

/** Values taken from 'http://www.webelements.com/' */
const double Atom::AtomicElementMass[NUMELEMENTS] = { 1.0,
    1.00794,   10.811,     12.0107,     14.0067,    15.9994,   18.9984032,
   30.973762,  32.065,     35.453,      79.904,     55.845,    40.078,
  126.90447,   24.3050,    63.546,       6.941,     39.0983,   85.4678,
  132.9054519, 65.38,      22.98976928, 26.9815386, 39.948,    74.92160,
  107.8682,   196.966569, 210,           9.012182, 137.327,   208.98040,
   51.9961,    58.933195, 112.411,     223,         69.723,    72.64,
    4.002602, 178.49,     200.59,      114.818,    192.217,    83.798,
   54.938045,  95.96,      20.1797,     58.6934,    92.90638, 190.23,
  106.42,     195.084,    207.2,       209,        101.07,    102.90550,
  186.207,    222,        226,          28.0855,    44.955912, 78.96,
   87.62,     118.710,    121.760,      47.867,     98,       127.60,
  180.94788,  204.3833,    50.9415,    183.84,     131.293,    91.224,
   88.90585,  174.9668,
    0.0
};  

// CONSTRUCTOR
Atom::Atom() : charge_(0.0), polar_(0.0), mass_(1.0), gb_radius_(0.0),
               gb_screen_(0.0), aname_(""), atype_(""), atype_index_(0),
               element_(UNKNOWN_ELEMENT), resnum_(0), mol_(0), chainID_(' ')
{}

/** If no 2 char element name provided, try to guess from name. */ 
Atom::Atom(NameType const& aname, char cid, const char* elt) :
  charge_(0.0), polar_(0.0), mass_(1.0), gb_radius_(0.0), gb_screen_(0.0),
  aname_(aname), atype_(""), atype_index_(0), element_(UNKNOWN_ELEMENT),
  resnum_(0), mol_(0), chainID_(cid)
{
  if (elt ==0 || (elt[0]==' ' && elt[1]==' '))
    SetElementFromName();
  else
    SetElementFromSymbol(elt[0], elt[1]);
  mass_ = AtomicElementMass[ element_ ];
}

/** Attempt to guess element from atom name. */ 
Atom::Atom( NameType const& aname, NameType const& atype, double q ) :
  charge_(q), polar_(0.0), mass_(1.0), gb_radius_(0.0), gb_screen_(0.0),
  aname_(aname), atype_(atype), atype_index_(0), element_(UNKNOWN_ELEMENT),
  resnum_(0), mol_(0), chainID_(' ')
{
  SetElementFromName();
  mass_ = AtomicElementMass[ element_ ];
}

/** Set type and type index. Attempt to guess element from atom name. */
Atom::Atom( NameType const& aname, NameType const& atype, int atidx ) :
  charge_(0.0), polar_(0.0), mass_(1.0), gb_radius_(0.0), gb_screen_(0.0),
  aname_(aname), atype_(atype), atype_index_(atidx), element_(UNKNOWN_ELEMENT),
  resnum_(0), mol_(0), chainID_(' ')
{
  SetElementFromName();
  mass_ = AtomicElementMass[ element_ ];
}

// Atom::DetermineElement()
/// Determine element from atomic number, mass, or atom name.
void Atom::DetermineElement(int atomicnum) {
  if (atomicnum>0) {
    // Determine element from atomic number
    for (int i = 1; i < (int)NUMELEMENTS; i++)
      if (AtomicElementNum[i] == atomicnum) {
        element_ = (AtomicElementType) i;
        break;
      }
  } else {
    // Determine element from mass. No mass means extra point.
    if ( mass_ == 0 )
      element_ = EXTRAPT;
    else
      SetElementFromMass();
  }
  // If element still unknown attempt to determine from name
  if (element_ == UNKNOWN_ELEMENT)
    SetElementFromName();
}

// CONSTRUCTOR
Atom::Atom(NameType const& aname, double charge, double mass, NameType const& atype) :
  charge_(charge), polar_(0.0), mass_(mass), gb_radius_(0.0), gb_screen_(0.0),
  aname_(aname), atype_(atype), atype_index_(0), element_(UNKNOWN_ELEMENT),
  resnum_(0), mol_(0), chainID_(' ')
{
  DetermineElement(0);
}

// CONSTRUCTOR
Atom::Atom( NameType const& name, double charge, double polar, int atomicnum, 
            double mass, int atidx, NameType const& type, double rad, double screen,
            char cID ) :
  charge_(charge), polar_(polar), mass_(mass), gb_radius_(rad), gb_screen_(screen),
  aname_(name), atype_(type), atype_index_(atidx), element_(UNKNOWN_ELEMENT),
  resnum_(0), mol_(0), chainID_(cID)
{
  DetermineElement(atomicnum);
}

// COPY CONSTRUCTOR
Atom::Atom(const Atom &rhs) :
  charge_(rhs.charge_),
  polar_(rhs.polar_),
  mass_(rhs.mass_),
  gb_radius_(rhs.gb_radius_),
  gb_screen_(rhs.gb_screen_),
  aname_(rhs.aname_),
  atype_(rhs.atype_),
  atype_index_(rhs.atype_index_),
  element_(rhs.element_),
  resnum_(rhs.resnum_),
  mol_(rhs.mol_),
  chainID_(rhs.chainID_),
  bonds_(rhs.bonds_),
  excluded_(rhs.excluded_)
{ }

// SWAP
void Atom::swap(Atom &first, Atom &second) {
  using std::swap;
  swap(first.charge_, second.charge_);
  swap(first.polar_, second.polar_);
  swap(first.mass_, second.mass_);
  swap(first.gb_radius_, second.gb_radius_);
  swap(first.gb_screen_, second.gb_screen_);
  swap(first.aname_, second.aname_);
  swap(first.atype_, second.atype_);
  swap(first.atype_index_, second.atype_index_);
  swap(first.element_, second.element_);
  swap(first.resnum_, second.resnum_);
  swap(first.mol_, second.mol_);
  swap(first.chainID_, second.chainID_);
  swap(first.bonds_, second.bonds_);
  swap(first.excluded_, second.excluded_);
}

// ASSIGNMENT via copy/swap idiom
Atom &Atom::operator=(Atom other) {
  swap(*this, other);
  return *this;
}

// Atom::AddBond()
void Atom::AddBond(int idxIn) {
  bonds_.push_back( idxIn );
}

// Atom::ClearBonds()
void Atom::ClearBonds() {
  bonds_.clear();
}

// Atom::SortBonds()
void Atom::SortBonds() {
  sort( bonds_.begin(), bonds_.end() );
}

// Atom::IsBondedTo()
bool Atom::IsBondedTo(int atomnum) const {
  for (std::vector<int>::const_iterator b = bonds_.begin(); b != bonds_.end(); ++b)
    if ( atomnum == *b ) return true;
  return false;
} 

// Atom::AddExclusionList()
void Atom::AddExclusionList(std::set<int> const& elist) {
  excluded_.clear();
  for (std::set<int>::const_iterator ei = elist.begin(); ei != elist.end(); ei++)
    excluded_.push_back( *ei );
}

// Atom::SetElementFromName()
/** If not already known, try to determine atomic element from atom name. 
  * Based on Amber standard atom names.
  */
void Atom::SetElementFromName() {
  if (element_!=UNKNOWN_ELEMENT) return;
  // Seek to first alpha
  const char* ptr = *aname_;
  while (*ptr != '\0' && !isalpha(*ptr)) ptr++;
  if (*ptr == '\0') return;
  char c1 = *ptr;
  char c2 = *(ptr+1);  
  switch (c1) {
    case 'H' : element_ = HYDROGEN; break;
    case 'C' :
      if (c2=='l' || c2=='L') element_ = CHLORINE;
      else if (c2=='0') element_ = CALCIUM;
      else if (c2=='s') element_ = CESIUM;
      else if (c2=='U') element_ = COPPER;
      else element_ = CARBON;
      break;
    case 'N' :
      if (c2=='a') element_ = SODIUM;
      else element_ = NITROGEN;
      break;
    case 'O' : element_ = OXYGEN; break;
    case 'F' :
      if (c2=='e' || c2=='E') element_ = IRON;
      else element_ = FLUORINE;
      break;
    case 'B' :
      if (c2=='r' || c2=='R') element_ = BROMINE;
      else element_ = BORON;
      break;
    case 'I' :
      if (c2=='M') element_ = CHLORINE;    // IM, Cl- in old FF94
      else if (c2=='P') element_ = SODIUM; // IP, Na+ in old FF94
      else element_ = IODINE;
      break;
    case 'P' : element_ = PHOSPHORUS; break;
    case 'S' : element_ = SULFUR; break;
    case 'M' :
      if (c2=='g' || c2=='G') element_ = MAGNESIUM;
      if (c2=='n' || c2 =='N') element_ = MANGANESE;
      break;
    case 'Z' :
      if (c2=='n' || c2 =='N') element_ = ZINC;
      break;
    case 'L' :
      if (c2=='i') element_ = LITHIUM;
      break;
    case 'K' : element_ = POTASSIUM; break;
    case 'R' :
      if (c2=='b') element_ = RUBIDIUM;
      break;
    default:
      SetElementFromSymbol(c1, c2);
      if (element_ == UNKNOWN_ELEMENT)
        mprintf("Warning: Could not determine atomic number from name [%s]\n",*aname_);
  }
}

// Atom::SetElementFromSymbol()
void Atom::SetElementFromSymbol(char c1, char c2) {
  if (element_ != UNKNOWN_ELEMENT) return;
  // Attempt to match up 1 or 2 char element name
  char en[2];
  bool oneChar = true;
  if (c1 != ' ' && c2 != ' ') {
    en[0] = toupper(c1);
    en[1] = toupper(c2);
    oneChar = false;
  } else if (c2 != ' ') {
    en[0] = toupper(c2);
    en[1] = ' ';
  } else if (c1 != ' ') {
    en[0] = toupper(c1);
    en[1] = ' ';
  } else
    return; // sanity check, both blank
  if (oneChar) {
    for (int i = 1; i < (int)NUMELEMENTS; i++)
      if (AtomicElementName[i][1]=='\0' && // 1 char
          en[0] == AtomicElementName[i][0])
      {
        element_ = (AtomicElementType)i;
        break;
      }
  } else {
    for (int i = 1; i < (int)NUMELEMENTS; i++)
      if (AtomicElementName[i][1]!='\0' && // 2 char
          en[0] == AtomicElementName[i][0] &&
          en[1] == AtomicElementName[i][1])
      {
        element_ = (AtomicElementType)i;
        break;
      }
  }
}

// Atom::SetElementFromMass()
/** Determine the atomic element from the first letter of the atom name
  * and the atomic mass. Based on the get_atomic_number() subroutine found
  * in '$AMBERHOME/AmberTools/src/sqm/qmmm_module'.
  */
void Atom::SetElementFromMass() {
  char c1 = aname_[0];
  switch (c1) {
    case 'a':
    case 'A':
      if (mass_ > 24.0 && mass_ <= 28.0)
        element_ = ALUMINUM; // 13 ! Aluminum
      else if (mass_ > 35.0 && mass_ <= 40.0) 
          element_ = ARGON; // 18 !Argon
      else if(mass_ > 73.0 && mass_ <= 77.0) 
          element_ = ARSENIC; // 33 !Arsenic
      else if(mass_ > 106.0 && mass_ <= 109.0) 
          element_ = SILVER; // 47 !Silver
      else if(mass_ > 195.0 && mass_ <= 199.0) 
          element_ = GOLD; // 79 !Gold
      else if(mass_ > 208.0 && mass_ <= 212.0) 
          element_ = ASTATINE; // 85 !Astatine
      break;
   
    case 'b':
    case 'B':
      if (mass_ > 8.0 && mass_ <= 10.0) 
          element_ = BERYLLIUM; // 4 !Beryllium
      else if(mass_ > 10.0 && mass_ <= 12.0) 
          element_ = BORON; // 5 !Boron
      else if(mass_ > 77.0 && mass_ <= 81.0) 
          element_ = BROMINE; // 35 !Bromine
      else if(mass_ > 135.0 && mass_ <= 139.0) 
          element_ = BARIUM; // 56 !Barium
      else if(mass_ > 207.0 && mass_ <= 211.0) 
          element_ = BISMUTH; // 83 !Bismuth
      break;

    case 'c':
    case 'C':
       if (mass_ > 10.0 && mass_ <= 14.0) 
          element_ = CARBON; //=  6 !Carbon
       else if(mass_ > 33.0 && mass_ <= 37.0) 
          element_ = CHLORINE; // 17 !Chlorine
       else if(mass_ > 38.0 && mass_ <= 42.0) 
          element_ = CALCIUM; // 20 !Calcium
       else if(mass_ > 50.0 && mass_ <= 54.0) 
          element_ = CHROMIUM; // 24 !Chromium
       else if(mass_ > 57.0 && mass_ <= 61.0) 
          element_ = COBALT; // 27 !Cobalt
       else if(mass_ > 61.0 && mass_ <= 65.0) 
          element_ = COPPER; // 29 !Copper
       else if(mass_ > 110.0 && mass_ <= 114.0) 
          element_ = CADMIUM; // 48 !Cadmium
       else if(mass_ > 131.0 && mass_ <= 135.0) 
          element_ = CESIUM; // 55 !Cesium
       break; 

    case 'f':
    case 'F':
       if (mass_ > 17.0 && mass_ <= 21.0) 
          element_ = FLUORINE; // 9 !Fluorine
       else if(mass_ > 54.0 && mass_ <= 58.0) 
          element_ = IRON; // 26 !Iron
       else if(mass_ > 218.0 && mass_ <= 228.0) 
          element_ = FRANCIUM; // 87 !Francium
       break;

    case 'g':
    case 'G':
       if (mass_ > 67.0 && mass_ <= 71.0) 
          element_ = GALLIUM; //  31 !Gallium
       else if(mass_ > 71.0 && mass_ <= 75.0) 
          element_ = GERMANIUM; // 32 !Germanium
       break;

    case 'h':
    case 'H':
       if(mass_ > 0.0 && mass_ <= 2.0) 
          element_ = HYDROGEN; // 1 !Hydrogen
       else if(mass_ > 3.0 && mass_ <= 5.0) 
          element_ = HELIUM; // 2 !Helium
       else if(mass_ > 176.0 && mass_ <= 180.0) 
          element_ = HAFNIUM; // 72 !Hafnium
       else if(mass_ > 198.0 && mass_ <= 202.0) 
          element_ = MERCURY; // 80 !Mercury
       break;

    case 'i':
    case 'I':
       if(mass_ > 112.0 && mass_ <= 116.0) 
          element_ = INDIUM; //49 !Indium
       else if(mass_ > 125.0 && mass_ <= 129.0) 
          element_ = IODINE; // 53 !Iodine
       else if(mass_ > 190.0 && mass_ <= 194.0) 
          element_ = IRIDIUM; // 77 !Iridium
       break;

    case 'k':
    case 'K':
       if(mass_ > 37.0 && mass_ <= 41.0) 
          element_ = POTASSIUM; //19 !Potassium
       else if(mass_ > 77.0 && mass_ <= 86.0) 
          element_ = KRYPTON; //36 !Krypton
       break;

    case 'l':
    case 'L':
       if(mass_ > 6.0 && mass_ <= 8.0) 
          element_ = LITHIUM; //3 !Lithium
       if (mass_ > 173.0 && mass_ < 177.0)
          element_ = LUTETIUM; // 71 !Lutetium
       break;

    case 'm':
    case 'M':
       if(mass_ > 22.0 && mass_ <= 26.0) 
          element_ = MAGNESIUM; //12 !Magnesium
       else if(mass_ > 53.0 && mass_ <= 57.0) 
          element_ = MANGANESE; //25 !Manganese
       else if(mass_ > 94.0 && mass_ <= 98.0) 
          element_ = MOLYBDENUM; //42 !Molybdenum
       break;

    case 'n':
    case 'N':
       if(mass_ > 13.0 && mass_ <= 15.0) 
          element_ = NITROGEN; //7 !Nitrogen
       else if(mass_ > 19.0 && mass_ <= 22.0) 
          element_ = NEON; //10 !Neon
       else if(mass_ > 22.1 && mass_ <= 23.0) 
          element_ = SODIUM; //11 !Sodium
       else if(mass_ > 57.0 && mass_ <= 61.0) 
          element_ = NICKEL; //28 !Nickel
       else if(mass_ > 95.0 && mass_ <= 99.0) 
          element_ = NIOBIUM; //41 !Niobium
       break;

    case 'o':
    case 'O':
       if(mass_ > 14.0 && mass_ <= 18.0) 
          element_ = OXYGEN; //8 !Oxygen
       else if(mass_ > 188.0 && mass_ <= 192.0) 
          element_ = OSMIUM; //76 !Osmium
       break;

    case 'p':
    case 'P':
       if(mass_ > 29.0 && mass_ <= 33.0) 
          element_ = PHOSPHORUS; //15 !Phosphorus
       else if(mass_ > 104.0 && mass_ <= 108.0) 
          element_ = PALLADIUM; //46 !Palladium
       else if(mass_ > 193.0 && mass_ <= 197.0) 
          element_ = PLATINUM; //78 !Platinum
       else if(mass_ > 205.0 && mass_ <= 208.0) 
          element_ = LEAD; //82 !Lead
       else if(mass_ > 208.0 && mass_ <= 212.0) 
          element_ = POLONIUM; //84 !Polonium
       break;

    case 'r':
    case 'R':
       if(mass_ > 84.0 && mass_ <= 88.0) 
          element_ = RUBIDIUM; //37 !Rubidium
       else if(mass_ > 99.0 && mass_ <= 102.0) 
          element_ = RUTHENIUM; //44 !Ruthenium
       else if(mass_ > 102.0 && mass_ <= 105.0) 
          element_ = RHODIUM; //45 !Rhodium
       else if(mass_ > 184.0 && mass_ <= 188.0) 
          element_ = RHENIUM; //75 !Rhenium
       else if(mass_ > 210.0 && mass_ <= 222.5) 
          element_ = RADON; //86 !Radon
       else if(mass_ > 223.0 && mass_ <= 229.0) 
          element_ = RADIUM; //88 !Radium
       break;

   case 's':
   case 'S':
       if(mass_ > 26.0 && mass_ <= 30.0) 
          element_ = SILICON; //14 !Silicon
       else if(mass_ > 30.0 && mass_ <= 34.0) 
          element_ = SULFUR; //16 !Sulfur
       else if(mass_ > 43.0 && mass_ <= 47.0) 
          element_ = SCANDIUM; //21 !Scandium
       else if(mass_ > 77.0 && mass_ <= 81.0) 
          element_ = SELENIUM; //34 !Selenium
       else if(mass_ > 86.0 && mass_ <= 89.0) 
          element_ = STRONTIUM; //38 !Strontium
       else if(mass_ > 116.0 && mass_ <= 120.0) 
          element_ = TIN; //50 !Tin
       else if(mass_ > 120.0 && mass_ <= 124.0) 
          element_ = ANTIMONY; //51 !Antimony
       break;

    case 't':
    case 'T':
       if(mass_ > 46.0 && mass_ <= 50.0) 
          element_ = TITANIUM; //22 !Titanium
       else if(mass_ > 96.0 && mass_ <= 100.0) 
          element_ = TECHNETIUM; //43 !Technetium
       else if(mass_ > 125.0 && mass_ <= 130.0) 
          element_ = TELLURIUM; //52 !Tellurium
       else if(mass_ > 179.0 && mass_ <= 183.0) 
          element_ = TANTALUM; //73 !Tantalum
       else if(mass_ > 201.0 && mass_ <= 206.0) 
          element_ = THALLIUM; //81 !Thallium
       break;

    case 'v':
    case 'V':
       if(mass_ > 49.0 && mass_ <= 53.0) 
          element_ = VANADIUM; //23 !Vanadium
       break;

    case 'w':
    case 'W':
       if(mass_ > 179.0 && mass_ <= 183.0) 
          element_ = TUNGSTEN; //74 !Tungsten
       break;

    case 'x':
    case 'X':
       if (mass_ > 127.0 && mass_ < 136.0)
          element_ = XENON; // 54 !Xenon 
       break;

    case 'y':
    case 'Y':
       if (mass_ > 87.0 && mass_ < 91.0)
          element_ = YTTRIUM; // 39 !Yttrium
      break;

    case 'z':
    case 'Z':
       if(mass_ > 61.0 && mass_ <= 69.0) 
          element_ = ZINC; //30 !Zinc
       else if(mass_ > 89.0 && mass_ <= 93.0) 
          element_ = ZIRCONIUM; //40 !Zirconium
       break;

    default:
      mprintf("Warning: Could not determine atomic number from mass (%lf) [%s]\n", 
              mass_, *aname_);
  }
}

// WarnBondLengthDefault()
void Atom::WarnBondLengthDefault(AtomicElementType atom1, AtomicElementType atom2, double cut) {
  mprintf("Warning: Bond length not found for %s - %s, using default= %f\n",
          AtomicElementName[atom1], AtomicElementName[atom2], cut);
}

/** Return optimal covalent bond distance based on the element types of atom1 
  * and atom2. When multiple hybridizations are possible the longest possible 
  * bond length is used.
  * Unless otherwise noted values taken from:
  * - Huheey, pps. A-21 to A-34; T.L. Cottrell, "The Strengths of Chemical Bonds," 
  *       2nd ed., Butterworths, London, 1958; 
  * - B. deB. Darwent, "National Standard Reference Data Series," National Bureau of Standards, 
  *       No. 31, Washington, DC, 1970; S.W. Benson, J. Chem. Educ., 42, 502 (1965).
  * Can be found on the web at:
  * - http://www.wiredchemist.com/chemistry/data/bond_energies_lengths.html
  */
// NOTE: Store cut^2 instead??
double Atom::GetBondLength(AtomicElementType atom1, AtomicElementType atom2) {
  // Default cutoff
  double cut = 1.60;
  if (atom1==atom2) {
    // Self
    switch (atom1) {
      case HYDROGEN  : cut=0.74; break;
      case CARBON    : cut=1.54; break;
      case NITROGEN  : cut=1.45; break;
      case OXYGEN    : cut=1.48; break;
      case PHOSPHORUS: cut=2.21; break;
      case SULFUR    : cut=2.05; break; // S-S gas-phase value; S=S is 1.49
      default: WarnBondLengthDefault(atom1,atom2,cut);
    }
  } else {
    AtomicElementType e1, e2;
    if (atom1 < atom2) {
      e1 = atom1;
      e2 = atom2;
    } else {
      e1 = atom2;
      e2 = atom1;
    }
    switch (e1) {
      case HYDROGEN: // Bonds to H
        switch (e2) {
          case CARBON    : cut=1.09; break;
          case NITROGEN  : cut=1.01; break;
          case OXYGEN    : cut=0.96; break;
          case PHOSPHORUS: cut=1.44; break;
          case SULFUR    : cut=1.34; break;
          default: WarnBondLengthDefault(e1,e2,cut);
        }
        break;
      case CARBON: // Bonds to C
        switch (e2) {
          case NITROGEN  : cut=1.47; break;
          case OXYGEN    : cut=1.43; break;
          case FLUORINE  : cut=1.35; break;
          case PHOSPHORUS: cut=1.84; break;
          case SULFUR    : cut=1.82; break;
          case CHLORINE  : cut=1.77; break;
          case BROMINE   : cut=1.94; break;
          default: WarnBondLengthDefault(e1,e2,cut);
        }
        break;
      case NITROGEN: // Bonds to N
        switch (e2) {
          case OXYGEN    : cut=1.40; break;
          case FLUORINE  : cut=1.36; break;
          case PHOSPHORUS: cut=1.71; // Avg over all nX-pX from gaff.dat
          case SULFUR    : cut=1.68; break; // Postma & Vos, Acta Cryst. (1973) B29, 915
          case CHLORINE  : cut=1.75; break;
          default: WarnBondLengthDefault(e1,e2,cut);
        }
        break;
      case OXYGEN: // Bonds to O
        switch (e2) {
          case FLUORINE  : cut=1.42; break;
          case PHOSPHORUS: cut=1.63; break;
          case SULFUR    : cut=1.48; break;
          default: WarnBondLengthDefault(e1,e2,cut);
        }
        break;
      case FLUORINE: // Bonds to F
        switch (e2) {
          case PHOSPHORUS: cut=1.54; break;
          case SULFUR    : cut=1.56; break;
          default: WarnBondLengthDefault(e1,e2,cut);
        }
        break;
      case PHOSPHORUS: // Bonds to P
        switch (e2) {
          case SULFUR  : cut=1.86; break;
          case CHLORINE: cut=2.03; break;
          default: WarnBondLengthDefault(e1,e2,cut);
        }
        break;
      case SULFUR: // Bonds to S
        switch (e2) {
          case CHLORINE: cut=2.07; break;
          default: WarnBondLengthDefault(e1,e2,cut);
        }
        break;
      default: WarnBondLengthDefault(e1,e2,cut);
    } // END switch(e1)
  }
  //mprintf("\t\tCUTOFF: [%s] -- [%s] = %lf\n",AtomicElementName[atom1],
  //        AtomicElementName[atom2],cut);
  return cut;
}

