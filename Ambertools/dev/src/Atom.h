#ifndef INC_ATOM_H
#define INC_ATOM_H
#include <vector>
#include <set> // For excluded 
#include "NameType.h"
// Class: AtomType
/// Hold information for an atom
class Atom {
  public:
    enum AtomicElementType { UNKNOWN_ELEMENT = 0,
      HYDROGEN,   BORON,      CARBON,   NITROGEN,  OXYGEN,     FLUORINE,
      PHOSPHORUS, SULFUR,     CHLORINE, BROMINE,   IRON,       CALCIUM,
      IODINE,     MAGNESIUM,  COPPER,   LITHIUM,   POTASSIUM,  RUBIDIUM,
      CESIUM,     ZINC,       SODIUM,   ALUMINUM,  ARGON,      ARSENIC,
      SILVER,     GOLD,       ASTATINE, BERYLLIUM, BARIUM,     BISMUTH,
      CHROMIUM,   COBALT,     CADMIUM,  FRANCIUM,  GALLIUM,    GERMANIUM,
      HELIUM,     HAFNIUM,    MERCURY,  INDIUM,    IRIDIUM,    KRYPTON,
      MANGANESE,  MOLYBDENUM, NEON,     NICKEL,    NIOBIUM,    OSMIUM,
      PALLADIUM,  PLATINUM,   LEAD,     POLONIUM,  RUTHENIUM,  RHODIUM,
      RHENIUM,    RADON,      RADIUM,   SILICON,   SCANDIUM,   SELENIUM,
      STRONTIUM,  TIN,        ANTIMONY, TITANIUM,  TECHNETIUM, TELLURIUM,
      TANTALUM,   THALLIUM,   VANADIUM, TUNGSTEN,  XENON,      ZIRCONIUM,
      YTTRIUM,    LUTETIUM,
      EXTRAPT 
    };
    // Constructors and assignment
    Atom();
    virtual ~Atom() {}
    /// Take atom name, chain ID, and (optional) 2 character element name.
    Atom(NameType const&, char, const char*);
    /// Take atom name, type name, and charge.
    Atom(NameType const&, NameType const&, double);
    /// Take atom name, type name, and type index.
    Atom(NameType const&, NameType const&, int);
    /// Take atom name, charge, mass, and type name
    Atom(NameType const&, double, double, NameType const&);
    /// Atom name, charge, polarizability atomic num, mass, type index, type name, gb radius and screen parameters, and chain ID.
    Atom(NameType const&, double, double, int, double, int, NameType const&, double, double, char);
    Atom(const Atom &);
    void swap(Atom &, Atom &);
    Atom &operator=(Atom);
    // Iterator over bonded atom #s
    typedef std::vector<int>::const_iterator bond_iterator;
    inline bond_iterator bondbegin() const { return bonds_.begin(); }
    inline bond_iterator bondend()   const { return bonds_.end();   }
    // Iterator over excluded atoms
    typedef std::vector<int>::const_iterator excluded_iterator;
    inline excluded_iterator excludedbegin() const { return excluded_.begin(); }
    inline excluded_iterator excludedend()   const { return excluded_.end();   }
    // Functions that set internal vars
    void SetResNum(int resnumIn)             { resnum_ = resnumIn;  }
    void SetMol(int molIn)                   { mol_ = molIn;        }
    void SetCharge(double qin)               { charge_ = qin;       }
    void SetGBradius(double rin)             { gb_radius_ = rin;    }
    void SetTypeIndex(int tin)               { atype_index_ = tin;  }
    // Inline functions returning internal vars
    inline bool NoMol()                const { return ( mol_ < 0 ); }
    inline const char *c_str()         const { return *aname_; }
    inline int ResNum()                const { return resnum_; }
    inline AtomicElementType Element() const { return element_; }
    inline int AtomicNumber()          const { return AtomicElementNum[element_];  }
    inline const char* ElementName()   const { return AtomicElementName[element_]; }
    inline const NameType& Name()      const { return aname_; }
    inline const NameType& Type()      const { return atype_; }
    inline int TypeIndex()             const { return atype_index_; }
    inline int MolNum()                const { return mol_; }
    inline char ChainID()              const { return chainID_; }
    inline int Nbonds()                const { return (int)bonds_.size(); }
    inline int Nexcluded()             const { return (int)excluded_.size(); }
    inline double Mass()               const { return mass_; }
    inline double Charge()             const { return charge_; }
    inline double Polar()              const { return polar_; }
    inline double GBRadius()           const { return gb_radius_; }
    inline double Screen()             const { return gb_screen_; }
    /// Add atom # to this atoms list of bonded atoms.
    void AddBond(int);
    void ClearBonds();
    void SortBonds();
    // TODO: Use this routine in AtomMap etc
    /// \return true if this atom is bonded to given atom number
    bool IsBondedTo(int) const;
    /// Create exclusion list from input set.
    void AddExclusionList(std::set<int> const&);
    /// \return Optimal bond length based on element types
    static double GetBondLength(AtomicElementType, AtomicElementType);
  protected:
    static const size_t NUMELEMENTS = 76;
  private:
    static const int AtomicElementNum[];
    static const char* AtomicElementName[];
    static const double AtomicElementMass[];
    double charge_;    ///< Charge in e-
    double polar_;     ///< Atomic polarizability in Ang^3
    double mass_;      ///< mass in amu
    double gb_radius_; ///< GB radius in Ang
    double gb_screen_; ///< GB screening parameter
    NameType aname_;   ///< Atom name
    NameType atype_;   ///< Atom type name
    int atype_index_;  ///< Atom type index for nonbond params
    AtomicElementType element_;
    int resnum_;
    int mol_;
    char chainID_;
    std::vector<int> bonds_;
    std::vector<int> excluded_;

    static void WarnBondLengthDefault(AtomicElementType, AtomicElementType, double);
    void DetermineElement(int);
    void SetElementFromName();
    void SetElementFromSymbol(char,char);
    void SetElementFromMass();
};
#endif
