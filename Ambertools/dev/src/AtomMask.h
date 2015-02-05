#ifndef INC_ATOMMASK_H
#define INC_ATOMMASK_H
#include <string>
#include <vector>
#include "MaskToken.h"
// Class: AtomMask
/// Hold info on selected atoms based on mask expression.
/** AtomMask is used to hold an array of integers that represent atom numbers
  * of atoms selected based on a mask string. 
  * First the mask string is set via SetMaskString, where it is converted into
  * mask tokens. Then the actual mask can be set up using a Topology in one of
  * two ways:
  * - Integer mask
  *    Although an array of ints becomes larger than a simple character mask 
  *    once more than 25% of the system is selected, it tends to be faster 
  *    than the character array up until about 80% of the system is selected, 
  *    at which point the speed is comparable. This is the default way to use
  *    AtomMask and is how most of the routines in the Frame class have been
  *    written to use AtomMask.
  * - Character mask
  *    This is the original way to use the atom mask, useful e.g. when 
  *    you need to know what atoms are not selected as well as what atoms
  *    are selected. Unlike the integer mask, the character mask is not
  *    directly accessible by outside routines, only by the AtomInCharMask
  *    routine.
  */
// *** NOTE ***
// invertMask, AddAtom, AddAtoms, and AddAtomRange currently only apply
// to the Selected array.
class AtomMask {
  public:
    AtomMask();
    AtomMask(std::string const&);
    ///< Create mask selecting atoms from begin to end.
    AtomMask(int,int);
    AtomMask(int);
    AtomMask(const AtomMask &);
    AtomMask& operator=(const AtomMask&);
    /// \return Internal selected atom array.
    std::vector<int> const& Selected()  const { return Selected_;            }
    /// AtomMask default iterator
    typedef std::vector<int>::const_iterator const_iterator;
    /// \return const iterator to the beginning of Selected
    const_iterator begin()              const { return Selected_.begin();    }
    /// \return const iterator at end of Selected
    const_iterator end()                const { return Selected_.end();      }
    /// \return last selected atom
    int back()                          const { return Selected_.back();     }
    /// \return number of selected atoms
    int Nselected()                     const { return nselected_;           }
    /// \return selected atom at idx
    const int& operator[](int idx)      const { return Selected_[idx];       }
    /// \return original mask expression as char*
    const char *MaskString()            const { return maskString_.c_str();  }
    /// \return original mask expression as std::string
    std::string const& MaskExpression() const { return maskString_;          }
    /// \return true if mask expression has been set.
    bool MaskStringSet()                const { return !maskString_.empty(); }
    /// \return true if no atoms selected. 
    bool None()                         const { return (nselected_==0);      }
    /// \return true is this is a character mask.
    bool IsCharMask()                   const { return (!CharMask_.empty()); }
    /// Reset atom mask
    void ResetMask();
    /// Clear any selected atoms in mask.
    void ClearSelected();
    /// Switch char used to denote selected atoms (T->F, F->T)
    void InvertMask();
    /// \return the number of atoms mask has in common with another mask
    int NumAtomsInCommon(AtomMask const&);
    /// Add atom to Selected array; assumes atoms will be in order.
    void AddSelectedAtom(int i) { Selected_.push_back( i ); nselected_=(int)Selected_.size(); }
    /// Add given atom to Selected array 
    void AddAtom(int);
    /// Add a list of atoms to mask
    void AddAtoms(std::vector<int> const&);
    /// Add minAtom <= atom < maxAtom to mask
    void AddAtomRange(int,int);
    /// Add atoms in given mask to this mask at positon, update position
    void AddMaskAtPosition(AtomMask const&, int);
    /// Print all mask atoms in to a line
    void PrintMaskAtoms(const char*) const;
    /// Set the mask string. If NULL, set * (all)
    int SetMaskString(const char*);
    /// Set the mask string.
    int SetMaskString( std::string const& );
    /// Set up Selected based on given char mask 
    void SetupIntMask(const char*,int,int);
    /// Set up CharMask based on given char mask 
    void SetupCharMask(const char*, int, int);
    /// \return true if given atom is T in CharMask
    bool AtomInCharMask(int) const;
    /// \return true if any atoms within given range are T in CharMask.
    bool AtomsInCharMask(int,int) const;
    /// Set number of atoms, needed for integer to char mask conversion.
    void SetNatom( int a) { Natom_ = a; }
    /// Convert from integer mask to char mask.
    int ConvertToCharMask();
    /// Convert from char mask to integer mask.
    int ConvertToIntMask();
    /// Print mask string and number of selected atoms.
    void MaskInfo() const;
    /// Print brief mask info
    void BriefMaskInfo() const;
    // Used by Topology to set up masks.
    typedef std::vector<MaskToken>::const_iterator token_iterator;
    inline token_iterator begintoken() const { return maskTokens_.begin(); }
    inline token_iterator endtoken()   const { return maskTokens_.end();   }
  private:
    int debug_;
    std::vector<char> CharMask_; ///< Char array of atoms, T if selected, F if not.
    char maskChar_;              ///< The character used to denote a selected atom (default 'T')
    std::string maskString_;     ///< String specifying atom selection
    /** Number of atoms mask was set-up with. Needed when converting from
      * integer mask to Character mask. */
    int Natom_;
    int nselected_;              ///< Number of selected atoms in mask
    std::vector<int> Selected_;  ///< Int array of selected atom numbers, 1 for each selected atom
    std::vector<MaskToken> maskTokens_;
    /// Convert mask expression to MaskTokens
    int Tokenize();
};
#endif
