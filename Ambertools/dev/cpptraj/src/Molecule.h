#ifndef INC_MOLECULE_H
#define INC_MOLECULE_H
// Class: Molecule
/// Hold information for a molecule
class Molecule {
  public:
    Molecule() : beginAtom_(0), endAtom_(0), isSolvent_(false) {}
    Molecule(int begin, int end) :
      beginAtom_(begin),
      endAtom_(end),
      isSolvent_(false)
    {}

    void SetFirst(int begin) { beginAtom_ = begin; }
    void SetLast(int last)   { endAtom_ = last;    }
    void SetSolvent()        { isSolvent_ = true;  }
    void SetNoSolvent()      { isSolvent_ = false; }

    inline int BeginAtom() const  { return beginAtom_;   }
    inline int EndAtom() const    { return endAtom_;     } 
    inline bool IsSolvent() const { return isSolvent_;   }
    inline int NumAtoms() const   { return (endAtom_ - beginAtom_); }

  private:
    int beginAtom_;
    int endAtom_;
    bool isSolvent_;
};
#endif
