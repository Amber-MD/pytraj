#ifndef INC_MASKTOKEN_H
#define INC_MASKTOKEN_H
#include "NameType.h"
/// Hold information use in mask selection. 
class MaskToken {
  public:
    enum MaskTokenType { 
      OP_NONE=0, ResNum, ResName, AtomNum, AtomName, AtomType, AtomElement, SelectAll,
      OP_AND, OP_OR, OP_NEG, OP_DIST
    };
    MaskToken();
    const char *TypeName() const;
    void Print() const;
    int SetToken( MaskTokenType, std::string const& );
    int SetDistance( std::string & );
    void SetOperator(MaskTokenType);

    inline MaskTokenType Type()   const { return type_;     }
    inline int Res1()             const { return res1_;     }
    inline int Res2()             const { return res2_;     }
    inline const NameType& Name() const { return name_;     }
    inline bool OnStack()         const { return onStack_;  }
    inline bool Within()          const { return d_within_; }
    inline bool ByAtom()          const { return d_atom_;   }
    inline double Distance()      const { return distance_; }

    void SetOnStack();
  private:
    static const char* MaskTypeString[];

    MaskTokenType type_;
    int res1_;
    int res2_;
    NameType name_;
    bool onStack_;
    // Distance criteria
    bool d_within_;
    bool d_atom_;
    double distance_;

    void MakeNameType();
};
#endif
