#ifndef INC_AXISTYPE_H
#define INC_AXISTYPE_H
#include "Topology.h"
/*! \file AxisType.h
    \brief Hold classes and functions used for NA structure analysis.
 */
/// Hold information for NA base.
class NA_Base {
  public:
    /// Type for each standard NA base.
    enum NAType { UNKNOWN_BASE, ADE, CYT, GUA, THY, URA };
    /// Type for phosphate/sugar atoms
    enum PSType { PHOS, O4p, C1p, C2p, C3p, C4p };
    NA_Base();
    NA_Base(const NA_Base&);
    NA_Base& operator=(const NA_Base&);
    static NAType ID_BaseFromName(NameType const&);
    NA_Base(Topology const&, int, NAType);
    void SetInputFrame(Frame const&);
    void PrintAtomNames() const;
    NAType Type()                  const { return type_;           }
    int ResNum()                   const { return rnum_;           }
    char BaseChar()                const { return bchar_;          }
    Frame const& Ref()             const { return Ref_;            }
    Frame const& Input()           const { return Inp_;            }
    AtomMask const& InputFitMask() const { return inpFitMask_;     }
    AtomMask const& RefFitMask()   const { return refFitMask_;     }
    const char* AtomName(int i)    const { return *(anames_[i]);   }
    bool HasPatom()                const { return atomIdx_[PHOS] != -1; }
    bool HasO4atom()               const { return atomIdx_[O4p] != -1;  }
    bool HasSugarAtoms()           const;
#   ifdef NASTRUCTDEBUG
    const char* ResName()       const { return *rname_;         }
    const char* RefName(int i)  const { return *(refnames_[i]); }
    int HBidx(int i)            const { return hbidx_[i];       }
#   endif
    const double* HBxyz(int i) const { return Inp_.XYZ(hbidx_[i]);      }
    const double* Pxyz()       const { return Inp_.XYZ(atomIdx_[PHOS]); }
    const double* O4xyz()      const { return Inp_.XYZ(atomIdx_[O4p]);  }
    const double* C1xyz()      const { return Inp_.XYZ(atomIdx_[C1p]);  }
    const double* C2xyz()      const { return Inp_.XYZ(atomIdx_[C2p]);  }
    const double* C3xyz()      const { return Inp_.XYZ(atomIdx_[C3p]);  }
    const double* C4xyz()      const { return Inp_.XYZ(atomIdx_[C4p]);  }
  private:
    int rnum_;                      ///< Original residue number
    char bchar_;                    ///< 1 char base name.
    NAType type_;                   ///< Base type.
    Frame Ref_;                     ///< Reference coords.
    std::vector<NameType> anames_;  ///< Atom names (Input)
#   ifdef NASTRUCTDEBUG
    NameType rname_;                 ///< Residue name
    std::vector<NameType> refnames_; ///< Atom names (Ref)
#   endif  
    Frame Inp_;                     ///< Input coords.
    int hbidx_[3];                  ///< Indices of h-bonding atoms
    int atomIdx_[6];                ///< Indices of phosphate/sugar agoms.
    AtomMask parmMask_;             ///< Mask corresponding to atoms in parm.
    AtomMask inpFitMask_;           ///< Mask of input atoms to be used in RMS fit.
    AtomMask refFitMask_;           ///< Mask of ref atoms to be used in RMS fit.

    int FindAtom(NameType const&) const;
};

/// Hold information for axis corresponding to base/base-pair.
class NA_Axis {
  public:
    NA_Axis();
    /// Used to set up base axis.
    void SetupBaseAxis(Matrix_3x3 const&, Vec3 const&, int);
    /// User to set up base pair axis.
    NA_Axis(int,int,bool);
    /// Used to set rotation matrix/origin for base pair axis
    void StoreRotMatrix(Matrix_3x3 const&, Vec3 const&);
    void PrintAxisInfo(const char*) const;
    void FlipYZ();
    void FlipXY();
    Matrix_3x3 const& Rot() const { return R_;      }
    Vec3 const& Oxyz()      const { return origin_; }
    Vec3 const& Rx()        const { return RX_;     }
    Vec3 const& Ry()        const { return RY_;     }
    Vec3 const& Rz()        const { return RZ_;     }
    int Res1()              const { return residue_number_; }
    int Res2()              const { return second_resnum_;  }
    bool IsAnti()           const { return isAnti_;         }
  private:
    Matrix_3x3 R_;       ///< Rotation matrix for this axis
    Vec3 origin_;        ///< Origin of this axis
    // Rotation X|Y|Z vecs are stored to avoid constant calls to R.
    Vec3 RX_;            ///< Rotation X vector (col 1)
    Vec3 RY_;            ///< Rotation Y vector (col 2)
    Vec3 RZ_;            ///< Rotation Z vector (col 3)
    int residue_number_; ///< Residue number
    int second_resnum_;  ///< Second residue if this is a base pair
    bool isAnti_;        ///< If basepair, true if pair is anti-parallel
};
#endif  
