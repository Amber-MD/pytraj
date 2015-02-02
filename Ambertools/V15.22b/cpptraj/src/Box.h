#ifndef INC_BOX_H
#define INC_BOX_H
#include "Matrix_3x3.h"
/// Hold box information; 3xlengths, 3xangles.
class Box {
  public:
    enum BoxType { NOBOX=0, ORTHO, TRUNCOCT, RHOMBIC, NONORTHO }; 

    Box();
    Box(const double*);
    Box(Matrix_3x3 const&);
    Box(const Box&);
    Box& operator=(const Box&);

    const char* TypeName() const; 

    void SetBetaLengths(double,double,double,double);
    void SetBox(const double*);
    void SetBox(Matrix_3x3 const&);
    void SetTruncOct();
    void SetNoBox();
    void SetMissingInfo(const Box&);

    double ToRecip(Matrix_3x3&, Matrix_3x3&) const;

    void SetX(double xin)     { box_[0] = xin; }
    void SetY(double yin)     { box_[1] = yin; }
    void SetZ(double zin)     { box_[2] = zin; }
    void SetAlpha(double ain) { box_[3] = ain; }
    void SetBeta(double bin)  { box_[4] = bin; }
    void SetGamma(double gin) { box_[5] = gin; }

    BoxType Type() const { return btype_;  }
    double BoxX()  const { return box_[0]; }
    double BoxY()  const { return box_[1]; }
    double BoxZ()  const { return box_[2]; }
    double Alpha() const { return box_[3]; }
    double Beta()  const { return box_[4]; }
    double Gamma() const { return box_[5]; }
    bool HasBox()  const { return (btype_ != NOBOX); }
    Vec3 Center()  const { return Vec3(box_[0]/2.0, box_[1]/2.0, box_[2]/2.0); }
    Vec3 Lengths() const { return Vec3(box_[0], box_[1], box_[2]);             }
    // For interfacing with file IO
    double* boxPtr()             { return box_; }
    const double* boxPtr() const { return box_; }

    double const& operator[](int idx) const { return box_[idx]; }
    double&       operator[](int idx)       { return box_[idx]; }
  private:
    static const double TRUNCOCTBETA;
    static const char* BoxNames[];
    //int debug_; // TODO: Replace with ifdefs or just comment out?
    BoxType btype_;
    double box_[6];

    void SetBoxType();
};
#endif
