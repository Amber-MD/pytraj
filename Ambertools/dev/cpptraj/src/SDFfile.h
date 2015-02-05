#ifndef INC_SDFFILE_H
#define INC_SDFFILE_H
#include "CpptrajFile.h"
#include "Atom.h"
/// Used to access SDF files.
class SDFfile : public CpptrajFile {
  public: 
    SDFfile();
    /// \return true if file is an SDF file
    static bool ID_SDF(CpptrajFile&);
    /// Read in and past header section of SDF file.
    bool ReadHeader();
    /// Read in the next SDF ATOM line. Get the X Y and Z coords.
    int SDF_XYZ(double *);
    /// Convert current line to Atom 
    Atom SDF_Atom();
    /// Read in the next SDF BOND line. Get the indices of the bonded atoms.
    int SDF_Bond(int &, int &);

    int SDF_Natoms()                    const { return Natoms_; }
    int SDF_Nbonds()                    const { return Nbonds_; }
    std::string const& SDF_Title()      const { return title_;  }
  private:
    int debug_;
    int Natoms_;
    int Nbonds_;
    char Name_[5]; 
    std::string title_;
};
#endif  
