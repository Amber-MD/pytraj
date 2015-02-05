#ifndef INC_TINKERFILE_H
#define INC_TINKERFILE_H
#include "BufferedLine.h"
#include "Box.h"
#include "Atom.h"
/// Use to access Tinker XYZ/ARC files.
class TinkerFile {
  public:
    TinkerFile();
    static bool ID_Tinker(CpptrajFile&);
    void SetTinkerName(std::string const& t) { tinkerName_ = t; }
    int TinkerNatom() const { return natom_; }
    Box const& TinkerBox() const { return box_; }
    std::string const& TinkerTitle() const { return title_; }

    int OpenTinker();
    int NextTinkerFrame();
    int ReadNextTinkerFrame(double*, double*);
    void Rewind() {
      file_.CloseFile();
      file_.OpenFileRead( tinkerName_ );
    }
    void CloseFile() { file_.CloseFile(); }
    FileName const& Filename() { return file_.Filename(); }
    std::vector<Atom> ReadTinkerAtoms(double*, std::vector<int>&);
  private:
    int CheckTitleLine();

    BufferedLine file_;
    int natom_; ///< Number of atoms in file.
    bool hasBox_; /// true if file has box coords.
    Box box_; ///< hold box info for first coords.
    std::string title_; ///< Title.
    std::string tinkerName_; ///< Full file name
};
#endif
