#ifndef INC_FILENAME_H
#define INC_FILENAME_H
#include <string>
/// Class to hold file name, extension, etc.
class FileName {
  public:
    FileName() {}
    FileName(const FileName&);
    FileName& operator=(const FileName&);
    /// Set file name and extensions; only name is known.
    int SetFileName(std::string const&);
    /// Set file name and extension; perform tilde expansion.
    int SetFileNameWithExpansion(std::string const&);
    /// Set file name and extensions; name and compressed status are known.
    int SetFileName(std::string const&, bool);
    /// Clear FileName
    void clear();
    /// \return true if string matches full or base file name.
    bool MatchFullOrBase(std::string const&) const;

    const std::string& Full()      const { return fullPathName_;         }
    const std::string& Base()      const { return baseName_;             }
    const char* full()             const { return fullPathName_.c_str(); }
    const char* base()             const { return baseName_.c_str();     }
    const std::string& Ext()       const { return extension_;            }
    const std::string& Compress()  const { return compressExt_;          }
    const std::string& DirPrefix() const { return dirPrefix_;            }
    bool empty()                   const { return fullPathName_.empty(); }
  private:
    enum CompressStatus { UNKNOWN=0, YES, NO };
    std::string fullPathName_;
    std::string baseName_;
    std::string extension_;
    std::string compressExt_;
    std::string dirPrefix_;

    int SetFileName(std::string const&, CompressStatus);
};
#endif
