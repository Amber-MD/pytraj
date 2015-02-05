#ifndef INC_CIFFILE_H
#define INC_CIFFILE_H
#include "Atom.h"
#include "BufferedLine.h"
#include <map>
/// Used to access CIF files
class CIFfile {
  public:
    class DataBlock;

    CIFfile() {}
    static bool ID_CIF( CpptrajFile& );
    int Read(std::string const&,int);
    /// \return const reference to the specified data block.
    DataBlock const& GetDataBlock(std::string const&) const;
    FileName const& CIFname() const { return file_.Filename(); }
  private:
    int AddDataBlock(DataBlock const&);

    enum mode { UNKNOWN = 0, SERIAL, LOOP };
    typedef std::vector<std::string> Sarray;
    BufferedLine file_;
    typedef std::map<std::string, DataBlock> CIF_DataType;
    CIF_DataType cifdata_;
    static const DataBlock emptyblock;
};
/// Used to hold CIF data blocks
class CIFfile::DataBlock {
  public:
    DataBlock() {}
    std::string const& Header() const { return dataHeader_;         }
    bool empty()                const { return dataHeader_.empty(); }
    int AddHeader(std::string const&);
    int AddSerialDataRecord(const char*, BufferedLine&);
    int AddLoopColumn(const char*);
    int AddLoopData(const char*, BufferedLine&);
    void ListData() const;
    int ColumnIndex(std::string const&) const;
    /// \return Serial data for given ID
    std::string Data(std::string const&) const;
    // Iterators
    typedef std::vector<Sarray>::const_iterator data_it;
    data_it begin() const { return columnData_.begin(); }
    data_it end()   const { return columnData_.end();   }
  private:
    static int ParseData(std::string const&, std::string&, std::string&);

    std::string dataHeader_; ///< The data header, e.g. '_atom_site'
    Sarray columnHeaders_;   ///< Column headers, e.g. 'label_atom_id'
    std::vector<Sarray> columnData_; ///< Array of column data, e.g.:
      /*
ATOM 1    N N    . SER A 1 1  ? -2.559 9.064   0.084   1.00 0.00 ? ? ? ? ? ? 1  SER A N    1
ATOM 2    C CA   . SER A 1 1  ? -3.245 8.118   0.982   1.00 0.00 ? ? ? ? ? ? 1  SER A CA   1
       */
};
#endif
