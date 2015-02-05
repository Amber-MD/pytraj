#include <cstring>
#include "CIFfile.h"
#include "ArgList.h"
#include "CpptrajStdio.h"

static inline int LineError(const char* msg, int num, const char* ptr) {
  mprinterr("Error: CIF line %i: %s\n", num, msg);
  if (ptr != 0) mprinterr("Error: '%s'\n", ptr);
  return 1;
}

/** Separate data header with format '_HEADER.ID' into _HEADER and ID. */
int CIFfile::DataBlock::ParseData( std::string const& sIn, std::string& Header, 
                                   std::string& ID)
{
  size_t found = sIn.find_first_of(".");
  if (found == std::string::npos) {
    mprinterr("Error: No '.' in data record: %s\n", sIn.c_str());
    return 1;
  }
  ID = sIn.substr(found+1);
  Header = sIn.substr(0, found);
  //mprintf("DEBUG:\t\t\tHeader=%s  ID=%s", Header.c_str(), ID.c_str());
  return 0;
}

/** If this is the first entry in a data block, make this the data blocks
  * header. If not, make sure this matches the previously set header.
  */
int CIFfile::DataBlock::AddHeader(std::string const& Header) {
  if (dataHeader_.empty()) { // First entry in this data block
    dataHeader_ = Header;
  } else if (dataHeader_ != Header) {
    mprinterr("Error: Data header in CIF file changes from %s to %s\n",
              dataHeader_.c_str(), Header.c_str());
    return 1;
  }
  return 0;
}

static inline bool IsQuoteChar(char qc) {
  return (qc == '\'' || qc == '"' || qc == ';');
}

/** Add entries to a serial data block. */
int CIFfile::DataBlock::AddSerialDataRecord( const char* ptr, BufferedLine& infile ) {
  if (ptr == 0) return 1;
  // Expect header.id data
  ArgList serialData( ptr, " " );
  std::string dataLine;
  if ( serialData.Nargs() < 2 ) {
    // Could be the data is spread across several lines. Expect that it
    // starts and ends with quote chars.
    const char* nextLine = infile.Line();
    if (nextLine == 0 || !IsQuoteChar(nextLine[0])) {
      mprinterr("Error: Line %i: Data '%s', expected quote character.\n",
                infile.LineNumber(), dataHeader_.c_str());
      return 1;
    }
    bool readMoreLines = true;
    while (readMoreLines) {
      if (nextLine == 0 || strncmp(nextLine, "loop_", 5) == 0) {
        mprinterr("Error: Line %i: Data '%s' record expected to have ID and data.\n",
                  infile.LineNumber(), dataHeader_.c_str());
        return 1;
      }
      dataLine.append( nextLine );
      const char* nextLine = infile.Line();
      if (nextLine!=0 && IsQuoteChar(nextLine[0])) // Terminal quote
        readMoreLines = false;
    }
  } else
    dataLine.assign( serialData[1]);
  std::string ID, Header;
  if (ParseData( serialData[0], Header, ID )) return 1;
  //mprintf("  Ndata=%i  Data=%s\n", serialData.Nargs(), serialData[1].c_str());
  if (AddHeader( Header )) return 1;

  columnHeaders_.push_back( ID );
  if (columnData_.empty()) columnData_.resize( 1 );
  columnData_[0].push_back( dataLine ); 
    
  return 0;
}

/** Add column label from loop section. */
int CIFfile::DataBlock::AddLoopColumn( const char* ptr ) {
  if (ptr == 0) return 1;
  // Expect header.id
  ArgList loopData( ptr, " " );
  if ( loopData.Nargs() > 1 ) {
    mprinterr("Error: Data record expected to have ID only.\n"
              "Error: '%s'\n", ptr);
    return 1;
  }
  std::string ID, Header;
  if (ParseData( loopData[0], Header, ID )) return 1;
  //mprintf("\n"); // DEBUG
  if (AddHeader( Header )) return 1;
  columnHeaders_.push_back( ID );

  return 0;
}  

/** Add loop data. */
int CIFfile::DataBlock::AddLoopData( const char* ptr, BufferedLine& infile ) {
  // Should be as much data as there are column headers
  ArgList loopData( ptr, " " );
  if ( loopData.Nargs() != (int)columnHeaders_.size()) {
    // Could be there are more lines of data. As long as we dont hit
    // another loop or data, add until we reach the correct number of
    // columns.
    int columns_read = loopData.Nargs();
    while (columns_read < (int)columnHeaders_.size()) {
      const char* nextLine = infile.Line();
      if (nextLine == 0 || nextLine[0] == '_' || strncmp(nextLine, "loop_", 5) == 0) {
        mprinterr("Error: Line %i: # of columns in loop data '%s' (%i) < # column headers (%zu)\n",
                  infile.LineNumber(), dataHeader_.c_str(),
                  loopData.Nargs(), columnHeaders_.size());
        return 1;
      }
      ArgList nextData( nextLine, " " );
      columns_read += nextData.Nargs();
      if (columns_read > (int)columnHeaders_.size()) {
        mprinterr("Error: Line %i: # of columns in loop data '%s' (%i) > # column headers (%zu)\n",
                  infile.LineNumber(), dataHeader_.c_str(),
                  loopData.Nargs(), columnHeaders_.size());
        return 1;
      }
      for (ArgList::const_iterator ia = nextData.begin(); ia != nextData.end(); ++ia)
        loopData.AddArg( *ia );
    }
  }
  columnData_.push_back( loopData.List() );
  return 0;
} 

/** List data currently stored in the block. */
void CIFfile::DataBlock::ListData() const {
  mprintf("DataBlock: %s\n", dataHeader_.c_str());
  for (Sarray::const_iterator colname = columnHeaders_.begin();
                              colname != columnHeaders_.end(); ++colname)
    mprintf("  Col %u name: %s\n", colname - columnHeaders_.begin(), (*colname).c_str());
  for (std::vector<Sarray>::const_iterator rec = columnData_.begin();
                                           rec != columnData_.end(); ++rec)
  {
    mprintf("    [%u] ", rec - columnData_.begin());
    for (Sarray::const_iterator col = (*rec).begin();
                                col != (*rec).end(); ++col)
      mprintf(" %s", (*col).c_str());
    mprintf("\n");
  }
}

/** \return the index of the specified column, -1 if not present. */
int CIFfile::DataBlock::ColumnIndex(std::string const& headerIn) const {
  for (Sarray::const_iterator col = columnHeaders_.begin();
                              col != columnHeaders_.end(); ++col)
    if (headerIn == *col)
      return (int)(col - columnHeaders_.begin());
  return -1;
}

std::string CIFfile::DataBlock::Data(std::string const& idIn) const {
  if (columnHeaders_.empty() || columnData_.empty()) return std::string("");
  int colnum = ColumnIndex( idIn );
  if (colnum == -1) return std::string("");
  return columnData_[0][colnum];
}

// -----------------------------------------------------------------------------
/// Used to return empty block for GetDataBlock
const CIFfile::DataBlock CIFfile::emptyblock = DataBlock();

/** Determine if fileIn is a CIF file. Look for entries beginning with 
  * an underscore (indicating data block), and a 'loop_' keyword or
  * '_entry.id' block.
  */
bool CIFfile::ID_CIF(CpptrajFile& fileIn) {
  // NOTE: ASSUME FILE SET UP FOR READ
  if (fileIn.OpenFile()) return false;
  int ndata = 0; // Number of '_XXX' entries seen
  bool foundLoop = false;
  bool foundEntryID = false;
  for (int i = 0; i < 10; i++) {
    std::string lineIn = fileIn.GetLine();
    if (lineIn[0] == '_') ndata++;
    if (lineIn.compare(0,5,"loop_")==0) foundLoop = true;
    if (lineIn.compare(0,9,"_entry.id")==0) foundEntryID = true;
  }
  fileIn.CloseFile();
  return ( ndata > 2 && (foundLoop || foundEntryID) );
}

// CIFfile::Read()
int CIFfile::Read(std::string const& fnameIn, int debugIn) {
  if (file_.OpenFileRead( fnameIn )) return 1;
  const char* ptr = file_.Line();
  mode currentMode = UNKNOWN;
  while (ptr != 0) {
    /// There are 3 places we can be; a data block, a looped data block,
    /// or unknown.
    if ( currentMode == UNKNOWN ) {
      // Are we at a data block yet?
      if (ptr[0] == '_')
        currentMode = SERIAL;
      else if ( strncmp(ptr, "loop_", 5) == 0 )
        currentMode = LOOP;
      else
        ptr = file_.Line();

    } else if ( currentMode == SERIAL ) {
      // SERIAL data block
      DataBlock serial;
      while ( ptr != 0 && ptr[0] == '_' ) {
        serial.AddSerialDataRecord(ptr, file_);
        ptr = file_.Line();
      }
      if (debugIn > 1) serial.ListData();
      currentMode = UNKNOWN;
      if (AddDataBlock( serial )) return 1;
    } else if ( currentMode == LOOP ) {
      DataBlock loop;
      ptr = file_.Line();
      if (ptr == 0 || ptr[0] != '_')
        return LineError("In CIF file, malformed loop.", file_.LineNumber(), ptr);
      while (ptr != 0 && ptr[0] == '_') {
        loop.AddLoopColumn(ptr);
        ptr = file_.Line();
      }
      // Should now be positioned at loop data
      if (ptr == 0)
        return LineError("In CIF file, no loop data.", file_.LineNumber(), ptr);
      while (ptr != 0 && ptr[0] != '_' && ptr[0] != '#') {
        loop.AddLoopData(ptr, file_);
        ptr = file_.Line();
      }
      if (debugIn > 1) loop.ListData();
      currentMode = UNKNOWN;
      if (AddDataBlock( loop )) return 1;
    }
  }
  if (debugIn > 0)    
    mprintf("\tCIF file '%s', %i lines.\n", file_.Filename().full(), file_.LineNumber());
  return 0;
}

// CIFfile::GetDataBlock()
CIFfile::DataBlock const& CIFfile::GetDataBlock(std::string const& header) const {
  CIF_DataType::const_iterator it = cifdata_.find( header );
  if (it == cifdata_.end()) {
    //mprinterr("Error: CIF data block '%s' not found.\n", header.c_str());
    return emptyblock;
  }
  return (*it).second;
}

// CIFfile::AddDataBlock()
int CIFfile::AddDataBlock( DataBlock const& block ) {
  if (block.Header().empty()) {
    mprinterr("Internal Error: Attempting to add empty CIF data block.\n");
    return 1;
  }
  CIF_DataType::const_iterator it = cifdata_.find( block.Header() );
  if (it != cifdata_.end()) {
    mprinterr("Error: Duplicate CIF block found: '%s'\n", block.Header().c_str());
    return 1;
  }
  cifdata_.insert( std::pair<std::string, DataBlock>(block.Header(), block) );
  return 0;
}
