#include <cstdlib> // atoi, atof
#include <cstring> // strchr
#include <cctype>  // isdigit, isalpha
#include "DataIO_Std.h"
#include "CpptrajStdio.h" 
#include "StringRoutines.h" // SetStringFormatString
#include "Array1D.h"
#include "BufferedLine.h"
#include "DataSet_MatrixDbl.h"
#include "DataSet_3D.h"

// CONSTRUCTOR
DataIO_Std::DataIO_Std() :
  DataIO(true, true, true), // Valid for 1D, 2D, 3D
  mode_(READ1D),
  isInverted_(false), 
  hasXcolumn_(true), 
  writeHeader_(true), 
  square2d_(false)
{}

static void PrintColumnError(int idx) {
  mprinterr("Error: Number of columns in file changes at line %i.\n", idx);
}

void DataIO_Std::ReadHelp() {
  mprintf("\tread1d:      Read data as 1D data sets (default).\n"
          "\tread2d:      Read data as 2D square matrix.\n"
          "\tvector:      Read data as vector: VX VY VZ [OX OY OZ]\n"
          "\tindex <col>: (1D) Use column # (starting from 1) as index column.\n");

}

const char* DataIO_Std::SEPARATORS = " ,\t"; // whitespace, comma, or tab-delimited

int DataIO_Std::processReadArgs(ArgList& argIn) {
  mode_ = READ1D;
  if (argIn.hasKey("read1d")) mode_ = READ1D;
  else if (argIn.hasKey("read2d")) mode_ = READ2D;
  else if (argIn.hasKey("vector")) mode_ = READVEC;
  indexcol_ = argIn.getKeyInt("index", -1);
  return 0;
}
  

// TODO: Set dimension labels
// DataIO_Std::ReadData()
int DataIO_Std::ReadData(std::string const& fname, 
                         DataSetList& dsl, std::string const& dsname)
{
  int err = 0;
  switch ( mode_ ) {
    case READ1D: err = Read_1D(fname, dsl, dsname); break;
    case READ2D: err = Read_2D(fname, dsl, dsname); break;
    case READVEC: err = Read_Vector(fname, dsl, dsname); break;
  }
  return err;
}

// DataIO_Std::Read_1D()
int DataIO_Std::Read_1D(std::string const& fname, 
                        DataSetList& datasetlist, std::string const& dsname)
{
  ArgList labels;
  bool hasLabels = false;
  // Column user args start from 1
  if (indexcol_ != -1)
    mprintf("\tUsing column %i as index column.\n", indexcol_--);

  // Buffer file
  BufferedLine buffer;
  if (buffer.OpenFileRead( fname )) return 1;

  // Read the first line. Attempt to determine the number of columns
  const char* linebuffer = buffer.Line();
  if (linebuffer == 0) return 1;
  int ntoken = buffer.TokenizeLine( SEPARATORS );
  if ( ntoken == 0 ) {
    mprinterr("Error: No columns detected in %s\n", buffer.Filename().full());
    return 1;
  }

  // Try to skip past any comments. If line begins with a '#', assume it
  // contains labels. 
  bool isCommentLine = true;
  const char* ptr = linebuffer;
  while (isCommentLine) {
    // Skip past any whitespace
    while ( *ptr != '\0' && isspace(*ptr) ) ++ptr;
    // Assume these are column labels until proven otherwise.
    if (*ptr == '#') {
      labels.SetList(ptr+1, SEPARATORS );
      if (!labels.empty()) {
        hasLabels = true;
        // If first label is Frame assume it is the index column
        if (labels[0] == "Frame" && indexcol_ == -1)
          indexcol_ = 0;
      }
      linebuffer = buffer.Line();
      ptr = linebuffer;
      if (ptr == 0) {
        mprinterr("Error: No data found in file.\n");
        return 1;
      }
    } else 
      // Not a recognized comment character, assume data.
      isCommentLine = false;
  }
  // Should be at first data line. Tokenize the line.
  ntoken = buffer.TokenizeLine( SEPARATORS );
  // If # of data columns does not match # labels, clear labels.
  if ( !labels.empty() && ntoken != labels.Nargs() ) {
    labels.ClearList();
    hasLabels = false;
  }
  if (indexcol_ != -1 && indexcol_ >= ntoken) {
    mprinterr("Error: Specified index column %i is out of range (%i columns).\n",
              indexcol_+1, ntoken);
    return 1;
  }
  // Determine the type of data stored in each column. Assume numbers should
  // be read with double precision.
  typedef std::vector<double> Darray;
  typedef std::vector<std::string> Sarray;
  std::vector<Darray> Dsets; // double data sets
  std::vector<Sarray> Ssets; // string data sets
  std::vector<int> SetIndices(ntoken, 0); // Indices into set arrays for each column
  for (int col = 0; col < ntoken; ++col) {
    const char* token = buffer.NextToken();
    // Determine data type
    if ( isdigit( token[0] )    || 
                  token[0]=='+' || 
                  token[0]=='-' ||
                  token[0]=='.'   )
    { // Number (double). Indices will start from 1.
      Dsets.push_back(Darray());
      SetIndices[col] = (int)Dsets.size();
    } else {
      // Assume string. STRING columns cannot be index columns
      if ( col == indexcol_ )
        mprintf("Warning: DataFile %s index column %i has string values and will be skipped.\n", 
                  buffer.Filename().full(), indexcol_+1);
      else {
        // Indices will start from -1.
        Ssets.push_back(Sarray());
        SetIndices[col] = -((int)Ssets.size());
      }
    }
  }
  //mprintf("DBG: SetIndices={");
  //for (std::vector<int>::const_iterator it = SetIndices.begin(); it != SetIndices.end(); ++it)
  //  mprintf(" %i", *it);
  //mprintf(" }\n");
  //mprintf("%zu double sets, %zu string sets.\n", Dsets.size(), Ssets.size());
  if (Dsets.empty() && Ssets.empty()) {
    mprinterr("Error: No data detected.\n");
    return 1;
  } 
  // Read in data.
  unsigned int Ndata = 0;
  do {
    if ( buffer.TokenizeLine( SEPARATORS ) != ntoken ) {
      PrintColumnError(buffer.LineNumber());
      break;
    }
    // Convert data in columns
    for (int i = 0; i < ntoken; ++i) {
      const char* token = buffer.NextToken();
      if (SetIndices[i] > 0)      // Double value
        Dsets[SetIndices[i]-1].push_back( atof(token) );
      else if (SetIndices[i] < 0) // String value
        Ssets[-SetIndices[i]-1].push_back( std::string(token) );
    }
    Ndata++;
  } while (buffer.Line() != 0); // Read in next line.
  buffer.CloseFile();
  mprintf("\tDataFile %s has %i columns, %i lines.\n", buffer.Filename().full(),
          ntoken, buffer.LineNumber());
  if (hasLabels) {
    mprintf("\tDataFile contains labels:\n");
    labels.PrintList();
  }
  // If no x column read in just use indices.
  Darray* Xptr = 0;
  Darray Xvals;
  if (indexcol_ == -1) {
    Xvals.reserve(Ndata);
    for (unsigned int i = 0; i < Ndata; i++)
      Xvals.push_back( i );
    Xptr = &Xvals;
  } else {
    mprintf("\tIndex column is %i\n", indexcol_ + 1);
    Xptr = &(Dsets[SetIndices[indexcol_]-1]);
    SetIndices[indexcol_] = 0;
  }
  for (int i = 0; i != ntoken; i++) {
    DataSet* ds = 0;
    if (SetIndices[i] != 0) {
      if (SetIndices[i] > 0)
        ds = datasetlist.AddOrAppendSet(dsname, i+1, "", *Xptr, Dsets[SetIndices[i]-1]);
      else //if (SetIndices[i] < 0)
        ds = datasetlist.AddOrAppendSet(dsname, i+1, "", Ssets[-SetIndices[i]-1]);
      if (ds == 0) return 1;
      if (hasLabels) ds->SetLegend( labels[i] );
    }
  }

  return 0;
}

// DataIO_Std::Read_2D()
int DataIO_Std::Read_2D(std::string const& fname, 
                        DataSetList& datasetlist, std::string const& dsname)
{
  // Buffer file
  BufferedLine buffer;
  if (buffer.OpenFileRead( fname )) return 1;
  mprintf("\tData will be read as a 2D square matrix.\n");
  // Skip comments
  const char* linebuffer = buffer.Line();
  while (linebuffer != 0 && linebuffer[0] == '#')
    linebuffer = buffer.Line();
  int ncols = -1;
  int nrows = 0;
  std::vector<double> matrixArray;
  while (linebuffer != 0) {
    int ntokens = buffer.TokenizeLine( SEPARATORS );
    if (ncols < 0) {
      ncols = ntokens;
      if (ntokens < 1) {
        mprinterr("Error: Could not tokenize line.\n");
        return 1;
      }
    } else if (ncols != ntokens) {
      mprinterr("Error: In 2D file, number of columns changes from %i to %i at line %i\n",
                ncols, ntokens, buffer.LineNumber());
      return 1;
    }
    for (int i = 0; i < ntokens; i++)
      matrixArray.push_back( atof( buffer.NextToken() ) );
    nrows++;
    linebuffer = buffer.Line();
  }
  if (ncols < 0) {
    mprinterr("Error: No data detected in %s\n", buffer.Filename().full());
    return 1;
  }
  DataSet* ds = datasetlist.AddSet(DataSet::MATRIX_DBL, dsname, "Mat");
  if (ds == 0) return 1;
  DataSet_MatrixDbl& Mat = static_cast<DataSet_MatrixDbl&>( *ds );
  Mat.SetTypeAndKind( DataSet_2D::DIST, DataSet_2D::FULL ); // TODO: FIXME 
  Mat.Allocate2D( ncols, nrows );
  std::copy( matrixArray.begin(), matrixArray.end(), Mat.begin() );

  return 0;
}

// DataIO_Std::Read_3D()
int DataIO_Std::Read_3D(std::string const& fname, 
                        DataSetList& datasetlist, std::string const& dsname)
{
  return 1;
}

// DataIO_Std::Read_Vector()
int DataIO_Std::Read_Vector(std::string const& fname, 
                            DataSetList& datasetlist, std::string const& dsname)
{
  // Buffer file
  BufferedLine buffer;
  if (buffer.OpenFileRead( fname )) return 1;
  mprintf("\tData will be read as a vector.\n");
  // Skip comments
  const char* linebuffer = buffer.Line();
  while (linebuffer != 0 && linebuffer[0] == '#')
    linebuffer = buffer.Line();
  double vecBuffer[6];
  std::fill(vecBuffer, vecBuffer+6, 0.0);
  size_t ndata = 0;
  int ncols = -1; // Should be 3, 6, or 9
  int nv = 0;     // Will be set to 3 or 6
  bool hasIndex = false;
  DataSet* vec = 0;
  while (linebuffer != 0) {
    int ntokens = buffer.TokenizeLine( SEPARATORS );
    if (ncols < 0) {
      ncols = ntokens;
      if (ntokens < 1) {
        mprinterr("Error: Could not tokenize line.\n");
        return 1;
      }
      if (ncols == 3 || ncols == 6 || ncols == 9)
        hasIndex = false;
      else if (ncols == 4 || ncols == 7 || ncols == 10) {
        hasIndex = true;
        mprintf("Warning: Not reading vector data indices.\n");
      } else {
        mprinterr("Error: Expected 3, 6, or 9 columns of vector data, got %i.\n", ncols);
        return 1;
      }
      if (ncols >= 6)
        nv = 6;
      else
        nv = 3;
      vec = datasetlist.AddSet(DataSet::VECTOR, dsname, "Vec");
    } else if (ncols != ntokens) {
      mprinterr("Error: In vector file, number of columns changes from %i to %i at line %i\n",
                ncols, ntokens, buffer.LineNumber());
      return 1;
    }
    if (hasIndex)
      buffer.NextToken(); // Skip index
    for (int i = 0; i < nv; i++)
      vecBuffer[i] = atof( buffer.NextToken() );
    vec->Add( ndata, vecBuffer ); 
    ndata++;
    linebuffer = buffer.Line();
  } 
  return 0;
} 

// -----------------------------------------------------------------------------
void DataIO_Std::WriteHelp() {
  mprintf("\tinvert:     Flip X/Y axes.\n"
          "\tnoxcol:     Do not print X (index) column.\n"
          "\tnoheader:   Do not print header line.\n"
          "\tsquare2d:   Write 2D data sets in matrix-like format.\n"
          "\tnosquare2d: Write 2D data sets as '<X> <Y> <Value>'.\n");
}

// DataIO_Std::processWriteArgs()
int DataIO_Std::processWriteArgs(ArgList &argIn) {
  isInverted_ = argIn.hasKey("invert");
  hasXcolumn_ = !argIn.hasKey("noxcol");
  writeHeader_ = !argIn.hasKey("noheader");
  square2d_ = argIn.hasKey("square2d");
  if (argIn.hasKey("nosquare2d")) square2d_ = false;
  return 0;
}

// WriteNameToBuffer()
void DataIO_Std::WriteNameToBuffer(CpptrajFile& fileIn, std::string const& label,
                                   int width,  bool leftAlign) 
{
  std::string temp_name = label;
  // If left aligning, add '#' to name; 
  if (leftAlign) {
    if (temp_name[0]!='#') {
      temp_name.insert(0,"#");
      // Ensure that name will not be larger than column width.
      if ((int)temp_name.size() > width)
        temp_name.resize( width );
    }
  }
  // Replace any spaces with underscores
  for (std::string::iterator tc = temp_name.begin(); tc != temp_name.end(); ++tc)
    if ( *tc == ' ' )
      *tc = '_';
  // Set up header format string
  std::string header_format = SetStringFormatString(width, leftAlign);
  fileIn.Printf(header_format.c_str(), temp_name.c_str());
}

// DataIO_Std::WriteData()
int DataIO_Std::WriteData(std::string const& fname, DataSetList const& SetList)
{
  int err = 0;
  // Open output file.
  CpptrajFile file;
  if (file.OpenWrite( fname )) return 1;
  if (isInverted_)
    err = WriteDataInverted(file, SetList);
  else
    err = WriteDataNormal(file, SetList);
  file.CloseFile();
  return err;
}

// DataIO_Std::WriteDataNormal()
int DataIO_Std::WriteDataNormal(CpptrajFile& file, DataSetList const& SetList) {
  std::string x_col_format;

  // Hold all 1D data sets.
  Array1D Sets( SetList );
  if (Sets.empty()) return 1;
  // For this output to work the X-dimension of all sets needs to match.
  // The most important things for output are min and step so just check that.
  // Use X dimension of set 0 for all set dimensions.
  Sets.CheckXDimension();
  // TODO: Check for empty dim.
  DataSet_1D const& Xdata = static_cast<DataSet_1D const&>( *Sets[0] );
  Dimension const& Xdim = static_cast<Dimension const&>( Xdata.Dim(0) );
  int xcol_width = Xdim.Label().size() + 1;
  if (xcol_width < 8) xcol_width = 8;
  int xcol_precision = 3;

  // Determine size of largest DataSet.
  size_t maxFrames = Sets.DetermineMax();

  // Set up X column.
  if (hasXcolumn_) {
    // Create format string for X column based on dimension in first data set.
    if (Xdim.Step() == 1.0) xcol_precision = 0;
    x_col_format = SetupCoordFormat( maxFrames, Xdim, xcol_width, xcol_precision ); 
  } else {
    // If not writing an X-column, set the format for the first dataset
    // to left-aligned.
    Sets[0]->SetDataSetFormat( true );
  }

  // Write header to buffer
  if (writeHeader_) {
    // If x-column present, write x-label
    if (hasXcolumn_)
      WriteNameToBuffer( file, Xdim.Label(), xcol_width, true );
    // To prevent truncation of DataSet legends, adjust the width of each
    // DataSet if necessary.
    bool labelLeftAligned = !hasXcolumn_;
    for (Array1D::const_iterator ds = Sets.begin(); ds != Sets.end(); ++ds) {
      int requiredColSize = (int)(*ds)->Legend().size();
      if (!labelLeftAligned || (ds == Sets.begin() && !hasXcolumn_))
        requiredColSize++;
      if ( requiredColSize > (*ds)->ColumnWidth() )
        (*ds)->SetWidth( (*ds)->Legend().size() );
      labelLeftAligned = false;
    }
    // Write dataset names to header, left-aligning first set if no X-column
    Array1D::const_iterator set = Sets.begin();
    if (!hasXcolumn_)
      WriteNameToBuffer( file, (*set)->Legend(), (*set)->ColumnWidth(), true  );
    else
      WriteNameToBuffer( file, (*set)->Legend(), (*set)->ColumnWidth(), false );
    ++set;
    for (; set != Sets.end(); ++set) 
      WriteNameToBuffer( file, (*set)->Legend(), (*set)->ColumnWidth(), false );
    file.Printf("\n"); 
  }

  // Write Data
  for (size_t frame=0L; frame < maxFrames; frame++) {
    // Output Frame for each set
    if (hasXcolumn_)
      file.Printf(x_col_format.c_str(), Xdata.Xcrd(frame));
    for (Array1D::const_iterator set = Sets.begin(); set != Sets.end(); ++set) 
      (*set)->WriteBuffer(file, frame);
    file.Printf("\n"); 
  }
  return 0;
}

// DataIO_Std::WriteDataInverted()
int DataIO_Std::WriteDataInverted(CpptrajFile& file, DataSetList const& SetList)
{
  // Hold all 1D data sets.
  Array1D Sets( SetList );
  if (Sets.empty()) return 1;
  // Determine size of largest DataSet.
  size_t maxFrames = Sets.DetermineMax();
  // Write each set to a line.
  for (Array1D::const_iterator set = Sets.begin(); set != Sets.end(); ++set) {
    // Write dataset name as first column.
    WriteNameToBuffer( file, (*set)->Legend(), (*set)->ColumnWidth(), false); 
    // Write each frame to subsequent columns
    for (size_t frame=0L; frame < maxFrames; frame++) 
      (*set)->WriteBuffer(file, frame);
    file.Printf("\n");
  }
  return 0;
}

// DataIO_Std::WriteData2D()
int DataIO_Std::WriteData2D( std::string const& fname, DataSetList const& setList) 
{
  // Open output file
  CpptrajFile file;
  if (file.OpenWrite( fname )) return 1;
  int err = 0;
  for (DataSetList::const_iterator set = setList.begin(); set != setList.end(); ++set)
  {
    if (set != setList.begin()) file.Printf("\n");
    err += WriteSet2D( *(*set), file );
  }
  file.CloseFile();
  return err;
}

// DataIO_Std::WriteSet2D()
int DataIO_Std::WriteSet2D( DataSet const& setIn, CpptrajFile& file ) {
  if (setIn.Ndim() != 2) {
    mprinterr("Internal Error: DataSet %s in DataFile %s has %zu dimensions, expected 2.\n",
              setIn.Legend().c_str(), file.Filename().full(), setIn.Ndim());
    return 1;
  }
  DataSet_2D const& set = static_cast<DataSet_2D const&>( setIn );
  int xcol_width = 8;
  int xcol_precision = 3;
  Dimension const& Xdim = static_cast<Dimension const&>(set.Dim(0));
  Dimension const& Ydim = static_cast<Dimension const&>(set.Dim(1));
  if (Xdim.Step() == 1.0) xcol_precision = 0;
  
  if (square2d_) {
    std::string ycoord_fmt;
    // Print XY values in a grid:
    // x0y0 x1y0 x2y0
    // x0y1 x1y1 x2y1
    // x0y2 x1y2 x2y2
    // If file has header, top-left value will be '#<Xlabel>-<Ylabel>',
    // followed by X coordinate values.
    if (writeHeader_) {
      ycoord_fmt = SetupCoordFormat( set.Nrows(), Ydim, xcol_width, xcol_precision );
      std::string header;
      if (Xdim.Label().empty() && Ydim.Label().empty())
        header = "#Frame";
      else
        header = "#" + Xdim.Label() + "-" + Ydim.Label();
      WriteNameToBuffer( file, header, xcol_width, true );
      std::string xcoord_fmt = SetupCoordFormat( set.Ncols(), Xdim, set.ColumnWidth(), 
                                                 xcol_precision );
      for (size_t ix = 0; ix < set.Ncols(); ix++)
        file.Printf(xcoord_fmt.c_str(), Xdim.Coord( ix ));
      file.Printf("\n");
    }
    for (size_t iy = 0; iy < set.Nrows(); iy++) {
      if (writeHeader_)
        file.Printf(ycoord_fmt.c_str(), Ydim.Coord( iy ));
      for (size_t ix = 0; ix < set.Ncols(); ix++)
        set.Write2D( file, ix, iy );
      file.Printf("\n");
    }
  } else {
    // Print X Y Values
    // x y val(x,y)
    if (writeHeader_)
      file.Printf("#%s %s %s\n", Xdim.Label().c_str(), 
                  Ydim.Label().c_str(), set.Legend().c_str());
    std::string col_fmt = SetupCoordFormat( set.Ncols(), Xdim, 8, 3 ) + " " +
                          SetupCoordFormat( set.Nrows(), Ydim, 8, 3 ); 
    for (size_t iy = 0; iy < set.Nrows(); ++iy) {
      for (size_t ix = 0; ix < set.Ncols(); ++ix) {
        file.Printf(col_fmt.c_str(), Xdim.Coord( ix ), Ydim.Coord( iy ));
        set.Write2D( file, ix, iy );
        file.Printf("\n");
      }
    }
  }
  return 0;
}

// DataIO_Std::WriteData3D()
int DataIO_Std::WriteData3D( std::string const& fname, DataSetList const& setList) 
{
  // Open output file
  CpptrajFile file;
  if (file.OpenWrite( fname )) return 1;
  int err = 0;
  for (DataSetList::const_iterator set = setList.begin(); set != setList.end(); ++set)
  {
    if (set != setList.begin()) file.Printf("\n");
    err += WriteSet3D( *(*set), file );
  }
  file.CloseFile();
  return err;

}

// DataIO_Std::WriteSet3D()
int DataIO_Std::WriteSet3D( DataSet const& setIn, CpptrajFile& file ) {
  if (setIn.Ndim() != 3) {
    mprinterr("Internal Error: DataSet %s in DataFile %s has %zu dimensions, expected 3.\n",
              setIn.Legend().c_str(), file.Filename().full(), setIn.Ndim());
    return 1;
  }
  DataSet_3D const& set = static_cast<DataSet_3D const&>( setIn );
  Dimension const& Xdim = static_cast<Dimension const&>(set.Dim(0));
  Dimension const& Ydim = static_cast<Dimension const&>(set.Dim(1));
  Dimension const& Zdim = static_cast<Dimension const&>(set.Dim(2));
  //if (Xdim.Step() == 1.0) xcol_precision = 0;
  
  // Print X Y Z Values
  // x y z val(x,y,z)
  if (writeHeader_)
    file.Printf("#%s %s %s %s\n", Xdim.Label().c_str(), 
                Ydim.Label().c_str(), Zdim.Label().c_str(), set.Legend().c_str());
  std::string col_fmt = SetupCoordFormat( set.NX(), Xdim, 8, 3 ) + " " +
                        SetupCoordFormat( set.NY(), Ydim, 8, 3 ) + " " +
                        SetupCoordFormat( set.NZ(), Zdim, 8, 3 );
  for (size_t iz = 0; iz < set.NZ(); ++iz) {
    for (size_t iy = 0; iy < set.NY(); ++iy) {
      for (size_t ix = 0; ix < set.NX(); ++ix) {
        file.Printf(col_fmt.c_str(), Xdim.Coord( ix ), Ydim.Coord( iy ), Zdim.Coord( iz ));
        set.Write3D( file, ix, iy, iz );
        file.Printf("\n");
      }
    }
  }
  return 0;
}
