#include <cstdio> // sscanf
#include "DataIO_Evecs.h"
#include "CpptrajStdio.h"
#include "BufferedFrame.h"
#include "DataSet_Modes.h"

// CONSTRUCTOR
DataIO_Evecs::DataIO_Evecs() : ibeg_(1), iend_(50), hasIend_(false) {
  SetValid( DataSet::MODES );
}

bool DataIO_Evecs::ID_DataFormat(CpptrajFile& infile) {
  if (infile.OpenFile()) return false;
  std::string firstLine = infile.GetLine();
  infile.CloseFile();
  return ( firstLine.compare(0, 18," Eigenvector file:") == 0 );
}

void DataIO_Evecs::ReadHelp() {
  mprintf("\tibeg <firstmode>: First mode to read in (default 1).\n"
          "\tiend <lastmode>:  Last mode to read in (default 50).\n");
}

int DataIO_Evecs::processReadArgs(ArgList& argIn) {
  ibeg_ = argIn.getKeyInt("ibeg",1);
  hasIend_ = argIn.Contains("iend");
  iend_ = argIn.getKeyInt("iend",50);
  if (iend_ < 1 || ibeg_ < 1) {
    mprinterr("Error: iend and ibeg must be > 0\n");
    return 1;
  }
  if (iend_ < ibeg_) {
    mprinterr("Error: iend cannot be less than ibeg\n");
    return 1;
  }
  return 0;
}

// DataIO_Evecs::ReadData()
int DataIO_Evecs::ReadData(std::string const& modesfile,
                           DataSetList& datasetlist, std::string const& dsname)
{
  // Process Arguments
  int modesToRead = iend_ - ibeg_ + 1;
  if (hasIend_)
    mprintf("\tAttempting to read %i modes (%i to %i) from %s\n", modesToRead,
            ibeg_, iend_, modesfile.c_str());
  else
    mprintf("\tReading modes from %s\n", modesfile.c_str());
  BufferedFrame infile;
  if (infile.OpenRead( modesfile)) return 1;
  // Read title line, convert to arg list
  const char* buffer = 0;
  if ( (buffer = infile.NextLine())==0 ) {
    mprinterr("Error: ReadEvecFile(): error while reading title (%s)\n",infile.Filename().full());
    return 1;
  }
  ArgList title(buffer);
  // Check if reduced
  bool reduced = title.hasKey("Reduced");
  // Allocate MODES dataset. No appending allowed.
  DataSet* mds = datasetlist.AddSet( DataSet::MODES, dsname, "Evecs" );
  if (mds == 0) return 1;
  DataSet_Modes& modesData = static_cast<DataSet_Modes&>( *mds );
  // Determine modes file type
  DataSet_2D::MatrixType modesType = DataSet_2D::NO_OP;
  for (int matidx = (int)DataSet_2D::NO_OP + 1;
           matidx != (int)DataSet_2D::NMAT; ++matidx)
  {
    if ( title.hasKey( DataSet_2D::MatrixOutputString((DataSet_2D::MatrixType)matidx) ))
    {
      modesType = (DataSet_2D::MatrixType)matidx;
      break;
    }
  }
  // For compatibility with quasih and nmode output
  if (modesType == DataSet_2D::NO_OP) {
    mprintf("Warning: ReadEvecFile(): Unrecognized type [%s]\n", title.ArgLine());
    mprintf("         Assuming MWCOVAR.\n");
    modesType = DataSet_2D::MWCOVAR;
  }
  modesData.SetType( modesType );
  // For newer modesfiles, get # of modes in file.
  int modesInFile = title.getKeyInt("nmodes",-1);
  if (modesInFile == -1) {
    modesInFile = modesToRead;
    mprintf("Warning: Older modes file, # of modes not known.\n"
            "Warning: Will try to read at least %i modes.\n", modesToRead);
  } else {
    mprintf("\tFile contains %i modes.\n", modesInFile);
    if (modesToRead > modesInFile) {
      mprintf("Warning: # modes to read (%i) > modes in file. Only reading %i modes.\n",
              modesToRead, modesInFile);
      modesToRead = modesInFile;
    } else if (!hasIend_ && modesToRead < modesInFile) {
      modesToRead = modesInFile;
      iend_ = modesInFile;
    }
  }
  // For newer modesfiles, get width of data elts
  int colwidth = title.getKeyInt("width", -1);
  if (colwidth == -1)
    colwidth = 11; // Default, 10 + 1 space
  modesData.SetPrecision(colwidth - 1, 5);
  // Read number of elements in avg coords and eigenvectors
  if ( (buffer = infile.NextLine())==0 ) {
    mprinterr("Error: ReadEvecFile(): error while reading number of atoms (%s)\n",
              infile.Filename().full());
    return 1;
  }
  int navgcrd = 0;
  int vecsize = 0;
  int nvals = sscanf(buffer, "%i %i", &navgcrd, &vecsize);
  if (nvals == 0) {
    mprinterr("Error: ReadEvecFile(): sscanf on coords failed (%s)\n",infile.Filename().full());
    return 1;
  } else if (nvals == 1) {
    mprintf("Warning: ReadEvecFile(): No value for eigenvector size found in %s,\n",
            infile.Filename().full());
    mprintf("         assuming it is equal to #average elements (%i)\n",navgcrd);
    vecsize = navgcrd;
  }
  // Allocate FrameBuffer
  int bufsize;
  if (navgcrd > vecsize)
    bufsize = navgcrd;
  else
    bufsize = vecsize;
  infile.SetupFrameBuffer( bufsize, colwidth, 7 );
  // Allocate memory for avg coords and read in
  modesData.AllocateAvgCoords( navgcrd );
  if (navgcrd > 0) {
    infile.ReadFrame( );
    infile.BufferToDouble( modesData.AvgFramePtr(), modesData.NavgCrd() );
    infile.BufferBegin(); // Reset buffer position
  }
  // Allocate buffer memory for eigenvalues and eigenvectors
  double* evalues = new double[ modesToRead ];
  double* evectors = 0;
  if (vecsize > 0)
    evectors = new double[ modesToRead * vecsize ];
  int nmodes = 0;      // Mode currently in memory
  int currentMode = 0; // Modes currently reading in from file
  bool firstRead = true;
  int error_status = 0;
  while ( (buffer = infile.NextLine())!=0 ) { // This should read in ' ****'
    if (buffer[0] != ' ' || buffer[1] != '*' || buffer[2] != '*') {
      mprinterr("Error: ReadEvecFile(): When reading eigenvector %i, expected ' ****',\n",
                currentMode+1);
      mprinterr("Error: got %s [%s]\n", buffer, infile.Filename().full());
      error_status = 1;
      break;
    }
    // Read number and eigenvalue 
    if ( (buffer = infile.NextLine())==0 ) {
      mprinterr("Error: ReadEvecFile(): error while reading number and eigenvalue (%s)\n",
                infile.Filename().full());
      error_status = 2;
      break;
    }
    if (sscanf(buffer, "%*i%lf", evalues + nmodes) != 1) {
      mprinterr("Error: ReadEvecFile(): error while scanning number and eigenvalue (%s)\n",
                infile.Filename().full());
      error_status = 3;
      break;
    }
    if (vecsize > 0) {
      // Read eigenvector
      // Older modesfiles could have vecsize > 0 but no eigenvectors, only 
      // blanks. If this is the case set vecsize to -1 to indicate a blank
      // read is needed after reading eigenvalues.
      double* Vec = evectors + (nmodes * vecsize);
      int vi = 0;
      while (vi < vecsize) {
        buffer = infile.NextLine();
        if (firstRead && (buffer[0] == '\n' || buffer[0] == '\r')) {
          mprintf("Warning: Old modes file with vecsize > 0 but no eigenvectors.\n");
          vecsize = -1;
          delete[] evectors;
          evectors = 0;
          break;
        }
        double tmpval[7];
        int nvals = sscanf(buffer, "%lf %lf %lf %lf %lf %lf %lf", tmpval,
                           tmpval+1, tmpval+2, tmpval+3, tmpval+4, tmpval+5, tmpval+6);
        for (int ti = 0; ti < nvals; ++ti)
          Vec[vi++] = tmpval[ti];
      }
      // Check if mode read was between ibeg and iend (which start from 1).
      // If so, increment number of modes.
      if (currentMode+1 >= ibeg_ && currentMode < iend_) ++nmodes;
      if (nmodes == modesToRead) break;
      ++currentMode;
    } else if (vecsize == -1) {
      // Blank read past empty eigenvector
      buffer = infile.NextLine();
    }
    firstRead = false;
  }
  infile.CloseFile();
  if (error_status == 0) {
    if (nmodes != modesToRead)
      mprintf("Warning: Number of read modes is %i, requested %i\n", nmodes, modesToRead);
    error_status = modesData.SetModes( reduced, nmodes, vecsize, evalues, evectors );
  }
  delete[] evalues;
  if (evectors != 0) delete[] evectors;
  return error_status;
}

// DataIO_Evecs::WriteData()
int DataIO_Evecs::WriteData(std::string const& fname, DataSetList const& SetList) {
  if (SetList.empty()) return 1;
  if (SetList.size() > 1)
    mprintf("Warning: Multiple sets not yet supported for Evecs write.\n");
  DataSet_Modes const& modesData = static_cast<DataSet_Modes const&>( *(*(SetList.begin())) );
  BufferedFrame outfile;
  if (outfile.OpenWrite( fname )) {
    mprinterr("Error: Could not open %s for writing.\n", fname.c_str());
    return 1;
  }
  if (modesData.IsReduced())
    outfile.Printf(" Reduced Eigenvector file: ");
  else
    outfile.Printf(" Eigenvector file: ");
  outfile.Printf("%s", DataSet_2D::MatrixOutputString(modesData.Type()));
  // Write out # of modes on title line to not break compat. with older modes files
  outfile.Printf(" nmodes %i", modesData.Nmodes());
  // Write out col width on title line to not break compat. with older modes files
  int colwidth = modesData.ColumnWidth();
  outfile.Printf(" width %i\n", colwidth);
  // First number is # avg coords, second is size of each vector
  outfile.Printf(" %4i %4i\n", modesData.NavgCrd(), modesData.VectorSize());
  // Set up framebuffer, default 7 columns
  int bufsize;
  if (modesData.NavgCrd() > modesData.VectorSize())
    bufsize = modesData.NavgCrd();
  else
    bufsize = modesData.VectorSize();
  outfile.SetupFrameBuffer( bufsize, colwidth, 7 );
  // Print average coords
  outfile.DoubleToBuffer( modesData.AvgFramePtr(), modesData.NavgCrd(),
                          modesData.DataFormat() );
  outfile.WriteFrame();
  // Eigenvectors and eigenvalues
  for (int mode = 0; mode < modesData.Nmodes(); ++mode) {
    outfile.Printf(" ****\n %4i ", mode+1);
    outfile.Printf(modesData.DataFormat(), modesData.Eigenvalue(mode));
    outfile.Printf("\n");
    if (modesData.Eigenvectors() != 0) {
      const double* Vec = modesData.Eigenvector(mode);
      outfile.BufferBegin();
      outfile.DoubleToBuffer( Vec, modesData.VectorSize(), modesData.DataFormat() );
      outfile.WriteFrame();
    }
  }
  outfile.CloseFile();
  return 0;
}
