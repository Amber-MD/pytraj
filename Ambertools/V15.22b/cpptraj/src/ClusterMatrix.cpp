#include <cfloat> // FLT_MAX
#include "ClusterMatrix.h"
#include "CpptrajFile.h"
#include "CpptrajStdio.h"

// NOTES:
//   Version 1: Add write of ignore array when reduced. Write nrows and
//              and nelements as 8 byte integers.
//   Version 2: Instead of nrows and nelements, write original nrows
//              and actual nrows to easily determine if this is a reduced
//              matrix. Also write sieve value.
//   Version 2 Update: Read/write sieve value as signed, negative
//                     value is random sieve. Variable is same #
//                     of bytes so should be backwards-compatible.
const unsigned char ClusterMatrix::Magic_[4] = {'C', 'T', 'M', 2};

// CONSTRUCTOR
/** Intended for use with cluster pairwise distance calculations
  * where frames may be sieved. The underlying TriangleMatrix will
  * only be set up to hold the actual number of frames based on
  * the sieve value. The Ignore array will be set up based on
  * the original number of frames.
  */
int ClusterMatrix::SetupWithSieve(size_t sizeIn, size_t sieveIn, int iseed)
{
  if (sievedFrames_.SetSieve( sieveIn, sizeIn, iseed )) return 1;
  // Sieved distances should be ignored.
  if (sievedFrames_.Type() != ClusterSieve::NONE) {
    // Set up the ignore array to ignore sieved frames
    ignore_.assign(sizeIn, true);
    size_t actual_nrows = 0;
    for (size_t frame = 0; frame < sizeIn; frame++)
      if (sievedFrames_.FrameToIdx(frame) != -1) {
        ignore_[frame] = false;
        ++actual_nrows;
      }
    // Set up underlying TriangleMatrix for sieved frames.
    Mat_.resize( 0L, actual_nrows );
    mprintf("\tPair-wise matrix set up with sieve, %zu frames, %zu sieved frames.\n",
            sizeIn, actual_nrows);
  } else {
    Mat_.resize( 0L, sizeIn );
    ignore_.assign(sizeIn, false);
    mprintf("\tPair-wise matrix set up, %zu frames\n", sizeIn);
  }
  return 0;
}

// COPY CONSTRUCTOR
ClusterMatrix::ClusterMatrix(const ClusterMatrix& rhs) :
  ignore_(rhs.ignore_),
  Mat_(rhs.Mat_)
{}

// ASSIGNMENT
ClusterMatrix& ClusterMatrix::operator=(const ClusterMatrix& rhs) {
  if (this == &rhs) return *this;
  ignore_ = rhs.ignore_;
  Mat_ = rhs.Mat_;
  return *this;
}

// ClusterMatrix::SaveFile()
/** Save the matrix to a binary file. Format is:
  *     Full: [4*char]     [uint_8] [uint_8]    [uint_8] [nelements*float]
  *  Reduced: [4*char]     [uint_8] [uint_8]    [uint_8] [nelements*float] [nrows*char]
  *     Vars: ['C''T''M'X] [nrows]  [nelements] [sieve]  [elements]        [ignore]
  */
int ClusterMatrix::SaveFile(std::string const& filename) const {
  CpptrajFile outfile;
  uint_8 ntemp;
  // No stdout write allowed.
  if (filename.empty()) {
    mprinterr("Internal Error: ClusterMatrix::SaveFile called with no filename.\n");
    return 1;
  }
  if (outfile.OpenWrite(filename)) {
    mprinterr("Error: ClusterMatrix::SaveFile: Could not open %s for write.\n", filename.c_str());
    return 1;
  }
  // Write magic byte
  outfile.Write( Magic_, 4 );
  // Write original nrows (size of ignore)
  ntemp = (uint_8)ignore_.size();
  outfile.Write( &ntemp, sizeof(uint_8) );
  // Write actual nrows
  ntemp = (uint_8)Mat_.Nrows();
  outfile.Write( &ntemp, sizeof(uint_8) );
  // Write out sieve value
  sint_8 stemp = (sint_8)sievedFrames_.Sieve();
  outfile.Write( &stemp, sizeof(sint_8) );
  // Write matrix elements
  outfile.Write( Mat_.Ptr(), Mat_.size()*sizeof(float) );
  // If this is a reduced matrix, write the ignore array as chars.
  if (sievedFrames_.Type() != ClusterSieve::NONE) {
    char* ignore_out = new char[ ignore_.size() ];
    int idx = 0;
    for (std::vector<bool>::const_iterator ig = ignore_.begin(); ig != ignore_.end(); ++ig) 
      if (*ig)
        ignore_out[idx++] = 'T';
      else
        ignore_out[idx++] = 'F';
    outfile.Write( ignore_out, ignore_.size()*sizeof(char) );
    delete[] ignore_out;
  }
  return 0;
}

// ClusterMatrix::LoadFile()
int ClusterMatrix::LoadFile(std::string const& filename, int sizeIn) {
  unsigned char magic[4];
  CpptrajFile infile;
  uint_8 ROWS, ELTS;
  sint_8 SIEVE;
  int sieve = 1;
  size_t actual_nrows = 0;
  // Open file for reading
  if (infile.OpenRead(filename)) {
    mprinterr("Error: ClusterMatrix::LoadFile: Could not open %s for read.\n", filename.c_str());
    return 1;
  }
  // Read and check magic byte
  infile.Read( magic, 4 );
  if ( magic[0]!=Magic_[0] || magic[1]!=Magic_[1] || magic[2]!=Magic_[2] ) {
    mprinterr("Error: ClusterMatrix::LoadFile: File %s is not a TriangleMatrix.\n",
              filename.c_str());
    return 1;
  }
  // Check version, read in nrows and nelements.
  if (magic[3] == 0) {
    int Ntemp = 0;
    infile.Read( &Ntemp, sizeof(int) );
    ROWS = (uint_8)Ntemp;
    actual_nrows = (size_t)ROWS;
    infile.Read( &Ntemp, sizeof(int) );
    ELTS = (uint_8)Ntemp;
  } else if (magic[3] == 1) {
    infile.Read( &ROWS, sizeof(uint_8) );
    actual_nrows = (size_t)ROWS;
    infile.Read( &ELTS, sizeof(uint_8) );
  } else if (magic[3] == 2) {
    infile.Read( &ROWS,  sizeof(uint_8) ); // V2: Original Nrows
    infile.Read( &ELTS,  sizeof(uint_8) ); // V2: Actual Nrows
    actual_nrows = (size_t)ELTS;
    infile.Read( &SIEVE, sizeof(sint_8) ); // V2: Sieve
    sieve = (int)SIEVE; 
  } else {
    mprinterr("Error: ClusterMatrix version %u is not recognized.\n", (unsigned int)magic[3]);
    return 1;
  }
  // If number of rows is not what was expected, abort
  if (ROWS != (uint_8)sizeIn) {
    mprinterr("Error: ClusterMatrix file %s has %lu rows, expected %i.\n",
              filename.c_str(), ROWS, sizeIn);
    return 1;
  }
  if (magic[3] == 0 || magic[3] == 1) {
    // Version 0/1: Actual # of rows is not known yet. Check that the # elements
    // in the file match the original # elements (i.e. matrix is not sieved).
    // If it is sieved this is not supported.
    uint_8 original_nelements = ( ROWS * (ROWS - 1UL) ) / 2UL;
    if ( original_nelements != ELTS ) {
      mprinterr("Error: Sieved data in ClusterMatrix file %s (version %u) not supported.\n",
                filename.c_str(), (unsigned int)magic[3]);
      return 1;
    }
    sieve = 1;
  }
  // Setup underlying TriangleMatrix for actual # of rows
  if ( Mat_.resize( 0L, actual_nrows ) ) return 1;
  // Set all ignore elements for original # rows to false.
  ignore_.assign(ROWS, false);
  // Read in matrix elements
  infile.Read( Mat_.Ptr(), Mat_.size()*sizeof(float) );
  // If sieved, read in the ignore array
  if (sieve != 1) {
    mprintf("Warning: ClusterMatrix %s contains sieved data.\n", filename.c_str());
    char* ignore_in = new char[ ROWS ]; // Original nrows
    infile.Read( ignore_in, ROWS*sizeof(char) );
    for (uint_8 row = 0; row < ROWS; ++row)
      if (ignore_in[row] == 'T')
        ignore_[row] = true;
    delete[] ignore_in;
  }
  // Setup sieve class
  if (sievedFrames_.SetSieve( sieve, ignore_ )) {
    mprinterr("Error: Could not set sieve from ClusterMatrix file.\n");
    return 1;
  }
  mprintf("\tLoaded %s: %u original rows, %u actual rows, %u elements, sieve=%i\n",
          filename.c_str(), ROWS, Mat_.Nrows(), Mat_.size(), sieve);
  return 0;
}

// ClusterMatrix::SetupMatrix()
int ClusterMatrix::SetupMatrix(size_t sizeIn) {
  if (Mat_.resize( 0L, sizeIn )) return 1;
  ignore_.assign( sizeIn, false );
  //sieve_ = 1;
  return 0;
}

// ClusterMatrix::FindMin()
/** Find the minimum; set corresponding row and column. Cannot currently
  * be used for sieved frames.
  */
double ClusterMatrix::FindMin(int& iOut, int& jOut) const {
  iOut = -1;
  jOut = -1;
  int iVal = 0;
  int jVal = 1;
  float min = FLT_MAX;
  for (size_t idx = 0UL; idx < Nelements(); ++idx) {
    if (ignore_[iVal] || ignore_[jVal]) {
      // If we dont care about this row/col, just increment
      jVal++;
      if (jVal >= (int)ignore_.size()) {
        iVal++;
        jVal = iVal + 1;
      }
    } else {
      // Otherwise search for minimum
      if ( Mat_[idx] < min ) {
        min = Mat_[idx];
        iOut = iVal;
        jOut = jVal;
      }
      // Increment indices
      jVal++;
      if (jVal >= (int)ignore_.size()) {
        iVal++;
        jVal = iVal + 1;
      }
    }
  }
  return (double)min;
}

// ClusterMatrix::PrintElements()
void ClusterMatrix::PrintElements() const {
  if (sievedFrames_.MaxFrames()==0) {
    // This is ClusterDistances matrix. Ignore size is # of sieved frames, sievedFrames is 0.
    unsigned int iVal = 0;
    unsigned int jVal = 1;
    for (size_t idx = 0UL; idx < Nelements(); ++idx) {
      if (!ignore_[iVal] && !ignore_[jVal])
        mprintf("\t%u %u %f\n",iVal,jVal,Mat_[idx]);
      // Increment indices
      jVal++;
      if (jVal >= ignore_.size()) {
        iVal++;
        jVal = iVal + 1;
      }
    }
  } else {
    // This is FrameDistances matrix. Ignore and sievedFrames size is # of original frames.
    for (unsigned int row = 0; row != Nframes(); row++)
      for (unsigned int col = row + 1; col != Nframes(); col++)
        if (!ignore_[row] && !ignore_[col])
          mprintf("\t%u %u %f\n", row+1, col+1, GetFdist(col, row));
  }
}

size_t ClusterMatrix::DataSize() const {
  return ( Mat_.DataSize() +
           (ignore_.capacity()*sizeof(bool) + sizeof(ignore_)) +
           sievedFrames_.DataSize() );
}
