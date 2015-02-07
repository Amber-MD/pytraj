// Traj_CharmmDcd
#include <cmath>   // for cos, acos
#include <cstring> // memset
#include "Traj_CharmmDcd.h"
#include "Constants.h"
#include "CpptrajStdio.h"
#include "ByteRoutines.h"

// CONSTRUCTOR
Traj_CharmmDcd::Traj_CharmmDcd() :
  dcdatom_(0),
  dcdframes_(0),
  isBigEndian_(false),
  is64bit_(false),
  blockSize_(4),
  dcd_dim_(3),
  boxBytes_(0),
  frame1Bytes_(0),
  frameNBytes_(0),
  coordinate_size_(0),
  nfixedat_(0),
  nfreeat_(0),
  freeat_(0),
  xcoord_(0),
  ycoord_(0),
  zcoord_(0)
{}

// DESTRUCTOR
Traj_CharmmDcd::~Traj_CharmmDcd() {
  if (freeat_!=0) delete[] freeat_;
  if (xcoord_!=0) delete[] xcoord_;
}

static inline bool CORD_32BIT(const unsigned char* buffer) {
  if (buffer[4] == 'C' && buffer[5] == 'O' && buffer[ 6] == 'R' && buffer[ 7] == 'D')
    return true; // 32 bit
  return false;
}

static inline bool CORD_64BIT(const unsigned char* buffer) {
  if (buffer[8] == 'C' && buffer[9] == 'O' && buffer[10] == 'R' && buffer[11] == 'D')
    return true; // 64 bit
  return false;
}

/** Charmm DCD is 32 or 64 bit int (84) followed by 'CORD'. Determine
  * endianness and bit size.
  */
bool Traj_CharmmDcd::ID_TrajFormat(CpptrajFile& fileIn) {
  byte8 firstbyte;
  unsigned char buffer[12];
  memset(buffer, ' ', 12);
  if (fileIn.OpenFile()) return false;
  if (fileIn.Read(buffer, 12) != 12) return false;
  fileIn.CloseFile();
  // If the second 4 or last 4 chars are C O R D, charmm DCD
  if (CORD_32BIT(buffer)) { 
    is64bit_ = false;
    firstbyte.i[1] = 0;
    blockSize_ = 4; 
  } else if (CORD_64BIT(buffer)) {
    is64bit_ = true;
    blockSize_ = 8;
  } else
    return false; // CORD not found.
  // Determine 32/64 bit
  memcpy(firstbyte.c, buffer, blockSize_ * sizeof(unsigned char));
  if (firstbyte.i[0] == 84)
    isBigEndian_ = false;
  else {
    // Swap bytes
    if (is64bit_)
      //endian_swap( firstbyte.i, 2 );
      endian_swap8( firstbyte.i, 1 );
    else
      endian_swap( firstbyte.i, 1 );
    if (firstbyte.i[0] == 84)
      isBigEndian_ = true;
    else
      return false;
  }
  return true; 
}

// Traj_CharmmDcd::openTrajin()
int Traj_CharmmDcd::openTrajin() {
  // Always read past the DCD header
  if (file_.OpenFile()) return 1;
  if (readDcdHeader()) return 1; 
  return 0;
}

// Traj_CharmmDcd::closeTraj()
/** Close the trajectory. If not reading, update the frame count
  * Header begins with header size, which could be a 4 bit or 8 bit number
  * depending on the size of an integer, followed by CORD. The number of 
  * frames is right after that so seek to size of int + 4.
  */
void Traj_CharmmDcd::closeTraj() {
  byte8 framecount;
  if (file_.IsOpen() && file_.Access() != CpptrajFile::READ) {
    file_.CloseFile();
    file_.OpenFile(CpptrajFile::UPDATE);
    file_.Seek( blockSize_+4 );
    framecount.i[1] = 0;
    framecount.i[0] = dcdframes_;
    if (debug_>0) mprintf("\tDEBUG: Updated DCD frame count is %i\n", dcdframes_);
    // NOTE: Here we are ensuring that ONLY 4 bytes are written. This could
    //       overflow for large # of frames.
    file_.Write(framecount.c,sizeof(unsigned char)*4); 
  }
  file_.CloseFile();
}

// ---------- Byte Swapping Routines -------------------------------------------
/*inline void endian_swap(unsigned int& x)
{
    x = (x>>24) | 
        ((x<<8) & 0x00FF0000) |
        ((x>>8) & 0x0000FF00) |
        (x<<24);
}*/
/*inline void endian_swap8(unsigned __int64& x)
{
    x = (x>>56) | 
        ((x<<40) & 0x00FF000000000000) |
        ((x<<24) & 0x0000FF0000000000) |
        ((x<<8)  & 0x000000FF00000000) |
        ((x>>8)  & 0x00000000FF000000) |
        ((x>>24) & 0x0000000000FF0000) |
        ((x>>40) & 0x000000000000FF00) |
        (x<<56);
}*/
// -----------------------------------------------------------------------------

// Traj_CharmmDcd::ReadBlock()
/** Read readByte number of bytes, convert to integer. If expected
  * is not -1, check that the integer matches expected.
  * Return the integer read on success, -1 on failure.
  */
int Traj_CharmmDcd::ReadBlock(int expected) {
  byte8 INbyte;
  // Read size of block
  INbyte.i[1] = 0;
  if (file_.Read(INbyte.c, sizeof(unsigned char)*blockSize_) < 1) {
    mprinterr("Error: Could not read block from DCD.\n");
    return -1;
  }
  // Swap endianness if necessary
  if (isBigEndian_) {
    if (is64bit_) 
      endian_swap( INbyte.i, 2 );
    else
      endian_swap( INbyte.i, 1 );
  }
  // Sum
  int val = INbyte.i[0] + INbyte.i[1];
  // If specified, check that this matches expected value
  if (expected != -1) {
    if ( val != expected ) {
      mprinterr("Error: Expected DCD block size of %i, got %i\n",expected,val);
      return -1;
    }
  }
  return val;
}

// Traj_CharmmDcd::WriteBlock()
/** Write given integer to charmm file. For now dont worry about the
  * endianness (use OS default).
  */
int Traj_CharmmDcd::WriteBlock(int blocksize) {
  byte8 OUTbyte;
  OUTbyte.i[1] = 0;
  OUTbyte.i[0] = blocksize;
  if (file_.Write(OUTbyte.c, sizeof(unsigned char)*blockSize_)) return 1;
  return 0;
}

/** Allocate space for coordinate reads/writes. */
void Traj_CharmmDcd::AllocateCoords() {
  coordinate_size_ = (size_t)dcdatom_ * sizeof(float);
  if (xcoord_ != 0) delete[] xcoord_;
  xcoord_ = new float[ 3 * dcdatom_ ];
  ycoord_ = xcoord_ + dcdatom_;
  zcoord_ = ycoord_ + dcdatom_;
}

// Traj_CharmmDcd::setupTrajin()
int Traj_CharmmDcd::setupTrajin(std::string const& fname, Topology* trajParm)
{
  if (file_.SetupRead( fname, debug_ )) return TRAJIN_ERR;
  // Call openTrajin, which will open and read past the DCD header.
  if ( openTrajin() ) return TRAJIN_ERR;
  // Warn if # atoms in parm does not match
  if (trajParm->Natom() != dcdatom_) {
    mprinterr("Error: # atoms in DCD file (%i) does not match # atoms in parm %s (%i)\n",
              dcdatom_, trajParm->c_str(), trajParm->Natom());
    return TRAJIN_ERR;
  }
  // Allocate space for coordinate reads.
  AllocateCoords();
  // DCD file may have less frames than is stored in the header.
  // Check the file size against the reported number of frames.
  size_t extraBytes;
  size_t file_size = (size_t)file_.UncompressedSize();
  if (file_size > 0) {
    size_t dimBytes = dcd_dim_ * sizeof(float);
    if (is64bit_)
      extraBytes = 4;
    else
      extraBytes = 2;
    frame1Bytes_ = (((size_t) dcdatom_ + extraBytes) * dimBytes) + boxBytes_;
    frameNBytes_ = (((size_t) nfreeat_ + extraBytes) * dimBytes) + boxBytes_;
    // Header size should be current position after open, which automatically
    // reads DCD header.
    headerBytes_ = (size_t)file_.Tell();
    file_size = file_size - headerBytes_ - frame1Bytes_;
    if ( (file_size % frameNBytes_) != 0 ) {
      mprintf("Warning: %s: Number of frames in DCD file could not be accurately determined.\n"
              "Warning:  File may be corrupted.\n", file_.Filename().base());
    }
    int nframes = (int)(file_size / frameNBytes_) + 1; // +1 for first frame
    if (nframes != dcdframes_) {
      mprintf("Warning: %s: Reported number of frames in DCD file is %i,\n",
              file_.Filename().base(), dcdframes_);
      mprintf("Warning:\tactual number of frames is %i. Only reading %i frames.\n",
              nframes, nframes);
      dcdframes_ = nframes;
    }
  } else {
    mprintf("Warning: Uncompressed size could not be determined. This is normal for\n");
    mprintf("Warning: bzip2 files. Cannot check # of frames. Will try to read %i\n",dcdframes_);
  }
  // Load box info so that it can be checked.
  double box[6];
  memset( box, 0, 6*sizeof(double));
  if (boxBytes_) {
    if (ReadBox( box )) return TRAJIN_ERR;
  }
  // Set traj info: No velocity, temperature, or time.
  SetCoordInfo( CoordinateInfo( Box(box), false, false, false ) );
  // If there are fixed atoms read the first frame now
  // TODO: Deal with fixed atoms
  closeTraj();
  return dcdframes_;
}

// Traj_CharmmDcd::readDcdHeader()
/** Read the header of a DCD file. Determine the endianness and bit size.
  * File must have already been opened. Return 1 on error, 0 on success.
  */
// TODO: Check if sizeof(int) is portable or if should be int32_t
int Traj_CharmmDcd::readDcdHeader() {
  // Endianness and 32/64 bit is determined in ID_TrajFormat. Read past
  // blockSize + 4 bytes (32/64 bit 84 + CORD).
  file_.Seek( blockSize_ + 4 );
  // ********** Step 1 - Read the rest of the first header block
  headerbyte buffer;
  if (file_.Read(buffer.c, sizeof(unsigned char)*80) < 1) {
    mprinterr("Error: Could not buffer DCD header.\n");
    return 1;
  }
  if (isBigEndian_) endian_swap(buffer.i, 20);
  // Print ICNTRL variables for debugging
  if (debug_>1) {
    for (int i=0; i < 20; i++) 
      mprintf("\ticntrl[%i]= %i\n",i,buffer.i[i]  );
  }
  // Make sure this is Charmm format; last integer in the header should not
  // be zero.
  if ( buffer.i[19] != 0 ) {
    if (debug_>0) mprintf("\tCharmm DCD\n");
    // Check for Charmm-specific flags
    //if ( buffer.i[10] != 0 ) dcdExtraBlock = true;
    if ( buffer.i[11] != 0 ) 
      dcd_dim_ = 4;
    else
      dcd_dim_ = 3;
  } else {
    mprinterr("\tNon-charmm DCD - currently unsupported.\n");
    return 1;
  }
  // Number of sets
  dcdframes_ = buffer.i[0];
  // Starting timestep
  //int istart = buffer.i[1];
  // Number of steps between frames
  //int nsavc  = buffer.i[2];
  // Number of fixed atoms
  nfixedat_  = buffer.i[8];
  // Box information
  if (buffer.i[10] != 0) 
    boxBytes_ = 48 + (2 * blockSize_); // 6(crds) * 8(double) + (hdr) + (end hdr)
  else
    boxBytes_ = 0;
  // Timestep
  float timestep = buffer.f[9];
  if (debug_>0) mprintf("\tTimestep is %f\n",timestep);
  // Read end size of first block, should also be 84
  if (ReadBlock(84)<0) return 1;
  // ********** Step 2 - Read title block
  // Read title block size
  int ntitle;
  std::string title;
  char dcdtitle[81];
  dcdtitle[80]='\0';
  int titleSize = ReadBlock(-1);
  if (titleSize < 0) return 1;
  // Read titles    
  if (debug_>1) mprintf("\tTitle block size %i\n",titleSize);
  if ( ((titleSize - 4) % 80) == 0 ) {
    // Read ntitle
    if (file_.Read(&ntitle,sizeof(int)) < 1) {
      mprintf("Error: DCD Reading ntitle.\n");
      return 1;
    }
    if (isBigEndian_) endian_swap(&ntitle,1);
    if (debug_>1) mprintf("\tNtitle %i\n",ntitle);
    for (int i=0; i < ntitle; i++) {
      file_.Read(dcdtitle,sizeof(char)*80);
      if (debug_>0) mprintf("\tTitle%i: [%s]\n",i+1,dcdtitle);
      title += dcdtitle;
    }
    SetTitle( title );
  }
  // Read title end block size
  if (ReadBlock(titleSize)<0) return 1;
  // ********** Step 3 - Read in natoms 
  // Read in next block size, should be 4
  if (ReadBlock(4)<0) return 1;
  // Read in number of atoms
  if (file_.Read(&dcdatom_,sizeof(int)) < 1) {
    mprintf("Error: DCD reading natom.\n");
    return 1;
  }
  if (isBigEndian_) endian_swap(&dcdatom_,1);
  if (debug_>0) mprintf("\tNatom %i\n",dcdatom_);
  // Read in end block size, should also be 4
  if (ReadBlock(4)<0) return 1;
  // ********** Step 4 - Read in free atom indices if necessary
  // Calculate number of free atoms (#total - #fixed)
  nfreeat_ = dcdatom_ - nfixedat_;
  // If number of fixed atoms not 0, need to read list of free atoms.
  if (nfixedat_ != 0) {
    mprintf("\tNumber of free atoms %i\n",nfreeat_);
    // Allocate space for nfreat atom indices
    if (freeat_ != 0) delete[] freeat_;
    freeat_ = new int[ nfreeat_ ];
    // Read index array size
    if (ReadBlock(nfreeat_ * 4) < 0) return 1;
    // Read index array
    if (file_.Read( freeat_, sizeof(int)*nfreeat_) < 1) {
      mprinterr("Error reading DCD free atom index array.\n");
      return 1;
    }
    if (isBigEndian_) endian_swap(freeat_, nfreeat_);
    // Read end index array size
    if (ReadBlock(nfreeat_ * 4) < 0) return 1;
  }
  return 0;
}

/** CHARMM and later versions of NAMD store the box angles as cos(angle), so
  * convert back to angle (degrees) if all of the angles are bounded between
  * -1 and 1. Special-case 90 degree angles for numerical stability of
  *  orthorhombic boxes
  */
static inline double BoxInDeg( double BoxInRad ) {
  if ( BoxInRad == 0 )
    return 90.0;
  return acos( BoxInRad ) * Constants::RADDEG;
}

// Traj_CharmmDcd::ReadBox()
int Traj_CharmmDcd::ReadBox(double* box) {
  double boxtmp[6];
  if ( ReadBlock(48) < 0) return 1;
  file_.Read(boxtmp, sizeof(double)*6);
  if (isBigEndian_) endian_swap8(boxtmp,6);
  if ( ReadBlock(-1) < 0) return 1;
  // Box lengths
  box[0] = boxtmp[0];
  box[1] = boxtmp[2];
  box[2] = boxtmp[5];
  if ( boxtmp[4] >= -1 && boxtmp[4] <= 1 &&
       boxtmp[3] >= -1 && boxtmp[3] <= 1 &&
       boxtmp[1] >= -1 && boxtmp[1] <= 1    ) 
  {
    // Angles are cos(angle)
    box[3] = BoxInDeg( boxtmp[4] );
    box[4] = BoxInDeg( boxtmp[3] );
    box[5] = BoxInDeg( boxtmp[1] );
  } else {
    // They are already in degrees.
    box[3] = boxtmp[4];
    box[4] = boxtmp[3];
    box[5] = boxtmp[1];
  }
  return 0;
}

// Traj_CharmmDcd::readFrame()
int Traj_CharmmDcd::readFrame(int set, Frame& frameIn) {
  if (set == 0)
    file_.Seek( headerBytes_ );
  else 
    file_.Seek( headerBytes_ + frame1Bytes_ + ((size_t)(set - 1) * frameNBytes_) );
  // Load box info
  if (boxBytes_ != 0) {
    if (ReadBox( frameIn.bAddress() )) return 1;
  }
  // NOTE: Only checking for file errors on X read attempt - is this OK?
  // Read X coordinates
  if (ReadBlock(-1) == -1) return 1;
  file_.Read(xcoord_, coordinate_size_);
  ReadBlock(-1);
  // Read Y coordinates
  ReadBlock(-1);
  file_.Read(ycoord_, coordinate_size_);
  ReadBlock(-1);
  // Read Z coordinates
  ReadBlock(-1);
  file_.Read(zcoord_, coordinate_size_);
  ReadBlock(-1);
  // Swap little->big endian if necessary
  if (isBigEndian_) 
    endian_swap(xcoord_, dcdatom_*3);
  // Put xyz values into coord array
  int x = 0;
  for (int n=0; n < dcdatom_; n++) {
    frameIn[x++] = (double) xcoord_[n];
    frameIn[x++] = (double) ycoord_[n];
    frameIn[x++] = (double) zcoord_[n];
  }
  return 0;
}

void Traj_CharmmDcd::WriteHelp() {
  mprintf("\tx64: Use 8 byte block size (default 4 bytes).\n");
}

// Traj_CharmmDcd::processWriteArgs()
int Traj_CharmmDcd::processWriteArgs(ArgList& argIn) {
  // Default is to write 32 bit
  is64bit_ = argIn.hasKey("x64");
  if (is64bit_)
    blockSize_ = 8;
  else
    blockSize_ = 4;
  // Determine 32 bit/64 bit block size 'if (sizeof(void*)!=4)'
  // TODO: Determine OS endianness
  isBigEndian_ = false;
  return 0;
}

// Traj_CharmmDcd::setupTrajout()
/** Set up the charmm dcd trajectory for writing. No writing is done here, the
  * actual write calls are in writeDcdHeader which is called from openTraj.
  * Set is64bit and isBigEndian, although they are not currently used during
  * writes; size and endianness will be OS default.
  */
// TODO: Check OS endianness!
int Traj_CharmmDcd::setupTrajout(std::string const& fname, Topology* trajParm,
                                 CoordinateInfo const& cInfoIn,
                                 int NframesToWrite, bool append)
{
  if (!append) {
    SetCoordInfo( cInfoIn );
    dcdatom_ = trajParm->Natom();
    // dcdframes = trajParm->parmFrames;
    dcdframes_ = 0;
    // Set up title
    if (Title().empty())
      SetTitle("Cpptraj generated dcd file.");
    // Allocate space for atom arrays
    AllocateCoords();
    // Calculate total dcd header size so we can seek past it on subsequent opens
    // (int + 4 + (int*20) + int) + (int + int + (ntitle*80) + int) + (int + int + int)
    // (22*int) + 4 + (3*int) + 80 + (3*int)
    // (28*int) + 84
    //dcdheadersize = (28*sizeof(int)) + 84;
    if (file_.SetupWrite( fname, debug_)) return 1;
    if (file_.OpenFile()) return 1;
    if (writeDcdHeader()) return 1;
  } else {
    if ( setupTrajin( fname, trajParm ) == TRAJIN_ERR ) return 1;
    mprintf("\tAppending to DCD file starting at frame %i\n", dcdframes_);
    // Re-open for appending
    if (file_.SetupAppend( fname, debug_ )) return 1;
    if (file_.OpenFile()) return 1;
  }
  return 0;
}

// Traj_CharmmDcd::writeDcdHeader()
/** Write the charmm dcd header. File should already be open.  All integers 
  * will be written at default OS size no matter what. For now alway write 
  * little-endian as well.
  */
int Traj_CharmmDcd::writeDcdHeader() {
  byte8 dcdkey;
  headerbyte buffer;
  // Write 84 - CORD + header size
  WriteBlock(84);
  // Write CORD header, 4 bytes only
  dcdkey.i[1] = 0;
  dcdkey.c[0]='C';
  dcdkey.c[1]='O';
  dcdkey.c[2]='R';
  dcdkey.c[3]='D';
  file_.Write(dcdkey.c, sizeof(unsigned char)*4);
  // Set up header information, 80 bytes
  memset(buffer.c, 0, 80*sizeof(unsigned char));
  // Frames
  //buffer.i[0] = trajParm->parmFrames;
  buffer.i[0] = 0;
  // Starting timestep
  buffer.i[1] = 1;
  // Number of steps between frames
  buffer.i[2] = 1;
  // Number of fixed atoms
  buffer.i[8] = 0;
  // Timestep
  buffer.f[9] = 0.001;
  // Charmm version - Should this just be set to 0?
  buffer.i[19] = 35;
  // Box information
  if (CoordInfo().HasBox()) {
    buffer.i[10] = 1;
    boxBytes_ = 48 + (2 * blockSize_);
  } else
    boxBytes_ = 0;
  // Write the header
  file_.Write(buffer.c, sizeof(unsigned char)*80);
  // Write endblock size
  WriteBlock(84);
  // Write title block - only 1 title for now
  // Title block size will be 4 + (NTITLE * 80)
  WriteBlock(84);
  // Write NTITLE
  dcdkey.i[0] = 1; 
  file_.Write(dcdkey.i, sizeof(int)*1);
  // If title is longer than 80 truncate it for now
  std::string title = Title();
  if (title.size() > 80) 
    mprintf("Warning: CharmmDCD: Title size is > 80 chars, truncating to 80.\n");
  title.resize(80);
  // Write title
  file_.Write((char*)title.c_str(), 80*sizeof(char));
  // Write title end block
  WriteBlock(84);
  // Write atom block - 4 bytes
  WriteBlock(4);
  dcdkey.i[0] = dcdatom_;
  file_.Write(dcdkey.i, sizeof(int));
  WriteBlock(4);

  return 0;
}

// Traj_CharmmDcd::writeFrame()
int Traj_CharmmDcd::writeFrame(int set, Frame const& frameOut) {
  // Box coords - 6 doubles, 48 bytes
  if (boxBytes_ != 0) {
    /* The format for the 'box' array used in cpptraj is not the same as the
     * one used for NAMD/CHARMM dcd files.  Refer to the reading routine above
     * for a description of the box info.
     */
    double boxtmp[6];
    boxtmp[0] = frameOut.BoxCrd().BoxX();
    boxtmp[2] = frameOut.BoxCrd().BoxY();
    boxtmp[5] = frameOut.BoxCrd().BoxZ();
    // The angles must be reported in cos(angle) format
    boxtmp[1] = cos(frameOut.BoxCrd().Gamma() * Constants::DEGRAD);
    boxtmp[3] = cos(frameOut.BoxCrd().Beta()  * Constants::DEGRAD);
    boxtmp[4] = cos(frameOut.BoxCrd().Alpha() * Constants::DEGRAD);
    WriteBlock(48);
    file_.Write(boxtmp, sizeof(double)*6);
    WriteBlock(48);
  }
  // Put X coords into xyz arrays
  int x = 0;
  for (int i = 0; i < dcdatom_; i++) {
    xcoord_[i] = (float)frameOut[x++];
    ycoord_[i] = (float)frameOut[x++];
    zcoord_[i] = (float)frameOut[x++];
  }
  // Write x coords
  WriteBlock(coordinate_size_);
  file_.Write(xcoord_, coordinate_size_);
  WriteBlock(coordinate_size_);
  // Write y coords
  WriteBlock(coordinate_size_);
  file_.Write(ycoord_, coordinate_size_);
  WriteBlock(coordinate_size_);
  // Write z coords 
  WriteBlock(coordinate_size_);
  file_.Write(zcoord_, coordinate_size_);
  WriteBlock(coordinate_size_);
  // Update frame count
  dcdframes_++;
  return 0;
}
 
// Traj_CharmmDcd::info()
void Traj_CharmmDcd::Info() {
  mprintf("is a CHARMM DCD file");
  if (isBigEndian_)
    mprintf(" Big Endian");
  else
    mprintf(" Little Endian");
  if (is64bit_) 
    mprintf(" 64 bit");
  else 
    mprintf(" 32 bit");
}
