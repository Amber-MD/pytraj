#include <cstdio>  // sscanf
#include <cstdlib> // atoi, atof
#include <cstring> // strncmp
#include "PDBfile.h"
#include "CpptrajStdio.h"

/// PDB record types
// NOTE: Must correspond with PDB_RECTYPE
const char* PDBfile::PDB_RECNAME[] = { "ATOM  ", "HETATM", "CRYST1", "TER   ",
  "END   ", "ANISOU", "EndRec", 0 };

// PDBfile::IsPDBkeyword()
bool PDBfile::IsPDBkeyword(std::string const& recname) {
  // Title Section
  if (recname.compare(0,6,"HEADER")==0) return true;
  if (recname.compare(0,6,"SOURCE")==0) return true;
  if (recname.compare(0,6,"AUTHOR")==0) return true;
  if (recname.compare(0,6,"OBSLTE")==0) return true;
  if (recname.compare(0,6,"KEYWDS")==0) return true;
  if (recname.compare(0,6,"REVDAT")==0) return true;
  if (recname.compare(0,6,"TITLE ")==0) return true;
  if (recname.compare(0,6,"EXPDTA")==0) return true;
  if (recname.compare(0,6,"SPRSDE")==0) return true;
  if (recname.compare(0,6,"SPLT  ")==0) return true;
  if (recname.compare(0,6,"NUMMDL")==0) return true;
  if (recname.compare(0,6,"JRNL  ")==0) return true;
  if (recname.compare(0,6,"CAVEAT")==0) return true;
  if (recname.compare(0,6,"MDLTYP")==0) return true;
  if (recname.compare(0,6,"REMARK")==0) return true;
  if (recname.compare(0,6,"COMPND")==0) return true;
  // Primary Structure Section
  if (recname.compare(0,5,"DBREF" )==0) return true; // DBREF[|1|2]
  if (recname.compare(0,6,"SEQADV")==0) return true;
  if (recname.compare(0,6,"MODRES")==0) return true;
  if (recname.compare(0,6,"SEQRES")==0) return true;
  // Coordinate Section
  if (recname.compare(0,6,"MODEL ")==0) return true;
  if (recname.compare(0,6,"ATOM  ")==0) return true;
  if (recname.compare(0,6,"ANISOU")==0) return true;
  if (recname.compare(0,3,"TER"   )==0) return true; // To recognize blank TER cards.
  if (recname.compare(0,6,"HETATM")==0) return true;
  // Crystallographic and Coordinate Transformation Section 
  if (recname.compare(0,6,"CRYST1")==0) return true;
  if (recname.compare(0,5,"SCALE" )==0) return true; // SCALEn
  if (recname.compare(0,5,"ORIGX" )==0) return true; // ORIGXn
  if (recname.compare(0,5,"MTRIX" )==0) return true; // MTRIXn
  return false;
}

// PDBfile::ID_PDB()
bool PDBfile::ID_PDB(CpptrajFile& fileIn) {
  // NOTE: ASSUME FILE SET UP FOR READ
  if (fileIn.OpenFile()) return false;
  std::string line1 = fileIn.GetLine();
  std::string line2 = fileIn.GetLine();
  fileIn.CloseFile();
  if (!IsPDBkeyword( line1 )) return false;
  if (!IsPDBkeyword( line2 )) return false;
  return true;
}

// PDBfile::NextRecord()
PDBfile::PDB_RECTYPE PDBfile::NextRecord() {
  if (NextLine() == 0) { // No more records to read
    recType_ = END_OF_FILE;
    return END_OF_FILE;
  }
  recType_ = UNKNOWN;
  // Try to ID the current record.
  if (strncmp(linebuffer_,"ATOM  ",6)==0 ||
      strncmp(linebuffer_,"HETATM",6)==0   )
    // Just return ATOM for ATOM/HETATM when reading.
    recType_ = ATOM;
  else if (strncmp(linebuffer_,"CRYST1",6)==0)
    recType_ = CRYST1;
  else if (linebuffer_[0]=='T' && linebuffer_[1]=='E' && linebuffer_[2]=='R')
    recType_ = TER;
  else if (linebuffer_[0]=='E' && linebuffer_[1]=='N' && linebuffer_[2]=='D')
    recType_ = END;
  return recType_;
}

Atom PDBfile::pdb_Atom() {
  // ATOM or HETATM keyword.
  // Check line length before any modification.
  size_t lineLength = strlen( linebuffer_ );
  // Atom number (6-11)
  // Atom name (12-16), Replace asterisks with single quotes.
  char savechar = linebuffer_[16];
  linebuffer_[16] = '\0';
  NameType aname(linebuffer_+12);
  aname.ReplaceAsterisk();
  linebuffer_[16] = savechar;
  // Chain ID (21)
  // Element (76-77), Protect against broken PDB files (lines too short).
  char eltString[2]; eltString[0] = ' '; eltString[1] = ' ';
  if (lineLength > 77) {
    eltString[0] = linebuffer_[76];
    eltString[1] = linebuffer_[77];
  } else if (!lineLengthWarning_) {
    mprintf("Warning: PDB line length is short (%zu chars, expected 80).\n", lineLength);
    lineLengthWarning_ = true;
  }
  // Set atom info.
  Atom outAtom(aname, linebuffer_[21], eltString);
  // NOTE: Additional values:
  //       10 chars between Bfactor and element: buffer[66] to buffer[75]
  //       Charge: buffer[78] and buffer[79]
  return outAtom;
}

void PDBfile::pdb_OccupanyAndBfactor(float& occ, float& bfac) {
  // Occupancy (54-59) | charge
  // B-factor (60-65) | radius
  sscanf(linebuffer_+54, "%f %f", &occ, &bfac);
}

// PDBfile::pdb_XYZ()
// NOTE: Should check for null Xout
void PDBfile::pdb_XYZ(double *Xout) {
  // X coord (30-38)
  char savechar = linebuffer_[38];
  linebuffer_[38] = '\0';
  Xout[0] = atof( linebuffer_+30 );
  linebuffer_[38] = savechar;
  // Y coord (38-46)
  savechar = linebuffer_[46];
  linebuffer_[46] = '\0';
  Xout[1] = atof( linebuffer_+38 );
  linebuffer_[46] = savechar;
  // Z coord (46-54)
  savechar = linebuffer_[54];
  linebuffer_[54] = '\0';
  Xout[2] = atof( linebuffer_+46 );
  linebuffer_[54] = savechar;
}

void PDBfile::pdb_Box(double* box) const {
  // CRYST1 keyword. RECORD A B C ALPHA BETA GAMMA SGROUP Z
  // A=6-15 B=15-24 C=24-33 alpha=33-40 beta=40-47 gamma=47-54
  sscanf(linebuffer_, "%*6s%9lf%9lf%9lf%7lf%7lf%7lf", box, box+1, box+2,
         box+3, box+4, box+5);
  mprintf("\tRead CRYST1 info from PDB: a=%g b=%g c=%g alpha=%g beta=%g gamma=%g\n",
          box[0], box[1], box[2], box[3], box[4], box[5]);
  // Warn if the box looks strange.
  if (box[0] == 1.0 && box[1] == 1.0 && box[2] == 1.0)
    mprintf("Warning: PDB cell lengths are all 1.0 Ang.;"
            " this usually indicates an invalid box.\n");
}

NameType PDBfile::pdb_ResName() {
  // Res name (16-20), Replace asterisks with single quotes.
  char savechar = linebuffer_[20];
  linebuffer_[20] = '\0';
  NameType resName(linebuffer_+16);
  linebuffer_[20] = savechar;
  resName.ReplaceAsterisk();
  return resName;
}

int PDBfile::pdb_ResNum() {
  // Res num (22-27)
  char savechar = linebuffer_[27];
  linebuffer_[27] = '\0';
  int resnum = atoi( linebuffer_+22 );
  linebuffer_[27] = savechar;
  return resnum;
}

// -----------------------------------------------------------------------------
// PDBfile::WriteRecordHeader()
void PDBfile::WriteRecordHeader(PDB_RECTYPE Record, int anum, NameType const& name,
                                NameType const& resnameIn, char chain, int resnum)
{
  char resName[5], atomName[5];

  resName[4]='\0';
  atomName[4]='\0';
  // Residue number in PDB format can only be 4 digits wide
  while (resnum>9999) resnum-=9999;
  // Atom number in PDB format can only be 5 digits wide
  while (anum>99999) anum-=99999;
  // Residue names in PDB format 3 chars long start at column 18. Residue 
  // names 2 chars long start at column 19. However in Amber residues can be 
  // 4 characters long padded with spaces at the end; adjust accordingly.
  resName[0] = resnameIn[0];
  resName[1] = resnameIn[1];
  if (resnameIn[2] == ' ')        // 2 chars
    resName[2] = '\0';
  else if (resnameIn[3] == ' ') { // 3 chars
    resName[2] = resnameIn[2];
    resName[3] = '\0';
  } else {                        // >= 4 chars
    resName[2] = resnameIn[2];
    resName[3] = resnameIn[3];
  }
  // Atom names in PDB format start from col 14 when <= 3 chars, 13 when 4 chars.
  if (name[3]!=' ') { // 4 chars
    atomName[0] = name[0];
    atomName[1] = name[1];
    atomName[2] = name[2];
    atomName[3] = name[3];
  } else {            // <= 3 chars
    atomName[0] = ' ';
    atomName[1] = name[0];
    atomName[2] = name[1];
    atomName[3] = name[2];
  }

  Printf("%-6s%5i %-4s%4s %c%4i",PDB_RECNAME[Record], anum, atomName,
               resName, chain, resnum);
}

// PDBfile::WriteTER()
void PDBfile::WriteTER(int anum, NameType const& resnameIn, char chain, int resnum)
{
  WriteRecordHeader(TER, anum, "", resnameIn, chain, resnum);
  Printf("\n"); 
}

// PDBfile::WriteHET()
void PDBfile::WriteHET(int res, double x, double y, double z) {
  WriteCoord(HETATM, anum_++, "XX", "XXX", ' ', 
             res, x, y, z, 0.0, 0.0, "", 0, false);
}

// PDBfile::WriteATOM()
void PDBfile::WriteATOM(int res, double x, double y, double z, 
                        const char* resnameIn, double Occ)
{
  WriteCoord(ATOM, anum_++, "XX", resnameIn, ' ',
             res, x, y, z, (float)Occ, 0.0, "", 0, false);
}

// PDBfile::WriteATOM()
void PDBfile::WriteATOM(const char* anameIn, int res, double x, double y, double z, 
                        const char* resnameIn, double Occ)
{
  WriteCoord(ATOM, anum_++, anameIn, resnameIn, ' ',
             res, x, y, z, (float)Occ, 0.0, "", 0, false);
}

// PDBfile::WriteCoord()
void PDBfile::WriteCoord(PDB_RECTYPE Record, int anum, NameType const& name,
                         NameType const& resnameIn, char chain, int resnum,
                         double X, double Y, double Z)
{
  WriteCoord(Record, anum, name, resnameIn, chain, 
             resnum, X, Y, Z, 0.0, 0.0, "", 0, false);
}

// PDBfile::WriteCoord()
void PDBfile::WriteCoord(PDB_RECTYPE Record, int anum, NameType const& name,
                         NameType const& resnameIn, char chain, int resnum,
                         double X, double Y, double Z, float Occ, float B, 
                         const char* Elt, int charge, bool highPrecision) 
{
  WriteRecordHeader(Record, anum, name, resnameIn, chain, resnum);
  if (highPrecision)
    Printf("    %8.3f%8.3f%8.3f%8.4f%8.4f      %2s%2s\n", X, Y, Z, Occ, B, Elt, "");
  else
    Printf("    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n", X, Y, Z, Occ, B, Elt, "");
}

// PDBfile::WriteANISOU()
void PDBfile::WriteANISOU(int anum, NameType const& name, 
                          NameType const& resnameIn, char chain, int resnum,
                          int u11, int u22, int u33, int u12, int u13, int u23,
                          const char* Elt, int charge)
{
  WriteRecordHeader(ANISOU, anum, name, resnameIn, chain, resnum);
  Printf("%c %7i%7i%7i%7i%7i%7i      %2s%2i\n", ' ', u11, u22, u33, 
         u12, u13, u23, Elt, charge);
}

// PDBfile::WriteTITLE()
void PDBfile::WriteTITLE(std::string const& titleIn) {
  std::string titleOut;
  titleOut.reserve(70);
  int lineNum = 1;
  for (std::string::const_iterator t = titleIn.begin(); t != titleIn.end(); ++t) {
    if (titleOut.empty()) { // Start of a new line
      if (lineNum > 1) // Need to write continuation
        Printf("TITLE   %2i ", lineNum);
      else
        Printf("TITLE      ");
    }
    titleOut.push_back( *t );
    if (titleOut.size() == 69) {
      // Flush line
      Printf("%-69s\n", titleOut.c_str());
      lineNum++;
      titleOut.clear();
    }
  }
  // Flush any unwritten chars.
  if (!titleOut.empty())
    Printf("%-69s\n", titleOut.c_str());
}

/** Expect x, y, z, alpha, beta, gamma */
void PDBfile::WriteCRYST1(const double* box, const char* space_group) {
  if (box==0) return;
  // RECROD A B C ALPHA BETA GAMMA SGROUP Z
  Printf("CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11s%4i\n",
         box[0], box[1], box[2], box[3], box[4], box[5], space_group, 1);
}

/** Write MODEL record: 1-6 MODEL, 11-14 model serial # */
void PDBfile::WriteMODEL(int model) {
  // Since num frames could be large, do not format the integer with width - OK?
  Printf("MODEL     %i\n", model);
}

void PDBfile::WriteENDMDL() { Printf("ENDMDL\n"); }

void PDBfile::WriteEND()    { Printf("END   \n"); }
