#include <cstdlib> // atoi
#include <cstdio>  // sscanf
#include <stdexcept> // std::runtime_error
#include "TinkerFile.h"
#include "ArgList.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h"

// CONSTRUCTOR
TinkerFile::TinkerFile() : natom_(0), hasBox_(false) {}

/// \return 1 if problem with or not a Tinker Atom/Title line.
static inline int SetNatomAndTitle(ArgList& lineIn, int& natom, std::string& title) {
  if (lineIn.Nargs() < 1) return 1;
  natom = lineIn.getNextInteger( -1 );
  if (natom < 1) return 1;
  std::string nextWord = lineIn.GetStringNext();
//if (nextWord.empty()) return 1;
  while (!nextWord.empty()) {
    if (!title.empty()) title += ' ';
    title.append( nextWord );
    nextWord = lineIn.GetStringNext();
  }
  return 0;
}

static inline bool IsAtomLine(ArgList& lineIn) {
  for (int i = 0; i < lineIn.Nargs(); i++) {
    std::string item = lineIn.GetStringNext();
    if (i == 0 || i >= 5) {
      try {
        convertToInteger( item );
      }
      catch (std::runtime_error e) {
        return false;
      }
    } else if (i >= 2 && i < 5) {
      try {
        convertToDouble( item );
      }
      catch (std::runtime_error e) {
        return false;
      }
    }
  }
  return true;
}

bool TinkerFile::ID_Tinker(CpptrajFile& fileIn) {
  // NOTE: ASSUME FILE SET UP FOR READ
  if (fileIn.OpenFile()) return false;
  ArgList firstLine( fileIn.NextLine() );
  ArgList secondLine( fileIn.NextLine() );
  ArgList thirdLine( fileIn.NextLine() );
  fileIn.CloseFile();
  // First line should have <natom> <title> only
  int natom = 0;
  std::string title;
  if ( SetNatomAndTitle(firstLine, natom, title) != 0 )
    return false;
  //mprinterr("Past SetNatomAndTitle\n");
  if (secondLine.Nargs() == 6) {
    bool isBoxLine = true;
    for (int i = 0; i < 6; i++) {
      // It is a box line if all 6 tokens are doubles
      try {
        convertToDouble( secondLine.GetStringNext() );
      } 
      catch (std::runtime_error e) {
        if (i != 1) return false;
        // We found a non-double on the second character -- it could be an atom
        // name. Check that the rest of the line matches an atom record
        isBoxLine = false;
        break;
      }
    }
    // If we are here it is not a box line, so make sure 
    if (!isBoxLine) {
      return IsAtomLine(secondLine);
    } else { // our second line WAS a box, now check the 3rd line
      return IsAtomLine(thirdLine);
    }
  }
  // There is no box, check that the second line is an atom line
  return IsAtomLine(secondLine);
}

/** Open tinker file. Read number of atoms and title from first frame. Set up
  * box coordinates if present.
  */
int TinkerFile::OpenTinker() {
  if (tinkerName_.empty()) {
    mprinterr("Internal Error: Tinker file name not set.\n");
    return 1;
  }
  if (file_.OpenFileRead( tinkerName_ )) return 1;
  ArgList line( file_.Line() );
  if ( SetNatomAndTitle(line, natom_, title_) ) {
    mprinterr("Error: Could not get # atoms / title from Tinker file.\n");
    return 1;
  }
  // Are box coords present? If not, next line read in should be 1 followed by
  // atom info and the line after that should be 2.
  hasBox_ = false;
  box_.SetNoBox();
  const char* secondptr = file_.Line();
  if (secondptr == 0) {
    mprinterr("Error: Could not get first atom line of Tinker file.\n");
    return 1;
  }
  // Third line could be zero if only 1 atom and no box.
  const char* thirdptr = file_.Line();
  if (natom_ == 1) {
    // If a third line was read, check if it is another title line. If so,
    // no box coordinates.
    if (thirdptr != 0) {
      line.SetList( std::string(thirdptr), " " );
      int natom2;
      std::string title2; // TODO: Check natom/title match?
      if (SetNatomAndTitle(line, natom2, title2))
        hasBox_ = true; // Not a title line, should be first atom so second should be box.
    } // else no third line read, no box.
  } else {
    if (thirdptr == 0) {
      mprinterr("Error: Could not get second atom line of Tinker file.\n");
      return 1;
    }
    // If the third line contains atom 1 there are box coords.
    file_.TokenizeLine(" ");
    int atomIdx = atoi( file_.NextToken() );
    if (atomIdx < 1) {
      mprinterr("Error: Third line contains invalid atom index.\n");
      mprinterr("Error: %s", thirdptr);
      return 1;
    }
    if (atomIdx == 1)
      hasBox_ = true;
  }
  // Set up box
  if (hasBox_) {
    double bp[6];
    if (secondptr == 0) return 1;
    if (sscanf(secondptr, "%lf %lf %lf %lf %lf %lf", bp, bp+1, bp+2, bp+3, bp+4, bp+5)!=6)
    {
      mprinterr("Error: Expected 6 box coordinates.\n");
      return 1;
    }
    box_.SetBox( bp );
  }
  // Close and reopen the file.
  file_.CloseFile();
  return file_.OpenFileRead( tinkerName_ );
}

/// \return 1 if number of atoms does not match what file was set up for.
int TinkerFile::CheckTitleLine() {
  file_.TokenizeLine(" ");
  int lineNatom = atoi( file_.NextToken() );
  if (lineNatom != natom_) {
    mprinterr("Error: Number of atoms in Tinker file changes from %i to %i\n",
              "Error: at line %i\n", natom_, lineNatom, file_.LineNumber());
    return 1;
  }
  return 0;
}

/** \return 0 if no more frames to read.
  * \return -1 if an error occurs.
  * \return 1 if more frames to read.
  */
int TinkerFile::NextTinkerFrame() {
  // Title line
  if (file_.Line() == 0) return 0;
  if (CheckTitleLine()) return -1;
  // Box line if necessary
  if (hasBox_) {
    if (file_.Line() == 0) {
      mprinterr("Error: Could not read Tinker box line (%i).\n", file_.LineNumber());
      return -1;
    }
  }
  for (int atidx = 0; atidx < natom_; atidx++)
    if (file_.Line() == 0) {
      mprinterr("Error: Could not read Tinker atom line (%i).\n", file_.LineNumber());
      return -1;
    }
  return 1;
}

/** \return 0 if no more frames to read.
  * \return -1 if an error occurs.
  * \return 1 if more frames to read.
  */
int TinkerFile::ReadNextTinkerFrame(double* Xptr, double* box) {
  // Title line
  if (file_.Line() == 0) return 0;
  if (CheckTitleLine()) return -1;
  // Box line
  if (hasBox_) {
    if (file_.Line() == 0) {
      mprinterr("Error: Could not read Tinker box line (%i).\n", file_.LineNumber());
      return -1;
    }
    int nbox = file_.TokenizeLine(" ");
    if (nbox != 6) {
      mprinterr("Error: In Tinker file line %i expected 6 box coords, got %i\n", 
                file_.LineNumber(), nbox);
      return -1;
    }
    for (int b = 0; b != nbox; b++)
      box[b] = atof( file_.NextToken() );
  }
  // Coords
  for (int atidx = 0; atidx < natom_; atidx++) {
    if (file_.Line() == 0) {
      mprinterr("Error: Could not read Tinker atom line (%i).\n", file_.LineNumber());
      return -1;
    }
    int ncol = file_.TokenizeLine(" ");
    if (ncol < 5) {
      mprinterr("Error: In Tinker file line %i expected at least 5 columns for atom, got %i\n",
                file_.LineNumber(), ncol);
      return -1;
    }
    file_.NextToken(); // Atom index
    file_.NextToken(); // Atom name
    Xptr[0] = atof( file_.NextToken() ); // X
    Xptr[1] = atof( file_.NextToken() ); // Y
    Xptr[2] = atof( file_.NextToken() ); // Z
    Xptr += 3;
  }
  return 1;
}

/* \return an array of Atoms from the current frame.
 * Assumes file has already been opened.
 */
std::vector<Atom> TinkerFile::ReadTinkerAtoms(double* XYZ, std::vector<int>& bonds)
{
  std::vector<Atom> atoms;
  if (XYZ == 0) {
    mprinterr("Internal Error: No space allocated for reading Tinker atom coordinates.\n");
    return atoms;
  }
  // Title line
  if (file_.Line() == 0) return atoms;
  if (CheckTitleLine()) return atoms;
  // Box line
  if (hasBox_) {
    if (file_.Line() == 0) return atoms;
  }
  // Read atoms
  atoms.reserve( natom_ );
  for (int atidx = 0; atidx < natom_; atidx++) {
    if (file_.Line() == 0) return std::vector<Atom>(0);
    int ncol = file_.TokenizeLine(" ");
    if (ncol < 6) {
      mprinterr("Error: In Tinker file line %i expected at least 5 columns for atom, got %i\n",
                file_.LineNumber(), ncol);
      return std::vector<Atom>(0);
    }
    file_.NextToken(); // Atom index
    NameType atom_name( file_.NextToken() );
    XYZ[0] = atof( file_.NextToken() ); // X
    XYZ[1] = atof( file_.NextToken() ); // Y
    XYZ[2] = atof( file_.NextToken() ); // Z
    XYZ += 3;
    const char* at_type_ptr = file_.NextToken(); // Atom Type Index
    int atom_type_index = atoi( at_type_ptr );
    NameType atom_type( at_type_ptr );
    // Read in any bonded partners.
    for (int col = 6; col != ncol; col++) {
      int bonded_atom = atoi(file_.NextToken()) - 1; // Tinker atoms start from 1
      if (atidx < bonded_atom) {
        bonds.push_back( atidx );
        bonds.push_back( bonded_atom );
      }
    }
    atoms.push_back( Atom(atom_name, atom_type, atom_type_index) );
  }
  return atoms;
}
