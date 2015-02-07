#include "Parm_SDF.h"
#include "SDFfile.h"
#include "CpptrajStdio.h"

// Parm_SDF::ID_ParmFormat() 
bool Parm_SDF::ID_ParmFormat(CpptrajFile& fileIn) {
  return SDFfile::ID_SDF(fileIn);
}
    
// Parm_SDF::ReadParm()
/** Read file as an SDF file. */
int Parm_SDF::ReadParm(std::string const& fname, Topology &parmOut) {
  SDFfile infile;

  if (infile.OpenRead(fname)) return 1;
  mprintf("    Reading SDF file %s as topology file.\n",infile.Filename().base());
  // Read header
  if (infile.ReadHeader()) return 1;
  parmOut.SetParmName( infile.SDF_Title(), infile.Filename() );
  // Read atoms
  double XYZ[3];
  for (int i = 0; i < infile.SDF_Natoms(); i++) {
    if ( infile.SDF_XYZ( XYZ ) ) {
      mprinterr("Error: Could not read atoms from SDF file.\n");
      return 1;
    }
    // Put everything in same residue
    parmOut.AddTopAtom( infile.SDF_Atom(), 0, NameType("LIG"), XYZ );
  }
  // Read bonds
  int at1, at2;
  for (int i = 0; i < infile.SDF_Nbonds(); i++) {
    if ( infile.SDF_Bond( at1, at2 ) ) {
      mprinterr("Error: Could not read bonds from SDF file.\n");
      return 1;
    }
    // SDF atom #s start from 1
    parmOut.AddBond( at1-1, at2-1 );
  }
  // No box
  parmOut.SetParmBox( Box() );

  mprintf("    SDF contains %i atoms, %zu bonds,\n", parmOut.Natom(), 
          parmOut.Bonds().size() + parmOut.BondsH().size());
  infile.CloseFile();

  return 0;
}

