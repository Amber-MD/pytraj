// Parm_Mol2.cpp
#include "Parm_Mol2.h"
#include "Mol2File.h"
#include "CpptrajStdio.h"

// Parm_Mol2::ID_ParmFormat() 
bool Parm_Mol2::ID_ParmFormat(CpptrajFile& fileIn) {
  return Mol2File::ID_Mol2(fileIn);
}
    
// Parm_Mol2::ReadParm()
/** Read file as a Tripos Mol2 file. */
int Parm_Mol2::ReadParm(std::string const& fname, Topology &parmOut) {
  Mol2File infile;
  int current_res = 0;
  if (infile.OpenRead(fname)) return 1;
  mprintf("    Reading Mol2 file %s as topology file.\n",infile.Filename().base());
  // Get @<TRIPOS>MOLECULE information
  if (infile.ReadMolecule()) return 1;
  parmOut.SetParmName( infile.Mol2Title(), infile.Filename() );

  // Get @<TRIPOS>ATOM information
  if (infile.ScanTo( Mol2File::ATOM)) return 1;
  double XYZ[3];
  for (int atom=0; atom < infile.Mol2Natoms(); atom++) {
    if ( infile.Mol2XYZ(XYZ) ) return 1;
    NameType mol2resname = infile.Mol2Residue(current_res);
    parmOut.AddTopAtom( infile.Mol2Atom(), current_res, mol2resname, XYZ );
  }

  // Get @<TRIPOS>BOND information [optional]
  int at1 = 0;
  int at2 = 0;
  if (infile.ScanTo(Mol2File::BOND)==0) {
    for (int bond=0; bond < infile.Mol2Nbonds(); bond++) {
      if (infile.Mol2Bond(at1, at2)) return 1;
      // mol2 atom #s start from 1
      parmOut.AddBond(at1-1, at2-1);
    }
  } else {
    mprintf("      Mol2 file does not contain bond information.\n");
  }

  // No box
  parmOut.SetParmBox( Box() );

  mprintf("    Mol2 contains %i atoms, %i residues,\n", parmOut.Natom(),parmOut.Nres());
  //mprintf("    %i bonds to H, %i other bonds.\n", parmOut.NbondsWithH,parmOut.NbondsWithoutH);

  infile.CloseFile();

  return 0;
}

