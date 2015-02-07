#include <cstdio> // sscanf
#include <cstring> // strlen, strncmp
#include <locale> // isspace
#include <cstdlib> // atoi, atof
#include "Parm_Amber.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // RemoveTrailingWhitespace, SetIntegerFormatString etc
#include "Constants.h" // ELECTOAMBER, AMBERTOELEC

// ---------- Constants and Enumerated types -----------------------------------
const int Parm_Amber::AMBERPOINTERS=31;
/// Enumerated type for FLAG_POINTERS section
/** These variables are part of the POINTERS section of the topology.
  NATOM;    total number of atoms in the system
  NTYPES;   number of AMBER atom types used, max is 60
  NBONH;    number of bonds containing hydrogen
  NBONA;    number of bonds without hydrogen
  NTHETH;   number of angles containing hydrogen
  NTHETA;   number of angles not containing hydrogen
  NPHIH;    number of dihedrals containing hydrogen
  NPHIA;    number of dihedrals not containing hydrogen
  NHPARM;   NOT USED
  NPARM;    1 for LES parm
  NNB;      total number of excluded atoms
  NTOTRS;   total number of residues
  MBONA;    NBONA + number of constraint bonds
  MTHETA;   NTHETA + number of constraint angles
  MPHIA;    NPHIA + number of constraint dihedral angles
  NUMBND;   total number of unique bond types
  NUMANG;   total number of unique angle types
  NPTRA;    total number of unique dihedral types
  NATYP;    number of atom types defined in parameter file
  NPHB;     number of types of hydrogen bonded pair interactions
  IFPERT;   =1 if perturbation info is to be read =0 otherwise
  NBPER;    number of bonds to be perturbed
  NGPER;    number of angles to be perturbed
  NDPER;    number of dihedrals to be perturbed
  MBPER;    num of pert bonds across boundary to non-pert groups
  MGPER;    num of pert angles across boundary to non-pert groups
  MDPER;    num of pert dihedrals across bndry to non-pert groups
  IFBOX;    >0 if periodic box info to be read =0 otherwise
  NMXRS;    number of atoms in the largest residue
  IFCAP;    =1 if CAP option was used in edit, =0 otherwise
  NUMEXTRA; number of extra points (aka lone pairs)
*/
enum topValues {
//0         1       2      3       4       5       6       7      8       9
  NATOM=0,  NTYPES, NBONH, NBONA,  NTHETH, NTHETA, NPHIH,  NPHIA, NHPARM, NPARM,
  NNB,      NRES,   MBONA, MTHETA, MPHIA,  NUMBND, NUMANG, NPTRA, NATYP,  NPHB,
  IFPERT,   NBPER,  NGPER, NDPER,  MBPER,  MGPER,  MDPER,  IFBOX, NMXRS,  IFCAP,
  NEXTRA
};
// FORTRAN format strings
static const char* F10I8 = "%FORMAT(10I8)";
static const char* F20a4 = "%FORMAT(20a4)";
static const char* F5E16 = "%FORMAT(5E16.8)";
static const char* F3I8  = "%FORMAT(3I8)";
static const char* F1a80 = "%FORMAT(1a80)";
static const char* F1I8  = "%FORMAT(1I8)";
/// Constant strings for Amber parm flags and fortran formats.
const Parm_Amber::ParmFlag Parm_Amber::FLAGS[] = {
  { "POINTERS",                   F10I8 }, ///< Described above in topValues
  { "ATOM_NAME",                  F20a4 }, ///< Atom names
  { "CHARGE",                     F5E16 }, ///< Atom charges
  { "MASS",                       F5E16 }, ///< Atom masses
  { "RESIDUE_LABEL",              F20a4 }, ///< Residue names
  { "RESIDUE_POINTER",            F10I8 }, ///< Residue boundaries (atoms)
  { "AMBER_ATOM_TYPE",            F20a4 }, ///< Atom types
  { "BONDS_INC_HYDROGEN",         F10I8 }, ///< Bonds to hydrogen
  { "BONDS_WITHOUT_HYDROGEN",     F10I8 }, ///< Bonds not including hydrogen
  { "SOLVENT_POINTERS",           F3I8  },
  { "ATOMS_PER_MOLECULE",         F10I8 },
  { "BOX_DIMENSIONS",             F5E16 },
  { "ATOM_TYPE_INDEX",            F10I8 },
  { "NUMBER_EXCLUDED_ATOMS",      F10I8 },
  { "NONBONDED_PARM_INDEX",       F10I8 },
  { "LENNARD_JONES_ACOEF",        F5E16 },
  { "LENNARD_JONES_BCOEF",        F5E16 },
  { "EXCLUDED_ATOMS_LIST",        F10I8 },
  { "RADII",                      F5E16 },
  { "SCREEN",                     F5E16 },
  { "BOND_FORCE_CONSTANT",        F5E16 },
  { "BOND_EQUIL_VALUE",           F5E16 },
  { "ANGLE_FORCE_CONSTANT",       F5E16 },
  { "ANGLE_EQUIL_VALUE",          F5E16 },
  { "DIHEDRAL_FORCE_CONSTANT",    F5E16 },
  { "DIHEDRAL_PERIODICITY",       F5E16 },
  { "DIHEDRAL_PHASE",             F5E16 },
  { "SCEE_SCALE_FACTOR",          F5E16 },
  { "SCNB_SCALE_FACTOR",          F5E16 },
  { "SOLTY",                      F5E16 },
  { "ANGLES_INC_HYDROGEN",        F10I8 },
  { "ANGLES_WITHOUT_HYDROGEN",    F10I8 },
  { "DIHEDRALS_INC_HYDROGEN",     F10I8 },
  { "DIHEDRALS_WITHOUT_HYDROGEN", F10I8 },
  { "HBOND_ACOEF",                F5E16 },
  { "HBOND_BCOEF",                F5E16 },
  { "HBCUT",                      F5E16 },
  { "TREE_CHAIN_CLASSIFICATION",  F20a4 },
  { "JOIN_ARRAY",                 F10I8 },
  { "IROTAT",                     F10I8 },
  { "ATOMIC_NUMBER",              F10I8 },
  { "TITLE",                      F20a4 },
  { "RADIUS_SET",                 F1a80 },
  { "LES_NTYP",                   F10I8 }, // Number of LES region types
  { "LES_TYPE",                   F10I8 }, // LES type for each atom
  { "LES_FAC",                    F5E16 }, // Scaling factor for typeA * typeB  
  { "LES_CNUM",                   F10I8 }, // Copy number for each atom; 0==in all
  { "LES_ID",                     F10I8 }, // LES region ID
  { "CAP_INFO",                   F10I8 },
  { "CAP_INFO2",                  F5E16 },
  { "IPOL",                       F1I8  }, // 0 for fixed charge, 1 for polarizable
  { "POLARIZABILITY",             F5E16 }, // Hold atom polarazabilities in Ang^3
  // CHAMBER parameters
  { "CTITLE",                             "%FORMAT(a80)" },
  { "CHARMM_UREY_BRADLEY_COUNT",          "%FORMAT(2I8)" }, // # UB terms and types
  { "CHARMM_UREY_BRADLEY",                F10I8          }, // UB: 2 atoms and param index
  { "CHARMM_UREY_BRADLEY_FORCE_CONSTANT", F5E16          }, 
  { "CHARMM_UREY_BRADLEY_EQUIL_VALUE",    F5E16          },
  { "CHARMM_NUM_IMPROPERS",               F10I8          }, // # improper terms
  { "CHARMM_IMPROPERS",                   F10I8          }, // Imp: 4 atoms and param index
  { "CHARMM_NUM_IMPR_TYPES",              "%FORMAT(I8)"  }, // # improper types
  { "CHARMM_IMPROPER_FORCE_CONSTANT",     F5E16          }, 
  { "CHARMM_IMPROPER_PHASE",              F5E16          },
  { "LENNARD_JONES_14_ACOEF",             "%FORMAT(3E24.16)" },
  { "LENNARD_JONES_14_BCOEF",             "%FORMAT(3E24.16)" },
  { "CHARMM_CMAP_COUNT",                  "%FORMAT(2I8)" }, // # CMAP terms, # unique CMAP params
  { "CHARMM_CMAP_RESOLUTION",             "%FORMAT(20I4)"}, // # steps along each Phi/Psi CMAP axis
  { "CHARMM_CMAP_PARAMETER_",             "%FORMAT(8(F9.5))"}, // CMAP grid
  { "CHARMM_CMAP_INDEX",                  "%FORMAT(6I8)" }, // Atom i,j,k,l,m of cross term and idx
  { "FORCE_FIELD_TYPE",                   "%FORMAT(i2,a78)"},// NOTE: Cannot use with SetFortranType
  // PDB extra info
  { "RESIDUE_NUMBER", "%FORMAT(20I4)" }, // PDB residue number
  { "RESIDUE_CHAINID", F20a4 }, // PDB chain ID
  { "RESIDUE_ICODE", F20a4 }, // PDB residue insertion code
  { "ATOM_ALTLOC", F20a4 } // PDB atom alt location indicator FIXME: format is guess
};

// -----------------------------------------------------------------------------
// CONSTRUCTOR
Parm_Amber::Parm_Amber() :
  debug_(0),
  nochamber_(false), 
  ptype_(OLDPARM),
  ftype_(UNKNOWN_FTYPE),
  fncols_(0),
  fprecision_(0),
  fwidth_(0),
  error_count_(0),
  buffer_(0),
  buffer_size_(0),
  buffer_max_size_(0)
{ }

// DESTRUCTOR
Parm_Amber::~Parm_Amber() {
  if (buffer_!=0) delete[] buffer_;
}

// Parm_Amber::ID_ParmFormat()
bool Parm_Amber::ID_ParmFormat(CpptrajFile& fileIn) {
  int iamber[12];
  char lineBuf[BUF_SIZE];
  // Assumes already set up for READ
  if (fileIn.OpenFile()) return false;
  fileIn.Gets(lineBuf, BUF_SIZE);
  // Check for %VERSION
  if (strncmp(lineBuf,"%VERSION",8)==0) {
    fileIn.Gets(lineBuf, BUF_SIZE);
    // Check for %FLAG
    if (strncmp(lineBuf,"%FLAG",5)==0) {
      if (debug_>0) mprintf("  AMBER TOPOLOGY file\n");
      ptype_ = NEWPARM;
      fileIn.CloseFile();
      return true;
    }
  } else {
    // Since %VERSION not present, If first line is 81 bytes and the second 
    // line has 12 numbers in 12I6 format, assume old-style Amber topology
    // NOTE: Could also be less than 81? Only look for 12 numbers?
    int line1size = (int)strlen(lineBuf);
    if (line1size == (81 + fileIn.IsDos())) {
      fileIn.Gets(lineBuf, BUF_SIZE);
      if ( sscanf(lineBuf,"%6i%6i%6i%6i%6i%6i%6i%6i%6i%6i%6i%6i", 
                  iamber,   iamber+1, iamber+2,  iamber+3, 
                  iamber+4, iamber+5, iamber+6,  iamber+7, 
                  iamber+8, iamber+9, iamber+10, iamber+11) == 12 )
      {
        if (debug_>0) mprintf("  AMBER TOPOLOGY, OLD FORMAT\n");
        ptype_ = OLDPARM;
        fileIn.CloseFile();
        return true;
      }
    }
  }
  fileIn.CloseFile();
  return false;
}

// Parm_Amber::ReadParm()
int Parm_Amber::ReadParm(std::string const& fname, Topology &TopIn ) {
  if (file_.OpenRead( fname )) return 1;
  int err = ReadAmberParm( TopIn );
  file_.CloseFile();
  return err;
}
// ----- AMBER WRITE ROUTINES --------------------------------------------------
// Parm_Amber::AmberIfbox()
/** \return Amber IFBOX type:
  *   0: No box
  *   1: Orthogonal
  *   2: Truncated octahedral
  *   3: General triclinic
  */
int Parm_Amber::AmberIfbox(const Box& boxIn) {
  if      (boxIn.Type() == Box::NOBOX   ) return 0;
  else if (boxIn.Type() == Box::ORTHO   ) return 1; 
  else if (boxIn.Type() == Box::TRUNCOCT) return 2;
  return 3;
}

/** Determine if name needs to be truncated. */
void Parm_Amber::CheckNameWidth(const char* typeIn, NameType const& nameIn) {
  if (nameIn[4] != '\0')
    mprintf("Warning: Parm_Amber: %s name (%s) is too large and will be truncated (4 chars max).\n",
            typeIn, *nameIn);
}

/** Convert internal bond type to integer array. */
static std::vector<int> BondArrayToIndex(BondArray const& bondIn, bool chamber) {
  std::vector<int> arrayOut;
  arrayOut.reserve( bondIn.size() * 3 );
  int amult = 3;
  int aoff = 0;
  if (chamber) {
    amult = 1;
    aoff = 1;
  }
  for (BondArray::const_iterator bnd = bondIn.begin(); bnd != bondIn.end(); ++bnd) {
    arrayOut.push_back( (bnd->A1()*amult)+aoff );
    arrayOut.push_back( (bnd->A2()*amult)+aoff );
    arrayOut.push_back( bnd->Idx() + 1 );
  }
  return arrayOut;
}

/** Convert internal angle type to integer array. */
static std::vector<int> AngleArrayToIndex(AngleArray const& angleIn) {
  std::vector<int> arrayOut;
  arrayOut.reserve( angleIn.size() * 4 );
  for (AngleArray::const_iterator ang = angleIn.begin(); ang != angleIn.end(); ++ang) {
    arrayOut.push_back( ang->A1()*3 );
    arrayOut.push_back( ang->A2()*3 );
    arrayOut.push_back( ang->A3()*3 );
    arrayOut.push_back( ang->Idx() + 1 );
  }
  return arrayOut;
}

/** Convert internal dihedral type to integer array. End/Improper dihedrals
  * have 3rd and 4th atoms respectively as negative numbers.
  */
static std::vector<int> DihedralArrayToIndex(DihedralArray const& dihedralIn, bool chamber) {
  std::vector<int> arrayOut;
  arrayOut.reserve( dihedralIn.size() * 5 );
  int amult = 3;
  int aoff = 0;
  if (chamber) {
    amult = 1;
    aoff = 1;
  }
  for (DihedralArray::const_iterator dih = dihedralIn.begin(); dih != dihedralIn.end(); ++dih) {
    arrayOut.push_back( (dih->A1()*amult)+aoff );
    arrayOut.push_back( (dih->A2()*amult)+aoff );
    if ( dih->Type() == DihedralType::BOTH || dih->Type() == DihedralType::END)
      arrayOut.push_back( -((dih->A3()*amult)+aoff) );
    else
      arrayOut.push_back( (dih->A3()*amult)+aoff );
    if ( dih->Type() == DihedralType::BOTH || dih->Type() == DihedralType::IMPROPER)
      arrayOut.push_back( -((dih->A4()*amult)+aoff) );
    else
      arrayOut.push_back( (dih->A4()*amult)+aoff );
    arrayOut.push_back( dih->Idx() + 1 );
  }
  return arrayOut;
} 

void Parm_Amber::WriteHelp() {
  mprintf("\tnochamber: Do not write CHAMBER information to topology (useful for e.g. using"
          " topology for visualization with VMD).\n");
}

int Parm_Amber::processWriteArgs(ArgList& argIn) {
  nochamber_ = argIn.hasKey("nochamber");
  return 0;
}

// Parm_Amber::WriteParm()
/** CHAMBER writes out topologies in slightly different order than LEaP,
  * namely the EXCLUDED and SOLVENT_POINTERS sections are in different
  * places. The LEaP order is used here since the ReadAmberParm routine
  * is optimized for that order.
  */
int Parm_Amber::WriteParm(std::string const& fname, Topology const& parmIn) {
  ptype_ = NEWPARM;
  // Create arrays of atom info
  std::vector<NameType> names, types, itree;
  std::vector<double> charge, mass, gb_radii, gb_screen, polar;
  std::vector<int> at_num, atype_index, numex, excluded, join, irotat;
  names.reserve(       parmIn.Natom() );
  types.reserve(       parmIn.Natom() );
  charge.reserve(      parmIn.Natom() );
  polar.reserve(       parmIn.Natom() );
  mass.reserve(        parmIn.Natom() );
  gb_radii.reserve(    parmIn.Natom() );
  gb_screen.reserve(   parmIn.Natom() );
  at_num.reserve(      parmIn.Natom() );
  atype_index.reserve( parmIn.Natom() );
  // Will not write GB params if all are zero.
  bool hasGB = false;
  for (Topology::atom_iterator atom = parmIn.begin(); atom != parmIn.end(); ++atom) 
  {
    names.push_back( atom->Name() );
    CheckNameWidth("Atom",names.back());
    charge.push_back( atom->Charge() * Constants::ELECTOAMBER );
    polar.push_back( atom->Polar() );
    at_num.push_back( atom->AtomicNumber() );
    mass.push_back( atom->Mass() );
    // TypeIndex needs to be shifted +1 for fortran
    atype_index.push_back( atom->TypeIndex() + 1 );
    types.push_back( atom->Type() );
    CheckNameWidth("Type",types.back());
    gb_radii.push_back( atom->GBRadius() );
    gb_screen.push_back( atom->Screen() );
    if (gb_radii.back() > 0.0 || gb_screen.back() > 0.0)
      hasGB = true;
    // Amber atom exclusion list prints a 0 placeholder for atoms with
    // no exclusions, so always print at least 1 for numex
    int nex = atom->Nexcluded();
    if (nex == 0) {
      numex.push_back( 1 );
      excluded.push_back( 0 );
    } else {
      numex.push_back( nex );
      for (Atom::excluded_iterator ex = atom->excludedbegin();
                                   ex != atom->excludedend(); ex++)
        // Amber atom #s start from 1
        excluded.push_back( (*ex) + 1 );
    }
  }

  // Create arrays of residue info
  std::vector<int> resnums;
  std::vector<NameType> resnames;
  resnums.reserve(  parmIn.Nres() );
  resnames.reserve( parmIn.Nres() );
  for (Topology::res_iterator res = parmIn.ResStart(); res != parmIn.ResEnd(); ++res)
  {
    resnames.push_back( res->Name() );
    // Amber atom #s start from 1
    resnums.push_back( res->FirstAtom()+1 );
  }

  // Create pointer array
  std::vector<int> values(AMBERPOINTERS, 0);
  values[NATOM] = parmIn.Natom();
  values[NTYPES] = parmIn.Nonbond().Ntypes();
  values[NBONH] = (int)parmIn.BondsH().size();
  values[NBONA] = (int)parmIn.Bonds().size(); 
  values[NTHETH] = (int)parmIn.AnglesH().size();
  values[NTHETA] = (int)parmIn.Angles().size();
  values[NPHIH] = (int)parmIn.DihedralsH().size(); 
  values[NPHIA] = (int)parmIn.Dihedrals().size();
  // FIXME: Currently LES info not 100% correct, in particular the excluded list
  if (parmIn.LES().HasLES()) {
    mprintf("Warning: Excluded atom list for LES info is not correct.\n");
    values[NPARM] = 1;
  }
  values[NNB] = (int)excluded.size();
  values[NRES] = parmIn.Nres();
  //   NOTE: Assuming MBONA == NBONA etc
  values[MBONA] = values[NBONA];
  values[MTHETA] = values[NTHETA];
  values[MPHIA] = values[NPHIA];
  values[NUMBND] = (int)parmIn.BondParm().size();
  values[NUMANG] = (int)parmIn.AngleParm().size();
  values[NPTRA] = (int)parmIn.DihedralParm().size();
  values[NATYP] = parmIn.NatomTypes(); // Only for SOLTY
  values[NPHB] = (int)parmIn.Nonbond().HBarray().size();
  values[IFBOX] = AmberIfbox( parmIn.ParmBox() );
  values[NMXRS] = parmIn.FindResidueMaxNatom();
  if (parmIn.Cap().NatCap() > 0)
    values[IFCAP] = 1;
  values[NEXTRA] = parmIn.NextraPts();

  // Determine if this is a CHAMBER topology
  AmberParmFlagType titleFlag = F_TITLE;
  if (parmIn.Chamber().HasChamber()) {
    if (nochamber_) 
      mprintf("\tnochamber: Removing CHAMBER info from topology.\n");
    else {
      titleFlag = F_CTITLE;
      ptype_ = CHAMBER;
    }
  }
 
  // Write parm
  if (file_.OpenWrite( fname )) return 1;
  // HEADER AND TITLE (4 lines, version, flag, format, title)
  file_.Printf("%-44s%s                  \n",
               "%VERSION  VERSION_STAMP = V0001.000  DATE = ",
               TimeString().c_str());
  std::string title = parmIn.ParmName();
  // Resize title to max 80 char
  if (title.size() > 80)
    title.resize(80);
  file_.Printf("%%FLAG %-74s\n%-80s\n%-80s\n", FLAGS[titleFlag].Flag, 
               FLAGS[titleFlag].Fmt, title.c_str());
  // POINTERS
  WriteInteger(F_POINTERS, values);
  // CHAMBER only - Version and FF type
  if (ptype_ == CHAMBER)
    file_.Printf("%%FLAG %-74s\n%-80s\n%2i%-78s\n",FLAGS[F_FF_TYPE].Flag,
                 FLAGS[F_FF_TYPE].Fmt, parmIn.Chamber().FF_Version(),
                 parmIn.Chamber().FF_Type().c_str());
  // NAMES, CHARGE, ATOMIC NUMBER, MASS, TYPE INDEX, EXCLUDED, NB INDEX
  WriteName(F_NAMES, names);
  WriteDouble(F_CHARGE, charge);
  WriteInteger(F_ATOMICNUM, at_num);
  WriteDouble(F_MASS, mass);
  WriteInteger(F_ATYPEIDX, atype_index);
  WriteInteger(F_NUMEX, numex);
  // NONBONDED INDICES - positive needs to be shifted by +1 for fortran
  std::vector<int> nbindex;
  nbindex.reserve( parmIn.Nonbond().NBindex().size() );
  for (std::vector<int>::const_iterator it = parmIn.Nonbond().NBindex().begin();
                                        it != parmIn.Nonbond().NBindex().end(); ++it)
    if (*it > -1)
      nbindex.push_back( *it + 1 );
    else
      nbindex.push_back( *it );
  WriteInteger(F_NB_INDEX, nbindex);
  // RESIDUE NAME, POINTER
  WriteName(F_RESNAMES, resnames);
  WriteInteger(F_RESNUMS, resnums);
  // BOND, ANGLE, and DIHEDRAL FORCE CONSTANT and EQUIL VALUES
  std::vector<double> Rk, Req;
  // Bond params
  Rk.reserve( parmIn.BondParm().size() );
  Req.reserve( parmIn.BondParm().size() );
  for (BondParmArray::const_iterator parm = parmIn.BondParm().begin(); 
                                     parm != parmIn.BondParm().end(); ++parm)
  {
    Rk.push_back( parm->Rk() );
    Req.push_back( parm->Req() );
  }
  WriteDouble(F_BONDRK,  Rk);
  WriteDouble(F_BONDREQ, Req);
  // TODO: Use resize(0) instead?
  Rk.clear();
  Req.clear();
  // Angle params
  Rk.reserve( parmIn.AngleParm().size() );
  Req.reserve( parmIn.AngleParm().size() );
  for (AngleParmArray::const_iterator parm = parmIn.AngleParm().begin(); 
                                      parm != parmIn.AngleParm().end(); ++parm)
  {
    Rk.push_back( parm->Tk() );
    Req.push_back( parm->Teq() );
  }
  WriteDouble(F_ANGLETK, Rk);
  WriteDouble(F_ANGLETEQ, Req);
  Rk.clear();
  Req.clear();
  // CHARMM only - Urey-Bradley
  if (ptype_ == CHAMBER) {
    std::vector<int> UBC(2);
    UBC[0] = parmIn.Chamber().UB().size();
    UBC[1] = parmIn.Chamber().UBparm().size();
    Rk.reserve(  UBC[1] );
    Req.reserve( UBC[1] );
    for (BondParmArray::const_iterator parm = parmIn.Chamber().UBparm().begin();
                                       parm != parmIn.Chamber().UBparm().end(); ++parm)
    {
      Rk.push_back( parm->Rk() );
      Req.push_back( parm->Req() );
    }
    WriteInteger(F_CHM_UBC, UBC);
    WriteInteger(F_CHM_UB, BondArrayToIndex(parmIn.Chamber().UB(), true));
    WriteDouble(F_CHM_UBFC, Rk);
    WriteDouble(F_CHM_UBEQ, Req);
    Rk.clear();
    Req.clear();
  }
  // Dihedral params
  Rk.reserve( parmIn.DihedralParm().size() );
  Req.reserve( parmIn.DihedralParm().size() );
  std::vector<double> phase, scee, scnb;
  phase.reserve( parmIn.DihedralParm().size() );
  scee.reserve( parmIn.DihedralParm().size() ); 
  scnb.reserve( parmIn.DihedralParm().size() ); 
  for (DihedralParmArray::const_iterator parm = parmIn.DihedralParm().begin(); 
                                         parm != parmIn.DihedralParm().end(); ++parm)
  {
    Rk.push_back( parm->Pk() );
    Req.push_back( parm->Pn() );
    phase.push_back( parm->Phase() );
    scee.push_back( parm->SCEE() );
    scnb.push_back( parm->SCNB() );
  }
  WriteDouble(F_DIHPK, Rk);
  WriteDouble(F_DIHPN, Req);
  WriteDouble(F_DIHPHASE, phase);
  WriteDouble(F_SCEE, scee);
  WriteDouble(F_SCNB, scnb);
  Rk.clear();
  Req.clear();
  phase.clear();
  scee.clear();
  scnb.clear();
  // CHAMBER only - Impropers
  if (ptype_ == CHAMBER) {
    std::vector<int> NIMP(1, parmIn.Chamber().Impropers().size());
    WriteInteger(F_CHM_NIMP, NIMP);
    WriteInteger(F_CHM_IMP, DihedralArrayToIndex(parmIn.Chamber().Impropers(),true));
    NIMP[0] = parmIn.Chamber().ImproperParm().size();
    Rk.reserve( NIMP[0] );
    phase.reserve( NIMP[0] );
    for (DihedralParmArray::const_iterator parm = parmIn.Chamber().ImproperParm().begin();
                                           parm != parmIn.Chamber().ImproperParm().end(); ++parm)
    {
      Rk.push_back( parm->Pk() );
      phase.push_back( parm->Phase() );
    }
    WriteInteger(F_CHM_NIMPT, NIMP);
    WriteDouble(F_CHM_IMPFC, Rk);
    WriteDouble(F_CHM_IMPP, phase);
    Rk.clear();
    phase.clear();
  }
  // SOLTY - Currently unused
  WriteDouble(F_SOLTY, std::vector<double>(values[NATYP], 0.0));
  // LJ params
  Rk.reserve( parmIn.Nonbond().NBarray().size() );
  Req.reserve( parmIn.Nonbond().NBarray().size() );
  for (NonbondArray::const_iterator nb = parmIn.Nonbond().NBarray().begin();
                                    nb != parmIn.Nonbond().NBarray().end(); ++nb)
  {
    Rk.push_back( nb->A() );
    Req.push_back( nb->B() );
  }
  WriteDouble(F_LJ_A, Rk);
  WriteDouble(F_LJ_B, Req);
  Rk.clear();
  Req.clear();
  // CHAMBER only - LJ 1-4: Same size as LJ arrays above, no need to reserve.
  if (ptype_ == CHAMBER) {
    for (NonbondArray::const_iterator nb = parmIn.Chamber().LJ14().begin();
                                      nb != parmIn.Chamber().LJ14().end(); ++nb)
    {
      Rk.push_back( nb->A() );
      Req.push_back( nb->B() );
    }
    WriteDouble(F_LJ14A, Rk);
    WriteDouble(F_LJ14B, Req);
    Rk.clear();
    Req.clear();
  } 
  // BONDS/ANGLES/DIHEDRAL INDICES WITH AND WITHOUT HYDROGEN
  WriteInteger(F_BONDSH,  BondArrayToIndex(parmIn.BondsH(), false)); 
  WriteInteger(F_BONDS,   BondArrayToIndex(parmIn.Bonds(), false));
  WriteInteger(F_ANGLESH, AngleArrayToIndex(parmIn.AnglesH()));
  WriteInteger(F_ANGLES,  AngleArrayToIndex(parmIn.Angles()));
  WriteInteger(F_DIHH,    DihedralArrayToIndex(parmIn.DihedralsH(),false));
  WriteInteger(F_DIH,     DihedralArrayToIndex(parmIn.Dihedrals(),false));
  // EXCLUDED ATOMS LIST
  WriteInteger(F_EXCLUDE, excluded);
  // HBOND
  Rk.reserve( parmIn.Nonbond().HBarray().size() );
  Req.reserve( parmIn.Nonbond().HBarray().size() );
  phase.reserve( parmIn.Nonbond().HBarray().size() );
  for (HB_ParmArray::const_iterator hb = parmIn.Nonbond().HBarray().begin();
                                    hb != parmIn.Nonbond().HBarray().end(); ++hb)
  {
    Rk.push_back( hb->Asol() );
    Req.push_back( hb->Bsol() );
    phase.push_back( hb->HBcut() );
  }
  WriteDouble(F_ASOL, Rk);
  WriteDouble(F_BSOL, Req);
  WriteDouble(F_HBCUT, phase);
  Rk.clear();
  Req.clear();
  phase.clear();
  // AMBER ATOM TYPE
  WriteName(F_TYPES, types);
  // TREE CHAIN CLASSIFICATION, JOIN, IROTAT
  // TODO: Generate automatically
  if (!parmIn.Extra().empty()) {
    itree.reserve(  parmIn.Natom() );
    join.reserve(   parmIn.Natom() );
    irotat.reserve( parmIn.Natom() );
    for (std::vector<AtomExtra>::const_iterator ex = parmIn.Extra().begin();
                                                ex != parmIn.Extra().end(); ++ex)
    {
      itree.push_back(  ex->Itree()  );
      join.push_back(   ex->Join()   );
      irotat.push_back( ex->Irotat() );
    }
  }
  WriteName(   F_ITREE,  itree);
  WriteInteger(F_JOIN,   join);
  WriteInteger(F_IROTAT, irotat);
  // Write solvent info if IFBOX>0
  if (values[IFBOX] > 0) {
    // Determine first solvent molecule 
    int firstSolventMol = -1;
    for (Topology::mol_iterator mol = parmIn.MolStart(); mol != parmIn.MolEnd(); ++mol) {
      if ( mol->IsSolvent() ) { 
        firstSolventMol = (int)(mol - parmIn.MolStart()); 
        break;
      }
    }
    // Determine final solute residue based on first solvent molecule.
    int finalSoluteRes = 0;
    if (firstSolventMol == -1)
      finalSoluteRes = parmIn.Nres(); // No solvent Molecules
    else if (firstSolventMol > 0) {
      int finalSoluteAtom = parmIn.Mol(firstSolventMol).BeginAtom() - 1;
      finalSoluteRes = parmIn[finalSoluteAtom].ResNum() + 1;
    }
    // If no solvent, just set to 1 beyond # of molecules
    if (firstSolventMol == -1)
      firstSolventMol = parmIn.Nmol();
    // Solvent Pointers
    std::vector<int> solvent_pointer(3);
    solvent_pointer[0] = finalSoluteRes; // Already +1
    solvent_pointer[1] = parmIn.Nmol();
    solvent_pointer[2] = firstSolventMol + 1;
    WriteInteger(F_SOLVENT_POINTER, solvent_pointer);
    // ATOMS PER MOLECULE
    std::vector<int> APM;
    APM.reserve( solvent_pointer[1] );
    for (Topology::mol_iterator mol = parmIn.MolStart(); mol != parmIn.MolEnd(); mol++)
      APM.push_back( (*mol).NumAtoms() );
    WriteInteger(F_ATOMSPERMOL, APM);
    // BOX DIMENSIONS
    std::vector<double> betaLengths(4); // Beta X Y Z
    betaLengths[0] = parmIn.ParmBox().Beta();
    betaLengths[1] = parmIn.ParmBox().BoxX();
    betaLengths[2] = parmIn.ParmBox().BoxY();
    betaLengths[3] = parmIn.ParmBox().BoxZ();
    WriteDouble(F_PARMBOX, betaLengths);
  }
  // CAP info
  if (values[IFCAP] == 1) {
    std::vector<int> CI(1, parmIn.Cap().NatCap()+1);
    WriteInteger(F_CAP_INFO, CI);
    std::vector<double> CD(4);
    CD[0] = parmIn.Cap().CutCap();
    CD[1] = parmIn.Cap().xCap();
    CD[2] = parmIn.Cap().yCap();
    CD[3] = parmIn.Cap().zCap();
    WriteDouble(F_CAP_INFO2, CD);
  }
  // GB RADIUS SET, RADII, SCREENING PARAMETERS
  if (hasGB) {
    std::string radius_set = parmIn.GBradiiSet();
    if (!radius_set.empty()) {
      WriteSetup(F_RADSET, 1);
      if (radius_set.size()>80)
        radius_set.resize(80);
      file_.Printf("%-80s\n",radius_set.c_str());
    }
    WriteDouble(F_RADII, gb_radii);
    WriteDouble(F_SCREEN, gb_screen);
  }
  // CHAMBER only - Write out CMAP parameters
  if (ptype_ == CHAMBER && parmIn.Chamber().HasCmap()) {
    std::vector<int> CMAP(2);
    CMAP[0] = parmIn.Chamber().Cmap().size(); // CMAP terms
    CMAP[1] = parmIn.Chamber().CmapGrid().size(); // CMAP grids
    WriteInteger(F_CHM_CMAPC, CMAP);
    std::vector<int> CMAP_RES;
    CMAP_RES.reserve( CMAP[1] );
    for (CmapGridArray::const_iterator grid = parmIn.Chamber().CmapGrid().begin();
                                       grid != parmIn.Chamber().CmapGrid().end(); ++grid)
      CMAP_RES.push_back( grid->Resolution() );
    WriteInteger(F_CHM_CMAPR, CMAP_RES);
    int ngrid = 1;
    for (CmapGridArray::const_iterator grid = parmIn.Chamber().CmapGrid().begin();
                                       grid != parmIn.Chamber().CmapGrid().end();
                                       ++grid, ++ngrid)
    {
      // Assign format string
      fformat_.assign( FLAGS[F_CHM_CMAPP].Fmt );
      std::string fflag(FLAGS[F_CHM_CMAPP].Flag);
      fflag.append( integerToString(ngrid, 2) );
      // Set type, cols, width, and precision from format string
      // Write FLAG and FORMAT lines
      if (WriteFlagAndFormat(fflag.c_str(), grid->Grid().size())) return 1;
      WriteDoubleArray(grid->Grid());
    }
    CMAP_RES.clear();
    CMAP_RES.reserve( parmIn.Chamber().Cmap().size()*6 );
    for (CmapArray::const_iterator c = parmIn.Chamber().Cmap().begin();
                                   c != parmIn.Chamber().Cmap().end(); ++c)
    {
      CMAP_RES.push_back( c->A1()+1  );
      CMAP_RES.push_back( c->A2()+1  );
      CMAP_RES.push_back( c->A3()+1  );
      CMAP_RES.push_back( c->A4()+1  );
      CMAP_RES.push_back( c->A5()+1  );
      CMAP_RES.push_back( c->Idx()+1 );
    }
    WriteInteger(F_CHM_CMAPI, CMAP_RES);
  }
  // Polarizability - only write if it needs to be there
  if (parmIn.Ipol() > 0) {
    WriteInteger(F_IPOL, std::vector<int>(1, parmIn.Ipol()));
    WriteDouble(F_POLAR, polar);
  }
  // LES parameters - FIXME: Not completely correct yet
  if (values[NPARM] == 1) {
    std::vector<int> LES_array(1, parmIn.LES().Ntypes());
    WriteInteger(F_LES_NTYP, LES_array);
    LES_array.clear();
    // Sanity check.
    if ( (int)parmIn.LES().Array().size() != parmIn.Natom() ) {
      mprinterr("Internal Error: # LES atoms (%zu) != # atoms in topology (%i).\n",
                parmIn.LES().Array().size(), parmIn.Natom());
      return 1;
    }
    LES_array.reserve( parmIn.LES().Array().size() );
    std::vector<int> LES_cnum, LES_id;
    LES_cnum.reserve( LES_array.size() );
    LES_id.reserve( LES_id.size() );
    for (LES_Array::const_iterator les = parmIn.LES().Array().begin();
                                   les != parmIn.LES().Array().end(); ++les)
    {
      LES_array.push_back( les->Type() );
      LES_cnum.push_back( les->Copy() );
      LES_id.push_back( les->ID() );
    }
    WriteInteger(F_LES_TYPE, LES_array);
    WriteDouble(F_LES_FAC, parmIn.LES().FAC());
    WriteInteger(F_LES_CNUM, LES_cnum);
    WriteInteger(F_LES_ID, LES_id);
  }
 
  file_.CloseFile();

  return 0;
}

// ---- AMBER READ ROUTINES ----------------------------------------------------
// BondIndexToArray()
static BondArray BondIndexToArray(std::vector<int> const& bondIdx, bool chamber, int& err) {
  BondArray bonds;
  if ( (bondIdx.size() % 3) != 0 ) {
    mprinterr("Internal Error: Size of Amber bonds array not divisible by 3.\n");
    err++;
    return bonds;
  }
  int adiv = 3;
  int aoff = 0;
  if (chamber) {
    adiv = 1;
    aoff = 1;
  }
  bonds.reserve(bondIdx.size() / 3);
  for (std::vector<int>::const_iterator it = bondIdx.begin(); it != bondIdx.end(); it += 3)
    bonds.push_back( BondType( (*it / adiv) - aoff, (*(it+1) / adiv) - aoff, *(it+2) - 1 ) );
  return bonds;
}
// BondParmToArray()
static BondParmArray BondParmToArray(std::vector<double> const& bondk,
                                     std::vector<double> const& bondeq, int& err)
{
  BondParmArray bp;
  if (bondk.size() != bondeq.size()) {
    mprinterr("Error: Size of bond parm arrays inconsistent.\n");
    err++;
    return bp;
  }
  if (bondk.empty()) return bp;
  bp.reserve( bondk.size() );
  for (unsigned int i = 0; i < bondk.size(); i++)
    bp.push_back( BondParmType( bondk[i], bondeq[i] ) );
  return bp;
}
// AngleIndexToArray()
static AngleArray AngleIndexToArray(std::vector<int> const& angleIdx, int& err) {
  AngleArray angles;
  if ( (angleIdx.size() % 4) != 0 ) {
    mprinterr("Internal Error: Size of Amber angles array not divisible by 4.\n");
    err++;
    return angles;
  }
  angles.reserve(angleIdx.size() / 4);
  for (std::vector<int>::const_iterator it = angleIdx.begin(); it != angleIdx.end(); it += 4)
    angles.push_back( AngleType( *it / 3, *(it+1) / 3, *(it+2) / 3, *(it+3) - 1 ) );
  return angles;
}
// AngleParmToArray()
static AngleParmArray AngleParmToArray(std::vector<double> const& anglek,
                                       std::vector<double> const& angleeq, int& err)
{
  AngleParmArray ap;
  if (anglek.size() != angleeq.size()) {
    mprinterr("Error: Size of angle parm arrays inconsistent.\n");
    err++;
    return ap;
  }
  if (anglek.empty()) return ap;
  ap.reserve( anglek.size() );
  for (unsigned int i = 0; i < anglek.size(); i++)
    ap.push_back( AngleParmType( anglek[i], angleeq[i] ) );
  return ap;
}
// DihedralIndexToArray()
static DihedralArray DihedralIndexToArray(std::vector<int> const& dihedralIdx,
                                          bool chamber, int& err)
{
  DihedralArray dihedrals;
  if ( (dihedralIdx.size() % 5) != 0 ) {
    mprinterr("Internal Error: Size of Amber dihedrals array not divisible by 5.\n");
    err++;
    return dihedrals;
  }
  int adiv = 3;
  int aoff = 0;
  if (chamber) {
    adiv = 1;
    aoff = 1;
  }
  dihedrals.reserve(dihedralIdx.size() / 5);
  for (std::vector<int>::const_iterator it = dihedralIdx.begin(); it != dihedralIdx.end(); it += 5)
    dihedrals.push_back( DihedralType( (*it/adiv)-aoff, (*(it+1)/adiv)-aoff, (*(it+2)/adiv)-aoff,
                                       (*(it+3)/adiv)-aoff, *(it+4) - 1 ) );
  return dihedrals;
}
// NonbondParmToArray()
static NonbondArray NonbondParmToArray(std::vector<double> const& LJ_A,
                                       std::vector<double> const& LJ_B)
{
  NonbondArray NBA;
  NBA.reserve( LJ_A.size() );
  for (unsigned int i = 0; i < LJ_A.size(); i++)
    NBA.push_back( NonbondType(LJ_A[i], LJ_B[i]) );
  return NBA;
}

// Parm_Amber::ReadAmberParm()
int Parm_Amber::ReadAmberParm( Topology &TopIn ) {
  Box parmbox;
  std::string title;
  int Npointers = AMBERPOINTERS;
  ChamberParmType chamberParm; 

  if (ptype_ == NEWPARM) {
    if (debug_>0) 
      mprintf("\tReading Amber Topology file %s\n",file_.Filename().base());
    // Title. If not found check for CTITLE (chamber)
    if (PositionFileAtFlag(F_TITLE)) {
      title = GetLine();
    } else {
      if (PositionFileAtFlag(F_CTITLE)) {
        title = GetLine();
        ptype_ = CHAMBER;
      } else {
        // No TITLE or CTITLE, weird, but dont crash out yet.
        mprintf("Warning: '%s' No TITLE in Amber Parm.\n",file_.Filename().base());
      }
    }
  } else {
    mprintf("\tReading old (<v7) Amber Topology file.\n");
    title = GetLine();
    // One less POINTERS (no NEXTRA)
    Npointers = 30;
  }
  if (debug_>0) mprintf("\tAmberParm Title: \"%s\"\n",title.c_str());
  TopIn.SetParmName( title, file_.Filename() );
  // POINTERS
  std::vector<int> values = GetFlagInteger(F_POINTERS, Npointers);
  if (values.empty()) {
    mprinterr("Error: '%s' Could not get POINTERS from Amber Topology.\n",file_.Filename().base());
    return 1;
  }
  // Warn about IFPERT
  if (values[IFPERT] > 0)
    mprintf("Warning: '%s' contains perturbation information.\n"
            "Warning:  Cpptraj currently does not read of write perturbation information.\n",
            file_.Filename().base());
  // CHAMBER only - read force field type
  if (ptype_ == CHAMBER) {
    if (PositionFileAtFlag(F_FF_TYPE)) {
      file_.Gets(lineBuffer_, BUF_SIZE);
      int ff_verno;
      sscanf(lineBuffer_, "%i", &ff_verno);
      std::string fftype(lineBuffer_+2);
      RemoveTrailingWhitespace(fftype);
      chamberParm.SetChamber( ff_verno, fftype );
      mprintf("\tCHAMBER topology: %i: %s\n", ff_verno, fftype.c_str());
    } else {
      mprinterr("Error: CHAMBER topology missing FORCE_FIELD_TYPE\n");
      return 1;
    }
  }
  // Read parm variables
  std::vector<NameType> names  = GetFlagName(F_NAMES, values[NATOM]);
  std::vector<double>   charge = GetFlagDouble(F_CHARGE, values[NATOM]);
  // New parm only: ATOMICNUM
  std::vector<int> at_num;
  if (ptype_ != OLDPARM)
    at_num = GetFlagInteger(F_ATOMICNUM, values[NATOM]);
  std::vector<double> mass = GetFlagDouble(F_MASS, values[NATOM]);
  std::vector<int> atype_index = GetFlagInteger(F_ATYPEIDX, values[NATOM]);
  // For old parm need to read past NUMEX
  if (ptype_ == OLDPARM)
    GetFlagInteger(F_NUMEX, values[NATOM]);
  std::vector<int> NB_index = GetFlagInteger(F_NB_INDEX, values[NTYPES]*values[NTYPES]);
  std::vector<NameType> resnames = GetFlagName(F_RESNAMES, values[NRES]);
  std::vector<int> resnums = GetFlagInteger(F_RESNUMS,values[NRES]);
  // Bond and angle parameters
  BondParmArray BPA = BondParmToArray(GetFlagDouble(F_BONDRK, values[NUMBND]),
                                      GetFlagDouble(F_BONDREQ,values[NUMBND]), error_count_);
  AngleParmArray APA = AngleParmToArray(GetFlagDouble(F_ANGLETK, values[NUMANG]),
                                        GetFlagDouble(F_ANGLETEQ,values[NUMANG]), error_count_);
  // CHAMBER only - read Urey-Bradley info
  if (ptype_ == CHAMBER) {
    std::vector<int> UBC = GetFlagInteger(F_CHM_UBC, 2); // terms, types
    if (UBC.size() != 2) {
      mprinterr("Error: Could not get CHARMM_UREY_BRADLEY_COUNT in CHAMBER topology.\n");
      return 1;
    }
    BondArray UB = BondIndexToArray( GetFlagInteger(F_CHM_UB, UBC[0]*3), true, error_count_ );
    std::vector<double> UBFC = GetFlagDouble(F_CHM_UBFC, UBC[1]);
    std::vector<double> UBEQ = GetFlagDouble(F_CHM_UBEQ, UBC[1]);
    BondParmArray UBparm = BondParmToArray(GetFlagDouble(F_CHM_UBFC, UBC[1]),
                                           GetFlagDouble(F_CHM_UBEQ, UBC[1]), error_count_);
    chamberParm.SetUB( UB, UBparm );
  }
  // Dihedral parameters
  std::vector<double> dihedral_pk = GetFlagDouble(F_DIHPK, values[NPTRA]);
  std::vector<double> dihedral_pn = GetFlagDouble(F_DIHPN, values[NPTRA]);
  std::vector<double> dihedral_phase = GetFlagDouble(F_DIHPHASE, values[NPTRA]);
  // New parm only: SCEE and SCNB
  std::vector<double> scee_scale;
  std::vector<double> scnb_scale;
  if (ptype_ != OLDPARM) {
    scee_scale = GetFlagDouble(F_SCEE,values[NPTRA]);
    scnb_scale = GetFlagDouble(F_SCNB,values[NPTRA]);
  }
  // CHAMBER only - read Impropers
  if (ptype_ == CHAMBER) {
    std::vector<int> NIMP = GetFlagInteger(F_CHM_NIMP, 1);
    if (NIMP.empty()) {
      mprinterr("Error: Could not get CHARMM_NUM_IMPROPERS in CHAMBER topology.\n");
      return 1;
    }
    DihedralArray IMP = DihedralIndexToArray(GetFlagInteger(F_CHM_IMP,NIMP[0]*5),true,error_count_);
    NIMP = GetFlagInteger(F_CHM_NIMPT, 1);
    if (NIMP.empty()) {
      mprinterr("Error: Could not get CHARMM_NUM_IMPR_TYPES in CHAMBER topology.\n");
      return 1;
    }
    std::vector<double> imp_fc    = GetFlagDouble(F_CHM_IMPFC,NIMP[0]);
    std::vector<double> imp_phase = GetFlagDouble(F_CHM_IMPP, NIMP[0]);
    DihedralParmArray IMPparm;
    IMPparm.reserve( NIMP[0] );
    for (unsigned int i = 0; i < imp_fc.size(); i++)
      IMPparm.push_back( DihedralParmType(imp_fc[i], imp_phase[i]) );
    chamberParm.SetImproper( IMP, IMPparm );
  }
  // SOLTY: currently unused
  std::vector<double> solty = GetFlagDouble(F_SOLTY, values[NATYP]);
  // Lennard-Jones A and B parameters
  int nlj = values[NTYPES] * (values[NTYPES]+1) / 2;
  NonbondArray NBA = NonbondParmToArray(GetFlagDouble(F_LJ_A, nlj),
                                        GetFlagDouble(F_LJ_B, nlj));
  // CHAMBER only - read LJ 1-4 terms
  if (ptype_ == CHAMBER)
    chamberParm.SetLJ14(NonbondParmToArray(GetFlagDouble(F_LJ14A, nlj),
                                           GetFlagDouble(F_LJ14B, nlj)));
  // Bonds, angles, dihedrals
  BondArray bondsh = BondIndexToArray(GetFlagInteger(F_BONDSH,values[NBONH]*3),false,error_count_);
  BondArray bonds  = BondIndexToArray(GetFlagInteger(F_BONDS, values[MBONA]*3),false,error_count_);
  AngleArray anglesh = AngleIndexToArray(GetFlagInteger(F_ANGLESH,values[NTHETH]*4), error_count_);
  AngleArray angles  = AngleIndexToArray(GetFlagInteger(F_ANGLES, values[MTHETA]*4), error_count_);
  DihedralArray dihedralsh = DihedralIndexToArray(GetFlagInteger(F_DIHH,values[NPHIH]*5),
                                                  false, error_count_);
  DihedralArray dihedrals  = DihedralIndexToArray(GetFlagInteger(F_DIH, values[MPHIA]*5),
                                                  false, error_count_);
  // For old parm need to read past EXCLUDE
  if (ptype_ == OLDPARM)
    GetFlagInteger(F_EXCLUDE, values[NNB]);
  std::vector<double> asol = GetFlagDouble(F_ASOL,values[NPHB]);
  std::vector<double> bsol = GetFlagDouble(F_BSOL,values[NPHB]);
  std::vector<double> hbcut = GetFlagDouble(F_HBCUT,values[NPHB]);
  std::vector<NameType> types = GetFlagName(F_TYPES,values[NATOM]);
  std::vector<NameType> itree = GetFlagName(F_ITREE,values[NATOM]);
  std::vector<int> join_array = GetFlagInteger(F_JOIN,values[NATOM]);
  std::vector<int> irotat = GetFlagInteger(F_IROTAT,values[NATOM]);
  // Get solvent info if IFBOX>0
  if (values[IFBOX]>0) {
    std::vector<int> solvent_pointer = GetFlagInteger(F_SOLVENT_POINTER,3);
    if (solvent_pointer.empty()) {
      mprinterr("Error: Could not get %s from Amber Topology file.\n",
                FLAGS[F_SOLVENT_POINTER].Flag);
      return 1;
    }
    //finalSoluteRes = solvent_pointer[0] - 1;
    int molecules = solvent_pointer[1];
    //firstSolvMol = solvent_pointer[2] - 1;
    std::vector<int> atomsPerMol = GetFlagInteger(F_ATOMSPERMOL,molecules);
    // boxFromParm = {OLDBETA, BOX(1), BOX(2), BOX(3)}
    std::vector<double> boxFromParm = GetFlagDouble(F_PARMBOX,4);
    // If no box information present in the parm (such as with Chamber prmtops)
    // set the box info if ifbox = 2, otherwise default is NOBOX; the box info 
    // will eventually be set by angles from the first trajectory associated 
    // with this parm.
    if (boxFromParm.empty()) {
      if (ptype_ != CHAMBER) mprintf("Warning: Prmtop missing Box information.\n");
      // ifbox 2: truncated octahedron for certain
      if (values[IFBOX] == 2) 
        parmbox.SetTruncOct(); 
    // Determine box type, set Box angles and lengths from beta (boxFromParm[0])
    } else {
      parmbox.SetBetaLengths( boxFromParm[0], boxFromParm[1], boxFromParm[2], boxFromParm[3] );
    }
    // Check for IFBOX/BoxType mismatch
    if (values[IFBOX]==2 && parmbox.Type() != Box::TRUNCOCT) {
      mprintf("Warning: Amber Parm Box should be Truncated Octahedron (ifbox==2)\n");
      mprintf("         but BOX_DIMENSIONS indicate %s - may cause imaging problems.\n",
              parmbox.TypeName());
    }
  }
  // Water Cap Info
  if (values[IFCAP]) {
    std::vector<int> CI = GetFlagInteger(F_CAP_INFO, 1);
    std::vector<double> CD = GetFlagDouble(F_CAP_INFO2, 4);
    if (CI.size() != 1 || CD.size() != 4) {
      mprinterr("Error: Could not read Cap info.\n");
      return 1;
    }
    TopIn.SetCap( CapParmType(CI[0]-1, CD[0], CD[1], CD[2], CD[3]) );
  }
  // New Parm only: GB parameters; radius set, radii, and screening parameters
  std::vector<double> gb_radii;
  std::vector<double> gb_screen;
  if (ptype_ != OLDPARM) {
    if (PositionFileAtFlag(F_RADSET)) {
      std::string radius_set = GetLine();
      if (debug_ > 0) mprintf("\tRadius Set: %s\n",radius_set.c_str());
      TopIn.SetGBradiiSet( radius_set );
    }
    gb_radii = GetFlagDouble(F_RADII,values[NATOM]);
    gb_screen = GetFlagDouble(F_SCREEN,values[NATOM]); 
  }
  // CHAMBER only - read CMAP info
  if (ptype_ == CHAMBER) {
    std::vector<int> CMAP = GetFlagInteger(F_CHM_CMAPC, 2); // cmap terms, cmap grids
    if (!CMAP.empty()) {
      // Get CMAP resolutions for each grid
      std::vector<int> CMAP_RES = GetFlagInteger(F_CHM_CMAPR, CMAP[1]);
      // Read CMAP grids.
      for (int i = 0; i < CMAP[1]; i++) {
        // Find flag for CMAP grid, set format line
        std::string fflag(FLAGS[F_CHM_CMAPP].Flag);
        fflag.append( integerToString(i+1, 2) );
        if (!PositionFileAtFlag(fflag.c_str())) {
          mprinterr("Error: Could not find CMAP grid %s\n", fflag.c_str());
          return 1;
        }
        // Set up cols, width, etc from format
        if (!SetFortranType()) return 1;
        // Add Grid
        chamberParm.AddCmapGrid( CmapGridType(CMAP_RES[i], 
                                              GetDouble(fwidth_,fncols_,CMAP_RES[i]*CMAP_RES[i])));
      }
      // Read CMAP parameters
      CMAP_RES = GetFlagInteger(F_CHM_CMAPI, CMAP[0]*6);
      for (std::vector<int>::const_iterator it = CMAP_RES.begin(); it != CMAP_RES.end(); it+=6)
        chamberParm.AddCmapTerm( CmapType(*it-1,     *(it+1)-1, *(it+2)-1,
                                          *(it+3)-1, *(it+4)-1, *(it+5)-1) );
    }
  }
  // Polarizability - new topology only
  std::vector<double> polar;
  if (ptype_ != OLDPARM) {
    std::vector<int> IPOL = GetFlagInteger(F_IPOL, 1);
    if (!IPOL.empty()) TopIn.SetIpol( IPOL[0] );
    if (TopIn.Ipol() > 0)
      polar = GetFlagDouble(F_POLAR, values[NATOM]);
  }
  // LES parameters
  if (values[NPARM] == 1) {
    std::vector<int> LES_array = GetFlagInteger(F_LES_NTYP, 1);
    if (LES_array.empty()) return 1;
    int nlestyp = LES_array[0];
    LES_array = GetFlagInteger(F_LES_TYPE, values[NATOM]);
    std::vector<double> LES_fac = GetFlagDouble(F_LES_FAC, nlestyp * nlestyp);
    std::vector<int> LES_cnum = GetFlagInteger(F_LES_CNUM, values[NATOM]);
    std::vector<int> LES_id = GetFlagInteger(F_LES_ID, values[NATOM]);
    LES_ParmType les_parameters( values[NATOM], nlestyp, LES_fac );
    for (int i = 0; i < values[NATOM]; i++)
      les_parameters.AddLES_Atom( LES_AtomType(LES_array[i], LES_cnum[i], LES_id[i]) );
    TopIn.SetLES( les_parameters );
  }
  // Extra info from PDB files.
  std::vector<int> pdb_resnum;
  std::vector<NameType> pdb_res_chainID;
  std::vector<NameType> pdb_res_icode;
  std::vector<NameType> pdb_atom_alt;
  if (ptype_ != OLDPARM) { // FIXME: No Chamber?
    pdb_resnum = GetFlagInteger(F_PDB_RES, values[NRES]);
    pdb_res_chainID = GetFlagName(F_PDB_CHAIN, values[NRES]);
    pdb_res_icode = GetFlagName(F_PDB_ICODE, values[NRES]);
    pdb_atom_alt = GetFlagName(F_PDB_ALT, values[NATOM]);
  }
  // Done reading. Set up topology if no errors encountered. 
  if (error_count_==0) {
    // Add total # atoms to resnums.
    resnums.push_back( values[NATOM]+1 ); 
    int ri = 0;
    // If any optional arrays are empty, fill with zeros.
    if (at_num.empty()) at_num.resize(values[NATOM], 0);
    if (gb_radii.empty()) gb_radii.resize(values[NATOM], 0);
    if (gb_screen.empty()) gb_screen.resize(values[NATOM], 0);
    if (polar.empty()) polar.resize(values[NATOM], 0);
    if (pdb_resnum.empty())
      for (int rnum = 1; rnum <= values[NRES]; rnum++)
        pdb_resnum.push_back( rnum );
    if (pdb_res_chainID.empty()) pdb_res_chainID.resize(values[NRES]);
    // Add atoms to topology.
    for (int ai = 0; ai < values[NATOM]; ai++) {
      if (ai + 1 == resnums[ri+1]) ++ri;
      // Add atom. Convert Amber charge to elec and shift type index by -1.
      // NOTE: If ever used, shift atom #s in excludedAtoms by -1 so they start from 0 
      TopIn.AddTopAtom( Atom(names[ai], charge[ai] * Constants::AMBERTOELEC, polar[ai],
                             at_num[ai], mass[ai], atype_index[ai] - 1, types[ai],
                             gb_radii[ai], gb_screen[ai], pdb_res_chainID[ri][0]),
                        pdb_resnum[ri], resnames[ri], 0 );
    }
    // NOTE: Shift indices in parm arrays by -1
    error_count_ += TopIn.SetBondInfo(bonds, bondsh, BPA);
    error_count_ += TopIn.SetAngleInfo(angles, anglesh, APA);
    if (dihedral_pk.size() != dihedral_pn.size() || 
        dihedral_pk.size() != dihedral_phase.size())
    {
      mprinterr("Error: Size of dihedral parm arrays inconsistent.\n");
      return 1;
    }
    DihedralParmArray DPA;
    DPA.reserve( dihedral_pk.size() );
    // Set default SCEE and SCNB if not already set
    if (scee_scale.empty()) scee_scale.resize(dihedral_pk.size(), 1.2);
    if (scnb_scale.empty()) scnb_scale.resize(dihedral_pk.size(), 2.0);
    for (unsigned int i = 0; i < dihedral_pk.size(); i++)
      DPA.push_back( DihedralParmType(dihedral_pk[i], dihedral_pn[i], dihedral_phase[i],
                                      scee_scale[i], scnb_scale[i]) );
    error_count_ += TopIn.SetDihedralInfo(dihedrals, dihedralsh, DPA);
    HB_ParmArray HBA;
    HBA.reserve( asol.size() );
    for (unsigned int i = 0; i < asol.size(); i++)
      HBA.push_back( HB_ParmType(asol[i], bsol[i], hbcut[i]) );
    // Shift positive indices in NONBONDED index array by -1.
    for (std::vector<int>::iterator it = NB_index.begin(); it != NB_index.end(); ++it)
      if (*it > 0)
        *it -= 1;
    error_count_ += TopIn.SetNonbondInfo( NonbondParmType(values[NTYPES], NB_index, NBA, HBA) );
    // Set up Amber extra atom info.
    std::vector<AtomExtra> extra;
    if (itree.size() != join_array.size() || itree.size() != irotat.size())
      mprintf("Warning: Size of Amber arrays ITREE/JOIN/IROTATE do not match. Setting to blank.\n");
    else {
      if (!itree.empty()) {
        extra.reserve( itree.size() );
        for (size_t n = 0; n != itree.size(); n++)
          extra.push_back( AtomExtra(itree[n], join_array[n], irotat[n], ' ') );
      }
    }
    error_count_ += TopIn.SetExtraAtomInfo(values[NATYP], extra);
    if (values[IFBOX]>0) 
      TopIn.SetParmBox( parmbox );
    TopIn.SetChamber( chamberParm );
  }

  return error_count_;
}

// -----------------------------------------------------------------------------
// Parm_Amber::GetLine()
std::string Parm_Amber::GetLine() {
  file_.Gets(lineBuffer_, BUF_SIZE);
  std::string line( lineBuffer_ );
  // Remove any trailing whitespace
  RemoveTrailingWhitespace( line );
  return line;
}

// Parm_Amber::GetInteger()
std::vector<int> Parm_Amber::GetInteger(int width, int ncols, int maxval) {
  std::vector<int> iarray;
  // Read prmtop section into buffer
  int err = AllocateAndRead( width, ncols, maxval );
  if (err == 0)
    return iarray;
  else if (err == -1) {
    mprinterr("Error in read of integer values from %s\n",file_.Filename().base());
    ++error_count_;
    return iarray;
  }
  // Reserve variable memory
  iarray.reserve( maxval );
  // Convert values in buffer to integer 
  char *ptrbegin = buffer_;
  char *ptrend = buffer_;
  for (int i = 0; i < maxval; i++) {
    // Advance past newlines / CR (dos)
    while (*ptrbegin=='\n' || *ptrbegin=='\r') ++ptrbegin;
    ptrend = ptrbegin + width;
    char lastchar = *ptrend;
    *ptrend = '\0';
    iarray.push_back( atoi(ptrbegin) );
    *ptrend = lastchar;
    ptrbegin = ptrend;
  }
  return iarray;
}

// Parm_Amber::GetDouble()
std::vector<double> Parm_Amber::GetDouble(int width, int ncols, int maxval) {
  std::vector<double> darray;
  // Read prmtop section into buffer
  int err = AllocateAndRead( width, ncols, maxval );
  if (err == 0)
    return darray;
  else if (err == -1) {
    mprinterr("Error in read of double values from %s\n",file_.Filename().base());
    ++error_count_;
    return darray;
  }
  // Reserve variable memory
  darray.reserve( maxval );
  // Convert values in buffer to integer 
  char *ptrbegin = buffer_;
  char *ptrend = buffer_;
  for (int i = 0; i < maxval; i++) {
    // Advance past newlines / CR (dos)
    while (*ptrbegin=='\n' || *ptrbegin=='\r') ++ptrbegin;
    ptrend = ptrbegin + width;
    char lastchar = *ptrend;
    *ptrend = '\0';
    darray.push_back( atof(ptrbegin) );
    *ptrend = lastchar;
    ptrbegin = ptrend;
  }
  return darray;
}

// Parm_Amber::GetName()
std::vector<NameType> Parm_Amber::GetName(int width, int ncols, int maxval) {
  std::vector<NameType> carray;
  // Read prmtop section into buffer
  int err = AllocateAndRead( width, ncols, maxval );
  if (err == 0)
    return carray;
  else if (err == -1) {
    mprinterr("Error in read of Name values from %s\n",file_.Filename().base());
    ++error_count_;
    return carray;
  }
  // Reserve variable memory
  carray.reserve( maxval );
  // Convert values in buffer to integer 
  char *ptrbegin = buffer_;
  char *ptrend = buffer_;
  for (int i = 0; i < maxval; i++) {
    // Advance past newlines / CR (dos)
    while (*ptrbegin=='\n' || *ptrbegin=='\r') ++ptrbegin;
    ptrend = ptrbegin + width;
    char lastchar = *ptrend;
    *ptrend = '\0';
    carray.push_back( ptrbegin );
    *ptrend = lastchar;
    ptrbegin = ptrend;
  }
  return carray;
}

// Parm_Amber::GetFlagInteger()
std::vector<int> Parm_Amber::GetFlagInteger(AmberParmFlagType fflag, int maxval) {
  std::vector<int> iarray;
  if (ptype_ != OLDPARM) {
    // Seek to flag and set up fncols, fwidth
    if (!SeekToFlag(fflag)) return iarray;
    // NOTE: Check that type matches?
    iarray = GetInteger( fwidth_, fncols_, maxval );
  } else {
    iarray = GetInteger( 6, 12, maxval );
  }
  if ((int)iarray.size() != maxval)
    mprinterr("Error: Reading INTEGER section %s\n",FLAGS[fflag].Flag);
  return iarray;
}

// Parm_Amber::GetFlagDouble()
std::vector<double> Parm_Amber::GetFlagDouble(AmberParmFlagType fflag, int maxval) {
  std::vector<double> darray;
  if (ptype_ != OLDPARM) {
    // Seek to flag and set up fncols, fwidth
    if (!SeekToFlag(fflag)) return darray;
    // NOTE: Check that type matches?
    darray = GetDouble( fwidth_, fncols_, maxval );
  } else {
    darray = GetDouble( 16, 5, maxval );
  }
  if ((int)darray.size() != maxval)
    mprinterr("Error: Reading DOUBLE section %s\n",FLAGS[fflag].Flag); 
  return darray;
}

// Parm_Amber::GetFlagName()
std::vector<NameType> Parm_Amber::GetFlagName(AmberParmFlagType fflag, int maxval) {
  std::vector<NameType> carray;
  if (ptype_ != OLDPARM) {
    // Seek to flag and set up fncols, fwidth
    if (!SeekToFlag(fflag)) return carray;
    // NOTE: Check that type matches?
    carray = GetName( fwidth_, fncols_, maxval );
  } else {
    carray = GetName( 4, 20, maxval );
  }
  if ((int)carray.size() != maxval)
    mprinterr("Error: Reading CHAR section %s\n",FLAGS[fflag].Flag);
  return carray;
}

// Parm_Amber::SeekToFlag()
bool Parm_Amber::SeekToFlag(AmberParmFlagType fflag) {
  // Find flag, set format line
  if (!PositionFileAtFlag(fflag)) return false;
  // Set up cols, width, etc from format
  if (!SetFortranType()) return false;
  return true;
}  

// Parm_Amber::AllocateAndRead()
int Parm_Amber::AllocateAndRead(int width, int ncols, int maxval) {
  char temp[3]; // Only for when maxval is 0, space for \n, \r, null
  // Sanity Check
  if (maxval < 0) {
    mprinterr("Internal Error: AllocateAndRead called with maxval < 0 (%i)\n",maxval);
    return -1;
  }
  // If # expected values is 0 there will still be a newline placeholder
  // in the parmtop. Read past that and return
  if (maxval==0) {
    file_.Gets(temp,2);
    return 0;
  }
  // Allocate buffer to read in entire section
  size_t BufferSize = GetFortranBufferSize(width, ncols, maxval);
  if (buffer_!=0) delete[] buffer_;
  buffer_ = new char[ BufferSize ];
  // Read section from file
  int err = file_.Read(buffer_,BufferSize);
  return err;
}

bool Parm_Amber::PositionFileAtFlag(AmberParmFlagType flag) {
  return PositionFileAtFlag( FLAGS[flag].Flag );
}

// Parm_Amber::PositionFileAtFlag()
bool Parm_Amber::PositionFileAtFlag(const char* Key) {
  char value[BUF_SIZE];
  bool searchFile = true;
  bool hasLooped = false;

  if (debug_ > 0) mprintf("Reading %s\n",Key);
  // Search for '%FLAG <Key>'
  while ( searchFile ) {
    while ( file_.Gets(lineBuffer_, BUF_SIZE) == 0 ) {
      if ( strncmp(lineBuffer_,"%FLAG",5)==0 ) {
        // Check flag Key
        sscanf(lineBuffer_, "%*s %s",value);
        if (strcmp(value,Key)==0) {
          if (debug_>1) mprintf("DEBUG: Found Flag Key [%s]\n",value);
          // Read next line; can be either a COMMENT or FORMAT. If COMMENT, 
          // read past until you get to the FORMAT line
          file_.Gets(lineBuffer_, BUF_SIZE); 
          while (strncmp(lineBuffer_, "%FORMAT",7)!=0)
            file_.Gets(lineBuffer_, BUF_SIZE);
          if (debug_>1) mprintf("DEBUG: Format line [%s]\n",lineBuffer_);
          // Set format
          // NOTE: Check against stored formats?
          fformat_.assign(lineBuffer_);
          return true;
        } // END found Key
      } // END found FLAG line
    } // END scan through file
    // If we havent yet tried to search from the beginning, try it now.
    // Otherwise the Key has not been found.
    if (!hasLooped) {
      file_.Rewind();
      hasLooped = true;
    } else
      searchFile = false;
  }

  // If we reached here Key was not found.
  if (debug_>0)
    mprintf("Warning: '%s' Could not find Key %s in file.\n",file_.Filename().base(),Key);
  fformat_.clear();
  return false;
}

// -----------------------------------------------------------------------------
int Parm_Amber::WriteFlagAndFormat(const char* flag, size_t Nelements)
{
  //mprintf("DEBUG: FlagFormat[%s], Nelements=%zu\n",fformat_.c_str(),Nelements);
  // Set type, cols, width, and precision from format string
  if (!SetFortranType()) return 1;
  // Write FLAG and FORMAT lines
  file_.Printf("%%FLAG %-74s\n%-80s\n", flag, fformat_.c_str());
  // If Nelements is 0 just print newline and exit
  if (Nelements == 0) {
    file_.Printf("\n");
    return 1;
  }
  // Allocate character buffer space for memory, +1 for null
  size_t bufsize = GetFortranBufferSize(fwidth_, fncols_, Nelements);
  //mprinterr("DEBUG: Current max=%zu  Current size=%zu  Requested size=%zu\n",
  //          buffer_max_size_, buffer_size_, bufsize);
  if (bufsize+1 > buffer_max_size_) {
    if (buffer_!=0) delete[] buffer_;
    buffer_ = new char[ bufsize+1 ];
    buffer_max_size_ = bufsize+1;
  }
  buffer_size_ = bufsize;
  return 0;
}

// Parm_Amber::WriteSetup()
int Parm_Amber::WriteSetup(AmberParmFlagType fflag, size_t Nelements) {
  // Assign format string
  fformat_.assign( FLAGS[fflag].Fmt );
  // For chamber, certain flags have different format (boo).
  if (ptype_ == CHAMBER) {
    if      (fflag == F_CHARGE  ) fformat_.assign("%FORMAT(3E24.16)");
    else if (fflag == F_ANGLETEQ) fformat_.assign("%FORMAT(3E25.17)");
    else if (fflag == F_LJ_A    ) fformat_.assign("%FORMAT(3E24.16)");
    else if (fflag == F_LJ_B    ) fformat_.assign("%FORMAT(3E24.16)");
  }
  return WriteFlagAndFormat(FLAGS[fflag].Flag, Nelements);
}

// Parm_Amber::WriteInteger()
int Parm_Amber::WriteInteger(AmberParmFlagType fflag, std::vector<int>const& iarray)
{
  std::string FS;
  if (WriteSetup(fflag, iarray.size())) return 0;
  // Set up printf format string
  FS = SetIntegerFormatString(fwidth_);
  const char *FORMAT = FS.c_str();
  char *ptr = buffer_;
  int col = 0;
  for (std::vector<int>::const_iterator it = iarray.begin(); it != iarray.end(); it++) {
    sprintf(ptr, FORMAT, *it);
    // If # chars written > width, this silently truncates
    ptr += fwidth_;
    ++col;
    if (col == fncols_) {
      sprintf(ptr,"\n");
      ++ptr;
      col = 0;
    }
  }
  //mprinterr("INT: Last col written = %i (%i)\n",col,fncols_);
  if (col != 0) sprintf(ptr,"\n");
  file_.Write(buffer_, buffer_size_);

  return 0;
}

// Parm_Amber::WriteDouble()
int Parm_Amber::WriteDouble(AmberParmFlagType fflag, std::vector<double>const& darray) {
  if (WriteSetup(fflag, darray.size())) return 0;
  return WriteDoubleArray( darray );
}

// NOTE: Needs to be separate to write CHARMM_CMAP_PARAMETER_
int Parm_Amber::WriteDoubleArray(std::vector<double>const& darray)
{
  // Set up printf format string
  std::string FS;
  if (ftype_ == FDOUBLE)
    FS = SetDoubleFormatString(fwidth_, fprecision_, 2);
  else
    FS = SetDoubleFormatString(fwidth_, fprecision_, 1);
  const char *FORMAT = FS.c_str();
  char *ptr = buffer_;
  int col = 0;
  for (std::vector<double>::const_iterator it = darray.begin(); it != darray.end(); it++) {
    sprintf(ptr, FORMAT, *it);
    // If # chars written > width, this silently truncates
    ptr += fwidth_;
    ++col;
    if (col == fncols_) {
      sprintf(ptr,"\n");
      ++ptr;
      col = 0;
    }
  }
  //mprinterr("DOUBLE: Last col written = %i (%i)\n",col,fncols_);
  if (col != 0) sprintf(ptr,"\n");
  file_.Write(buffer_, buffer_size_);

  return 0;
}

// Parm_Amber::WriteName()
int Parm_Amber::WriteName(AmberParmFlagType fflag, std::vector<NameType>const& carray)
{
  std::string FS;
  if (WriteSetup(fflag, carray.size())) return 0;
  // Set up printf format string; true == no leading space
  FS = SetStringFormatString(fwidth_, true);
  const char *FORMAT = FS.c_str();
  char *ptr = buffer_;
  int col = 0;
  for (std::vector<NameType>::const_iterator it = carray.begin(); it != carray.end(); it++) {
    sprintf(ptr, FORMAT, *(*it));
    // If # chars written > width, this silently truncates
    ptr += fwidth_;
    ++col;
    if (col == fncols_) {
      sprintf(ptr,"\n");
      ++ptr;
      col = 0;
    }
  }
  //mprinterr("NAME: Last col written = %i (%i)\n",col,fncols_);
  if (col != 0) sprintf(ptr,"\n");
  file_.Write(buffer_, buffer_size_);

  return 0;
}

// -----------------------------------------------------------------------------
// Parm_Amber::GetFortranBufferSize()
/** Given number of columns and the width of each column, return the 
  * necessary char buffer size for N data elements.
  */
size_t Parm_Amber::GetFortranBufferSize(int width, int ncols, int N) {
  size_t bufferLines=0;
  size_t BufferSize=0;

  BufferSize = N * width;
  bufferLines = N / ncols;
  if ((N % ncols)!=0) ++bufferLines;
  // If DOS file there are CRs before Newlines
  if (file_.IsDos()) bufferLines *= 2;
  BufferSize += bufferLines;
  //if (debug_>0) 
  //  fprintf(stdout,"*** Buffer size is %i including %i newlines.\n",BufferSize,bufferLines);
  return BufferSize;
}

// Parm_Amber::SetFortranType()
/** Given a fortran-type format string, set the corresponding fortran
  * type. Set fncols (if present), fwidth, and fprecision (if present).
  */
// 01234567
// %FORMAT([<cols>][(]<type><width>[<precision>][)])
bool Parm_Amber::SetFortranType() {
  std::locale loc;
  std::string arg;

  if ( fformat_.empty() ) return false;
  // Make sure characters are upper case.
  for (std::string::iterator p = fformat_.begin(); p != fformat_.end(); p++)
    toupper(*p, loc);
  // Advance past left parentheses
  std::string::iterator ptr = fformat_.begin() + 7;
  while (*ptr=='(') ++ptr;
  // If digit, have number of data columns. Min # is 1
  fncols_ = 1;
  if (isdigit(*ptr, loc)) {
    while (ptr!=fformat_.end() && isdigit(*ptr, loc)) {
      arg += *ptr;
      ++ptr;
    }
    fncols_ = atoi( arg.c_str() );
  }
  // Advance past any more left parentheses
  while (ptr!=fformat_.end() && *ptr=='(') ++ptr;
  // Type
  if (ptr==fformat_.end()) {
    mprinterr("Error: Malformed fortran format string (%s) in Amber Topology %s\n",
              fformat_.c_str(), file_.Filename().base());
    return false;
  }
  switch (*ptr) {
    case 'I' : ftype_ = FINT;    break;
    case 'E' : ftype_ = FDOUBLE; break;
    case 'A' : ftype_ = FCHAR;   break;
    case 'F' : ftype_ = FFLOAT;  break;
    default  : ftype_ = UNKNOWN_FTYPE;
  }
  ++ptr;
  // Width
  fwidth_ = 0;
  arg.clear(); 
  while (isdigit(*ptr,loc)) {
    arg += *ptr;
    ++ptr;
  }
  fwidth_ = atoi( arg.c_str() );
  // Precision
  fprecision_ = 0;
  if (*ptr == '.') {
    ++ptr;
    arg.clear();
    while (isdigit(*ptr,loc)) {
      arg += *ptr;
      ++ptr;
    }
    fprecision_ = atoi( arg.c_str() );
  }
  if (debug_ > 2)
    mprintf("[%s]: cols=%i type=%i width=%i precision=%i\n",fformat_.c_str(),
            fncols_,(int)ftype_,fwidth_,fprecision_);

  return true;
}
