// Traj_PDBfile
#include "Traj_PDBfile.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Traj_PDBfile::Traj_PDBfile() :
  pdbAtom_(0),
  currentSet_(0),
  ter_num_(0),
  pdbWriteMode_(NONE),
  dumpq_(false),
  dumpr_(false),
  pdbres_(false),
  pdbatom_(false),
  write_cryst1_(false),
  pdbTop_(0),
  chainchar_(' ')
{}

//------------------------------------------------------------------------
bool Traj_PDBfile::ID_TrajFormat(CpptrajFile& fileIn) {
  return PDBfile::ID_PDB( fileIn );
}

// Traj_PDBfile::closeTraj()
/** If not writing one PDB per frame write the END record. */
void Traj_PDBfile::closeTraj() {
  if ( (pdbWriteMode_ == SINGLE || pdbWriteMode_ == MODEL) &&
        file_.IsOpen() )
    file_.WriteEND();
  if (pdbWriteMode_ != MULTI)
    file_.CloseFile();
}

// Traj_PDBfile::openTrajin()
int Traj_PDBfile::openTrajin() {
  currentSet_ = 0;
  return file_.OpenFile();
}

// Traj_PDBfile::setupTrajin()
/** Scan PDB file to determine number of frames (models). The first frame will 
  * also be checked to ensure that the atom names match those in the parm file
  * in TrajectoryFile.
  */
int Traj_PDBfile::setupTrajin(std::string const& fname, Topology* trajParm)
{
  int atom;
  pdbWriteMode_ = NONE;
  if (file_.SetupRead( fname, debug_ )) return TRAJIN_ERR;
  if (file_.OpenFile()) return TRAJIN_ERR;
  // Two strats - check for MODEL keywords or see how many times natom ATOM
  // records can be read. Currently employing the latter.
  int Frames = 0;
  int numMismatch = 0;
  bool scanPDB = true;
  Box boxInfo;
  while (scanPDB) {
    atom = 0;
    while ( atom < trajParm->Natom() ) {
      //fprintf(stdout,"DEBUG: PDB Read atom %i\n",atom);
      if ( file_.NextRecord() == PDBfile::END_OF_FILE ) {
        scanPDB = false;
        break;
      }
      //fprintf(stdout,"DEBUG: PDB buffer %i [%s]\n",atom,buffer);
      if (file_.RecType() ==  PDBfile::CRYST1) {
        // Read in box information
        double box_crd[6];
        file_.pdb_Box( box_crd );
        boxInfo.SetBox( box_crd );
      } 
      // Skip non-ATOM records
      if (file_.RecType() != PDBfile::ATOM) continue;
      // If still on first frame, check pdb atom name against the name in the 
      // associated parm file.
      if (Frames==0) {
        Atom pdbAtom = file_.pdb_Atom();
        if ( pdbAtom.Name() != (*trajParm)[atom].Name() ) {
          if (debug_>1) 
            mprintf("Warning: %s: PDB atom %i name [%s] does not match parm atom name [%s]\n",
                    file_.Filename().base(), atom+1, *(pdbAtom.Name()), 
                    *((*trajParm)[atom].Name()));
          ++numMismatch;
        }
      }
      ++atom;
    }
    if (Frames==0) {
      // First frame #atoms 
      pdbAtom_ = atom;
    } else {
      // Check that # atoms read in this frame match the first frame
      if (atom>0 && pdbAtom_!=atom) {
        mprintf("Warning: PDB %s: Reading frame %i, got %i atoms, expected %i.\n",
                file_.Filename().base(), Frames+1, atom, pdbAtom_);
        mprintf("Warning: Only using frames 1-%i\n", Frames);
        scanPDB = false;
        break;
      }
    }  
    if (scanPDB) ++Frames;
  }
  file_.CloseFile(); 
  if (Frames<1) {
    mprinterr("Error: PDB %s: No frames read. atom=%i expected %i.\n",
              file_.Filename().base(), atom, trajParm->Natom());
    return TRAJIN_ERR;
  }
  if (debug_>0) mprintf("Traj_PDBfile: %s has %i atoms, %i frames.\n",
                        file_.Filename().base(), pdbAtom_, Frames);
  // Report mismatches of pdb atom names against parm names
  if (numMismatch > 0)
    mprintf("Warning: In PDB file %s: %i name mismatches with parm %s.\n",
            file_.Filename().base(), numMismatch, trajParm->c_str());
  // Set traj info - no velocity, temperature, time
  SetCoordInfo( CoordinateInfo(boxInfo, false, false, false) );
  return Frames;
}

// Traj_PDBfile::readFrame()
/** Read frame (model) from PDB file. */
int Traj_PDBfile::readFrame(int set, Frame& frameIn)
{
  int atom;
  if (set < currentSet_) {
    file_.Rewind();
    currentSet_ = 0;
  }
  // Position file at group of ATOM keywords for specified set
  while (currentSet_ < set) {
    atom = 0;
    while (atom < pdbAtom_) {
      if ( file_.NextRecord() == PDBfile::END_OF_FILE ) return 1;
      if ( file_.RecType() == PDBfile::ATOM ) ++atom;
    }
    currentSet_++;
  }
  atom = 0;
  double *Xptr = frameIn.xAddress(); 
  while (atom < pdbAtom_) {
    if ( file_.NextRecord() == PDBfile::END_OF_FILE ) return 1;
    // Skip non-ATOM records
    if ( file_.RecType() == PDBfile::CRYST1 )
      file_.pdb_Box( frameIn.bAddress() );
    else if ( file_.RecType() == PDBfile::ATOM ) {
      // Read current PDB record XYZ into Frame
      file_.pdb_XYZ( Xptr );
      ++atom; 
      Xptr += 3;
    }
  }
  currentSet_++;

  return 0;
}

void Traj_PDBfile::WriteHelp() {
  mprintf("\tdumpq:       Write atom charge/radius in occupancy/B-factor columns (PQR format).\n"
          "\tpdbres:      Use PDB V3 residue names.\n"
          "\tpdbatom:     Use PDB V3 atom names.\n"
          "\tpdbv3:       Use PDB V3 residue/atom names.\n"
          "\tteradvance:  Increment record (atom) # for TER records (default no).\n"
          "\tmodel:       Write to single file separated by MODEL records.\n"
          "\tmulti:       Write each frame to separate files.\n"
          "\tchainid <c>: Write character 'c' in chain ID column.\n"
          "\tsg <group>:  Space group for CRYST1 record, only if box coordinates written.\n");
}

// Traj_PDBfile::processWriteArgs()
int Traj_PDBfile::processWriteArgs(ArgList& argIn) {
  pdbWriteMode_ = SINGLE;
  if (argIn.hasKey("dumpq")) {
   dumpq_ = true; 
   dumpr_ = true;
  }
  pdbres_ = argIn.hasKey("pdbres");
  pdbatom_ = argIn.hasKey("pdbatom");
  if (argIn.hasKey("pdbv3")) {
    pdbres_ = true;
    pdbatom_ = true;
  }
  if (argIn.hasKey("teradvance")) ter_num_ = 1;
  if (argIn.hasKey("model")) pdbWriteMode_ = MODEL;
  if (argIn.hasKey("multi")) pdbWriteMode_ = MULTI;
  space_group_ = argIn.GetStringKey("sg");
  std::string temp = argIn.GetStringKey("chainid");
  if (!temp.empty()) chainchar_ = temp[0];
  return 0;
}

// Traj_PDBfile::setupTrajout()
/** Set parm information needed for write, and check write mode against
  * number of frames to be written.
  */ 
int Traj_PDBfile::setupTrajout(std::string const& fname, Topology* trajParm,
                               CoordinateInfo const& cInfoIn,
                               int NframesToWrite, bool append)
{
  if (trajParm==0) return 1;
  SetCoordInfo( cInfoIn );
  pdbTop_ = trajParm;
  pdbAtom_ = pdbTop_->Natom();
  // Set up file
  if (append && pdbWriteMode_ != MULTI) {
    if (file_.SetupAppend( fname, debug_)) return 1;
  } else {
    if (append && pdbWriteMode_ == MULTI)
      mprintf("Warning: 'append' not compatible with 'multi' pdb write.\n");
    if (file_.SetupWrite( fname, debug_ )) return 1;
  }
  // Set a chainID for each atom
  // TODO: Set different chain ID for solute mols and solvent
  chainID_.clear();
  if (chainchar_ == ' ') {
    chainID_.reserve( pdbAtom_ );
    for (Topology::atom_iterator atom = trajParm->begin(); atom != trajParm->end(); ++atom)
      chainID_.push_back( (*atom).ChainID() );
  } else
    chainID_.resize(pdbAtom_, chainchar_);
        
  // Save residue names. If pdbres specified convert to PDBV3 residue names.
  resNames_.clear();
  resNames_.reserve( trajParm->Nres() );
  if (pdbres_) {
    for (Topology::res_iterator res = trajParm->ResStart();
                                res != trajParm->ResEnd(); ++res) {
      NameType rname = (*res).Name();
      // convert protein residue names back to more like PDBV3 format:
      if (rname == "HID " || rname == "HIE " ||
          rname == "HIP " || rname == "HIC "   )
        rname = "HIS ";
      else if (rname == "CYX " || rname == "CYM ")
        rname = "CYS ";
      else if (rname == "MEM ") 
        rname = "MET ";
      else if (rname == "ASH ")
        rname = "ASP ";
      else if (rname == "GLH ")
        rname = "GLU ";
      // also for nucleic acid names:
      else if ( rname[2] == ' ' && rname[3] == ' ' ) {
        // RNA names
        if      ( rname[0] == 'G' ) rname="  G ";
        else if ( rname[0] == 'C' ) rname="  C ";
        else if ( rname[0] == 'A' ) rname="  A ";
        else if ( rname[0] == 'U' ) rname="  U ";
      } else if ( rname[0] == 'D' ) {
        // DNA names
        if      ( rname[1] == 'G' ) rname=" DG ";
        else if ( rname[1] == 'C' ) rname=" DC ";
        else if ( rname[1] == 'A' ) rname=" DA ";
        else if ( rname[1] == 'T' ) rname=" DT ";
      }
      resNames_.push_back( rname );
    }
  } else {
    for (Topology::res_iterator res = trajParm->ResStart();
                                res != trajParm->ResEnd(); ++res)
      resNames_.push_back( (*res).Name() );
  }
  // If number of frames to write > 1 and not doing 1 pdb file per frame,
  // set write mode to MODEL
  if (append || (pdbWriteMode_==SINGLE && NframesToWrite>1)) 
    pdbWriteMode_ = MODEL;
  // TODO: Setup title
  // Open here if writing to single file
  if (pdbWriteMode_ != MULTI) {
    if ( file_.OpenFile() ) return 1;
    if (!Title().empty()) file_.WriteTITLE( Title() );
  }
  write_cryst1_ = (CoordInfo().TrajBox().Type() != Box::NOBOX);
  if (write_cryst1_) {
    if (pdbWriteMode_==MODEL)
      mprintf("Warning: For PDB with MODEL, box coords for first frame only will be written to CRYST1.\n");
    if (space_group_.empty())
      mprintf("Warning: No PDB space group specified.\n");
  }
  return 0;
}

// Traj_PDBfile::writeFrame()
/** Write the frame (model) to PDB file. */
int Traj_PDBfile::writeFrame(int set, Frame const& frameOut) {
  if (pdbWriteMode_==MULTI) {
    // If writing 1 pdb per frame set up output filename and open
    if (file_.OpenWriteNumbered( set + 1 )) return 1;
    if (!Title().empty()) 
      file_.WriteTITLE( Title() );
    if (write_cryst1_)
      file_.WriteCRYST1( frameOut.BoxCrd().boxPtr(), space_group_.c_str() );
  } else {
    // Write box coords, first frame only.
    if (write_cryst1_) {
      file_.WriteCRYST1( frameOut.BoxCrd().boxPtr(), space_group_.c_str() );
      write_cryst1_ = false;
    }
  }
  // If specified, write MODEL keyword
  if (pdbWriteMode_==MODEL)
    file_.WriteMODEL(set + 1); 

  float Occ = 1.0; 
  float B = 0.0;
  int anum = 1; // Actual PDB atom number
  int aidx = 0; // Atom index in topology
  Topology::mol_iterator mol = pdbTop_->MolStart();
  int lastAtomInMol;
  if (pdbTop_->Nmol() > 0)
    lastAtomInMol = (*mol).EndAtom();
  else
    lastAtomInMol = -1;
  const double *Xptr = frameOut.xAddress();
  for (Topology::atom_iterator atom = pdbTop_->begin(); atom != pdbTop_->end(); ++atom, ++aidx) {
    int res = atom->ResNum();
    // If this atom belongs to a new molecule print a TER card
    // Use res instead of res+1 since this TER belongs to last mol/res
    if (aidx == lastAtomInMol) {
      file_.WriteTER( anum, resNames_[res-1], chainID_[aidx-1], pdbTop_->Res(res-1).OriginalResNum() );
      anum += ter_num_;
      ++mol;
      lastAtomInMol = mol->EndAtom();
    }
    if (!pdbTop_->Extra().empty()) {
      Occ = pdbTop_->Extra()[aidx].Occupancy();
      B   = pdbTop_->Extra()[aidx].Bfactor();
    }
    if (dumpq_) Occ = (float) atom->Charge();
    if (dumpr_) B = (float) atom->GBRadius();
    // If pdbatom change amber atom names to pdb v3
    NameType atomName = atom->Name();
    if (pdbatom_) {
      if      (atomName == "H5'1") atomName = "H5'";
      else if (atomName == "H5'2") atomName = "H5''";
      else if (atomName == "H2'1") atomName = "H2'";
      else if (atomName == "H2'2") atomName = "H2''";
      else if (atomName == "O1P ") atomName = "OP1";
      else if (atomName == "O2P ") atomName = "OP2";
      else if (atomName == "H5T ") atomName = "HO5'";
      else if (atomName == "H3T ") atomName = "HO3'";
      else if (atomName == "HO'2") atomName = "HO2'";
    }
    file_.WriteCoord(PDBfile::ATOM, anum++, atomName, resNames_[res],
                     chainID_[aidx], pdbTop_->Res(res).OriginalResNum(), 
                     Xptr[0], Xptr[1], Xptr[2], Occ, B, 
                     atom->ElementName(), 0, dumpq_);
    Xptr += 3;
  }
  if (pdbWriteMode_==MULTI) {
    // If writing 1 pdb per frame, close output file
    file_.WriteEND();
    file_.CloseFile();
  } else if (pdbWriteMode_==MODEL) {
    // If MODEL keyword was written, write corresponding ENDMDL record
    file_.WriteENDMDL();
  }

  return 0;
}

// Traj_PDBfile::Info()
void Traj_PDBfile::Info() {
  mprintf("is a PDB file");
  if (pdbWriteMode_ != NONE) {
    if (pdbWriteMode_==MULTI)
      mprintf(" (1 file per frame)");
    else if (pdbWriteMode_==MODEL)
      mprintf(" (1 MODEL per frame)");
    if (dumpq_ && !dumpr_) 
      mprintf(", writing charges to occupancy column");
    else if (dumpr_ && !dumpq_) 
      mprintf(", writing GB radii to B-factor column");
    else if (dumpr_ && dumpq_)
      mprintf(", writing charges/GB radii to occupancy/B-factor columns");
    if (pdbres_ && pdbatom_)
      mprintf(", using PDB V3 res/atom names");
    else if (pdbres_)
      mprintf(", using PDB V3 residue names");
    else if (pdbatom_)
      mprintf(", using PDB V3 atom names");
  }
}
