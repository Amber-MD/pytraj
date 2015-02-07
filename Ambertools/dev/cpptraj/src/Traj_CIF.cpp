// Traj_CIF
#include "Traj_CIF.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // convertToInteger, convertToDouble

bool Traj_CIF::ID_TrajFormat(CpptrajFile& fileIn) {
  return CIFfile::ID_CIF( fileIn );
}

// Traj_CIF::openTrajin()
// FIXME: Currently CIF files are always read in and stored in memory.
//        No write. Everything handled by setupTrajin and readFrame.
int Traj_CIF::openTrajin() {
  return 0;
}

// Traj_CIF::setupTrajin()
/** Read in entire CIF file. */
int Traj_CIF::setupTrajin(std::string const& fname, Topology* trajParm)
{
  if (file_.Read( fname, debug_ )) return TRAJIN_ERR;
  CIFfile::DataBlock const& block = file_.GetDataBlock("_atom_site");
  if (block.empty()) return TRAJIN_ERR;
  // Get coordinate x/y/z columns
  Cartn_x_col_ = block.ColumnIndex("Cartn_x");
  Cartn_y_col_ = block.ColumnIndex("Cartn_y");
  Cartn_z_col_ = block.ColumnIndex("Cartn_z");
  if (Cartn_x_col_ == -1 || Cartn_y_col_ == -1 || Cartn_z_col_ == -1) {
    mprinterr("Error: Could not find Cartn_x|y|z columns in CIF file.\n");
    return TRAJIN_ERR;
  }
  // Determine # atoms and # models
  Nmodels_ = 0;
  int model_col = block.ColumnIndex("pdbx_PDB_model_num");
  int id_col    = block.ColumnIndex("id");
  if (id_col == -1) {
    mprinterr("Error: No ID column found in _atom_site block.\n");
    return TRAJIN_ERR;
  }
  CIFfile::DataBlock::data_it line = block.end();
  --line; // Go to last _atom_site line
  // totalAtoms will == #models * natom
  int totalAtoms = convertToInteger( (*line)[id_col] );
  if (model_col == -1) {
    // No model # column; assume 1 model
    Nmodels_ = 1;
  } else {
    Nmodels_ = convertToInteger( (*line)[model_col] );
  }
  if ( (totalAtoms % Nmodels_) != 0 ) {
    mprintf("Warning: Total number of atoms in CIF (%i) is not divisible by\n"
            "Warning:  number of models (%i). This indicates the number of atoms\n"
            "Warning:  in each model is not the same. Only reading %i atoms of\n"
            "Warning:  the first model.\n",
            totalAtoms, Nmodels_, trajParm->Natom());
    Nmodels_ = 1;
    Natoms_ = trajParm->Natom();
  } else {
    Natoms_ = totalAtoms / Nmodels_;
    if (Natoms_ != trajParm->Natom()) {
      mprinterr("Error: Number of atoms in CIF (%i) does not equal number of atoms\n"
                "Error: in associated topology '%s' (%i)\n", Natoms_,
                trajParm->c_str(), trajParm->Natom());
      return TRAJIN_ERR;
    }
  }
  mprintf("\t%i atoms, %i models.\n", Natoms_, Nmodels_);
  // Get unit cell parameters if present.
  boxInfo_.SetNoBox();
  CIFfile::DataBlock const& cellblock = file_.GetDataBlock("_cell");
  if (!cellblock.empty()) {
    double cif_box[6];
    cif_box[0] = convertToDouble( cellblock.Data("length_a") );
    cif_box[1] = convertToDouble( cellblock.Data("length_b") );
    cif_box[2] = convertToDouble( cellblock.Data("length_c") );
    cif_box[3] = convertToDouble( cellblock.Data("angle_alpha") );
    cif_box[4] = convertToDouble( cellblock.Data("angle_beta" ) );
    cif_box[5] = convertToDouble( cellblock.Data("angle_gamma") );
    mprintf("\tRead cell info from CIF: a=%g b=%g c=%g alpha=%g beta=%g gamma=%g\n",
              cif_box[0], cif_box[1], cif_box[2], cif_box[3], cif_box[4], cif_box[5]);
    boxInfo_.SetBox( cif_box);
  }
  // Set traj info - No velocity, temperature, time.
  SetCoordInfo( CoordinateInfo( boxInfo_, false, false, false ) );
  // Get title. 
  CIFfile::DataBlock const& entryblock = file_.GetDataBlock("_entry");
  if (!entryblock.empty())
    SetTitle( entryblock.Data("id") );

  return Nmodels_;
}

// Traj_CIF::readFrame()
int Traj_CIF::readFrame(int set, Frame& frameIn) {
  //if (set >= Nmodels_) return 1;
  // FIXME: Shouldnt have to always search for the block
  CIFfile::DataBlock const& block = file_.GetDataBlock("_atom_site");
  CIFfile::DataBlock::data_it line = block.begin() + (set * Natoms_);
  CIFfile::DataBlock::data_it end  = line + Natoms_;
  double *Xptr = frameIn.xAddress(); 
  for (; line != end; ++line) {
    *(Xptr++) = convertToDouble( (*line)[Cartn_x_col_] );
    *(Xptr++) = convertToDouble( (*line)[Cartn_y_col_] );
    *(Xptr++) = convertToDouble( (*line)[Cartn_z_col_] );
  }
  std::copy( boxInfo_.boxPtr(), boxInfo_.boxPtr() + 6, frameIn.bAddress() );
  return 0;
}

// Traj_CIF::info()
void Traj_CIF::Info() {
  mprintf("is a CIF file");
}
