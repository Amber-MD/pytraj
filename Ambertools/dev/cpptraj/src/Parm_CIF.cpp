#include "Parm_CIF.h"
#include "CIFfile.h"
#include "StringRoutines.h"
#include "CpptrajStdio.h"

// NOTE: MUST correspond to EntryType!
const char* Parm_CIF::Entries[] = {
  "label_atom_id", "label_comp_id", "Cartn_x", "Cartn_y", "Cartn_z", 
  "label_seq_id", "label_asym_id"
};

// Parm_CIF::ReadParm()
int Parm_CIF::ReadParm(std::string const& fname, Topology &TopIn) {
  CIFfile infile;
  CIFfile::DataBlock::data_it line;

  if (infile.Read( fname, debug_ )) return 1;
  CIFfile::DataBlock const& block = infile.GetDataBlock("_atom_site");
  if (block.empty()) {
    mprinterr("Error: CIF data block '_atom_site' not found.\n");
    return 1;
  }
  // Does this CIF contain multiple models?
  int Nmodels = 0;
  int model_col = block.ColumnIndex("pdbx_PDB_model_num");
  if (model_col != -1) {
    line = block.end();
    --line;
    Nmodels = convertToInteger( (*line)[model_col] );
    if (Nmodels > 1)
      mprintf("Warning: CIF '%s' contains %i models. Using first model for topology.\n", 
              fname.c_str(), Nmodels);
  }
  // Get essential columns
  int COL[NENTRY];
  for (int i = 0; i < (int)NENTRY; i++) {
    COL[i] = block.ColumnIndex(Entries[i]);
    if (COL[i] == -1) {
      mprinterr("Error: In CIF file '%s' could not find entry '%s' in block '%s'\n",
                fname.c_str(), Entries[i], block.Header().c_str());
      return 1;
    }
    if (debug_>0) mprintf("DEBUG: '%s' column = %i\n", Entries[i], COL[i]);
  }
  // Get optional columns
  int occ_col = block.ColumnIndex("occupancy");
  int bfac_col = block.ColumnIndex("B_iso_or_equiv");
  std::vector<AtomExtra> extra;

  // Loop over all atom sites
  int current_res = 0;
  double XYZ[3];
  double occupancy = 1.0;
  double bfactor = 0.0;
  for (line = block.begin(); line != block.end(); ++line) {
    // If more than 1 model check if we are done.
    if (Nmodels > 1) {
      if ( convertToInteger( (*line)[model_col] ) > 1 )
        break;
    }
    if (occ_col != -1) occupancy = convertToDouble( (*line)[ occ_col ] );
    if (bfac_col != -1) bfactor = convertToDouble( (*line)[ bfac_col ] );
    extra.push_back( AtomExtra(occupancy, bfactor) );
    XYZ[0] = convertToDouble( (*line)[ COL[X] ] );
    XYZ[1] = convertToDouble( (*line)[ COL[Y] ] );
    XYZ[2] = convertToDouble( (*line)[ COL[Z] ] );
    NameType currentResName( (*line)[ COL[RNAME] ] );
    // It seems that in some CIF files, there doesnt have to be a residue
    // number. Check if residue name has changed.
    if ( (*line)[ COL[RNUM] ][0] == '.' ) {
      Topology::res_iterator lastResidue = TopIn.ResEnd();
      --lastResidue;
      if ( currentResName != (*lastResidue).Name() )
        current_res = TopIn.Nres() + 1;
    } else
      current_res = convertToInteger( (*line)[ COL[RNUM] ] );
    TopIn.AddTopAtom( Atom((*line)[ COL[ANAME] ], (*line)[ COL[CHAINID] ][0], "  "),
                      current_res, currentResName, XYZ );
  }
  TopIn.SetExtraAtomInfo( 0, extra );
  // Get title. 
  CIFfile::DataBlock const& entryblock = infile.GetDataBlock("_entry");
  std::string ciftitle;
  if (!entryblock.empty())
    ciftitle = entryblock.Data("id");
  TopIn.SetParmName( ciftitle, infile.CIFname() );
  // Get unit cell parameters if present.
  CIFfile::DataBlock const& cellblock = infile.GetDataBlock("_cell");
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
    TopIn.SetParmBox( Box(cif_box) ); 
  }
  
  return 0;
}

// Parm_CIF::ID_ParmFormat()
bool Parm_CIF::ID_ParmFormat(CpptrajFile& fileIn) {
  return CIFfile::ID_CIF( fileIn );
}
