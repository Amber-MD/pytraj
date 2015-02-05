// Action_Molsurf
#include <cstring> // strcpy, memset
#include "Action_Molsurf.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_Molsurf::Action_Molsurf() {
  sasa_=0;
  atom_ = 0;
  probe_rad_ = 1.4;
  rad_offset_ = 0;

  upper_neighbors=0;
  neighbors=0;
  toruslist=0;
  probelist=0;

  concave_face=0;
  saddle_face=0;
  convex_face=0;
  cone_face=0;
  broken_concave_face=0;
  concave_cycle=0;

  vertexlist=0;
  concave_edge_list=0;
  convex_edge_list=0;
  convex_circle_list=0;
  concave_circle_list=0;

  cyclelist=0;
  low_torus=0;
  cusp_edge=0;
  cusp_pair=0;
} 

// DESTRUCTOR
Action_Molsurf::~Action_Molsurf() {
  ClearMemory();
  if (atom_!=0) delete[] atom_;
}

void Action_Molsurf::Help() {
  mprintf("\t[<name>] [<mask1>] [out filename] [probe <probe_rad>] [offset <rad_offset>]\n"
          "  Calculate Connolly surface area of atoms in <mask1>\n");
}

// MolSurf::ClearMemory()
/// Clear mem used by molsurf data structures
void Action_Molsurf::ClearMemory() {
  if (upper_neighbors!=0) delete[] upper_neighbors;
  if (neighbors!=0) delete[] neighbors;
  if (probelist!=0) delete[] probelist;
  if (toruslist!=0) delete[] toruslist;
  if (convex_circle_list!=0) delete[] convex_circle_list;
  if (concave_circle_list!=0) delete[] concave_circle_list;
  if (concave_face!=0) delete[] concave_face;
  if (convex_face!=0) delete[] convex_face;
  if (saddle_face!=0) delete[] saddle_face;
  if (cone_face!=0) delete[] cone_face;
  if (broken_concave_face!=0) delete[] broken_concave_face;
  if (concave_cycle!=0) delete[] concave_cycle;
  if (cyclelist!=0) delete[] cyclelist;
  if (vertexlist!=0) delete[] vertexlist;
  if (concave_edge_list!=0) delete[] concave_edge_list;
  if (convex_edge_list!=0) delete[] convex_edge_list;
  if (low_torus!=0) delete[] low_torus;
  if (cusp_edge!=0) delete[] cusp_edge;
  if (cusp_pair!=0) delete[] cusp_pair;
}

// Action_Molsurf::AllocateMemory()
/// Allocate mem used by molsurf internal data structures
int Action_Molsurf::AllocateMemory() {
  int error_status = 0;
  int natm_sel = Mask1_.Nselected();
  // ---------- Allocate Memory For molsurf routines --------------------
  upper_neighbors = new NEIGHBOR_TORUS[ NUM_NEIGHBOR * natm_sel ];
  neighbors = new NEIGHBOR [ NUM_NEIGHBOR * natm_sel ]; 
  probelist = new PROBE[ NUM_PROBE * natm_sel ];
  toruslist = new TORUS[ NUM_TORUS * natm_sel ];
  convex_circle_list = new CIRCLE[ NUM_CIRCLE * natm_sel ];
  concave_circle_list = new CIRCLE[ NUM_CIRCLE * natm_sel ];
  concave_face = new CONCAVE_FACE[ NUM_FACE * natm_sel ];
  convex_face = new CONVEX_FACE[ NUM_FACE * natm_sel ];
  saddle_face = new SADDLE_FACE[ NUM_FACE * natm_sel ];
  cone_face = new CONE_FACE[ NUM_FACE * natm_sel ];  
  broken_concave_face = new BROKEN_CONCAVE_FACE[ NUM_FACE * natm_sel ];
  concave_cycle = new CONCAVE_CYCLE[ NUM_CYCLE * natm_sel ];
  cyclelist = new CYCLE[ NUM_CYCLE * natm_sel ];
  vertexlist = new VERTEX[ NUM_VERTEX * natm_sel ];
  concave_edge_list = new EDGE[ NUM_EDGE * natm_sel ];
  convex_edge_list = new EDGE[ NUM_EDGE * natm_sel ];
  low_torus = new LOW_TORUS[ NUM_TORUS * natm_sel ];
  // NOTE: cusp_edge must be initialized. Handled in action before 1st call
  cusp_edge = new CUSP_EDGE[ NUM_EDGE * natm_sel ];
  cusp_pair = new CUSP_PAIR[ NUM_CUSP * natm_sel ]; 
  // ------------------------------------------------------------
  return error_status;
}

// Action_Molsurf::Init()
Action::RetType Action_Molsurf::Init(ArgList& actionArgs, TopologyList* PFL, DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  // Get keywords
  DataFile* outfile = DFL->AddDataFile( actionArgs.GetStringKey("out"), actionArgs );
  probe_rad_ = actionArgs.getKeyDouble("probe",1.4);
  rad_offset_ = actionArgs.getKeyDouble("offset",0.0);
  std::string submaskString = actionArgs.GetStringKey("submask");
  while (!submaskString.empty()) {
    SubMasks_.push_back( AtomMask(submaskString) );
    submaskString = actionArgs.GetStringKey("submask");
  }

  // Get Masks
  Mask1_.SetMaskString( actionArgs.GetMaskNext() );

  // Dataset to store angles
  sasa_ = DSL->AddSet(DataSet::DOUBLE, actionArgs.GetStringNext(),"MSURF");
  if (sasa_==0) return Action::ERR;
  // Add dataset to data file list
  if (outfile != 0) outfile->AddSet(sasa_);
  // Submask string data sets
  for (Marray::const_iterator mask = SubMasks_.begin(); mask != SubMasks_.end(); ++mask) {
    DataSet* ds = DSL->AddSetIdxAspect(DataSet::FLOAT, sasa_->Name(), 
                                       mask-SubMasks_.begin(), "submask", mask->MaskExpression());
    if (ds == 0) return Action::ERR;
    if (outfile != 0) outfile->AddSet(ds);
    SubData_.push_back( ds );
  }

  mprintf("    MOLSURF: [%s] Probe Radius=%.3lf\n",Mask1_.MaskString(),probe_rad_);
  if (rad_offset_>0)
    mprintf("\tRadii will be incremented by %.3lf\n",rad_offset_);
  if (!SubMasks_.empty())
    mprintf("\tThe contribution to the total area for %zu sub-masks will be calculated.\n",
            SubMasks_.size());

  return Action::OK;
}

// Action_Molsurf::Setup()
/** Set mask up for this parmtop. Allocate the ATOM structure array used 
  * by the molsurf C routines and set everything but the atom coords.
  */
Action::RetType Action_Molsurf::Setup(Topology* currentParm, Topology** parmAddress) {
  if ( currentParm->SetupIntegerMask(Mask1_) ) return Action::ERR;
  if (Mask1_.None()) {
    mprintf("Warning: Mask contains 0 atoms.\n");
    return Action::ERR;
  }

  mprintf("\tCalculating surface area for %i atoms.\n",Mask1_.Nselected());
  // NOTE: If Mask is * dont include any solvent?
  // Set up submasks
  if (!SubMasks_.empty()) {
    // Set up an array so we can index from parm atoms to Mask1 indices.
    mask1idx_.assign( currentParm->Natom(), -1 );
    int idx = 0;
    for (int at = 0; at != currentParm->Natom(); at++) {
      if (at == Mask1_[idx])
        mask1idx_[at] = idx++;
      mprintf("DBG: mask1idx_[%i]= %i\n", at, mask1idx_[at]);
      if (idx == Mask1_.Nselected()) break;
    }
    for (Marray::iterator mask = SubMasks_.begin(); mask != SubMasks_.end(); ++mask)
    {
      if (currentParm->SetupIntegerMask(*mask)) return Action::ERR;
      if (mask->None())
        mprintf("Warning: No atoms selected for mask '%s'\n", mask->MaskString());
      else {
        mask->MaskInfo();
        // Check that submask atoms are all in common with main mask.
        // TODO: Use mask1idx?
        if (Mask1_.NumAtomsInCommon(*mask) != mask->Nselected()) {
          mprinterr("Error: Sub-mask '%s' atoms are not a subset of main mask '%s'\n",
                    mask->MaskString(), Mask1_.MaskString());
          return Action::ERR;
        }
      }
    }
  }
  // The ATOM structure is how molsurf organizes atomic data. Allocate
  // here and fill in parm info. Coords will be filled in during action. 
  if (atom_!=0) delete[] atom_;
  atom_ = new ATOM[ Mask1_.Nselected() ];
  if (atom_==0) {
    mprinterr("Error: Molsurf::Setup Could not allocate memory for ATOMs.\n");
    return Action::ERR;
  }
  // Set up parm info for atoms in mask
  if ( (*currentParm)[0].GBRadius() == 0 ) {
    mprinterr("Error: Molsurf::Setup: Molsurf requires radii, but no radii in %s\n",
              currentParm->c_str());
    return Action::ERR;
  }
  ATOM *atm_ptr = atom_;
  for (AtomMask::const_iterator parmatom = Mask1_.begin();
                                parmatom != Mask1_.end();
                                ++parmatom, ++atm_ptr)
  {
    atm_ptr->anum = *parmatom + 1; // anum is for debug output only, atoms start from 1
    const Atom patom = (*currentParm)[*parmatom];
    int nres = patom.ResNum();
    atm_ptr->rnum = nres+1; // for debug output only, residues start from 1
    patom.Name().ToBuffer( atm_ptr->anam );
    strcpy(atm_ptr->rnam, currentParm->Res(nres).c_str());
    atm_ptr->pos[0] = 0;
    atm_ptr->pos[1] = 0;
    atm_ptr->pos[2] = 0;
    atm_ptr->q = patom.Charge();
    atm_ptr->rad = patom.GBRadius() + rad_offset_;
    atm_ptr->area = 0.0;
  }

  // De-allocate memory first since # atoms may have changed
  ClearMemory();
  if (AllocateMemory()) return Action::ERR;

  //if (debug>0) memory_usage();

  return Action::OK;  
}

// Action_Molsurf::DoAction()
Action::RetType Action_Molsurf::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {
  // Set up coordinates for atoms in mask
  ATOM *atm_ptr = atom_;
  for (AtomMask::const_iterator maskatom = Mask1_.begin(); maskatom != Mask1_.end();
       ++maskatom, ++atm_ptr)
  {
    const double* XYZ = currentFrame->XYZ( *maskatom );
    std::copy( XYZ, XYZ+3, atm_ptr->pos );
    atm_ptr->area = 0.0;
  }

  // NOTE: cusp_edge is the only data structure that requires initialization 
  memset( cusp_edge, 0, NUM_EDGE * Mask1_.Nselected() * sizeof(CUSP_EDGE));
  double molsurf_sasa = molsurf( probe_rad_, atom_, Mask1_.Nselected(),
                                 upper_neighbors, neighbors,
                                 toruslist, probelist, concave_face,
                                 saddle_face, convex_face,
                                 cone_face, broken_concave_face,
                                 concave_cycle, vertexlist,
                                 concave_edge_list, convex_edge_list,
                                 convex_circle_list, concave_circle_list,
                                 cyclelist, low_torus, cusp_edge,
                                 cusp_pair); 

  sasa_->Add(frameNum, &molsurf_sasa);

  // Get any submask areas
  DSarray::const_iterator ds = SubData_.begin(); 
  for (Marray::const_iterator mask = SubMasks_.begin(); mask != SubMasks_.end(); ++mask, ++ds)
  {
    double areaSum = 0.0;
    for (AtomMask::const_iterator maskatom = mask->begin(); maskatom != mask->end(); ++maskatom)
      areaSum += atom_[ mask1idx_[*maskatom] ].area;
    float fval = (float)areaSum;
    (*ds)->Add(frameNum, &fval);
  }

  return Action::OK;
} 
