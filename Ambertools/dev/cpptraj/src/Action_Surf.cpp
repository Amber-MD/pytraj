// Action_Surf 
#include <cmath>  // sqrt
#include <cctype> // toupper
#include "Action_Surf.h"
#include "Constants.h" // For FOURPI, TWOPI
#include "CpptrajStdio.h"
#include "DistRoutines.h"
#ifdef _OPENMP
#  include "omp.h"
#endif
// CONSTRUCTOR
Action_Surf::Action_Surf() : surf_(0) {} 

void Action_Surf::Help() {
  mprintf("\t<name> <mask1> [out filename]\n"
          "  Calculate LCPO surface area of atoms in <mask1>\n");
}

// Action_Surf::Init()
Action::RetType Action_Surf::Init(ArgList& actionArgs, TopologyList* PFL, DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  // Get keywords
  DataFile* outfile = DFL->AddDataFile( actionArgs.GetStringKey("out"), actionArgs);

  // Get Masks
  Mask1_.SetMaskString( actionArgs.GetMaskNext() );

  // Dataset to store surface area 
  surf_ = DSL->AddSet(DataSet::DOUBLE, actionArgs.GetStringNext(), "SA");
  if (surf_==0) return Action::ERR;
  // Add dataset to data file list
  if (outfile != 0) outfile->AddSet( surf_ );

  mprintf("    SURF: Calculating surface area for atoms in mask [%s]\n",Mask1_.MaskString());
  mprintf("#Citation: Weiser, J.; Shenkin, P. S.; Still, W. C.; \"Approximate atomic\n"
          "#          surfaces from linear combinations of pairwise overlaps (LCPO).\"\n"
          "#          J. Comp. Chem. (1999), V.20, pp.217-230.\n");

  return Action::OK;
}

// Action_Surf::Setup()
/** Set LCPO surface area calc parameters for this parmtop if not already set. 
  * Get the mask, and check that the atoms in mask belong to solute. 
  */
Action::RetType Action_Surf::Setup(Topology* currentParm, Topology** parmAddress) {
  SurfInfo SI;

  if (currentParm->SetupIntegerMask( Mask1_ )) return Action::ERR;
  if (Mask1_.None()) {
    mprintf("Warning: Mask '%s' corresponds to 0 atoms.\n", Mask1_.MaskString());
    return Action::ERR;
  }
  mprintf("\tLCPO surface area will be calculated for %i atoms.\n",Mask1_.Nselected());

  // Setup surface area calc for this parm.
  // Check that each atom in Mask1 is part of solute
  // Create a separate mask for building the atom neighbor list that only
  // includes atoms with vdw radius > 2.5. Consider all atoms for icosa, only
  // non-H's for LCPO
  atomi_neighborMask_.ResetMask();
  atomi_noNeighborMask_.ResetMask();
  atomj_neighborMask_.ResetMask();
  SurfaceInfo_neighbor_.clear();
  SurfaceInfo_noNeighbor_.clear();
  int soluteAtoms = 0;
  for (AtomMask::const_iterator atomi = Mask1_.begin(); atomi!=Mask1_.end(); atomi++) {
    int molNum = (*currentParm)[ *atomi ].MolNum();
    if (currentParm->Mol( molNum ).IsSolvent()) {
      mprinterr("Error: Atom %i in mask %s does not belong to solute.\n",
                *atomi+1, Mask1_.MaskString());
      return Action::ERR;
    }
    ++soluteAtoms;
    SetAtomLCPO( *currentParm, *atomi, &SI ); 
    if (SI.vdwradii > 2.5) {
      atomi_neighborMask_.AddAtom(*atomi);
      SurfaceInfo_neighbor_.push_back( SI );
    } else {
      atomi_noNeighborMask_.AddAtom(*atomi);
      SurfaceInfo_noNeighbor_.push_back( SI );
    }
  }
  mprintf("\t%i solute atoms.\n",soluteAtoms);
  if (soluteAtoms <= 0) {
    mprinterr("Error: No solute atoms in %s.\n",currentParm->c_str());
    return Action::ERR;
  }
  // From all solute atoms, create a second mask for building atom 
  // neighbor list that only includes atoms with vdw radius > 2.5.
  VDW_.clear();
  VDW_.reserve( soluteAtoms );
  if (currentParm->Nmol() < 1) {
    mprinterr("Error: Topology %s has no molecule information, LCPO surface area\n"
              "Error:   cannot be calculated. Try using 'fixatomorder' prior to 'surf' command.\n",
              currentParm->c_str());
    return Action::ERR;
  }
  for (Topology::mol_iterator mol = currentParm->MolStart(); 
                              mol != currentParm->MolEnd(); ++mol)
  {
    if (!mol->IsSolvent()) {
      for (int atomj=mol->BeginAtom(); atomj != mol->EndAtom(); atomj++) {
        SetAtomLCPO( *currentParm, atomj, &SI );
        VDW_.push_back( SI.vdwradii );
        if (SI.vdwradii > 2.5)
          atomj_neighborMask_.AddAtom(atomj);
      }
    }
  }
  return Action::OK;  
}

// Action_Surf::DoAction()
/** Calculate surface area. */
Action::RetType Action_Surf::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {
  double SA;
  int atomi, idx;
  AtomMask::const_iterator atomj; 
  std::vector<int> ineighbor; // TODO: Class var
  std::vector<double> Distances_i_j; // TODO: Class var
  int max_atomi_neighbormask = atomi_neighborMask_.Nselected();

  // Set up neighbor list for each atom in mask and calc its surface 
  // area contribution. Sum these up to get the total surface area.
  SA = 0.0;
#ifdef _OPENMP
#pragma omp parallel private(atomi,idx,atomj,ineighbor,Distances_i_j) reduction(+: SA)
{
#pragma omp for 
#endif
  for (idx = 0; idx < max_atomi_neighbormask; idx++) {
    atomi = atomi_neighborMask_[idx];
    // Vdw of atom i
    double vdwi = VDW_[atomi];
    // Set up neighbor list for atom i
    ineighbor.clear();
    Distances_i_j.clear();
    for (atomj = atomj_neighborMask_.begin(); atomj != atomj_neighborMask_.end(); atomj++)
    {
      if (atomi != *atomj) {
        double dij = sqrt( DIST2_NoImage(currentFrame->XYZ(atomi), currentFrame->XYZ(*atomj)) );
        // Count atoms as neighbors if their VDW radii touch
        if ( (vdwi + VDW_[*atomj]) > dij ) {
          ineighbor.push_back(*atomj);
          Distances_i_j.push_back(dij);
        }
        //mprintf("SURF_NEIG:  %i %i %lf\n",atomi,*atomj,dij);
      }
    }
    // Calculate surface area
    // -------------------------------------------------------------------------
    // Calculate LCPO surface area Ai for atomi:
    // Ai = P1*S1 + P2*Sum(Aij) + P3*Sum(Ajk) + P4*Sum(Aij * Sum(Ajk))
    double sumaij = 0.0;
    double sumajk = 0.0;
    double sumaijajk = 0.0;

    // DEBUG - print neighbor list
    //mprintf("SURF: Neighbors for atom %i:",atomi);
    //for (std::vector<int>::iterator jt = ineighbor.begin(); jt != ineighbor.end(); jt++) {
    //  mprintf(" %i",*jt);
    //}
    //mprintf("\n");

    // Calculate surface area of atom i
    double vdwi2 = vdwi * vdwi;
    double Si = vdwi2 * Constants::FOURPI;

    // Loop over all neighbors of atomi (j)
    // NOTE: Factor through the 2 in aij and ajk?
    std::vector<double>::const_iterator Dij = Distances_i_j.begin();
    for (std::vector<int>::const_iterator jt = ineighbor.begin(); 
                                          jt != ineighbor.end(); 
                                          ++jt, ++Dij) 
    {
      //printf("i,j %i %i\n",atomi + 1,(*jt)+1);
      double vdwj = VDW_[*jt];
      double vdwj2 = vdwj * vdwj;
      double dij = *Dij;
      double tmpaij = vdwi - (dij * 0.5) - ( (vdwi2 - vdwj2)/(2.0 * dij) );
      double aij = Constants::TWOPI * vdwi * tmpaij;
      sumaij += aij;

      // Find which neighbors of atom i (j and k) are themselves neighbors
      double sumajk_2 = 0.0;
      for (std::vector<int>::const_iterator kt = ineighbor.begin(); 
                                            kt != ineighbor.end(); 
                                            kt++) 
      {
        if ( (*kt) == (*jt) ) continue;
        //printf("i,j,k %i %i %i\n",atomi + 1,(*jt)+1,(*kt)+1);
        double vdwk = VDW_[*kt];
        double djk = sqrt(DIST2_NoImage(currentFrame->XYZ(*jt), currentFrame->XYZ(*kt)));
        //printf("%4s%6i%6i%12.8lf\n","DJK ",(*jt)+1,(*kt)+1,djk);
        //printf("%6s%6.2lf%6.2lf\n","AVD ",vdwj,vdwk);
        if ( (vdwj + vdwk) > djk ) {
          double vdw2dif = vdwj2 - (vdwk * vdwk);
          double tmpajk = (2.0*vdwj) - djk - (vdw2dif / djk);
          double ajk = Constants::PI*vdwj*tmpajk;
          //tmpajk = vdwj - (djk *0.5) - ( (vdwj2 - (vdwk * vdwk))/(2.0 * djk) );
          //ajk = 2.0 * PI * vdwi * tmpajk;
          //printf("%4s%6i%6i%12.8lf%12.8lf%12.8lf\n","AJK ",(*jt)+1,(*kt)+1,ajk,vdw2dif,tmpajk);
          sumajk += ajk;
          sumajk_2 += ajk;
        }
      } // END loop over neighbor-neighbor pairs of atom i (kt)

      sumaijajk += (aij * sumajk_2);

      // DEBUG
      //printf("%4s%20.8lf %20.8lf %20.8lf\n","AJK ",aij,sumajk,sumaijajk);

    } // END Loop over neighbors of atom i (jt)
    SA += ( (SurfaceInfo_neighbor_[idx].P1 * Si       ) +
            (SurfaceInfo_neighbor_[idx].P2 * sumaij   ) +
            (SurfaceInfo_neighbor_[idx].P3 * sumajk   ) +
            (SurfaceInfo_neighbor_[idx].P4 * sumaijajk) 
          );
    // -------------------------------------------------------------------------
  } // END Loop over atoms in mask (atomi) 
#ifdef _OPENMP
} // END pragma omp parallel
#endif
  // Second loop over atoms with no neighbors (vdw <= 2.5)
  int idxj=0;
  for (atomj = atomi_noNeighborMask_.begin(); atomj != atomi_noNeighborMask_.end(); atomj++)
  {
    // Vdw of atom i
    double vdwi = VDW_[*atomj];
    // Calculate surface area of atom i
    double vdwi2 = vdwi * vdwi;
    double Si = vdwi2 * Constants::FOURPI; 
    SA += (SurfaceInfo_noNeighbor_[idxj++].P1 * Si);
  }

  surf_->Add(frameNum, &SA);

  return Action::OK;
} 

// -----------------------------------------------------------------------------
// Action_Surf::AssignLCPO()
/** Assign parameters for LCPO method. All radii are incremented by 1.4 Ang.
  */
void Action_Surf::AssignLCPO(SurfInfo *S, double vdwradii, double P1, double P2,
                      double P3, double P4) 
{
  S->vdwradii = vdwradii + 1.4;
  S->P1 = P1;
  S->P2 = P2;
  S->P3 = P3;
  S->P4 = P4;
}

// WarnLCPO()
/// Called when the number of bonds to the atom of type atype is not usual.
static void WarnLCPO(NameType const& atype, int atom, int numBonds) {
  mprintf("Warning: Unusual number of bonds for atom %i (%i), type %s.\n",
          atom, numBonds, *atype);
  mprintf("Using default atom parameters.\n");
}

// Action_Surf::SetAtomLCPO()
/** Set up parameters only used in surface area calcs.
  * Adapted from gbsa=1 method in SANDER, mdread.F90
  * \param currentParm The Topology containing atom information. 
  * \param atidx The atom number to set up parameters for.
  * \param SIptr Address to store the SI parameters.
  */
void Action_Surf::SetAtomLCPO(Topology const& currentParm, int atidx, SurfInfo* SIptr) 
{
  const Atom& atom = currentParm[atidx];
  // Get the number of non-H bonded neighbors to this atom
  int numBonds = 0;
  for (Atom::bond_iterator batom = atom.bondbegin(); batom != atom.bondend(); batom++)
    if ( currentParm[ *batom ].Element() != Atom::HYDROGEN )
      ++numBonds;
  char atype0 = toupper(atom.Type()[0]);
  char atype1 = toupper(atom.Type()[1]);
  // TODO: Only set parameters for solute atoms?
  // Set vdw radii and LCPO parameters for this atom
  switch (atom.Element()) {
    case Atom::CARBON:
      if (atom.Nbonds() == 4) {
        switch ( numBonds ) {
          case 1: AssignLCPO(SIptr, 1.70, 0.77887, -0.28063, -0.0012968, 0.00039328); break;
          case 2: AssignLCPO(SIptr, 1.70, 0.56482, -0.19608, -0.0010219, 0.0002658);  break;
          case 3: AssignLCPO(SIptr, 1.70, 0.23348, -0.072627, -0.00020079, 0.00007967); break;
          case 4: AssignLCPO(SIptr, 1.70, 0.00000, 0.00000, 0.00000, 0.00000); break;
          default: WarnLCPO(atom.Type(),atidx + 1,numBonds);
                   AssignLCPO(SIptr, 1.70, 0.77887, -0.28063, -0.0012968, 0.00039328);
        }
      } else {
        switch ( numBonds ) {
          case 2: AssignLCPO(SIptr, 1.70, 0.51245, -0.15966, -0.00019781, 0.00016392); break;
          case 3: AssignLCPO(SIptr, 1.70, 0.070344, -0.019015, -0.000022009, 0.000016875); break;
          default: WarnLCPO(atom.Type(),atidx + 1,numBonds);
                   AssignLCPO(SIptr, 1.70, 0.77887, -0.28063, -0.0012968, 0.00039328);
        }
      }
      break;
    case Atom::OXYGEN:
      if (atype0=='O' && atype1==' ') 
        AssignLCPO(SIptr, 1.60, 0.68563, -0.1868, -0.00135573, 0.00023743);
      else if (atype0=='O' && atype1=='2')
        AssignLCPO(SIptr, 1.60, 0.88857, -0.33421, -0.0018683, 0.00049372);
      else {
        switch ( numBonds ) {
          case 1: AssignLCPO(SIptr, 1.60, 0.77914, -0.25262, -0.0016056, 0.00035071); break;
          case 2: AssignLCPO(SIptr, 1.60, 0.49392, -0.16038, -0.00015512, 0.00016453); break;
          default: WarnLCPO(atom.Type(),atidx + 1,numBonds);
                   AssignLCPO(SIptr, 1.60, 0.77914, -0.25262, -0.0016056, 0.00035071);
        }
      }
      break;
    case Atom::NITROGEN:
      if (atype0=='N' && atype1=='3') {
        switch ( numBonds ) {
          case 1: AssignLCPO(SIptr, 1.65, 0.078602, -0.29198, -0.0006537, 0.00036247); break;
          case 2: AssignLCPO(SIptr, 1.65, 0.22599, -0.036648, -0.0012297, 0.000080038); break;
          case 3: AssignLCPO(SIptr, 1.65, 0.051481, -0.012603, -0.00032006, 0.000024774); break;
          default: WarnLCPO(atom.Type(),atidx + 1,numBonds);
                   AssignLCPO(SIptr, 1.65, 0.078602, -0.29198, -0.0006537, 0.00036247);
        }
      } else {
        switch ( numBonds ) {
          case 1: AssignLCPO(SIptr, 1.65, 0.73511, -0.22116, -0.00089148, 0.0002523); break;
          case 2: AssignLCPO(SIptr, 1.65, 0.41102, -0.12254, -0.000075448, 0.00011804); break;
          case 3: AssignLCPO(SIptr, 1.65, 0.062577, -0.017874, -0.00008312, 0.000019849); break;
          default: WarnLCPO(atom.Type(),atidx + 1,numBonds);
                   AssignLCPO(SIptr, 1.65, 0.078602, -0.29198, -0.0006537, 0.00036247);
        }
      }
      break;
    case Atom::SULFUR:
      if (atype0=='S' && atype1=='H') 
        AssignLCPO(SIptr, 1.90, 0.7722, -0.26393, 0.0010629, 0.0002179);
      else 
        AssignLCPO(SIptr, 1.90, 0.54581, -0.19477, -0.0012873, 0.00029247);
      break;
    case Atom::PHOSPHORUS:
      switch ( numBonds ) {
        case 3: AssignLCPO(SIptr, 1.90, 0.3865, -0.18249, -0.0036598, 0.0004264); break;
        case 4: AssignLCPO(SIptr, 1.90, 0.03873, -0.0089339, 0.0000083582, 0.0000030381); break;
        default: WarnLCPO(atom.Type(),atidx + 1,numBonds);
                 AssignLCPO(SIptr, 1.90, 0.3865, -0.18249, -0.0036598, 0.0004264);
      }
      break;
    case Atom::HYDROGEN:
      AssignLCPO(SIptr, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000); 
      break;
    default:
      if (atype0=='Z') 
        AssignLCPO(SIptr, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000);
      else if (atype0=='M' && atype1=='G') 
        //  Mg radius = 0.99A: ref. 21 in J. Chem. Phys. 1997, 107, 5422
        //  Mg radius = 1.18A: ref. 30 in J. Chem. Phys. 1997, 107, 5422
        //  Mg radius = 1.45A: Aqvist 1992
        //  The P1-4 values were taken from O.sp3 with two bonded 
        //  neighbors -> O has the smallest van der Waals radius 
        //  compared to all other elements which had been parametrized
        AssignLCPO(SIptr, 1.18, 0.49392, -0.16038, -0.00015512, 0.00016453);
      else if (atype0=='F')
        AssignLCPO(SIptr, 1.47, 0.68563, -0.1868, -0.00135573, 0.00023743);
      else {
        mprintf("Warning: Using carbon SA parms for unknown atom %i type %s\n",
                atidx + 1, *(atom.Type()));
        AssignLCPO(SIptr, 1.70, 0.51245, -0.15966, -0.00019781, 0.00016392);
      }
  } // END switch atom.Element() 
}
