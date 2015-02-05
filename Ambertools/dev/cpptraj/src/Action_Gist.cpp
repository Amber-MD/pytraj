// Gist 
#include <cmath>
//#include <ctime>
using namespace std;
#include "Action_Gist.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" 
#include "Constants.h" // GASK_KCAL and SMALL

// CONSTRUCTOR
Action_Gist::Action_Gist() :
  CurrentParm_(0),
//  watermodel_(false),
//  useTIP3P_(false),
//  useTIP4P_(false),
//  useTIP4PEW_(false),
//  useTIP5P_(false),
//  useTIP3PFW_(false),
//  useSPCE_(false),
//  useSPCFW_(false),
  doOrder_(false),
  doEij_(false)
{
  BULK_DENS_ = 0;
  gridcntr_[0] = -1;
  gridcntr_[1] = -1;
  gridcntr_[2] = -1;
    
  gridorig_[0] = -1;
  gridorig_[1] = -1;
  gridorig_[2] = -1;
  
  gridspacn_ = 0;
 } 


void Action_Gist::Help() {
//  mprintf("<watermodel>[{tip3p|tip4p|tip4pew}] [doorder] [doeij] [gridcntr <xval> <yval> <zval>] [griddim <xval> <yval> <zval>] [gridspacn <spaceval>] [out <filename>] \n");
  mprintf("\t[doorder] [doeij] [refdens <rdval>] [gridcntr <xval> <yval> <zval>]\n"
          "\t[griddim <xval> <yval> <zval>] [gridspacn <spaceval>]\n"
          "\t[out <filename>]\n");
/*  mprintf("\tGIST needs the specification of the water model being used. Supported water models are: \n");
  mprintf("\ta) TIP3P specified as tip3p. \n");
  mprintf("\tb) TIP4P specified as tip4p. \n");
  mprintf("\tc) TIP4PEW specified as tip4pew. \n");
  mprintf("\td) TIP5P specified as tip5p. \n");
  mprintf("\te) TIP3PFW specified as tip3pfw. \n");
  mprintf("\tf) SPCE specified as spce. \n");
  mprintf("\tg) SPCFW specified as spcfw. \n");
  mprintf("  Calculate GIST between water molecules in selected region \n");
*/}

// Action_Gist::Init()
Action::RetType Action_Gist::Init(ArgList& actionArgs, TopologyList* PFL, DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  if (DSL->EnsembleNum() > -1) {
    mprinterr("Error: GIST currently cannot be used in ensemble mode.\n");
    return Action::ERR;
  }
  gist_init_.Start();
  // Get keywords
  // Dataset to store gist results
  datafile_ = actionArgs.GetStringKey("out");
  // Generate the data set name, and hold onto the master data set list
  /*string ds_name = actionArgs.GetStringKey("name");
  ds_name = myDSL_.GenerateDefaultName("GIST");
  // We have 4?? data sets Add them here
  // Now add all of the data sets
  for (int i = 0; i < 4; i++) {
    myDSL_.AddSetAspect(DataSet::DOUBLE, ds_name,
      integerToString(i+1).c_str());
      }
  //  myDSL_.AddSet(DataSet::DOUBLE, ds_name, NULL);
  */  
/*  useTIP3P_ = actionArgs.hasKey("tip3p");
  useTIP4P_ = actionArgs.hasKey("tip4p");
  useTIP4PEW_ = actionArgs.hasKey("tip4pew");
  useTIP5P_ = actionArgs.hasKey("tip5p");
  useTIP3PFW_ = actionArgs.hasKey("tip3pfw");
  useSPCE_ = actionArgs.hasKey("spce");
  useSPCFW_ = actionArgs.hasKey("spcfw");
  if (!useTIP3P_ && !useTIP4P_ && !useTIP4PEW_ && !useTIP5P_ && !useTIP3PFW_ && !useSPCE_ && !useSPCFW_) {
    mprinterr("Init Error: Only water models supported are TIP3P, TIP4P, TIP4PEW, TIP5P, TIP3P/FW, SPC/E, SPC/FW\n");
    return Action::ERR;
  }
*/
  doOrder_ = actionArgs.hasKey("doorder");
  doEij_ = actionArgs.hasKey("doeij");
  gridspacn_ = actionArgs.getKeyDouble("gridspacn", 0.50);
  // Set Bulk Energy based on water model
/*  if (useTIP3P_) BULK_E_ = -19.0653;
  if (useTIP4PEW_) BULK_E_ = -22.071;
  if (useTIP4P_) BULK_E_ = -19.71152;
  if (useTIP5P_) BULK_E_ = -19.19174;
  if (useTIP3PFW_) BULK_E_ = -22.7374;
  if (useSPCE_) BULK_E_ = -22.2574;
  if (useSPCFW_) BULK_E_ = -23.7458;
//  if (usePOL3_) BULK_E_ = -22.071;
  mprintf("\tGIST bulk energy: %10.5f\n", BULK_E_);
*/  
  // Set Bulk Density 55.5M
  BULK_DENS_ = actionArgs.getKeyDouble("refdens", -1);
  // Grid center
  if ( actionArgs.hasKey("gridcntr") ) {
    gridcntr_[0] = actionArgs.getNextDouble(-1);
    gridcntr_[1] = actionArgs.getNextDouble(-1);
    gridcntr_[2] = actionArgs.getNextDouble(-1);
  } else {
    mprintf("Warning: No grid center values specified, using default\n");
    gridcntr_[0] = 0.0;
    gridcntr_[1] = 0.0;
    gridcntr_[2] = 0.0;
  }
  // Grid dimensions
  griddim_.clear();
  griddim_.resize( 3 );
  if ( actionArgs.hasKey("griddim") ) {
    griddim_[0] = actionArgs.getNextInteger(-1);
    griddim_[1] = actionArgs.getNextInteger(-1);
    griddim_[2] = actionArgs.getNextInteger(-1);
  } else {
    mprintf("Warning: No grid dimension values specified, using default\n");
    griddim_[0] = 40;
    griddim_[1] = 40;
    griddim_[2] = 40;
  }

  InitImaging(true); // Always image

  mprintf("    GIST:\n");
  if(doOrder_)
    mprintf("\tDo Order calculation\n");
  else
    mprintf("\tSkip Order calculation\n");
  if(doEij_)
    mprintf("\tCompute and print water-water Eij matrix\n");
  else
    mprintf("\tSkip water-water Eij matrix\n");
  if (BULK_DENS_ < 0) {
    BULK_DENS_ = 0.0334;
  //BULK_DENS_ = 0.033422885325;
    mprintf("\tNo water reference density specified, using default: %6.4f, equivalent to 1g/cc\n",
            BULK_DENS_);
  } else {
    mprintf("\tWater reference density: %6.4f\n", BULK_DENS_);
    if ( BULK_DENS_ > (0.0334*1.2) )
      mprintf("Warning: water reference density is high, consider using 0.0334 for 1g/cc water density\n");
    else if ( BULK_DENS_ < (0.0334*0.8) )
      mprintf("Warning: water reference density is low, consider using 0.0334 for 1g/cc water density\n");
  }
  mprintf("\tGIST grid center: %5.3f %5.3f %5.3f\n", gridcntr_[0],gridcntr_[1],gridcntr_[2]);
  mprintf("\tGIST grid dimension: %d %d %d \n", griddim_[0],griddim_[1],griddim_[2]);
  mprintf("\tGIST grid spacing: %5.3f A^3\n", gridspacn_);
  mprintf("\t#Please cite these papers if you use GIST results in a publication:\n"
          "\t#    Crystal Nguyen, Michael K. Gilson, and Tom Young, arXiv:1108.4876v1 (2011)\n"
          "\t#    Crystal N. Nguyen, Tom Kurtzman Young, and Michael K. Gilson,\n"
          "\t#      J. Chem. Phys. 137, 044101 (2012)\n"
          "\t#    Lazaridis, J. Phys. Chem. B 102, 3531–3541 (1998)\n");
  gist_init_.Stop();
  return Action::OK;
}

// Action_Gist::Setup()
/** Set Gist up for this parmtop. Get masks etc.
  */
Action::RetType Action_Gist::Setup(Topology* currentParm, Topology** parmAddress) {
  gist_setup_.Start();
  CurrentParm_ = currentParm;      
  NFRAME_ = 0;
  max_nwat_ = 0;

  MAX_GRID_PT_ = griddim_[0] * griddim_[1] * griddim_[2];
  Vvox_ = gridspacn_*gridspacn_*gridspacn_;
  G_max_x_ = griddim_[0] * gridspacn_ + 1.5 ;
  G_max_y_ = griddim_[1] * gridspacn_ + 1.5 ;
  G_max_z_ = griddim_[2] * gridspacn_ + 1.5 ;
  
  //mprintf("\tGIST Setup: %d %d %d %d %f \n", griddim_[0], griddim_[1], 
  //        griddim_[2], MAX_GRID_PT_, Vvox_);
  mprintf("\tGIST number of voxels: %d, voxel volume: %f A^3\n",  MAX_GRID_PT_, Vvox_);

  // Set up grid origin
  gridorig_[0] = gridcntr_[0] - 0.5*griddim_[0]*gridspacn_;
  gridorig_[1] = gridcntr_[1] - 0.5*griddim_[1]*gridspacn_;
  gridorig_[2] = gridcntr_[2] - 0.5*griddim_[2]*gridspacn_;
  mprintf("\tGIST grid origin: %5.3f %5.3f %5.3f\n", 
          gridorig_[0], gridorig_[1], gridorig_[2]);

  // Set up cumulative energy arrays
  /*  x_.clear();
  x_.resize(5, 0.0);
  y_.clear();
  y_.resize(5, 0.0);
  z_.clear();
  z_.resize(5, 0.0);*/
  wh_evdw_.clear();
  wh_evdw_.resize(MAX_GRID_PT_, 0.0);
  wh_eelec_.clear();
  wh_eelec_.resize(MAX_GRID_PT_, 0.0);
  ww_evdw_.clear();
  ww_evdw_.resize(MAX_GRID_PT_, 0.0);
  ww_eelec_.clear();
  ww_eelec_.resize(MAX_GRID_PT_, 0.0);

  //voxel coords
  grid_x_.clear();    
  grid_x_.resize(MAX_GRID_PT_, 0.0); 
  grid_y_.clear();          
  grid_y_.resize(MAX_GRID_PT_, 0.0);
  grid_z_.clear();          
  grid_z_.resize(MAX_GRID_PT_, 0.0); 


  // get the actual voxel coordinates
  voxel_ = 0;
  for (int i = 0; i < griddim_[0]; ++i) {
    for (int j = 0; j < griddim_[1]; ++j) {
      for (int k = 0; k < griddim_[2]; ++k) {
        grid_x_[voxel_] = Xcrd(i);
        grid_y_[voxel_] = Ycrd(j);
        grid_z_[voxel_] = Zcrd(k);
        voxel_++;
      }
    }
  }

  Esw_dens_.clear();
  Esw_dens_.resize(MAX_GRID_PT_, 0.0);
  Esw_norm_.clear();
  Esw_norm_.resize(MAX_GRID_PT_, 0.0);
  Eww_norm_.clear();
  Eww_norm_.resize(MAX_GRID_PT_, 0.0);
  Eww_dens_.clear();
  Eww_dens_.resize(MAX_GRID_PT_, 0.0);

  if(doEij_) {
    ww_Eij_.clear();
    ww_Eij_.resize(MAX_GRID_PT_);
    for(int i = 1; i < MAX_GRID_PT_; i++) ww_Eij_[i].resize(i);
    
    //CN: need to initialize ww_Eij_ to 0.0 but not Euler angles
    for (int a=1; a<MAX_GRID_PT_; a++)
      for (int l=0; l<a; l++) ww_Eij_[a][l]=0.0;  
  }
  the_vox_.clear();
  the_vox_.resize(MAX_GRID_PT_);
  phi_vox_.clear();
  phi_vox_.resize(MAX_GRID_PT_);
  psi_vox_.clear();
  psi_vox_.resize(MAX_GRID_PT_);

  dTStrans_dens_.clear();
  dTStrans_dens_.resize(MAX_GRID_PT_, 0.0);
  dTStrans_norm_.clear();
  dTStrans_norm_.resize(MAX_GRID_PT_, 0.0); 
  dTSorient_dens_.clear();
  dTSorient_dens_.resize(MAX_GRID_PT_, 0.0);
  dTSorient_norm_.clear();
  dTSorient_norm_.resize(MAX_GRID_PT_, 0.0);

  nwat_.clear();
  nwat_.resize(MAX_GRID_PT_, 0);
  nH_.clear();
  nH_.resize(MAX_GRID_PT_, 0);
  nw_angle_.clear();
  nw_angle_.resize(MAX_GRID_PT_, 0);
  dens_.clear();
  dens_.resize(MAX_GRID_PT_, 0.0);
  g_.clear();
  g_.resize(MAX_GRID_PT_, 0.0);
  gH_.clear();
  gH_.resize(MAX_GRID_PT_, 0.0);
  dipolex_.clear();
  dipolex_.resize(MAX_GRID_PT_, 0.0);
  dipoley_.clear();
  dipoley_.resize(MAX_GRID_PT_, 0.0);
  dipolez_.clear();
  dipolez_.resize(MAX_GRID_PT_, 0.0);
  neighbor_.clear();
  neighbor_.resize(MAX_GRID_PT_, 0.0);
  neighbor_dens_.clear();
  neighbor_dens_.resize(MAX_GRID_PT_, 0.0);
  neighbor_norm_.clear();
  neighbor_norm_.resize(MAX_GRID_PT_, 0.0);
  qtet_.clear();
  qtet_.resize(MAX_GRID_PT_, 0.0);
  pol_.clear();
  pol_.resize(MAX_GRID_PT_, 0.0);

  gridwat_.clear();
  gridwat_.resize( currentParm->Nsolvent() );

  // We need box info
  if (currentParm->BoxType() == Box::NOBOX) {
    mprinterr("Error: Gist: Must have explicit solvent with periodic boundaries!");
    return Action::ERR;
  }
  SetupImaging( currentParm->BoxType() );

  resnum_ = 0;
  voxel_ = 0;
  gist_setup_.Stop();
  return Action::OK;  
}

// Action_Gist::DoAction()
Action::RetType Action_Gist::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {
  NFRAME_ ++;
//  if (NFRAME_==1) mprintf("GIST Action \n");

  // Simulation box length - assign here because it can vary for npt simulation
  //Lx_ = currentFrame->BoxCrd().BoxX();
  //Ly_ = currentFrame->BoxCrd().BoxY();
  //Lz_ = currentFrame->BoxCrd().BoxZ();
//  if (NFRAME_==1) mprintf("GIST Action box length: %f %f %f \n", Lx_, Ly_, Lz_);
  
  int solventMolecules = CurrentParm_->Nsolvent();
  resnum_ = 0;
  voxel_ = 0;
  resindex1_ = 0;
  for (solvmol_ = CurrentParm_->MolStart();
       solvmol_ != CurrentParm_->MolEnd(); ++solvmol_)
  {
    resindex1_++;
    if (!solvmol_->IsSolvent()) continue;
    gist_grid_.Start();
    Grid( currentFrame );
    gist_grid_.Stop();
    voxel_ = gridwat_[resnum_];
    resnum_++;
    gist_nonbond_.Start();
    NonbondEnergy( currentFrame );
    gist_nonbond_.Stop();
    if (voxel_ >= MAX_GRID_PT_) continue;
    gist_euler_.Start();
    EulerAngle( currentFrame );
    gist_euler_.Stop();
    gist_dipole_.Start();
    Dipole( currentFrame );
    gist_dipole_.Stop();
  }
  if(doOrder_) Order( currentFrame );
  
  //Debug
//  if (NFRAME_==1) mprintf("GIST  DoAction:  Found %d solvent residues \n", resnum_);
  if (solventMolecules != resnum_) {
    mprinterr("GIST  DoAction Error: Number of solvent molecules don't match %d %d\n", solventMolecules, resnum_);
  }
  
  return Action::OK;
}

// Action_Gist::NonbondEnergy()
void Action_Gist::NonbondEnergy(Frame *currentFrame) {
  double rij2, rij, r2, r6, r12, f12, f6, e_vdw, e_elec;
  int satom, satom2, atom1, atom2;
  
  int  voxel2 = 0;
  double q1, q2;
  
  // Setup imaging info
  Matrix_3x3 ucell, recip;
  if (ImagingEnabled())
    currentFrame->BoxCrd().ToRecip(ucell, recip);

  // Inner loop has both solute and solvent
  resnum2_=0;
  resindex2_ = 1;
  // skip if water2 has index larger than water1 so that every pair is only evaluated once
  solvmol2_ = CurrentParm_->MolStart();
  for (resindex2_=1; resindex2_<resindex1_; resindex2_++)
  {    
    if (!solvmol2_->IsSolvent()) {
      // Outer loop is not water, break inner loop if water 1 is outside the grid
      if (voxel_ >= MAX_GRID_PT_) {
        ++solvmol2_;
        continue;
      }
    } else { 
      // Inner loop is water
      voxel2 = gridwat_[resnum2_];
      resnum2_++;
      // skip if both waters are outside the grid
      if (voxel_>=MAX_GRID_PT_ && voxel2>=MAX_GRID_PT_) {
        ++solvmol2_;
        continue;
      }
    }
      
    // Loop over all solvent atoms of water 1
    atom1=0;
    for (satom = solvmol_->BeginAtom(); satom < solvmol_->EndAtom(); ++satom)
    {
      // Set up coord index for this atom
      const double* XYZ =  currentFrame->XYZ( satom );
      atom2=0;
      for (satom2 = solvmol2_->BeginAtom(); satom2 < solvmol2_->EndAtom(); ++satom2)
      {    
        // Set up coord index for this atom
        const double* XYZ2 = currentFrame->XYZ( satom2 );
        // Calculate the vector pointing from atom2 to atom1
        rij2 = DIST2(XYZ, XYZ2, ImageType(), currentFrame->BoxCrd(), ucell, recip);
        rij = sqrt(rij2);
        // LJ energy
        NonbondType const& LJ = CurrentParm_->GetLJparam(satom, satom2); 
        r2    = 1 / rij2;
        r6    = r2 * r2 * r2;
        r12   = r6 * r6;
        f12   = LJ.A() * r12;  // A/r^12
        f6    = LJ.B() * r6;   // B/r^6
        e_vdw = f12 - f6;     // (A/r^12)-(B/r^6)
        // LJ Force 
        // Coulomb energy 
        q1 = (*CurrentParm_)[satom].Charge() * Constants::ELECTOAMBER;
        q2 = (*CurrentParm_)[satom2].Charge() * Constants::ELECTOAMBER;
        e_elec = (q1*q2/rij);
        if (!solvmol2_->IsSolvent()) {
          // solute-solvent interaction
          wh_evdw_[voxel_] +=  e_vdw;
          wh_eelec_[voxel_] += e_elec;
        } else {
          // solvent-solvent interaction, need to compute for all waters,
          // even those outside the grid but only one water needs to be 
          // inside the grid. 
          if (voxel_<MAX_GRID_PT_) {
            ww_evdw_[voxel_] +=  e_vdw;
            ww_eelec_[voxel_] += e_elec;
            // Store the water neighbor using only O-O distance
            if (atom2==0 && atom1==0 && rij<3.5)
              neighbor_[voxel_] += 1.0;
          }
          // CN: only store Eij[voxel1][voxel2] if both voxels lie on the grid.
          if (voxel2<MAX_GRID_PT_) {
            ww_evdw_[voxel2] +=  e_vdw;
            ww_eelec_[voxel2] += e_elec;
            // Store the water neighbor using only O-O distance
            if (atom2==0 && atom1==0 && rij<3.5)
              neighbor_[voxel2] += 1.0;
            if (doEij_ && (voxel_<MAX_GRID_PT_)) {
              if (voxel_>voxel2) {
                ww_Eij_[voxel_][voxel2] += e_vdw;
                ww_Eij_[voxel_][voxel2] += e_elec;
              } else {
                ww_Eij_[voxel2][voxel_] += e_vdw;
                ww_Eij_[voxel2][voxel_] += e_elec;
              }
            } //print Eij && voxel<MAX_GRID_PT_
          }
        }//IF is solvent
        atom2++;
      } // END Inner loop ALL atoms
      atom1++;
    } // END Outer loop solvent atoms
    ++solvmol2_;
  }  // END Inner loop ALL molecules
}

// Action_Gist::Grid()
void Action_Gist::Grid(Frame *frameIn) {
  int  i, gridindex[3], nH;
  Vec3 comp,  atom_coord;
  i = solvmol_->BeginAtom();

  gridwat_[resnum_] = MAX_GRID_PT_ + 1;
  atom_coord = Vec3(frameIn->XYZ(i));
  // get the components of the water vector
  comp = Vec3(atom_coord) - Vec3(gridorig_);
  nH=0;
  //If Oxygen is far from grid, 1.5A or more in any durection, skip calculation
  if (comp[0] <= G_max_x_ && comp[1] <= G_max_y_ && comp[2] <= G_max_z_ && 
      comp[0] >= -1.5     && comp[1] >= -1.5     && comp[2] >= -1.5 )
  {
    //if (comp[0]<= G_max_x || comp[1]<= G_max_y || comp[2]<= G_max_z ||
    //    comp[0]>= -1.5 || comp[1]>= -1.5 || comp[2]>= -1.5 ) {
    //Water is at most 1.5A away from grid, so we need to check for H even if O is outside grid
    nH=2;
    
    //O is inside grid only if comp is >=0
    if (comp[0]>=0 && comp[1]>=0 && comp[2]>=0 ){
      comp /= gridspacn_;
      gridindex[0] = (int) comp[0];
      gridindex[1] = (int) comp[1];
      gridindex[2] = (int) comp[2];
      
      if ((gridindex[0]<griddim_[0]) && (gridindex[1]<griddim_[1]) && (gridindex[2]<griddim_[2]))
      {
        // this water belongs to grid point gridindex[0], gridindex[1], gridindex[2]
        voxel_ = (gridindex[0]*griddim_[1] + gridindex[1])*griddim_[2] + gridindex[2];
        gridwat_[resnum_] = voxel_;
        nwat_[voxel_]++;
        if (max_nwat_ < nwat_[voxel_]) max_nwat_ = nwat_[voxel_];
      }
    }
    
    // evaluate hydrogen atoms
    for (int a=1; a<=nH; a++) {
      atom_coord = Vec3(frameIn->XYZ(i+a));
      comp = Vec3(atom_coord) - Vec3(gridorig_);
      if (comp[0]<0 || comp[1]<0 || comp[2]<0) continue;
      comp /= gridspacn_;
      gridindex[0] = (int) comp[0];
      gridindex[1] = (int) comp[1];
      gridindex[2] = (int) comp[2];
      if ((gridindex[0]<griddim_[0]) && 
          (gridindex[1]<griddim_[1]) && 
          (gridindex[2]<griddim_[2]))
      {
        voxel_ = (gridindex[0]*griddim_[1] + gridindex[1])*griddim_[2] + gridindex[2];
        nH_[voxel_]++;
      }
    } 
  }
}

// Action_Gist::EulerAngle()
void Action_Gist::EulerAngle(Frame *frameIn) {
  //if (NFRAME_==1) mprintf("GIST Euler Angles \n");
  Vec3 x_lab, y_lab, z_lab, O_wat, H1_wat, H2_wat, x_wat, y_wat, z_wat, node, v;
  double dp;

  int i = solvmol_->BeginAtom();
  O_wat = Vec3(frameIn->XYZ(i));
  H1_wat = Vec3(frameIn->XYZ(i+1)) - O_wat;
  H2_wat = Vec3(frameIn->XYZ(i+2)) - O_wat;
  
  // make sure the first three atoms are oxygen followed by two hydrogen
  if ((*CurrentParm_)[i].Element() != Atom::OXYGEN) {
    mprintf("Warning: GIST: First coordinates do not belong to oxygen atom (%s)\n",
            (*CurrentParm_)[i].ElementName());
  }
  if ((*CurrentParm_)[i+1].Element() != Atom::HYDROGEN || 
      (*CurrentParm_)[i+2].Element() != Atom::HYDROGEN)
  {
    mprintf("Warning: GIST: second and third coordinates do not belong to hydrogen atoms (%s, %s)\n",
            (*CurrentParm_)[i+1].ElementName(), (*CurrentParm_)[i+2].ElementName());
  } 
  
  // Define lab frame of reference
  x_lab[0]=1.0; x_lab[1]=0;   x_lab[2]=0;
  y_lab[0]=0;   y_lab[1]=1.0; y_lab[2]=0;
  z_lab[0]=0;   z_lab[1]=0;   z_lab[2]=1.0;     
  
  // Define the water frame of reference - all axes must be normalized
  // make h1 the water x-axis (but first need to normalized)
  x_wat = H1_wat;
  x_wat.Normalize();
  // the normalized z-axis is the cross product of h1 and h2 
  z_wat = x_wat.Cross( H2_wat );
  z_wat.Normalize();
  // make y-axis as the cross product of h1 and z-axis
  y_wat = z_wat.Cross( x_wat );
  y_wat.Normalize();
  
  // Find the X-convention Z-X'-Z'' Euler angles between the water frame and the lab/host frame
  // First, theta = angle between the water z-axis of the two frames
  dp = z_lab*( z_wat);
  theta_ = acos(dp);
  //  if (theta>0 && theta<PI) {
  if (theta_>1E-5 && theta_<Constants::PI-1E-5) {
    // phi = angle between the projection of the water x-axis and the node
    // line of node is where the two xy planes meet = must be perpendicular to both z axes
    // direction of the lines of node = cross product of two normals (z axes)
    // acos of x always gives the angle between 0 and pi, which is okay for theta since theta ranges from 0 to pi
    node = z_lab.Cross( z_wat );
    node.Normalize();
    
    // Second, find the angle phi, which is between x_lab and the node
    dp = node*( x_lab );
    if (dp <= -1.0)     phi_ = Constants::PI;
    else if (dp >= 1.0) phi_ = Constants::PI;
    else                phi_ = acos(dp);
    // check angle phi
    if (phi_>0 && phi_<(Constants::TWOPI)) {
      // method 2
      v = x_lab.Cross( node );
      dp = v*( z_lab );
      if (dp<0) phi_ = Constants::TWOPI - phi_;
    }
    
    // Third, rotate the node to x_wat about the z_wat axis by an angle psi
    // psi = angle between x_wat and the node 
    dp = x_wat*( node );
    if (dp<=-1.0)     psi_ = Constants::PI;
    else if (dp>=1.0) psi_ = 0;
    else              psi_ = acos(dp);
    // check angle psi
    if (psi_>0 && psi_<(Constants::TWOPI)) {
      // method 2
      Vec3 v = node.Cross( x_wat );
      dp = v*( z_wat );
      if (dp<0) psi_ = Constants::TWOPI - psi_;
    }
    
    if (!(theta_<=Constants::PI && theta_>=0 && 
          phi_<=Constants::TWOPI && phi_>=0 && psi_<=Constants::TWOPI && psi_>=0))
    {
//      mprintf("GIST: angles: %f %f %f\n", theta_, phi_, psi_);
//      H1_wat.Print("H1_wat");
//      H2_wat.Print("H2_wat");
      mprinterr("Error: Euler: angles don't fall into range.\n");
      //break; 
    }
    
    the_vox_[voxel_].push_back(theta_);
    phi_vox_[voxel_].push_back(phi_);
    psi_vox_[voxel_].push_back(psi_);
    nw_angle_[voxel_]++;
  }
  //else mprintf("%i: gimbal lock problem, two z_wat paralell\n", resnum-1);
} 

// Action_Gist::Dipole()
void Action_Gist::Dipole(Frame *frameIn) {
  
  //if (NFRAME_==1) mprintf("GIST Dipole \n");
  double dipolar_vector[3], charge;
  int satom;

  dipolar_vector[0] = 0.0;
  dipolar_vector[1] = 0.0;
  dipolar_vector[2] = 0.0;
  // Loop over solvent atoms
  for (satom = solvmol_->BeginAtom(); satom < solvmol_->EndAtom(); ++satom)
  {
    const double* XYZ = frameIn->XYZ( satom );
    // Calculate dipole vector. The oxygen of the solvent is used to 
    // assign the voxel index to the water.
    // NOTE: the total charge on the solvent should be neutral for this 
    //       to have any meaning
    charge = (*CurrentParm_)[satom].Charge();
    //      mprintf("%i %f %f %f %f\n", resnum-1, charge, XYZ[0], XYZ[1], XYZ[2]);
    dipolar_vector[0] += (charge * XYZ[0]);
    dipolar_vector[1] += (charge * XYZ[1]);
    dipolar_vector[2] += (charge * XYZ[2]);
  }
  dipolex_[voxel_] += dipolar_vector[0];
  dipoley_[voxel_] += dipolar_vector[1];
  dipolez_[voxel_] += dipolar_vector[2];
}

// Action_Gist::Order() 
void Action_Gist::Order(Frame *frameIn) {
//  if (NFRAME_==1) mprintf("GIST Order Parameter \n");
  int i;
  double cos, sum, r1, r2, r3, r4, rij2, x[5], y[5], z[5];
  Vec3 neighbor1(0.0), neighbor2(0.0), neighbor3(0.0), neighbor4(0.0);
  Vec3 O_wat1, O_wat2, O_wat3, v1, v2;
  resnum_=0;

  for (solvmol_ = CurrentParm_->MolStart();
       solvmol_ != CurrentParm_->MolEnd(); ++solvmol_)
  {
    if (!solvmol_->IsSolvent()) continue;

    // obtain 4 closest neighbors for every water
    resnum_++;
    voxel_ = gridwat_[resnum_-1];
    if (voxel_>=MAX_GRID_PT_) continue;
    // assume that oxygen is the first atom
    i = solvmol_->BeginAtom();
    O_wat1 = Vec3(frameIn->XYZ( i ));

    r1=1000; r2=1000; r3=1000; r4=1000; resnum2_=0;
    // Can't make into triangular matrix
    for (solvmol2_ = CurrentParm_->MolStart();
         solvmol2_ != CurrentParm_->MolEnd(); ++solvmol2_)
    {
      if (!solvmol2_->IsSolvent()) continue;
      resnum2_++;
      if (resnum_ == resnum2_) continue;
      i = solvmol2_->BeginAtom();
      O_wat2 = Vec3(frameIn->XYZ( i ));      
      rij2 = DIST2_NoImage(O_wat1, O_wat2);
      if (rij2<r1) {
        r4 = r3;
        r3 = r2;
        r2 = r1;
        r1 = rij2;
        neighbor4 = neighbor3;
        neighbor3 = neighbor2;
        neighbor2 = neighbor1;
        neighbor1 = O_wat2;
      }
      else if (rij2<r2) {
        r4 = r3;
        r3 = r2;
        r2 = rij2;
        neighbor4 = neighbor3;
        neighbor3 = neighbor2;
        neighbor2 = O_wat2;
      }
      else if (rij2<r3) {
        r4 = r3;
        r3 = rij2;
        neighbor4 = neighbor3;
        neighbor3 = O_wat2;       
      }
      else if (rij2<r4) {
        r4 = rij2;
        neighbor4 = O_wat2;
     }        
    }
    x[1]=neighbor1[0]; y[1]=neighbor1[1]; z[1]=neighbor1[2];
    x[2]=neighbor2[0]; y[2]=neighbor2[1]; z[2]=neighbor2[2];
    x[3]=neighbor3[0]; y[3]=neighbor3[1]; z[3]=neighbor3[2];
    x[4]=neighbor4[0]; y[4]=neighbor4[1]; z[4]=neighbor4[2];    
    // Compute the tetrahedral order parameter
    sum=0;
    for (int mol1=1; mol1<=3; mol1++) {
      for (int mol2=mol1+1; mol2<=4; mol2++) {
        O_wat2[0] = x[mol1];
        O_wat2[1] = y[mol1];
        O_wat2[2] = z[mol1];
        O_wat3[0] = x[mol2];
        O_wat3[1] = y[mol2];
        O_wat3[2] = z[mol2];
        v1 = O_wat2 - O_wat1;
        v2 = O_wat3 - O_wat1;    
        r1 = v1.Magnitude2();
        r2 = v2.Magnitude2();
        cos = (v1*( v2))/sqrt(r1*r2);
        sum += (cos + 1.0/3)*(cos + 1.0/3);
      }
    }
    qtet_[voxel_] += (1.0 - (3.0/8)*sum);
/*    double dbl = (1.0 - (3.0/8)*sum)
 *    if (dbl<-3.0 || dbl>1.0) {
      mprintf("BAD! voxel=%d, q=%9.5f\n", voxel_, dbl);
    }*/
  }
}

void Action_Gist::Print() {
  gist_print_.Start();
  // Implement NN to compute orientational entropy for each voxel
  double NNr, rx, ry, rz, rR, dbl;
  dTSorienttot_=0;
  for (int gr_pt=0; gr_pt<MAX_GRID_PT_; gr_pt++) {
    dTSorient_dens_[gr_pt]=0; dTSorient_norm_[gr_pt]=0;  
    int nwtot = nw_angle_[gr_pt];
    if (nwtot<=1) continue;
    for (int n=0; n<nwtot; n++) {
      NNr=10000;
      for (int l=0; l<nwtot; l++) {
        if (l==n) continue;
        rx = cos(the_vox_[gr_pt][l]) - cos(the_vox_[gr_pt][n]);
        ry = phi_vox_[gr_pt][l] - phi_vox_[gr_pt][n];
        rz = psi_vox_[gr_pt][l] - psi_vox_[gr_pt][n];
        if      (ry>Constants::PI) ry = Constants::TWOPI-ry;
        else if (ry<-Constants::PI) ry = Constants::TWOPI+ry;
        if      (rz>Constants::PI) rz = Constants::TWOPI-rz;
        else if (rz<-Constants::PI) rz = Constants::TWOPI+rz;
        rR = sqrt(rx*rx + ry*ry + rz*rz);
        if (rR>0 && rR<NNr) NNr = rR;
      }
      if (NNr<9999 && NNr>0) {
        dbl = log(NNr*NNr*NNr*nwtot/(3.0*Constants::TWOPI));
        dTSorient_norm_[gr_pt] += dbl;
      }  
    }
    dTSorient_norm_[gr_pt] = Constants::GASK_KCAL*300*(dTSorient_norm_[gr_pt]/nwtot+0.5772);
    dTSorient_dens_[gr_pt] = dTSorient_norm_[gr_pt]*nwat_[gr_pt]/(NFRAME_*Vvox_);
    dTSorienttot_ += dTSorient_dens_[gr_pt];
  }
  dTSorienttot_ *= Vvox_;
  mprintf("Maximum number of waters found in one voxel for %d frames = %d\n", NFRAME_, max_nwat_);
  mprintf("Total referenced orientational entropy of the grid: dTSorient = %9.5f kcal/mol, Nf=%d\n",
          dTSorienttot_, NFRAME_);
  
  // Compute translational entropy for each voxel
  dTStranstot_=0;
  for (int a=0; a<MAX_GRID_PT_; a++) {
    dens_[a] = 1.0*nwat_[a]/(NFRAME_*Vvox_);
    g_[a] = dens_[a]/BULK_DENS_;
    gH_[a] = 1.0*nH_[a]/(NFRAME_*Vvox_*2*BULK_DENS_);
    if (nwat_[a]>1) {
       dTStrans_dens_[a] = -Constants::GASK_KCAL*BULK_DENS_*300*g_[a]*log(g_[a]);
       dTStrans_norm_[a] = dTStrans_dens_[a]/dens_[a];
       dTStranstot_ += dTStrans_dens_[a];
    } else {
       dTStrans_dens_[a]=0; dTStrans_norm_[a]=0;
    }
  }
  dTStranstot_ *= Vvox_;
  mprintf("Total referenced translational entropy of the grid: dTStrans = %9.5f kcal/mol, Nf=%d\n",
          dTStranstot_, NFRAME_);

  // Compute average voxel energy
  double Eswtot = 0.0;
  double Ewwtot = 0.0;
  for (int a=0; a<MAX_GRID_PT_; a++) {
    if (nwat_[a]>1) {
       Esw_dens_[a] = (wh_evdw_[a]+wh_eelec_[a])/(NFRAME_*Vvox_);
       Esw_norm_[a] = (wh_evdw_[a]+wh_eelec_[a])/nwat_[a];
       Eww_dens_[a] = (ww_evdw_[a]+ww_eelec_[a])/(2*NFRAME_*Vvox_);
       Eww_norm_[a] = (ww_evdw_[a]+ww_eelec_[a])/(2*nwat_[a]); 
       Eswtot += Esw_dens_[a];
       Ewwtot += Eww_dens_[a];
    } else {
       Esw_dens_[a]=0; Esw_norm_[a]=0; Eww_norm_[a]=0; Eww_dens_[a]=0;
    }
    // Compute the average number of water neighbor, average order parameter, and average dipole density 
    if (nwat_[a]>0) {
      qtet_[a] /= nwat_[a];
      neighbor_norm_[a] = 1.0*neighbor_[a]/nwat_[a];
    }
    neighbor_dens_[a] = 1.0*neighbor_[a]/(NFRAME_*Vvox_);
    dipolex_[a] /= (0.20822678*NFRAME_*Vvox_);
    dipoley_[a] /= (0.20822678*NFRAME_*Vvox_);
    dipolez_[a] /= (0.20822678*NFRAME_*Vvox_);
    pol_[a] = sqrt(dipolex_[a]*dipolex_[a] + dipoley_[a]*dipoley_[a] + dipolez_[a]*dipolez_[a]);
  }
  Eswtot *= Vvox_;
  Ewwtot *= Vvox_;
  mprintf("Total water-solute energy of the grid: Esw = %9.5f kcal/mol\n", Eswtot);
  mprintf("Total unreferenced water-water energy of the grid: Eww = %9.5f kcal/mol\n", Ewwtot);

  // Print the gist info file
  // Print the energy data
  /*if (!datafile_.empty()) {
  // Now write the data file with all of the GIST energies
  DataFile dfl;
  ArgList dummy;
  dfl.SetupDatafile(datafile_, dummy, 0);
  for (int i = 0; i < myDSL_.size(); i++) {
  dfl.AddSet(myDSL_[i]);
  }
  
  dfl.Write();
  
  }*/
  //stored as float
  PrintDX("gist-gO.dx", g_);
  PrintDX("gist-gH.dx", gH_);
  PrintDX("gist-Esw-dens.dx", Esw_dens_);
  PrintDX("gist-Eww-dens.dx", Eww_dens_);
  PrintDX("gist-dTStrans-dens.dx", dTStrans_dens_);
  PrintDX("gist-dTSorient-dens.dx", dTSorient_dens_);
  PrintDX("gist-neighbor-norm.dx", neighbor_norm_); 
  PrintDX("gist-dipole-dens.dx", pol_);
  //stored as doubles
  PrintDX_double("gist-order-norm.dx", qtet_);
  PrintDX_double("gist-dipolex-dens.dx", dipolex_);
  PrintDX_double("gist-dipoley-dens.dx", dipoley_);
  PrintDX_double("gist-dipolez-dens.dx", dipolez_);

  if (!datafile_.empty())
    PrintOutput(datafile_);
  else
    PrintOutput("gist-output.dat");
  gist_print_.Stop();
  double total = gist_grid_.Total() + gist_nonbond_.Total() + 
                 gist_euler_.Total() + gist_dipole_.Total() +
                 gist_init_.Total() + gist_setup_.Total() + 
                 gist_print_.Total();
  mprintf("\tGIST timings:\n");
  gist_init_.WriteTiming(1,    "Init: ", total);
  gist_setup_.WriteTiming(1,   "Setup:", total);
  gist_grid_.WriteTiming(2,    "Grid:   ", total);
  gist_nonbond_.WriteTiming(2, "Nonbond:", total);
  gist_euler_.WriteTiming(2,   "Euler:  ", total);
  gist_dipole_.WriteTiming(2,  "Dipole: ", total);
  gist_print_.WriteTiming(1,   "Print:", total);
  mprintf("TIME:\tTotal: %.4f s\n", total);
}

// Print GIST data in dx format
void Action_Gist::PrintDX_double(string const& filename, std::vector<double>& data)
{
  CpptrajFile outfile;
  if (outfile.OpenWrite(filename)) {
    mprinterr("Print Error: Could not open OpenDX output file.\n");
    return;
  }
  // Print the OpenDX header
  outfile.Printf("object 1 class gridpositions counts %d %d %d\n",
                 griddim_[0], griddim_[1], griddim_[2]);
  outfile.Printf("origin %lg %lg %lg\n", gridorig_[0], gridorig_[1], gridorig_[2]);
  outfile.Printf("delta %lg 0 0\n", gridspacn_);
  outfile.Printf("delta 0 %lg 0\n", gridspacn_);
  outfile.Printf("delta 0 0 %lg\n", gridspacn_);
  outfile.Printf("object 2 class gridconnections counts %d %d %d\n",
                 griddim_[0], griddim_[1], griddim_[2]);
  outfile.Printf(
    "object 3 class array type double rank 0 items %d data follows\n",
    MAX_GRID_PT_);

  // Now print out the data. It is already in row-major form (z-axis changes
  // fastest), so no need to do any kind of data adjustment
  for (int i = 0; i < MAX_GRID_PT_ - 2; i += 3)
    outfile.Printf("%g %g %g\n", data[i], data[i+1], data[i+2]);
  // Print out any points we may have missed
  switch (MAX_GRID_PT_ % 3) {
    case 2: outfile.Printf("%g %g\n", data[MAX_GRID_PT_-2], data[MAX_GRID_PT_-1]); break;
    case 1: outfile.Printf("%g\n", data[MAX_GRID_PT_-1]); break;
  }

  outfile.CloseFile();
}

// Print GIST data in dx format
void Action_Gist::PrintDX(string const& filename, std::vector<float>& data)
{
  CpptrajFile outfile;
  if (outfile.OpenWrite(filename)) {
    mprinterr("Print Error: Could not open OpenDX output file.\n");
    return;
  }
  // Print the OpenDX header
  outfile.Printf("object 1 class gridpositions counts %d %d %d\n",
                 griddim_[0], griddim_[1], griddim_[2]);
  outfile.Printf("origin %lg %lg %lg\n", gridorig_[0], gridorig_[1], gridorig_[2]);
  outfile.Printf("delta %lg 0 0\n", gridspacn_);
  outfile.Printf("delta 0 %lg 0\n", gridspacn_);
  outfile.Printf("delta 0 0 %lg\n", gridspacn_);
  outfile.Printf("object 2 class gridconnections counts %d %d %d\n",
                 griddim_[0], griddim_[1], griddim_[2]);
  outfile.Printf(
    "object 3 class array type float rank 0 items %d data follows\n",
    MAX_GRID_PT_);

  // Now print out the data. It is already in row-major form (z-axis changes
  // fastest), so no need to do any kind of data adjustment
  for (int i = 0; i < MAX_GRID_PT_ - 2; i += 3)
    outfile.Printf("%g %g %g\n", data[i], data[i+1], data[i+2]);
  // Print out any points we may have missed
  switch (MAX_GRID_PT_ % 3) {
    case 2: outfile.Printf("%g %g\n", data[MAX_GRID_PT_-2], data[MAX_GRID_PT_-1]); break;
    case 1: outfile.Printf("%g\n", data[MAX_GRID_PT_-1]); break;
  }

  outfile.CloseFile();
}

// Print GIST data in dx format
void Action_Gist::PrintOutput(string const& filename)
{
  CpptrajFile outfile;
  if (outfile.OpenWrite(filename)) {
    mprinterr("Print Error: Could not open GISToutput file.\n");
    return;
  }
  //ofstream myfile;
  //  myfile.open("gist-output.dat");  
  //myfile.open(filename);  
  outfile.Printf("GIST Output, information printed per voxel\n");
  outfile.Printf("voxel xcoord ycoord zcoord population g_O g_H ");
  outfile.Printf("dTStrans-dens(kcal/mol/A^3) dTStrans-norm(kcal/mol) dTSorient-dens(kcal/mol/A^3) dTSorient-norm(kcal/mol) ");
  outfile.Printf("Esw-dens(kcal/mol/A^3) Esw-norm(kcal/mol) ");
  outfile.Printf("Eww-dens(kcal/mol/A^3) Eww-norm-unref(kcal/mol) ");
  outfile.Printf("Dipole_x-dens(D/A^3) Dipole_y-dens(D/A^3) Dipole_z-dens(D/A^3) Dipole-dens(D/A^3) neighbor-dens(1/A^3) neighbor-norm order-norm\n");
  // Now print out the data. 
  for (int i=0; i<MAX_GRID_PT_; i++){
    outfile.Printf( "%d %g %g %g %d %g %g ",i , grid_x_[i] , grid_y_[i], grid_z_[i], nwat_[i] , g_[i], gH_[i] );
    outfile.Printf( "%g %g %g %g ",dTStrans_dens_[i], dTStrans_norm_[i], dTSorient_dens_[i] , dTSorient_norm_[i]);
    outfile.Printf( "%g %g ",Esw_dens_[i], Esw_norm_[i] );
    outfile.Printf( "%g %g ",Eww_dens_[i] , Eww_norm_[i] );
    outfile.Printf( "%g %g %g %g ",dipolex_[i] , dipoley_[i] , dipolez_[i] , pol_[i] );
    outfile.Printf( "%g %g %g \n",neighbor_dens_[i] , neighbor_norm_[i] , qtet_[i]);
    }
  outfile.CloseFile();

  //double Eijtot=0;
  if(doEij_) {
    if (outfile.OpenWrite("Eww_ij.dat") == 0) {
      double dbl;
      for (int a=1; a < MAX_GRID_PT_; a++) {
        for (int l=0; l<a; l++) {
          dbl = ww_Eij_[a][l];
          if (dbl != 0) {
            dbl /= (NFRAME_*2);
           outfile.Printf("%10d %10d %12.5E\n", a, l, dbl);
           //Eijtot += dbl;
           }
        }
      }
      outfile.CloseFile();
    } else
      mprinterr("Error: Could not open 'Eww_ij.dat' for writing.\n");
  }
  //Eijtot *= 2.0;
  // Debug Eww_ij
  //mprintf("\tTotal grid energy: %9.5f\n", Eijtot);
}
