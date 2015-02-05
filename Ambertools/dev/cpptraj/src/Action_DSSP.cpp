#include <cmath> // sqrt
#include "Action_DSSP.h"
#include "CpptrajStdio.h"
#include "DistRoutines.h"
/// Hbond energy calc prefactor for kcal/mol: q1*q2*E, 0.42*0.20*332
const double Action_DSSP::DSSP_fac = 27.888;

const int Action_DSSP::NSSTYPE = 8;

// CONSTRUCTOR
Action_DSSP::Action_DSSP() :
  ensembleNum_(-1),
  debug_(0),
  outfile_(0),
  dsspFile_(0),
  Nres_(0),
  Nselected_(0),
  Nframe_(0),
  printString_(false),
  masterDSL_(0),
  BB_N_("N"),
  BB_H_("H"),
  BB_C_("C"),
  BB_O_("O"),
  BB_CA_("CA")
{}

void Action_DSSP::Help() {
  mprintf("\t[out <filename>] [<mask>] [sumout <filename>] [assignout <filename>]\n"
          "\t[ptrajformat] [namen <N name>] [nameh <H name>] [nameca <CA name>]\n"
          "\t[namec <C name>] [nameo <O name>] [totalout <filename>]\n"
          "  Calculate secondary structure content for residues in <mask>.\n"
          "  If sumout not specified, the filename specified by out is used with .sum suffix.\n");
}

const char  Action_DSSP::dssp_char[] = { ' ', 'E', 'B', 'G', 'H', 'I', 'T', 'S' };
const char* Action_DSSP::SSchar[]    = { "0", "b", "B", "G", "H", "I", "T", "S" };
const char* Action_DSSP::SSname[]={"None", "Para", "Anti", "3-10", "Alpha", "Pi", "Turn", "Bend"};

// Action_DSSP::Init()
Action::RetType Action_DSSP::Init(ArgList& actionArgs, TopologyList* PFL, DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  ensembleNum_ = DSL->EnsembleNum();
  debug_ = debugIn;
  Nframe_ = 0;
  // Get keywords
  outfile_ = DFL->AddDataFile( actionArgs.GetStringKey("out"), actionArgs );
  std::string temp = actionArgs.GetStringKey("sumout");
  if (temp.empty() && outfile_ != 0) 
    temp = outfile_->DataFilename().Full() + ".sum";
  dsspFile_ = DFL->AddDataFile( temp );
  DataFile* totalout = DFL->AddDataFile( actionArgs.GetStringKey("totalout"), actionArgs );
  assignout_ = actionArgs.GetStringKey("assignout");
  printString_ = actionArgs.hasKey("ptrajformat");
  temp = actionArgs.GetStringKey("namen");
  if (!temp.empty()) BB_N_ = temp;
  temp = actionArgs.GetStringKey("nameh");
  if (!temp.empty()) BB_H_ = temp;
  temp = actionArgs.GetStringKey("namec");
  if (!temp.empty()) BB_C_ = temp;
  temp = actionArgs.GetStringKey("nameo");
  if (!temp.empty()) BB_O_ = temp;
  temp = actionArgs.GetStringKey("nameca");
  if (!temp.empty()) BB_CA_ = temp;
  // Get masks
  Mask_.SetMaskString( actionArgs.GetMaskNext() );

  // Set up the DSSP data set
  dsetname_ = actionArgs.GetStringNext();
  // Set default name if none specified
  if (dsetname_.empty())
    dsetname_ = DSL->GenerateDefaultName("DSSP");
  // Set up Z labels
  if (outfile_ != 0)
    outfile_->ProcessArgs("zlabels None,Para,Anti,3-10,Alpha,Pi,Turn,Bend");
  // Create data sets for total fraction SS vs time.
  for (int i = 0; i < NSSTYPE; i++) {
    totalDS_[i] = DSL->AddSetAspect(DataSet::FLOAT, dsetname_, SSname[i]);
    if (totalDS_[i] == 0) {
      mprinterr("Error: Could not create DSSP total frac v time data set.\n");
      return Action::ERR;
    }
    // For now dont add 'None' so colors match up.
    if (totalout != 0 && i > 0) totalout->AddSet( totalDS_[i] );
  }

  mprintf( "    SECSTRUCT: Calculating secondary structure using mask [%s]\n",Mask_.MaskString());
  if (outfile_ != 0) 
    mprintf("\tDumping results to %s\n", outfile_->DataFilename().full());
  if (dsspFile_ != 0)
    mprintf("\tSum results to %s\n", dsspFile_->DataFilename().full());
  if (printString_) { 
    mprintf("\tSS data for each residue will be stored as a string.\n");
    for (int i = 0; i < NSSTYPE; i++)
      mprintf("\t\t%s = %s\n", SSchar[i], SSname[i]);
  } else {
    mprintf("\tSS data for each residue will be stored as integers.\n");
    for (int i = 0; i < NSSTYPE; i++)
      mprintf("\t\t%i = %s\n", i, SSname[i]);
  }
  if (!assignout_.empty())
    mprintf("\tOverall assigned SS will be written to %s\n", assignout_.c_str());
  mprintf("\tBackbone Atom Names: N=[%s]  H=[%s]  C=[%s]  O=[%s]  CA=[%s]\n",
          *BB_N_, *BB_H_, *BB_C_, *BB_O_, *BB_CA_ );
  mprintf("# Citation: Kabsch, W.; Sander, C.; \"Dictionary of Protein Secondary Structure:\n"
          "#           Pattern Recognition of Hydrogen-Bonded and Geometrical Features.\"\n"
          "#           Biopolymers (1983), V.22, pp.2577-2637.\n" );
  DSL->SetDataSetsPending(true);
  masterDSL_ = DSL;
  return Action::OK;
}

// Action_DSSP::Setup()
/** Set up secondary structure calculation for all residues selected by the
  * mask expression. A residue is selected if at least one of the following 
  * atoms named "C   ", "O   ", "N   ", or "H   " (i.e. standard atom protein 
  * BB atom names) is selected. The residue is only initialized if it has not
  * been previously selected and set up by a prior call to setup.
  */
// NOTE: Currently relatively memory-intensive. Eventually set up so that SecStruct and
// CO_HN_Hbond members exist only for selected residues? Use Map?
Action::RetType Action_DSSP::Setup(Topology* currentParm, Topology** parmAddress) {
  // Set up mask for this parm
  if ( currentParm->SetupIntegerMask( Mask_ ) ) return Action::ERR;
  if ( Mask_.None() ) {
    mprintf("Warning: DSSP: Mask has no atoms.\n");
    return Action::ERR;
  }

  // Initially mark all residues already set up as not selected and 
  // reset all atom coordinate indices.
  for (unsigned int res = 0; res != SecStruct_.size(); res++) {
    SecStruct_[res].isSelected = false;
    SecStruct_[res].C  = -1;
    SecStruct_[res].H  = -1;
    SecStruct_[res].N  = -1;
    SecStruct_[res].O  = -1;
    SecStruct_[res].CA = -1;
  }

  // Set up SecStruct for each solute residue
  Range soluteRes = currentParm->SoluteResidues();
  Nres_ = soluteRes.Back() + 1;
  if (debug_>0) mprintf("\tDSSP: Setting up for %i residues.\n", Nres_);

  // Set up for each residue of the current Parm if not already set-up.
  SSres RES;
  RES.sstype = NONE;
  RES.isSelected = false;
  RES.C  = -1;
  RES.O  = -1;
  RES.N  = -1;
  RES.H  = -1;
  RES.CA = -1;
  RES.CO_HN_Hbond.assign( Nres_, 0 );
  std::fill( RES.SSprob, RES.SSprob + NSSTYPE, 0 );
  RES.resDataSet = 0;
  // Only resize SecStruct if current # residues > previous # residues
  if (Nres_ > (int)SecStruct_.size())
    SecStruct_.resize(Nres_, RES);

  // Go through all atoms in mask. Determine which residues have their C,
  // O, N, or H atoms selected. Store the actual coordinate index 
  // (i.e. atom# * 3) instead of atom # for slight speed gain. 
  for (AtomMask::const_iterator atom = Mask_.begin(); atom!=Mask_.end(); ++atom) {
    int res = (*currentParm)[*atom].ResNum();
    // If residue is out of bounds skip it
    if ( res < Nres_ ) {
      //fprintf(stdout,"DEBUG: Atom %i Res %i [%s]\n",*atom,res,P->names[*atom]);
      SecStruct_[res].isSelected = true;
      if (      (*currentParm)[*atom].Name() == BB_C_)
        SecStruct_[res].C = (*atom) * 3;
      else if ( (*currentParm)[*atom].Name() == BB_O_)
        SecStruct_[res].O = (*atom) * 3;
      else if ( (*currentParm)[*atom].Name() == BB_N_)
        SecStruct_[res].N = (*atom) * 3;
      else if ( (*currentParm)[*atom].Name() == BB_H_)
        SecStruct_[res].H = (*atom) * 3;
      else if ( (*currentParm)[*atom].Name() == BB_CA_)
        SecStruct_[res].CA = (*atom) * 3;
    }
  }

  // For each residue selected in the mask, check if residue is missing atoms. 
  // Set up DataSet if necessary. 
  Nselected_ = 0;
  std::vector<std::string> missingResidues;
  for (int res = 0; res < Nres_; ++res) {
    if (SecStruct_[res].isSelected) {
      // Check if C-O/N-H selected.
      SecStruct_[res].hasCO = (SecStruct_[res].C != -1 &&
                               SecStruct_[res].O != -1);
      SecStruct_[res].hasNH = (SecStruct_[res].N != -1 &&
                               SecStruct_[res].H != -1);
      if (!SecStruct_[res].hasCO || !SecStruct_[res].hasNH || SecStruct_[res].CA == -1)
      {
        missingResidues.push_back( currentParm->TruncResNameNum( res ) );
        if (debug_ > 0) {
          mprintf("Warning: Not all BB atoms found for res %s:", missingResidues.back().c_str());
          if (SecStruct_[res].N==-1) mprintf(" N");
          if (SecStruct_[res].H==-1) mprintf(" H");
          if (SecStruct_[res].C==-1) mprintf(" C");
          if (SecStruct_[res].O==-1) mprintf(" O");
          if (SecStruct_[res].CA==-1) mprintf(" CA");
          mprintf("\n");
        }
      }
      // Set up dataset if necessary 
      if (SecStruct_[res].resDataSet == 0) {
        // Setup dataset for this residue
        if (printString_)
          SecStruct_[res].resDataSet =
            masterDSL_->AddSetIdxAspect( DataSet::STRING, dsetname_, res+1, "res");
        else
          SecStruct_[res].resDataSet = 
            masterDSL_->AddSetIdxAspect( DataSet::INTEGER, dsetname_, res+1, "res");
        if (SecStruct_[res].resDataSet == 0) {
          mprinterr("Error: Could not allocate DSSP data set for residue %i\n", res+1);
          return Action::ERR;
        }
        if (outfile_ != 0) outfile_->AddSet(SecStruct_[res].resDataSet);
        SecStruct_[res].resDataSet->SetLegend( currentParm->TruncResNameNum(res) );
      }
      ++Nselected_;
    }
  }
  if (!missingResidues.empty()) {
    mprintf("Warning: Not all BB atoms found for %u residues:", missingResidues.size());
    for (std::vector<std::string>::const_iterator mr = missingResidues.begin();
                                                  mr != missingResidues.end(); ++mr)
      mprintf(" %s", mr->c_str());
    mprintf("\nInfo: This is expected for Proline and terminal/non-standard residues.\n"
              "Info: Expected BB atom names: N=[%s]  H=[%s]  C=[%s]  O=[%s]  CA=[%s]\n",
            *BB_N_, *BB_H_, *BB_C_, *BB_O_, *BB_CA_ );
    mprintf("Info: Re-run with action debug level >= 1 to see which residues are missing atoms.\n");
  }

  // Count number of selected residues
  mprintf("\tMask [%s] corresponds to %u residues.\n", Mask_.MaskString(), Nselected_);
# ifdef DSSPDEBUG
  // DEBUG - Print atom nums for each residue set up
  for (int res=0; res < Nres_; res++) {
    if (SecStruct_[res].isSelected) {
      mprintf("DEBUG: Res %i", res + 1);
      if (SecStruct_[res].hasCO)
        mprintf(" C=%s O=%s",currentParm->AtomMaskName(SecStruct_[res].C/3).c_str(),
                             currentParm->AtomMaskName(SecStruct_[res].O/3).c_str());
      if (SecStruct_[res].hasNH)
        mprintf(" N=%s H=%s",currentParm->AtomMaskName(SecStruct_[res].N/3).c_str(),
                             currentParm->AtomMaskName(SecStruct_[res].H/3).c_str());
      if (SecStruct_[res].CA != -1)
        mprintf(" CA=%s",currentParm->AtomMaskName(SecStruct_[res].CA/3).c_str());
      mprintf("\n");
    }
  }
# endif
  return Action::OK;
}

// Action_DSSP::isBonded()
/** Return 1 if residue 1 CO bonded to residue 2 NH.
  * Ensure residue numbers are valid and residues are selected.
  */
int Action_DSSP::isBonded(int res1, int res2) {
  if (res1<0 || res2<0 || res1>=Nres_ || res2>=Nres_) return 0;
  return SecStruct_[res1].CO_HN_Hbond[res2];
}

/// \return true if type1 has priority over type2.
bool Action_DSSP::HasPriority( SStype type1, SStype type2 ) {
  switch (type1) {
    case ALPHA: return true;
    case ANTI:
      if (type2 != ALPHA) return true;
      break;
    case PARA:
      if (type2 != ALPHA && type2 != ANTI) return true;
      break;
    case H3_10:
      if (type2 != ALPHA && type2 != ANTI && type2 != PARA) return true;
      break;
    case HPI:
      if (type2 != ALPHA && type2 != ANTI && type2 != PARA &&
          type2 != H3_10) return true;
      break;
    case TURN:
      if (type2 != ALPHA && type2 != ANTI && type2 != PARA &&
          type2 != H3_10 && type2 != HPI) return true;
      break;
    case BEND:
      if (type2 != ALPHA && type2 != ANTI && type2 != PARA &&
          type2 != H3_10 && type2 != HPI && type2 != TURN) return true;
      break;
    case NONE: break;
  }
  return false;
}

// Action_DSSP::SSassign()
/** Assign all residues from res1 to res2-1 the secondary structure sstype
  * only if it has priority. 
  * Assumes given residue range is valid.
  */
void Action_DSSP::SSassign(int res1, int res2, SStype typeIn, bool force) {
# ifdef DSSPDEBUG
  mprintf("DEBUG:\tCalling SSassign from %i to %i, %s:", res1+1, res2, SSname[typeIn]);
# endif
  for (int res = res1; res < res2; res++) {
    if (res==Nres_) break;
    if ( HasPriority(typeIn, SecStruct_[res].sstype) || force ) {
#     ifdef DSSPDEBUG
      mprintf(" %i", res+1); // DEBUG
#     endif
      SecStruct_[res].sstype = typeIn;
    }
  }
# ifdef DSSPDEBUG
  mprintf("\n"); // DEBUG
# endif
}
 
// Action_DSSP::DoAction()
/** Determine secondary structure by hydrogen bonding pattern. */    
Action::RetType Action_DSSP::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {
  int resi, resj;
  const double *C, *O, *H, *N;
  double rON, rCH, rOH, rCN, E;
  // Determine C=O to H-N hydrogen bonds for each residue to each other residue
#ifdef _OPENMP
#pragma omp parallel private(resi,resj,C,O,H,N,rON, rCH, rOH, rCN, E)
{
#pragma omp for
#endif
  for (resi = 0; resi < Nres_; resi++) {
    if (SecStruct_[resi].isSelected) {
      // Reset previous SS assignment/Hbond status
      SecStruct_[resi].sstype = NONE;
      SecStruct_[resi].CO_HN_Hbond.assign( Nres_, 0 );
      if (SecStruct_[resi].hasCO) {
        C = currentFrame->CRD(SecStruct_[resi].C);
        O = currentFrame->CRD(SecStruct_[resi].O);
        for (resj = 0; resj < Nres_; resj++) {
          if (SecStruct_[resj].isSelected && resi != resj && SecStruct_[resj].hasNH)
          {
            N = currentFrame->CRD(SecStruct_[resj].N);
            H = currentFrame->CRD(SecStruct_[resj].H);

            rON = 1.0/sqrt(DIST2_NoImage(O, N));
            rCH = 1.0/sqrt(DIST2_NoImage(C, H));
            rOH = 1.0/sqrt(DIST2_NoImage(O, H));
            rCN = 1.0/sqrt(DIST2_NoImage(C, N));

            E = DSSP_fac * (rON + rCH - rOH - rCN);
            if (E < -0.5) {
#             ifdef DSSPDEBUG
              mprintf("DEBUG: %i-CO --> %i-NH  E= %g\n", resi+1, resj+1, E);
#             endif
              SecStruct_[resi].CO_HN_Hbond[resj] = 1;
            }
          }
        }
      }
    }
  }
#ifdef _OPENMP
} // END pragma omp parallel
#endif

  // Determine Secondary Structure based on Hbonding pattern.
  // In case of structural overlap, priority is given to the structure first 
  // in this list (see p. 2587 & 2595 in the Kabsch & Sander paper):
  //   H = 4-helix (alpha helix)
  //   B = residue in isolated beta bridge (anti)
  //   E = extended strand, participates in beta-ladder (para)
  //   G = 3-helix (310 helix)
  //   I = 5-helix (pi helix)
  //   T = H-bonded turn
  //   S = bend
  for (resi=0; resi < Nres_; resi++) {
    if (SecStruct_[resi].isSelected) {
#     ifdef DSSPDEBUG 
      mprintf("DEBUG: Residue %i -----\n", resi+1);
#     endif
      // Alpha helices
      if ( isBonded( resi - 1, resi+3 ) && isBonded( resi, resi + 4) )
        SSassign(resi, resi+4, ALPHA, false);

      // Beta sheets - only needed if SS not already assigned to alpha
      if ( SecStruct_[resi].sstype != ALPHA ) {
        for (resj=0; resj < Nres_; resj++) {
          if (SecStruct_[resj].isSelected) {
            // Only consider residues spaced more than 2 apart
            int abs_resi_resj = resi - resj;
            if (abs_resi_resj<0) abs_resi_resj = -abs_resi_resj;
            if (abs_resi_resj > 2) {
              if ( (isBonded(resi-1, resj) && isBonded(resj, resi+1)) ||
                   (isBonded(resj-1, resi) && isBonded(resi, resj+1)) )
              {
                // Parallel
                // NOTE: Not checking if ANTI since not geometrically possible
                //       to be anti-parallel and parallel at the same time.
#               ifdef DSSPDEBUG
                mprintf("DEBUG:\tAssigning %i to parallel beta.\n", resi+1);
#               endif
                SecStruct_[resi].sstype = PARA;
                break;
              } else if ( (isBonded(resi-1, resj+1) && isBonded(resj-1, resi+1)) ||
                          (isBonded(resi,   resj  ) && isBonded(resj,   resi  )) )
              {
                // Anti-parallel
#               ifdef DSSPDEBUG
                mprintf("DEBUG:\tAssigning %i to anti-parallel beta.\n", resi+1);
#               endif
                SecStruct_[resi].sstype = ANTI;
                break;
              }
            }
          }
        }
      }

      // 3-10 helix
      if ( isBonded( resi - 1, resi+2 ) && isBonded( resi, resi + 3) )
        SSassign(resi, resi+3, H3_10, false);

      // Pi helix
      if ( isBonded( resi - 1, resi+4 ) && isBonded( resi, resi + 5) )
        SSassign(resi, resi+5, HPI, false);
    
      // n-Turn, n=3,4,5
      for (int n=5; n > 2; n--) {
        if ( isBonded(resi, resi + n) ) {
          SSassign(resi+1, resi + n, TURN, false); // FIXME: Should this be resi?
          break;
        }
      }

      // Bend - has lowest priority, so only do if no assignment.
      if (SecStruct_[resi].sstype == NONE && resi > 1 && resi < Nres_ - 2)
      {
        if (SecStruct_[resi-2].CA != -1 &&
            SecStruct_[resi  ].CA != -1 &&
            SecStruct_[resi+2].CA != -1)
        {
          const double* CAm2 = currentFrame->CRD(SecStruct_[resi-2].CA);
          const double* CA0  = currentFrame->CRD(SecStruct_[resi  ].CA);
          const double* CAp2 = currentFrame->CRD(SecStruct_[resi+2].CA);
          Vec3 CA1( CA0[0]-CAm2[0], CA0[1]-CAm2[1], CA0[2]-CAm2[2] );
          Vec3 CA2( CAp2[0]-CA0[0], CAp2[1]-CA0[1], CAp2[2]-CA0[2] );
          CA1.Normalize();
          CA2.Normalize();
          // 1.221730476 rad = 70 degrees
          if (CA1.Angle(CA2) > 1.221730476) {
#           ifdef DSSPDEBUG
            mprintf("DEBUG: Bend calc %i-%i-%i: %g rad.\n", resi-1, resi+1, resi+3, CA1.Angle(CA2));
#           endif
            SecStruct_[resi].sstype = BEND;
          }
        }
      }
    }
  } // End Initial SS assignment over all residues

  // Change 3-10 and Pi helices that are less than minimal size to Turn
  SStype lastType = NONE;
  int resStart = -1;
  for (resi = 0; resi < Nres_; resi++) {
    if (SecStruct_[resi].isSelected) {
      if (lastType != SecStruct_[resi].sstype) {
        // Secondary structure type has changed.
        if (lastType == H3_10) {
#         ifdef DSSPDEBUG
          mprintf("DEBUG: 3-10 helix length is %i\n", resi - resStart);
#         endif
          if (resi - resStart < 3)
            SSassign(resStart, resi, TURN, true);
        } else if (lastType == HPI) {
#         ifdef DSSPDEBUG
          mprintf("DEBUG: PI helix length is %i\n", resi - resStart);
#         endif
          if (resi - resStart < 5)
            SSassign(resStart, resi, TURN, true);
        }
        resStart = resi;
#       ifdef DSSPDEBUG
        mprintf("DEBUG: ResStart=%i for type %s\n", resi+1, SSname[SecStruct_[resi].sstype]);
#       endif
      }
      lastType = SecStruct_[resi].sstype;
    }
  }

  // Store data for each residue. Get statistics.
  int totalSS[NSSTYPE];
  std::fill( totalSS, totalSS + NSSTYPE, 0 ); 
  for (resi=0; resi < Nres_; resi++) {
    if (SecStruct_[resi].isSelected) {
      totalSS[SecStruct_[resi].sstype]++; 
      SecStruct_[resi].SSprob[SecStruct_[resi].sstype]++;
      if (printString_)
        SecStruct_[resi].resDataSet->Add(frameNum, SSchar[SecStruct_[resi].sstype]); 
      else
        SecStruct_[resi].resDataSet->Add(frameNum, &(SecStruct_[resi].sstype));
    }
  }
  for (int i = 0; i < NSSTYPE; i++) {
    float fvar = (float)totalSS[i];
    fvar /= (float)Nselected_;
    totalDS_[i]->Add(frameNum, &fvar);
  }
  ++Nframe_;

  return Action::OK;
}

// Action_DSSP::Print()
void Action_DSSP::Print() {
  if (dsetname_.empty()) return;
  // Try not to print empty residues. Find the minimum and maximum residue
  // for which there is data. Output res nums start from 1.
  int min_res = -1;
  int max_res = -1;
  for (int resi = 0; resi != (int)SecStruct_.size(); resi++) {
    if (SecStruct_[resi].resDataSet != 0) {
      if (min_res < 0) min_res = resi;
      if (resi > max_res) max_res = resi;
    }
  }
  if (min_res < 0 || max_res < min_res) {
    mprinterr("Error: No residues have SS data.\n");
    return;
  }
  // Calculate average of each SS type across all residues.
  if (dsspFile_ != 0) {
    std::vector<DataSet*> dsspData_(NSSTYPE);
    Dimension Xdim( min_res + 1, 1, max_res - min_res + 1, "Residue" );
    // Set up a dataset for each SS type. TODO: NONE type?
    for (int ss = 1; ss < NSSTYPE; ss++) {
      dsspData_[ss] = masterDSL_->AddSetIdxAspect(DataSet::DOUBLE, dsetname_, ss, "avgss");
      dsspData_[ss]->SetLegend( SSname[ss] );
      dsspData_[ss]->SetDim(Dimension::X, Xdim);
      dsspFile_->AddSet( dsspData_[ss] ); 
    }
    
    // Calc the avg SS type for each residue that has data.
    int idx = 0; 
    for (int resi = min_res; resi < max_res+1; resi++) {
      if (SecStruct_[resi].resDataSet != 0) {
        for (int ss = 1; ss < NSSTYPE; ss++) {
          double avg = (double)SecStruct_[resi].SSprob[ss];
          avg /= (double)Nframe_;
          dsspData_[ss]->Add(idx, &avg);
        }
        ++idx;
      }
    }
  }
  // Print out SS assignment like PDB
  if (!assignout_.empty()) {
    CpptrajFile outfile;
    if (outfile.OpenEnsembleWrite(assignout_, ensembleNum_) == 0) {
      int total = 0;
      int startRes = -1;
      std::string resLine, ssLine;
      for (int resi = min_res; resi < max_res+1; resi++) {
        if (startRes == -1) startRes = resi;
        // Convert residue name.
        resLine += Residue::ConvertResName( SecStruct_[resi].resDataSet->Legend() );
        // Figure out which SS element is dominant for res if selected
        if (SecStruct_[resi].resDataSet != 0) {
          int dominantType = 0;
          int ssmax = 0;
          for (int ss = 0; ss < NSSTYPE; ss++) {
            if ( SecStruct_[resi].SSprob[ss] > ssmax ) {
              ssmax = SecStruct_[resi].SSprob[ss];
              dominantType = ss;
            }
          }
          ssLine += dssp_char[dominantType];
        } else
          ssLine += '-';
        total++;
        if ((total % 50) == 0 || resi == max_res) {
          outfile.Printf("%-8i %s\n", startRes+1, resLine.c_str());
          outfile.Printf("%8s %s\n\n", " ", ssLine.c_str());
          startRes = -1;
          resLine.clear();
          ssLine.clear();
        } else if ((total % 10) == 0) {
          resLine += ' '; 
          ssLine += ' ';
        }
      }
    }
  }
}
