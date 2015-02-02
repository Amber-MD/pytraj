// Action_Jcoupling
#include <cstdlib> //getenv
#include <cstdio> //sscanf
#include <cstring> //strcpy, strlen
#include <cmath> //cos
#include "Action_Jcoupling.h"
#include "CpptrajStdio.h"
#include "Constants.h" // DEGRAD, RADDEG
#include "TorsionRoutines.h"

// CONSTRUCTOR
Action_Jcoupling::Action_Jcoupling() :
  debug_(0),
  Nconstants_(0),
  CurrentParm_(0),
  masterDSL_(0),
  setcount_(1)
{} 

void Action_Jcoupling::Help() {
  mprintf("\t<mask1> [outfile <filename>] [kfile <param file>] [out <filename>]\n"
          "\t[name <dsname>]\n"
          "  Calculate J-coupling values for all dihedrals found in <mask1>.\n");
}

// DESTRUCTOR
Action_Jcoupling::~Action_Jcoupling() {
  //fprintf(stderr,"Jcoupling Destructor.\n");
  for (karplusConstantMap::iterator reslist = KarplusConstants_.begin();
                                    reslist != KarplusConstants_.end(); ++reslist)
  {
    karplusConstantList* currentList = (*reslist).second;
    delete currentList;
  }
  // Close output
  outputfile_.CloseFile();
}

// Action_Jcoupling::loadKarplus()
/** Load Karplus parameters from input file.
  * Expected format:
  * - {type}<+|-| ><a[4]><+|-| ><b[4]><+|-| ><c[4]><+|-| ><d[4]><A[6]><B[6]><C[6]>{<D[6]>}
  *   <reslabel[4]>* 
  * \return 0 on success, 1 on error
  */
int Action_Jcoupling::loadKarplus(std::string filename) {
  char buffer[512],residue[5];
  char *end, *ptr;
  int i;
  CpptrajFile KarplusFile;
  karplusConstant KC;
  karplusConstantList* currentResList=0;
  std::string CurrentRes;
  karplusConstantMap::iterator reslist;

  if (filename.empty()) {
    mprinterr("Error: jcoupling: Could not find Karplus parameter file.\n");
    return 1;
  }
  if (KarplusFile.OpenRead( filename )) {
    mprinterr("Error: jcoupling: Could not read Karplus parameter file %s\n",
              filename.c_str());
    mprinterr("Error: Ensure the file exists and is readable.\n");
    return 1;
  }
  // residue is only for reading in 4 chars for residue names
  residue[4]='\0'; 
  // Read through all lines of the file
  while (KarplusFile.Gets(buffer,512)==0) {
    // Skip newlines and comments
    if (buffer[0]=='\n' || buffer[0]=='#') continue;
    ptr=buffer;
    // First char is optional type. If optional type is C, then the Karplus 
    // function specified in Perez et al. JACS (2001) 123 will be used, and 
    // A, B, and C will be taken as C0, C1, and C2.
    if(ptr[0]=='C') {
      KC.type=1;
      ptr++;
    } else {
      KC.type=0;
    }
    // Read atom names with optional preceding character (+, -)
    for (i=0; i<4; i++) {
      if      (*ptr=='+') KC.offset[i]=1;
      else if (*ptr=='-') KC.offset[i]=-1;
      else                     KC.offset[i]=0;
      ++ptr;
      char *endchar = ptr + 4;
      char savechar = *endchar;
      *endchar = '\0';
      KC.atomName[i] = ptr;
      *endchar = savechar;
      ptr += 4;
      //mprintf("DEBUG:\tAtomName %i [%s]\n",i,KC.atomName[i]);
    }
    // Read parameters
    // NOTE: Using sscanf here instead of atof since the 4th parameter is
    //       optional, behavior is undefined for accessing uninitialized
    //       portion of buffer.
    i = sscanf(ptr, "%6lf%6lf%6lf%6lf",KC.C,KC.C+1,KC.C+2,KC.C+3);
    if (i<3) {
      mprintf("Error: jcoupling: Expected at least 3 Karplus parameters, got %i\n",i);
      mprintf("       Line: [%s]\n",buffer);
      return 1;
    } else if (i==3) KC.C[3]=0.0;
    KC.C[3]*=Constants::DEGRAD;
    // Place the read-in karplus constants in a map indexed by residue name 
    // so that all karplus constants for a given residue are in one place. 
    KarplusFile.Gets(buffer,512);
    // end will hold the end of the read-in buffer string
    end = buffer + strlen(buffer);
    for (ptr = buffer; ptr < end; ptr+=4) {
      if (*ptr=='\n') continue;
      residue[0] = ptr[0];
      residue[1] = ptr[1];
      residue[2] = ptr[2];
      residue[3] = ptr[3];
      CurrentRes.assign(residue);
      //mprintf("DEBUG:\t[%s]\n",CurrentRes.c_str());
      reslist = KarplusConstants_.find(CurrentRes);
      if (reslist == KarplusConstants_.end() ) {
        // List does not exist for residue yet, create it.
        currentResList = new karplusConstantList;
        KarplusConstants_.insert( reslist, 
                                  std::pair<std::string,karplusConstantList*>(
                                    CurrentRes,currentResList) );
      } else
        // Retrieve list for residue.
        currentResList = (*reslist).second;

      currentResList->push_back(KC);
      ++Nconstants_;
    } // END loop over residues in residue line 
  } // END Gets over input file
  KarplusFile.CloseFile();
  // DEBUG - Print out all parameters
  if (debug_>0) {
      mprintf("    KARPLUS PARAMETERS:\n");
      for (reslist=KarplusConstants_.begin(); reslist!=KarplusConstants_.end(); ++reslist) 
      {
        mprintf("\t[%4s]\n",(*reslist).first.c_str());
        for (karplusConstantList::iterator kc = currentResList->begin();
                                           kc != currentResList->end(); ++kc) 
        {
          mprintf("\t\t%1i",(*kc).type);
          mprintf(" %4s",*((*kc).atomName[0]));
          mprintf(" %4s",*((*kc).atomName[1]));
          mprintf(" %4s",*((*kc).atomName[2]));
          mprintf(" %4s",*((*kc).atomName[3]));
          mprintf(" %i %i %i %i",(*kc).offset[0],(*kc).offset[1],(*kc).offset[2],(*kc).offset[3]);
          mprintf(" %6.2lf %6.2lf %6.2lf %6.2lf\n",(*kc).C[0],(*kc).C[1],(*kc).C[2],(*kc).C[3]);
        }
      }
  }
  return 0;
}

// -----------------------------------------------------------------------------
// Action_Jcoupling::Init()
Action::RetType Action_Jcoupling::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  debug_ = debugIn;
  outfile_ = 0;
  // Get Keywords
  std::string outfilename = actionArgs.GetStringKey("outfile");
  outfile_ = DFL->AddDataFile( actionArgs.GetStringKey("out"), actionArgs );
  std::string karpluspath = actionArgs.GetStringKey("kfile");
  setname_ = actionArgs.GetStringKey("name");
  // Get Masks
  Mask1_.SetMaskString( actionArgs.GetMaskNext() );

  // If no Karplus params specified check environment vars. 
  if (karpluspath.empty()) {
    // Check if the KARPLUS env var is set.
    const char* env = getenv("KARPLUS");
    if (env != 0) {
      mprintf("Info: Using parameter file defined by $KARPLUS environment variable.\n");
      karpluspath.assign(env);
    } else {
      // If KARPLUS not set check for $AMBERHOME/dat/Karplus.txt
      env = getenv("AMBERHOME");
      if (env == 0) {
        mprinterr("Error: Either AMBERHOME must be set or KARPLUS must point\n"
                  "Error:   to the file containing Karplus parameters.\n");
        return Action::ERR;
      }
      mprintf("Info: Using parameter file in '$AMBERHOME/dat/'.\n");
      karpluspath.assign(env);
      karpluspath += "/dat/Karplus.txt";
    }
  }
  // Load Karplus parameters
  if (loadKarplus(karpluspath)) 
    return Action::ERR;

  mprintf("    J-COUPLING: Searching for dihedrals in mask [%s].\n"
          "\tUsing Karplus parameters in \"%s\"\n"
          "\t%i parameters found for %zu residues.\n",
          Mask1_.MaskString(), karpluspath.c_str(), Nconstants_, KarplusConstants_.size());
  if (outfile_ != 0)
    mprintf("\tDataSets will be written to %s\n", outfile_->DataFilename().full());
  if (!outfilename.empty())
    mprintf("                Writing fixed-format output to %s\n",outfilename.c_str());
  mprintf("# Citations: Chou et al. JACS (2003) 125 p.8959-8966\n"
          "#            Perez et al. JACS (2001) 123 p.7081-7093\n");
  // Open output
  if (!outfilename.empty()) {
    if ( outputfile_.OpenEnsembleWrite( outfilename, DSL->EnsembleNum() ) ) return Action::ERR;
  }
  DSL->SetDataSetsPending(true);
  masterDSL_ = DSL;
  return Action::OK;
}

// Action_Jcoupling::Setup()
/** Set up a j-coupling calculation for dihedrals defined by atoms within
  * the mask.
  */
Action::RetType Action_Jcoupling::Setup(Topology* currentParm, Topology** parmAddress) {
  std::string resName;
  karplusConstantList* currentResList=0;
  int MaxResidues;
  jcouplingInfo JC;

  if ( currentParm->SetupCharMask(Mask1_) ) return Action::ERR;
  if (Mask1_.None()) {
    mprinterr("Error: Mask specifies no atoms.\n");
    return Action::ERR;
  }
  // If JcouplingInfo has already been set up, print a warning and reset for
  // new parm.
  if (!JcouplingInfo_.empty()) {
    mprintf("Warning: Jcoupling has been set up for another parm.\n"
            "Warning:   Resetting jcoupling info for new parm %s\n",currentParm->c_str());
    //JcouplingInfo_.clear();
  }

  // For each residue, set up 1 jcoupling calc for each parameter defined in
  // KarplusConstants for this residue. Only set up the Jcoupling calc if all
  // atoms involved are present in the mask.
  Range resRange = currentParm->SoluteResidues();
  for (Range::const_iterator residue = resRange.begin();
                             residue != resRange.end(); ++residue)
  {
    // Skip residue if no atoms within residue are selected.
    if (!Mask1_.AtomsInCharMask(currentParm->Res(*residue).FirstAtom(),
                                currentParm->Res(*residue).LastAtom())) continue;
    resName.assign(currentParm->Res(*residue).c_str());
    karplusConstantMap::iterator reslist = KarplusConstants_.find(resName);
    // If list does not exist for residue, skip it.
    if (reslist == KarplusConstants_.end() ) {
      mprintf("Warning: Karplus parameters not found for residue [%i:%s]\n",
              *residue+1, resName.c_str());
      continue;
    }
    currentResList = (*reslist).second;
    // For each parameter set in the list find the corresponding atoms.
    for (karplusConstantList::iterator kc = currentResList->begin();
                                       kc != currentResList->end(); ++kc) 
    {
      // Init jcoupling info. Constants will point inside KarplusConstants.
      JC.residue = *residue;
      JC.atom[0] = -1;
      JC.atom[1] = -1;
      JC.atom[2] = -1;
      JC.atom[3] = -1;
      JC.C = kc->C;
      JC.type = kc->type;
      // For each atom in the dihedral specified in this Karplus constant, find
      // corresponding atoms in parm.
      bool allAtomsFound = true;
      for (int idx=0; idx < 4; idx++) {
        JC.atom[idx] = currentParm->FindAtomInResidue(*residue + kc->offset[idx],
                                                      kc->atomName[idx]       );
        if (JC.atom[idx] == -1) {
          mprintf("Warning: Atom '%s' at position %i not found for residue %i\n",
                    *(kc->atomName[idx]), idx, *residue+kc->offset[idx]+1);
          allAtomsFound = false;
        }
      }
      if (allAtomsFound) { 
        // Check that all the atoms involved in this Jcouple dihedral are
        // in the atom mask. If so, add jcoupling info to the list.
        if (Mask1_.AtomInCharMask(JC.atom[0]) && Mask1_.AtomInCharMask(JC.atom[1]) &&
            Mask1_.AtomInCharMask(JC.atom[2]) && Mask1_.AtomInCharMask(JC.atom[3]))
        {
          // TODO: Look for previously set up matching data set
          if (setname_.empty())
            setname_ = masterDSL_->GenerateDefaultName("JC");
          JC.data_ = masterDSL_->AddSetIdx( DataSet::FLOAT, setname_, setcount_++ );
          if ( JC.data_ != 0 ) {
            JC.data_->SetLegend( currentParm->TruncResNameNum(JC.residue) + "_" +
                                 (*currentParm)[JC.atom[0]].Name().Truncated() + "-" +
                                 (*currentParm)[JC.atom[1]].Name().Truncated() + "-" +
                                 (*currentParm)[JC.atom[2]].Name().Truncated() + "-" +
                                 (*currentParm)[JC.atom[3]].Name().Truncated()  );
            if (outfile_ != 0)
              outfile_->AddSet( JC.data_ ); 
            JcouplingInfo_.push_back(JC);
          } else {
            mprinterr("Internal Error: Could not set up Jcoupling data set for res %i\n",
                      JC.residue+1);
          }
        }
      }
    } // END loop over karplus parameters for this residue
  } // END loop over all residues

  // Print info for this parm
  mprintf("    J-COUPLING: [%s] Will calculate J-coupling for %zu dihedrals.\n",
          Mask1_.MaskString(), JcouplingInfo_.size());
  if (JcouplingInfo_.empty()) {
    mprintf("Warning: No dihedrals found for J-coupling calculation!\n"
            "Warning:   Check that all atoms of dihedrals are included in mask [%s]\n"
            "Warning:   and/or that dihedrals are defined in Karplus parameter file.\n",
            Mask1_.MaskString());
    return Action::ERR;
  }
  // DEBUG
  if (debug_>0) {
    MaxResidues=1;
    for (std::vector<jcouplingInfo>::iterator jc = JcouplingInfo_.begin();
                                              jc != JcouplingInfo_.end(); ++jc) 
    {
      mprintf("%8i [%i:%4s]",MaxResidues,jc->residue, currentParm->Res(jc->residue).c_str());
      mprintf(" %6i:%-4s",jc->atom[0],(*currentParm)[jc->atom[0]].c_str());
      mprintf(" %6i:%-4s",jc->atom[1],(*currentParm)[jc->atom[1]].c_str());
      mprintf(" %6i:%-4s",jc->atom[2],(*currentParm)[jc->atom[2]].c_str());
      mprintf(" %6i:%-4s",jc->atom[3],(*currentParm)[jc->atom[3]].c_str());
      mprintf(" %6.2lf%6.2lf%6.2lf%6.2lf %i\n",jc->C[0],jc->C[1],jc->C[2],jc->C[3],
              jc->type);
      MaxResidues++;
    }
  }
  CurrentParm_ = currentParm;    
  return Action::OK;  
}

// Action_Jcoupling::DoAction()
/** For each dihedral defined in JcouplingInfo, perform the dihedral and
  * Jcoupling calculation.
  */
Action::RetType Action_Jcoupling::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {
  double Jval;

  if (outputfile_.IsOpen())
    outputfile_.Printf("#Frame %i\n",frameNum+1);

  for (std::vector<jcouplingInfo>::iterator jc = JcouplingInfo_.begin();
                                            jc !=JcouplingInfo_.end(); ++jc)
  {
    double phi = Torsion(currentFrame->XYZ(jc->atom[0]),
                         currentFrame->XYZ(jc->atom[1]),
                         currentFrame->XYZ(jc->atom[2]),
                         currentFrame->XYZ(jc->atom[3]) );
    if (jc->type==1) {
      //phitemp = phi + jc->C[3]; // Only necessary if offsets become used in perez-type calc
      Jval = jc->C[0] + (jc->C[1] * cos(phi)) + (jc->C[2] * cos(phi * 2.0)); 
    } else {
      double phitemp = cos( phi + jc->C[3] );
      Jval = (jc->C[0] * phitemp * phitemp) + (jc->C[1] * phitemp) + jc->C[2];
    }
    float fval = (float)Jval;
    jc->data_->Add(frameNum, &fval);

    int residue = jc->residue;
    // Output
    if (outputfile_.IsOpen())
      outputfile_.Printf("%5i %4s%4s%4s%4s%4s%12f%12f\n",
                         residue+1, CurrentParm_->Res(residue).c_str(),
                         (*CurrentParm_)[jc->atom[0]].c_str(), 
                         (*CurrentParm_)[jc->atom[1]].c_str(),
                         (*CurrentParm_)[jc->atom[2]].c_str(), 
                         (*CurrentParm_)[jc->atom[3]].c_str(),
                         phi*Constants::RADDEG, Jval);
  }

  return Action::OK;
} 
