#include "Action_Contacts.h"
#include "CpptrajStdio.h"
#include "DistRoutines.h"
// TODO: Separate byResidue stuff

// CONSTRUCTOR
Action_Contacts::Action_Contacts() :
  byResidue_(false),
  distance_(49.0), // 7.0^2
  dt_(1.0),
  first_(false),
  CurrentParm_(0)
{ }

void Action_Contacts::Help() {
  mprintf("\t[ first | reference | ref <ref> | refindex <#> ] [byresidue]\n"
          "\t[out <filename>] [time <interval>] [distance <cutoff>] [<mask>]\n"
          "  Calculate contacts for each frame based on a reference.\n"
          "    byresidue: calculate number of contacts and save results per residue\n");
}

// DESTRUCTOR
Action_Contacts::~Action_Contacts() {
  outfile_.CloseFile();
  if (byResidue_)
    outfile2_.CloseFile();
}

/** Set up native contacts based on reference structure. */
int Action_Contacts::SetupContacts(Frame const& refframe, Topology const& refparm) {
  // Determine which pairs of atoms satisfy cutoff, build contact list
  AtomMask::const_iterator atom1end = Mask_.end() - 1;
  for (AtomMask::const_iterator atom1 = Mask_.begin();
                                atom1 != atom1end; ++atom1)
  {
    for (AtomMask::const_iterator atom2 = atom1 + 1;
                                  atom2 != Mask_.end(); ++atom2)
    {
      double d2 = DIST2_NoImage(refframe.XYZ(*atom1), refframe.XYZ(*atom2));
      if (d2 < distance_)
        nativecontacts_.insert( contactType(*atom1, *atom2) );
    }
  }
  // DEBUG - Print contacts
  mprintf("\tSetup %zu native contacts:\n", nativecontacts_.size());
  for (contactListType::iterator contact = nativecontacts_.begin();
                                 contact != nativecontacts_.end(); ++contact)
  {
    int a1 = (*contact).first;
    int a2 = (*contact).second;
    mprintf("\t\tAtom %i[%s] to %i[%s]\n", a1+1, refparm[a1].c_str(),
            a2+1, refparm[a2].c_str());
  }
  return 0;
}

// Action_Contacts::Init()
Action::RetType Action_Contacts::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{

  byResidue_ = actionArgs.hasKey("byresidue");
  double dist = actionArgs.getKeyDouble("distance", 7.0);
  dt_ = actionArgs.getKeyDouble("time", 1.0);
  // Square the cutoff
  distance_ = dist * dist;
  first_ = actionArgs.hasKey("first");
  // Get reference
  ReferenceFrame REF = FL->GetFrameFromArgs( actionArgs );
  if (REF.error()) return Action::ERR;
  std::string outfilename = actionArgs.GetStringKey("out"); 
  if (outfile_.OpenEnsembleWrite(outfilename, DSL->EnsembleNum()))
    return Action::ERR;
  if (byResidue_) {
    if (outfilename.empty()) {
      mprinterr("Error: Contacts 'byresidue' requires output filename.\n");
      return Action::ERR;
    }
    if (outfile2_.OpenEnsembleWrite( outfilename + ".native", DSL->EnsembleNum() ))
      return Action::ERR;
  }

  // Get Mask
  std::string mask0 = actionArgs.GetMaskNext();
  if (mask0.empty() && byResidue_)
    Mask_.SetMaskString("@CA");
  else
    Mask_.SetMaskString( mask0 );
  
  // Initialize reference. If no reference mask is given mask0 will be used.
  // First arg 'nofit' set to true, no fitting with contacts. Allows last arg
  // 'RefTrans' to be null.
  //if (RefInit(true, false, Mask_.MaskString(), actionArgs, FL, PFL, 0)!=0)
  //  return 1;
  if (!first_ && REF.empty()) {
    mprintf("\tNo reference structure specified. Defaulting to first.\n");
    first_ = true;
  }

  if (!first_) {
    // TODO: Convert FrameList to return frame reference?
    // Set up atom mask for reference frame
    if (REF.Parm().SetupIntegerMask(Mask_, REF.Coord())) return Action::ERR;
    // Set up reference contacts 
    SetupContacts(REF.Coord(), REF.Parm());
  }

  // Output file header - only if not byresidue
  if (!byResidue_) {
    outfile_.Printf("#time\tContacts\tnative Contacts ");
    if (!first_)
      outfile_.Printf("(number of natives: %zu)", nativecontacts_.size());
    outfile_.Printf("\n");
  }

  mprintf("    CONTACTS: [%s] Calculating current contacts and comparing results to",
          Mask_.MaskString());
  if (first_)
    mprintf(" first frame.\n");
  else
    mprintf(" reference structure.\n");
  mprintf("                   Distance cutoff is %lf angstroms.\n", dist);
  if (outfilename.empty())
    mprintf("              Results will be written to stdout");
  else
    mprintf("              Writing results to %s", outfilename.c_str());
  mprintf("\n");
  if (byResidue_)
    mprintf("              Results are output on a per-residue basis.\n");

  return Action::OK;
}

Action::RetType Action_Contacts::Setup(Topology* currentParm, Topology** parmAddress) {
  //if (first_) 
  //  RefParm_ = currentParm;
  // Set up atom mask 
  if (currentParm->SetupIntegerMask(Mask_)) return Action::ERR;

  // Determine which residues are active based on the mask
  activeResidues_.clear();
  for (AtomMask::const_iterator atom = Mask_.begin();
                                atom != Mask_.end(); ++atom)
  {
    int resnum = (*currentParm)[*atom].ResNum();
    activeResidues_.insert( resnum );
  }

  // byresidue header - only on first time through
  if (residueContacts_.empty() && byResidue_) {
    outfile_.Printf("#time");
    outfile2_.Printf("#time");
    for (std::set<int>::iterator res = activeResidues_.begin();
                                 res != activeResidues_.end(); ++res)
    {
      outfile_.Printf("\tresidue %i", *res);
      outfile2_.Printf("\tresidue %i", *res);
    }
    outfile_.Printf("\tTotal\n");
    outfile2_.Printf("\tTotal\n");
  }

  // Reserve space for residue contact counts
  residueContacts_.reserve( currentParm->Nres() );
  residueNative_.reserve( currentParm->Nres() );
  
  CurrentParm_ = currentParm;
  return Action::OK;
}

Action::RetType Action_Contacts::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {
  if (first_) {
    SetupContacts( *currentFrame, *CurrentParm_ );
    first_ = false;
  }

  // Determine how many contacts and how many native contacts are present.
  //contactListType::iterator native = nativecontacts_.begin();
  residueContacts_.assign( CurrentParm_->Nres(), 0 );
  residueNative_.assign( CurrentParm_->Nres(), 0 );

  // Determine which pairs of atoms satisfy cutoff
  int numcontacts = 0;
  int numnative = 0;
  AtomMask::const_iterator atom1end = Mask_.end() - 1;
  for (AtomMask::const_iterator atom1 = Mask_.begin();
                                atom1 != atom1end; ++atom1)
  {
    int res1 = (*CurrentParm_)[*atom1].ResNum();
    for (AtomMask::const_iterator atom2 = atom1 + 1;
                                  atom2 != Mask_.end(); ++atom2)
    {
      double d2 = DIST2_NoImage(currentFrame->XYZ(*atom1), currentFrame->XYZ(*atom2));
      // Contact?
      if (d2 < distance_) {
        ++numcontacts;
        int res2 = (*CurrentParm_)[*atom2].ResNum();
        mprintf("CONTACT: %i res %i to %i res %i [%i]",*atom1+1,res1+1,*atom2+1,res2+1,numcontacts);
        ++residueContacts_[res1];
        ++residueContacts_[res2];
        // Is this a native contact?
        contactListType::iterator nativebegin = nativecontacts_.lower_bound( *atom1 );
        if (nativebegin != nativecontacts_.end()) {
          contactListType::iterator nativeend = nativecontacts_.upper_bound( *atom1 );
          for (contactListType::iterator native = nativebegin; native!=nativeend; ++native) {
            if ( *atom2 == (*native).second ) {
              ++numnative;
              mprintf(" NATIVE [%i]",numnative);
              ++residueNative_[res1];
              ++residueNative_[res2];
              break;
            }
          }
        }
        mprintf("\n");
      }
/*
      if (native != nativecontacts_.end()) {
        if ( *atom1 == (*native).first && *atom2 == (*native).second ) {
          if (d2 < distance_)
            ++numnative;
          ++native;
        }
      }
*/
    }
  }

  // The total # of contacts is multiplied by 2 since each contact
  // is bidirectional, i.e. atom1->atom2 and atom2->atom1 each
  // count as a separate contact.
  if (!byResidue_) {
    outfile_.Printf("%10.2f\t%i\t%i\n", (double)(frameNum+1) * dt_,
                    numcontacts*2, numnative*2);
  } else {
    outfile_.Printf("%10.2f", (double)(frameNum+1) * dt_);
    outfile2_.Printf("%10.2f", (double)(frameNum+1) * dt_);
    for (std::set<int>::iterator res = activeResidues_.begin();
                                 res != activeResidues_.end(); ++res)
    {
      outfile_.Printf("\t%i", residueContacts_[ *res ]);
      outfile2_.Printf("\t%i", residueNative_[ *res ]);
    }
    outfile_.Printf("\t%i\n", numcontacts*2);
    outfile2_.Printf("\t%i\n", numnative*2);
  }

  return Action::OK;
}
     
