#include "Action_AutoImage.h"
#include "CpptrajStdio.h"
#include "DistRoutines.h"
#include "ImageRoutines.h"

// CONSTRUCTOR
Action_AutoImage::Action_AutoImage() :
  origin_(false),
  ortho_(false),
  usecom_(true),
  truncoct_(false),
  useMass_(false),
  triclinic_(OFF)
{}

void Action_AutoImage::Help() {
  mprintf("\t[<mask> | anchor <mask> [fixed <fmask>] [mobile <mmask>]]\n"
          "\t[origin] [firstatom] [familiar | triclinic]\n"
          "  Automatically center and image periodic trajectory.\n"
          "  The 'anchor' molecule (default the first molecule) will be centered;\n"
          "  all 'fixed' molecules will be imaged only if imaging brings them closer\n"
          "  to the 'anchor' molecule; default for 'fixed' molecules is all\n"
          "  non-solvent non-ion molecules. All other molecules (referred to as\n"
          "  'mobile') will be imaged freely.\n");
}

// Action_AutoImage::Init()
Action::RetType Action_AutoImage::Init(ArgList& actionArgs, TopologyList* PFL, DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  // Get keywords
  origin_ = actionArgs.hasKey("origin");
  usecom_ = !actionArgs.hasKey("firstatom");
  if (actionArgs.hasKey("familiar")) triclinic_ = FAMILIAR;
  if (actionArgs.hasKey("triclinic")) triclinic_ = FORCE;
  anchor_ = actionArgs.GetStringKey("anchor");
  fixed_  = actionArgs.GetStringKey("fixed");
  mobile_ = actionArgs.GetStringKey("mobile");
  // Get mask expression for anchor if none yet specified
  if (anchor_.empty())  
    anchor_ = actionArgs.GetMaskNext();

  mprintf("    AUTOIMAGE: To");
  if (origin_)
    mprintf(" origin");
  else
    mprintf(" box center");
  mprintf(" based on");
  if (usecom_)
    mprintf(" center of mass");
  else
    mprintf(" first atom position");
  if (!anchor_.empty())
    mprintf(", anchor mask is [%s]\n", anchor_.c_str());
  else
    mprintf(", anchor is first molecule.\n");
  if (!fixed_.empty())
    mprintf("\tAtoms in mask [%s] will be fixed to anchor region.\n", fixed_.c_str());
  if (!mobile_.empty())
    mprintf("\tAtoms in mask [%s] will be imaged independently of anchor region.\n",
            mobile_.c_str());

  return Action::OK;
}

// Action_AutoImage::SetupAtomRanges()
/** Based on the given atom mask expression determine what molecules are
  * selected by the mask.
  * \return A list of atom pairs that mark the beginning and end of each
  *         selected molecule.
  */
Action_AutoImage::pairList Action_AutoImage::SetupAtomRanges( Topology* currentParm, 
                                             std::string const& maskexpr )
{
  pairList imageList;
  AtomMask Mask1( maskexpr.c_str() );

  if (currentParm->SetupCharMask( Mask1 )) return imageList;
  if (Mask1.None()) return imageList;
  for (Topology::mol_iterator mol = currentParm->MolStart();
                              mol != currentParm->MolEnd(); mol++)
  {
    int firstAtom = (*mol).BeginAtom();
    int lastAtom = (*mol).EndAtom();
    // Check that each atom in the range is in Mask1
    bool rangeIsValid = true;
    for (int atom = firstAtom; atom < lastAtom; ++atom) {
      if (!Mask1.AtomInCharMask(atom)) {
        rangeIsValid = false;
        break;
      }
    }
    if (rangeIsValid) {
      imageList.push_back( firstAtom );
      imageList.push_back( lastAtom );
    }
  }
  mprintf("\tMask [%s] corresponds to %zu molecules\n", Mask1.MaskString(), imageList.size()/2);
  return imageList;
}

// Action_AutoImage::Setup()
Action::RetType Action_AutoImage::Setup(Topology* currentParm, Topology** parmAddress) {
  bool fixedauto = false;
  bool mobileauto = false;

  if (currentParm->Nmol() < 1) {
    mprintf("Warning: autoimage: Parm %s does not contain molecule information\n",
            currentParm->c_str());
    return Action::ERR;
  }
  // Determine Box info
  if (currentParm->BoxType()==Box::NOBOX) {
    mprintf("Warning: autoimage: Parm %s does not contain box information.\n",
            currentParm->c_str());
    return Action::ERR;
  }
  ortho_ = false;
  if (currentParm->BoxType()==Box::ORTHO && triclinic_==OFF) ortho_=true;
  // If box is originally truncated oct and not forcing triclinic, 
  // turn familiar on.
  if (currentParm->BoxType()==Box::TRUNCOCT && triclinic_!=FORCE && triclinic_!=FAMILIAR) {
    mprintf("\tOriginal box is truncated octahedron, turning on 'familiar'.\n");
    triclinic_=FAMILIAR;
  }

  // Set up anchor region
  if (!anchor_.empty()) {
    anchorList_ = SetupAtomRanges( currentParm, anchor_ );
  } else {
    anchorList_.clear();
    anchorList_.push_back( currentParm->Mol(0).BeginAtom() );
    anchorList_.push_back( currentParm->Mol(0).EndAtom() );
  }
  if (anchorList_.empty() || anchorList_.size() > 2) {
    mprinterr("Error: Anchor mask [%s] corresponds to %zu mols, should only be 1.\n",
              anchor_.c_str(), anchorList_.size() / 2);
    return Action::ERR;
  }
  // Set up mask for centering anchor
  anchorMask_.ResetMask();
  anchorMask_.AddAtomRange( anchorList_[0], anchorList_[1] );
  int anchormolnum = (*currentParm)[ anchorList_[0] ].MolNum();
  mprintf("\tAnchor molecule is %i\n", anchormolnum+1);
  // Set up fixed region
  if (!fixed_.empty()) 
    fixedList_ = SetupAtomRanges( currentParm, fixed_ );
  else { 
    fixedauto = true;
    fixedList_.clear();
  }
  // Set up mobile region
  if (!mobile_.empty())
    mobileList_ = SetupAtomRanges( currentParm, mobile_ );
  else {
    mobileauto = true;
    mobileList_.clear();
  }
  // Automatic search through molecules for fixed/mobile
  if (fixedauto || mobileauto) {
    int molnum = 0;
    for (Topology::mol_iterator mol = currentParm->MolStart();
                                mol != currentParm->MolEnd(); mol++)
    {
      // Skip the anchor molecule
      if (molnum != anchormolnum) { 
        // Solvent and 1 atom molecules (prob. ions) go in mobile list,
        // everything else into fixed list.
        if ( (*mol).IsSolvent() || (*mol).NumAtoms() == 1 ) {
          if (mobileauto) {
            mobileList_.push_back( (*mol).BeginAtom() );
            mobileList_.push_back( (*mol).EndAtom()   );
          }
        } else {
          if (fixedauto) {
            fixedList_.push_back( (*mol).BeginAtom() );
            fixedList_.push_back( (*mol).EndAtom()   );
          }
        }
      }
      ++molnum;
    }
  }
  // DEBUG: Print fixed and mobile lists
  if (!fixedList_.empty()) {
    mprintf("\tThe following molecules are fixed to anchor:");
    for (pairList::iterator atom = fixedList_.begin(); 
                            atom != fixedList_.end(); atom += 2)
      mprintf(" %i", (*currentParm)[ *atom ].MolNum()+1 );
    mprintf("\n");
  }
  mprintf("\t%zu molecules are mobile.\n", mobileList_.size() / 2 );
  //mprintf("\tThe following molecules are mobile:\n");
  //for (pairList::iterator atom = mobileList_.begin(); 
  //                        atom != mobileList_.end(); atom += 2)
  //  mprintf("\t\t%i\n", (*currentParm)[ *atom ].MolNum()+1 );

  truncoct_ = (triclinic_==FAMILIAR);

  return Action::OK;
}

// Action_AutoImage::DoAction()
Action::RetType Action_AutoImage::DoAction(int frameNum, Frame* currentFrame, Frame** frameAddress) {
  Matrix_3x3 ucell, recip;
  Vec3 fcom;
  Vec3 bp, bm, offset(0.0);
  Vec3 Trans, framecenter, imagedcenter, anchorcenter;

  if (!ortho_) currentFrame->BoxCrd().ToRecip(ucell, recip);
  // Center w.r.t. anchor
  if (useMass_)
    fcom = currentFrame->VCenterOfMass( anchorMask_ );
  else
    fcom = currentFrame->VGeometricCenter( anchorMask_ );
  if (origin_) {
    fcom.Neg(); // Shift to coordinate origin (0,0,0)
    anchorcenter.Zero();
  } else {
    if (ortho_ || truncoct_) // Center is box xyz over 2
      anchorcenter = currentFrame->BoxCrd().Center();
    else                     // Center in frac coords is (0.5,0.5,0.5)
      anchorcenter = ucell.TransposeMult(Vec3(0.5));
    fcom = anchorcenter - fcom;
  }
  currentFrame->Translate(fcom);

  // Setup imaging, and image everything in currentFrame 
  // according to mobileList. 
  if (ortho_) {
    if (Image::SetupOrtho(currentFrame->BoxCrd(), bp, bm, origin_)) {
      mprintf("Warning: autoimage: Frame %i imaging failed, box lengths are zero.\n",frameNum+1);
      // TODO: Return OK for now so next frame is tried; eventually indicate SKIP?
      return Action::OK;
    }
    Image::Ortho(*currentFrame, bp, bm, offset, usecom_, useMass_, mobileList_);
  } else {
    if (truncoct_)
      fcom = Image::SetupTruncoct( *currentFrame, 0, useMass_, origin_ );
    Image::Nonortho(*currentFrame, origin_, fcom, offset, ucell, recip, truncoct_,
                    usecom_, useMass_, mobileList_);
  }  

  // For each molecule defined by atom pairs in fixedList, determine if the
  // imaged position is closer to anchor center than the current position.
  // Always use molecule center when imaging fixedList.
  for (pairList::iterator atom1 = fixedList_.begin();
                          atom1 != fixedList_.end(); ++atom1)
  {
    int firstAtom = *atom1;
    ++atom1;
    int lastAtom = *atom1;
    if (useMass_) 
      framecenter = currentFrame->VCenterOfMass(firstAtom, lastAtom);
    else
      framecenter = currentFrame->VGeometricCenter(firstAtom, lastAtom);
    // NOTE: imaging routines will modify input coords.
    //imagedcenter[0] = framecenter[0];
    //imagedcenter[1] = framecenter[1];
    //imagedcenter[2] = framecenter[2];
    if (ortho_)
      Trans = Image::Ortho(framecenter, bp, bm, currentFrame->BoxCrd());
    else
      Trans = Image::Nonortho(framecenter, truncoct_, origin_, ucell, recip, fcom, -1.0);
    // If molecule was imaged, determine whether imaged position is closer to anchor.
    if (Trans[0] != 0 || Trans[1] != 0 || Trans[2] != 0) {
      imagedcenter = framecenter;
      imagedcenter += Trans;
      double framedist2 = DIST2_NoImage( anchorcenter, framecenter );
      double imageddist2 = DIST2_NoImage( anchorcenter, imagedcenter );
      //mprintf("DBG: [%5i] Fixed @%i-%i frame dist2=%lf, imaged dist2=%lf\n", frameNum,
      //        firstAtom+1, lastAtom+1,
      //        framedist2, imageddist2);
      if (imageddist2 < framedist2) {
        // Imaging these atoms moved them closer to anchor. Update coords in currentFrame.
        currentFrame->Translate(Trans, firstAtom, lastAtom);
        //for (int idx = firstAtom*3; idx < lastAtom*3; ++idx)
        //  (*currentFrame)[idx] = fixedFrame[idx];
      }
    }
  }
    
  return Action::OK;
}

