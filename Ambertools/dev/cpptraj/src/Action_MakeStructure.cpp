#include <cctype> // islpha
#include "Action_MakeStructure.h"
#include "CpptrajStdio.h"
#include "Constants.h" // DEGRAD
#include "TorsionRoutines.h"
#include "StringRoutines.h" // convertToDouble, convertToInteger

// CONSTRUCTOR
Action_MakeStructure::Action_MakeStructure() : debug_(0) {
  // Initially known structure types. 
  SS.push_back(SS_TYPE(  -57.8,  -47.0,    0.0,   0.0, 0, "alpha"    ));
  SS.push_back(SS_TYPE(   57.8,   47.0,    0.0,   0.0, 0, "left"     ));
  SS.push_back(SS_TYPE(  -75.0,  145.0,    0.0,   0.0, 0, "pp2"      ));
  SS.push_back(SS_TYPE( -100.0,  130.0,    0.0,   0.0, 0, "hairpin"  ));
  SS.push_back(SS_TYPE( -150.0,  155.0,    0.0,   0.0, 0, "extended" ));
  SS.push_back(SS_TYPE(  -60.0,  -30.0,  -90.0,   0.0, 1, "typeI"    )); 	
  SS.push_back(SS_TYPE(  -60.0,  120.0,   80.0,   0.0, 1, "typeII"   ));	
  SS.push_back(SS_TYPE(  -60.0,  -30.0, -120.0, 120.0, 1, "typeVIII" )); 	
  SS.push_back(SS_TYPE(   60.0,   30.0,   90.0,   0.0, 1, "typeI'"   )); 	
  SS.push_back(SS_TYPE(   60.0, -120.0,  -80.0,   0.0, 1, "typeII'"  )); 	
  SS.push_back(SS_TYPE(  -60.0,  120.0,  -90.0,   0.0, 1, "typeVIa1" )); // 2nd res must be cis-PRO
  SS.push_back(SS_TYPE( -120.0,  120.0,  -60.0,   0.0, 1, "typeVIa2" )); // 2nd res must be cis-PRO
  SS.push_back(SS_TYPE( -135.0,  135.0,  -75.0, 160.0, 1, "typeVIb"  )); // 2nd res must be cis-PRO
}

/// Used to indicate type not found
#define SS_EMPTY -1
// Action_MakeStructure::FindSStype()
int Action_MakeStructure::FindSStype(std::string const& typeIn) 
{
  for (unsigned int sstype = 0; sstype < SS.size(); ++sstype)
    if (typeIn == SS[sstype].type_arg)
      return (int)sstype;
  return SS_EMPTY;
} 

void Action_MakeStructure::Help() {
  mprintf("\t<List of Args>\n"
          "  Apply dihedrals to specified residues using arguments found in <List of Args>,\n"
          "  where an argument is 1 or more of the following arg types:\n"
          "\t'<sstype>:<res range>' Apply SS type (phi/psi) to residue range.\n"
          "\t\t<sstype> standard = alpha, left, pp2, hairpin, extended\n"
          "\t\t<sstype> turn = typeI, typeII, typeVIII, typeI', typeII,\n"
          "\t\t                typeVIa1, typeVIa2, typeVIb\n"
          "\t\tTurns are applied to 2 residues at a time, so resrange must be divisible by 4.\n"
          "\t'<custom ss>:<res range>:<phi>:<psi>' Apply custom <phi>/<psi> to residue range.\n"
          "\t'<custom turn>:<res range>:<phi1>:<psi1>:<phi2>:<psi2>' Apply custom turn <phi>/<psi> pair to residue range.\n"
          "\t'<custom dih>:<res range>:<dih type>:<angle>' Apply <angle> to dihedrals in range.\n"
          "\t\t<dih type> =");
  DihedralSearch::ListKnownTypes();
  mprintf("\t'<custom dih>:<res range>:<at0>:<at1>:<at2>:<at3>:<angle>[:<offset>]' Apply <angle> to dihedral defined by atoms <at1>, <at2>, <at3>, and <at4>.\n");
  DihedralSearch::OffsetHelp();
  mprintf("\t'ref:<range>:<refname>[:<ref range>]' Apply dihedrals from reference <refname>.\n");
}

// Action_MakeStructure::Init()
Action::RetType Action_MakeStructure::Init(ArgList& actionArgs, TopologyList* PFL, DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  debug_ = debugIn;
  secstruct_.clear();
  // Get all arguments 
  std::string ss_expr = actionArgs.GetStringNext();
  while ( !ss_expr.empty() ) {
    ArgList ss_arg(ss_expr, ":");
    if (ss_arg.Nargs() < 2) {
      mprinterr("Error: Malformed SS arg.\n");
      Help();
      return Action::ERR;
    }
    // Type is 1st arg, range is 2nd arg.
    SecStructHolder ss_holder(ss_arg[1], FindSStype(ss_arg[0]));

    if (ss_arg.Nargs() == 2) {
    // Find SS type: <ss type>:<range>
      if (ss_holder.sstype_idx == SS_EMPTY) {
        mprinterr("Error: SS type %s not found.\n", ss_arg[0].c_str());
        return Action::ERR;
      } 
      ss_holder.dihSearch_.SearchFor(DihedralSearch::PHI);
      ss_holder.dihSearch_.SearchFor(DihedralSearch::PSI);
      secstruct_.push_back( ss_holder );

    } else if (ss_arg[0] == "ref") {
    // Use dihedrals from reference structure
      if (ss_arg.Nargs() < 3) {
        mprinterr("Error: Invalid 'ref' arg. Requires 'ref:<range>:<refname>[:<ref range>]'\n");
        return Action::ERR;
      }
      ss_arg.MarkArg(0);
      ss_arg.MarkArg(1);
      // Sanity check: Currently only unique args of this type are allowed
      if (ss_holder.sstype_idx != SS_EMPTY) {
        mprinterr("Error: Ref backbone types must be unique [%s]\n", ss_arg[0].c_str());
        return Action::ERR;
      }
      // Use backbone phi/psi from reference structure
      ss_holder.dihSearch_.SearchFor(DihedralSearch::PHI);
      ss_holder.dihSearch_.SearchFor(DihedralSearch::PSI);
      // Get reference structure
      DataSet_Coords_REF* REF = (DataSet_Coords_REF*)
                                DSL->GetReferenceFrame(ss_arg.GetStringNext()); // ss_arg[2]
      if (REF == 0) {
        mprinterr("Error: Could not get reference structure [%s]\n", ss_arg[2].c_str());
        return Action::ERR;
      }
      // Get reference residue range, or use resRange
      Range refRange(ss_arg.GetStringNext(), -1); // ss_arg[3]
      if (!refRange.Empty()) {
        if (ss_holder.resRange.Size() != refRange.Size()) {
          mprinterr("Error: Reference range [%s] must match residue range [%s]\n",
                    refRange.RangeArg(), ss_holder.resRange.RangeArg());
          return Action::ERR;
        }
      } else
        refRange = ss_holder.resRange;
      // Look for phi/psi only in reference
      DihedralSearch refSearch;
      refSearch.SearchFor(DihedralSearch::PHI);
      refSearch.SearchFor(DihedralSearch::PSI);
      if (refSearch.FindDihedrals( REF->Top(), refRange )) return Action::ERR;
      // For each found dihedral, set theta 
      for (DihedralSearch::mask_it dih = refSearch.begin(); dih != refSearch.end(); ++dih)
      {
        double torsion = Torsion( REF->RefFrame().XYZ(dih->A0()),
                                  REF->RefFrame().XYZ(dih->A1()),
                                  REF->RefFrame().XYZ(dih->A2()),
                                  REF->RefFrame().XYZ(dih->A3()) );
        ss_holder.thetas_.push_back( (float)torsion );
      }
      secstruct_.push_back( ss_holder );

    } else if (ss_arg.Nargs() == 4 && isalpha(ss_arg[2][0])) {
    // Single dihedral type: <name>:<range>:<dih type>:<angle>
      DihedralSearch::DihedralType dtype = DihedralSearch::GetType(ss_arg[2]);
      if (ss_holder.sstype_idx == SS_EMPTY) {
        // Type not yet defined. Create new type. 
        if (dtype == DihedralSearch::NDIHTYPE) {
          mprinterr("Error: Dihedral type %s not found.\n", ss_arg[2].c_str());
          return Action::ERR;
        }
        if (!validDouble(ss_arg[3])) {
          mprinterr("Error: 4th arg (angle) is not a valid number.\n");
          return Action::ERR;
        }
        SS.push_back( SS_TYPE(convertToDouble(ss_arg[3]), 0.0, 0.0, 0.0, 2, ss_arg[0]) );
        ss_holder.sstype_idx = (int)(SS.size() - 1);
      }
      ss_holder.dihSearch_.SearchFor( dtype ); 
      secstruct_.push_back( ss_holder );

    } else if (ss_arg.Nargs() == 7 || ss_arg.Nargs() == 8) {
    // Single custom dihedral type: <name>:<range>:<at0>:<at1>:<at2>:<at3>:<angle>[:<offset>]
      if (ss_holder.sstype_idx == SS_EMPTY) {
        // Type not yet defined. Create new type.
        if (!validDouble(ss_arg[6])) {
          mprinterr("Error: 7th arg (angle) is not a valid number.\n");
          return Action::ERR;
        }
        SS.push_back( SS_TYPE(convertToDouble(ss_arg[6]), 0.0, 0.0, 0.0, 2, ss_arg[0]) );
        ss_holder.sstype_idx = (int)(SS.size() - 1);
      }
      int offset = 0;
      if (ss_arg.Nargs() == 8) {
        if (!validInteger(ss_arg[7])) {
          mprinterr("Error: 8th arg (offset) is not a valid number.\n");
          return Action::ERR;
        }
        offset = convertToInteger(ss_arg[7]);
      }
      ss_holder.dihSearch_.SearchForNewType(offset,ss_arg[2],ss_arg[3],ss_arg[4],ss_arg[5],
                                            ss_arg[0]);
      secstruct_.push_back( ss_holder );

    } else if (ss_arg.Nargs() == 4 || ss_arg.Nargs() == 6) {
    // Custom SS/turn type: <name>:<range>:<phi1>:<psi1>[:<phi2>:<psi2>]
      if (ss_holder.sstype_idx == SS_EMPTY) {
        // Type not yet defined. Create new type.
        if (!validDouble(ss_arg[2]) || !validDouble(ss_arg[3])) {
          mprinterr("Error: 3rd or 4th arg (phi1/psi1) is not a valid number.\n");
          return Action::ERR;
        }
        double phi1 = convertToDouble(ss_arg[2]);
        double psi1 = convertToDouble(ss_arg[3]);
        int isTurn = 0;
        double phi2 = 0.0;
        double psi2 = 0.0;
        if (ss_arg.Nargs() == 6) {
          isTurn = 1;
          if (!validDouble(ss_arg[4]) || !validDouble(ss_arg[5])) {
            mprinterr("Error: 5th or 6th arg (phi2/psi2) is not a valid number.\n");
            return Action::ERR;
          }
          phi2 = convertToDouble(ss_arg[4]);
          psi2 = convertToDouble(ss_arg[5]);
        }
        SS.push_back(SS_TYPE(phi1, psi1, phi2, psi2, isTurn, ss_arg[0] ));
        ss_holder.sstype_idx = (int)(SS.size() - 1);
      }
      ss_holder.dihSearch_.SearchFor(DihedralSearch::PHI);
      ss_holder.dihSearch_.SearchFor(DihedralSearch::PSI);
      secstruct_.push_back( ss_holder );

    } else {
      mprinterr("Error: SS arg type [%s] not recognized.\n", ss_arg[0].c_str());
      return Action::ERR;
    }
    ss_expr = actionArgs.GetStringNext();
  } // End loop over args
  if (secstruct_.empty()) {
    mprinterr("Error: No SS types defined.\n");
    return Action::ERR;
  }
  mprintf("    MAKESTRUCTURE:\n");
  for (std::vector<SecStructHolder>::iterator ss = secstruct_.begin();
                                              ss != secstruct_.end(); ++ss)
  {
    if ((*ss).sstype_idx != SS_EMPTY) {
      const SS_TYPE& myType = SS[(*ss).sstype_idx];
      switch ( myType.isTurn ) {
        case 0:
          mprintf("\tSS type %s will be applied to residue(s) %s\n",
                 myType.type_arg.c_str(), (*ss).resRange.RangeArg());
          break;
        case 1:
          mprintf("\tTurn type %s will be applied to residue(s) %s\n",
                  myType.type_arg.c_str(), (*ss).resRange.RangeArg());
          break;
        case 2:
          mprintf("\tDihedral value of %.2f will be applied to %s dihedrals in residue(s) %s\n",
                  myType.phi, myType.type_arg.c_str(), (*ss).resRange.RangeArg());
      }
    } else 
      mprintf("\tBackbone angles from reference will be applied to residue(s) %s\n",
              (*ss).resRange.RangeArg());
  }
  return Action::OK;
}

// Action_MakeStructure::Setup()
Action::RetType Action_MakeStructure::Setup(Topology* currentParm, Topology** parmAddress) {
  // Set up each SS type
  for (std::vector<SecStructHolder>::iterator ss = secstruct_.begin();
                                              ss != secstruct_.end(); ++ss)
  {
    if ((*ss).dihSearch_.FindDihedrals(*currentParm, (*ss).resRange))
      return Action::ERR;
    if ( (*ss).sstype_idx != SS_EMPTY) {
      SS_TYPE const& myType = SS[(*ss).sstype_idx];
      mprintf("\tResRange=[%s] Type=%s, %i dihedrals", (*ss).resRange.RangeArg(), 
              myType.type_arg.c_str(), (*ss).dihSearch_.Ndihedrals());
      // Set up found dihedrals 
      // TODO: Check that # dihedrals is multiple of 2?
      (*ss).Rmasks_.clear();
      (*ss).thetas_.clear();
      if (myType.isTurn == 1) {
        // Turn: Require phi/psi residue pairs; must be multiple of 4
        if ( ((*ss).dihSearch_.Ndihedrals() % 4) != 0) {
          mprintf("Error: Assigning turn SS requires residue phi/psi pairs.\n");
          return Action::ERR;
        }
        // Turns require that each pair of residues is consecutive
        for (DihedralSearch::mask_it dih = (*ss).dihSearch_.begin();
                                     dih != (*ss).dihSearch_.end(); ++dih)
        {
          // First has to be phi
          if (debug_>0) mprintf(" %i:%s", (*dih).ResNum()+1, (*dih).Name().c_str());
          if ( (*dih).Name() != "phi" ) {
            mprinterr("Error: Assigning turn SS requires 1st dihedral be phi.\n");
            return Action::ERR;
          }
          int res1num = (*dih).ResNum();
          (*ss).thetas_.push_back((float)(myType.phi * Constants::DEGRAD));
          (*ss).Rmasks_.push_back( DihedralSearch::MovingAtoms(*currentParm, 
                                   (*dih).A1(), (*dih).A2()) );
          ++dih;
          // Second has to be psi and +0
          if (debug_>0) mprintf("-%i:%s", (*dih).ResNum()+1, (*dih).Name().c_str());
          if ((*dih).Name() != "psi" || ((*dih).ResNum() != res1num)) {
            mprinterr("Error: Assigning turn SS requires 2nd dihedral be psi and consecutive.\n");
            return Action::ERR;
          }
          (*ss).thetas_.push_back((float)(myType.psi * Constants::DEGRAD));
          (*ss).Rmasks_.push_back( DihedralSearch::MovingAtoms(*currentParm, 
                                   (*dih).A1(), (*dih).A2()) );
          ++dih;
          // Third has to be phi and +1
          if (debug_>0) mprintf("-%i:%s", (*dih).ResNum()+1, (*dih).Name().c_str());
          if ((*dih).Name() != "phi" || ((*dih).ResNum() - res1num) != 1) {
            mprinterr("Error: Assigning turn SS requires 3rd dihedral be phi and consecutive.\n");
            return Action::ERR;
          }
          (*ss).thetas_.push_back((float)(myType.phi2 * Constants::DEGRAD));
          (*ss).Rmasks_.push_back( DihedralSearch::MovingAtoms(*currentParm, 
                                   (*dih).A1(), (*dih).A2()) );
          ++dih;
          // Fourth has to be phi and +1
          if (debug_>0) mprintf("-%i:%s", (*dih).ResNum()+1, (*dih).Name().c_str());
          if ((*dih).Name() != "psi" || ((*dih).ResNum() - res1num) != 1) {
            mprinterr("Error: Assigning turn SS requires 4th dihedral be psi and consecutive.\n");
            return Action::ERR;
          }
          (*ss).thetas_.push_back((float)(myType.psi2 * Constants::DEGRAD));
          (*ss).Rmasks_.push_back( DihedralSearch::MovingAtoms(*currentParm, 
                                   (*dih).A1(), (*dih).A2()) );
        }
      } else if (myType.isTurn == 2) {
        // Dihedrals of a single type
        for (DihedralSearch::mask_it dih = (*ss).dihSearch_.begin();
                                     dih != (*ss).dihSearch_.end(); ++dih)
        {
          if (debug_>0) mprintf(" %i:%s", (*dih).ResNum()+1, (*dih).Name().c_str());
          (*ss).thetas_.push_back((float)(myType.phi * Constants::DEGRAD));
          (*ss).Rmasks_.push_back( DihedralSearch::MovingAtoms(*currentParm, 
                                   (*dih).A1(), (*dih).A2()) );
        }
      } else if (myType.isTurn == 0) {
        // Assign SS.
        for (DihedralSearch::mask_it dih = (*ss).dihSearch_.begin();
                                     dih != (*ss).dihSearch_.end(); ++dih)
        {
          if (debug_>0) mprintf(" %i:%s", (*dih).ResNum()+1, (*dih).Name().c_str());
          if ((*dih).Name() == "phi") {
            (*ss).thetas_.push_back((float)(myType.phi * Constants::DEGRAD));
            (*ss).Rmasks_.push_back( DihedralSearch::MovingAtoms(*currentParm, 
                                     (*dih).A1(), (*dih).A2()) );
          } else {
            (*ss).thetas_.push_back((float)(myType.psi * Constants::DEGRAD));
            (*ss).Rmasks_.push_back( DihedralSearch::MovingAtoms(*currentParm, 
                                     (*dih).A1(), (*dih).A2()) );
          }
        }
      }
    } else {
      // type was empty. Make sure the number of found dihedrals equals
      // number of reference dihedrals.
      if ( (*ss).dihSearch_.Ndihedrals() != (int)(*ss).thetas_.size() ) {
        mprinterr("Error: Number of found dihedrals (%i) != number reference dihedrals (%u)\n",
                  (*ss).dihSearch_.Ndihedrals(), (*ss).thetas_.size());
        return Action::ERR;
      }
      std::vector<float>::const_iterator theta = (*ss).thetas_.begin();
      for (DihedralSearch::mask_it dih = (*ss).dihSearch_.begin();
                                   dih != (*ss).dihSearch_.end(); ++dih)
      {
        mprintf("\tDihedral in residue %i = %f\n", (*dih).ResNum()+1, *(theta++)*Constants::RADDEG);
        (*ss).Rmasks_.push_back( DihedralSearch::MovingAtoms(*currentParm,
                                 (*dih).A1(), (*dih).A2()) );
      }
    }
    mprintf("\n");
  } // END loop over SS types
  //CurrentParm_ = currentParm; // DEBUG
  return Action::OK;
}

// Action_MakeStructure::DoAction()
Action::RetType Action_MakeStructure::DoAction(int frameNum, Frame* currentFrame, 
                                               Frame** frameAddress) 
{
  Matrix_3x3 rotationMatrix;
  for (std::vector<SecStructHolder>::iterator ss = secstruct_.begin();
                                              ss != secstruct_.end(); ++ss)
  {
    std::vector<float>::iterator theta = (*ss).thetas_.begin();
    std::vector<AtomMask>::iterator Rmask = (*ss).Rmasks_.begin();
    for (DihedralSearch::mask_it dih = (*ss).dihSearch_.begin();
                                 dih != (*ss).dihSearch_.end(); ++dih, ++theta, ++Rmask)
    {
      double theta_in_radians = (double)*theta;
      // Calculate current value of dihedral
      double torsion = Torsion( currentFrame->XYZ( (*dih).A0() ),
                                currentFrame->XYZ( (*dih).A1() ),
                                currentFrame->XYZ( (*dih).A2() ),
                                currentFrame->XYZ( (*dih).A3() ) );
      // Calculate delta needed to get to theta
      double delta = theta_in_radians - torsion;
      // Set axis of rotation
      Vec3 axisOfRotation = currentFrame->SetAxisOfRotation((*dih).A1(), (*dih).A2());
      // Calculate rotation matrix for delta 
      rotationMatrix.CalcRotationMatrix(axisOfRotation, delta);
      if (debug_ > 0) {
        std::string a0name = CurrentParm_->TruncResAtomName( (*dih).A0() );
        std::string a1name = CurrentParm_->TruncResAtomName( (*dih).A1() );
        std::string a2name = CurrentParm_->TruncResAtomName( (*dih).A2() );
        std::string a3name = CurrentParm_->TruncResAtomName( (*dih).A3() );
          mprintf("\tRotating Dih %i:%s (%i-%i-%i-%i) (@%.2f) by %.2f deg to get to %.2f.\n",
                  (*dih).ResNum()+1, (*dih).Name().c_str(),
                  (*dih).A0() + 1, (*dih).A1() + 1, (*dih).A2() + 1, (*dih).A3() + 1, 
                  torsion*Constants::RADDEG, delta*Constants::RADDEG, theta_in_radians*Constants::RADDEG);
      }
      // Rotate around axis
      currentFrame->Rotate(rotationMatrix, *Rmask);
    }
  }
  return Action::OK;
}
#undef SS_EMPTY
