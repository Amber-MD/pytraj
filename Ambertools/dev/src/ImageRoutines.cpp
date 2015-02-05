#include <cmath> // floor
#include "ImageRoutines.h"
#include "DistRoutines.h"

/** Check that at least 1 atom in the range is in Mask1 */
static inline void CheckRange(Image::PairType& atomPairs, AtomMask const& MaskIn, 
                              int firstAtom, int lastAtom)
{
  bool rangeIsValid = false;
  for (int atom = firstAtom; atom < lastAtom; ++atom) {
    if (MaskIn.AtomInCharMask(atom)) {
      rangeIsValid = true;
      break;
    }
  }
  if (rangeIsValid) {
    atomPairs.push_back( firstAtom );
    atomPairs.push_back( lastAtom );
  }
}

// Image::CreatePairList() 
/** Create an atom pair list by molecule, residue, or atom.
  * NOTE: The mask is passed in by VALUE, not REFERENCE so that it can
  *       be set up and used here.
  */
Image::PairType Image::CreatePairList(Topology const& Parm, Mode modeIn, AtomMask maskIn) {
  PairType atomPairs;
  // Set up mask based on desired imaging mode.
  if ( modeIn == BYMOL || modeIn == BYRES ) {
    if ( Parm.SetupCharMask( maskIn ) ) return atomPairs;
  } else { // BYATOM
    if ( Parm.SetupIntegerMask( maskIn ) ) return atomPairs;
  }
  if (maskIn.None()) {
//    mprintf("Warning: Mask '%s' selects no atoms for topology '%s'.\n",
//            maskIn.MaskString(), Parm.c_str());
    return atomPairs;
  }
  // Set up atom range for each entity to be imaged.
  switch (modeIn) {
    case BYMOL:
      atomPairs.reserve( Parm.Nmol()*2 );
      for (Topology::mol_iterator mol = Parm.MolStart();
                                  mol != Parm.MolEnd(); ++mol)
        CheckRange( atomPairs, maskIn, (*mol).BeginAtom(), (*mol).EndAtom());
     break;
    case BYRES:
      atomPairs.reserve( Parm.Nres()*2 );
      for (Topology::res_iterator residue = Parm.ResStart();
                                  residue != Parm.ResEnd(); ++residue)
        CheckRange( atomPairs, maskIn, (*residue).FirstAtom(), (*residue).LastAtom() );
      break;
    case BYATOM:
      atomPairs.reserve( Parm.Natom()*2 );
      for (AtomMask::const_iterator atom = maskIn.begin(); atom != maskIn.end(); ++atom) {
        atomPairs.push_back(  *atom    );
        atomPairs.push_back( (*atom)+1 );
      }
      break;
  }
//  mprintf("\tNumber of %ss to be imaged is %zu based on mask '%s'\n",
//           ModeString[modeIn], atomPairs.size()/2, maskIn.MaskString());
  return atomPairs;
}

// -----------------------------------------------------------------------------
// Image::SetupTruncoct()
/** Set up centering if putting nonortho cell into familiar trunc. oct. shape.
  * \param frameIn Frame to set up for.
  * \param ComMask If not null center is calcd w.r.t. center of atoms in mask.
  * \param useMass If true calculate COM, otherwise calc geometric center.
  * \param origin If true and ComMask is null use origin, otherwise use box center.
  * \return Coordinates of center.
  */
Vec3 Image::SetupTruncoct( Frame const& frameIn, AtomMask* ComMask, bool useMass, bool origin)
{
  if (ComMask!=0) {
    // Use center of atoms in mask
    if (useMass)
      return frameIn.VCenterOfMass( *ComMask );
    else
      return frameIn.VGeometricCenter( *ComMask );
  } else if (!origin) {
    // Use box center
    return frameIn.BoxCrd().Center(); 
  }
  //fprintf(stdout,"DEBUG: fcom = %lf %lf %lf\n",fcom[0],fcom[1],fcom[2]);
  return Vec3(0.0, 0.0, 0.0); // Default is origin {0,0,0}
}

// Image::Nonortho()
/** \param frameIn Frame to image.
  * \param origin If true image w.r.t. coordinate origin.
  * \param fcom If truncoct is true, calc distance w.r.t. this coordinate.
  * \param ucell Unit cell matrix.
  * \param recip Reciprocal coordinates matrix.
  * \param truncoct If true imaging will occur using truncated octahedron shape.
  * \param center If true image w.r.t. center coords, otherwise use first atom coords.
  * \param useMass If true use COM, otherwise geometric center.
  * \param AtomPairs Atom pairs to image.
  */
void Image::Nonortho(Frame& frameIn, bool origin, Vec3 const& fcom, Vec3 const& offIn, 
                     Matrix_3x3 const& ucell, Matrix_3x3 const& recip,
                     bool truncoct, bool center,
                     bool useMass, PairType const& AtomPairs)
{
  Vec3 Coord;
  Vec3 offset = ucell.TransposeMult( offIn );
  double min = -1.0;

  if (truncoct)
    min = 100.0 * (frameIn.BoxCrd().BoxX()*frameIn.BoxCrd().BoxX()+
                   frameIn.BoxCrd().BoxY()*frameIn.BoxCrd().BoxY()+
                   frameIn.BoxCrd().BoxZ()*frameIn.BoxCrd().BoxZ());

  // Loop over atom pairs
  for (PairType::const_iterator atom = AtomPairs.begin();
                                atom != AtomPairs.end(); ++atom)
  {
    int firstAtom = *atom;
    ++atom;
    int lastAtom = *atom;
    //if (debug>2)
    //  mprintf( "  IMAGE processing atoms %i to %i\n", firstAtom+1, lastAtom);
    // Set up Coord with position to check for imaging based on first atom or 
    // center of mass of atoms first to last.
    if (center) {
      if (useMass)
        Coord = frameIn.VCenterOfMass(firstAtom,lastAtom);
      else
        Coord = frameIn.VGeometricCenter(firstAtom,lastAtom);
    } else 
      Coord = frameIn.XYZ( firstAtom );

    // boxTrans will hold calculated translation needed to move atoms back into box
    Vec3 boxTrans = Nonortho(Coord, truncoct, origin, ucell, recip, fcom, min) + offset;

    frameIn.Translate(boxTrans, firstAtom, lastAtom);

  } // END loop over atom pairs
}

// Image::Nonortho()
/** \param Coord Coordinate to image.
  * \param truncoct If true, image in truncated octahedral shape.
  * \param origin If true, image w.r.t. coordinate origin.
  * \param ucell Unit cell matrix.
  * \param recip Reciprocal coordinates matrix.
  * \param fcom If truncoct, image translated coordinate w.r.t. this coord.
  * \return Vector containing image translation.
  */
Vec3 Image::Nonortho(Vec3 const& Coord, bool truncoct, 
                     bool origin, Matrix_3x3 const& ucell, Matrix_3x3 const& recip, 
                     Vec3 const& fcom, double min)
{
  int ixyz[3];

  Vec3 fc = recip * Coord;

  if ( origin )
    fc += 0.5; 

  Vec3 boxTransOut = ucell.TransposeMult( Vec3(floor(fc[0]), floor(fc[1]), floor(fc[2])) );
  boxTransOut.Neg();

  // Put into familiar trunc. oct. shape
  if (truncoct) {
    Vec3 TransCoord = recip * (Coord + boxTransOut);
    Vec3 f2 = recip * fcom;

    if (origin) {
      TransCoord += 0.5;
      f2 += 0.5;
    }

    DIST2_ImageNonOrthoRecip(TransCoord, f2, min, ixyz, ucell);
    if (ixyz[0] != 0 || ixyz[1] != 0 || ixyz[2] != 0) {
      boxTransOut += ucell.TransposeMult( ixyz );
      //if (debug > 2)
      //  mprintf( "  IMAGING, FAMILIAR OFFSETS ARE %i %i %i\n", 
      //          ixyz[0], ixyz[1], ixyz[2]);
    }
  }
  return boxTransOut;
}

// -----------------------------------------------------------------------------
// Image::SetupOrtho()
/** \param boxIn Box coordinates of Frame to image.
  * \param bp Output: Box + boundary.
  * \param bm Output: Box - boundary.
  * \param origin If true, image w.r.t. coordinate origin, otherwise box center.
  * \return 1 if box lengths are zero, 0 if setup completed successfully.
  */
int Image::SetupOrtho(Box const& boxIn, Vec3& bp, Vec3& bm, bool origin) {
  // Set up boundary information for orthorhombic cell
  if (origin) {
    bp = boxIn.Center();
    bm.SetVec( -bp[0], -bp[1], -bp[2] );
  } else {
    bp.SetVec( boxIn.BoxX(), boxIn.BoxY(), boxIn.BoxZ()  );
    bm.Zero();
  }
  if (bp.IsZero()) return 1;
  return 0;
}

// Image::Ortho()
/** \param frameIn Frame to image.
  * \param bp Box + boundary.
  * \param bm Box - boundary.
  * \param center If true image w.r.t. center of atoms, otherwise first atom.
  * \param useMass If true calc center of mass, otherwise geometric center.
  */
void Image::Ortho(Frame& frameIn, Vec3 const& bp, Vec3 const& bm, Vec3 const& offIn,
                  bool center, bool useMass, PairType const& AtomPairs)
{
  Vec3 Coord;
  Vec3 offset(offIn[0] * frameIn.BoxCrd()[0],
              offIn[1] * frameIn.BoxCrd()[1],
              offIn[2] * frameIn.BoxCrd()[2]);
  // Loop over atom pairs
  for (PairType::const_iterator atom = AtomPairs.begin();
                                atom != AtomPairs.end(); atom++)
  {
    int firstAtom = *atom;
    ++atom;
    int lastAtom = *atom;
    //if (debug>2)
    //  mprintf( "  IMAGE processing atoms %i to %i\n", firstAtom+1, lastAtom);
    // Set up Coord with position to check for imaging based on first atom or 
    // center of mass of atoms first to last.
    if (center) {
      if (useMass)
        Coord = frameIn.VCenterOfMass(firstAtom,lastAtom);
      else
        Coord = frameIn.VGeometricCenter(firstAtom,lastAtom);
    } else 
      Coord = frameIn.XYZ( firstAtom );

    // boxTrans will hold calculated translation needed to move atoms back into box
    Vec3 boxTrans = Ortho(Coord, bp, bm, frameIn.BoxCrd()) + offset;

    // Translate atoms according to Coord
    frameIn.Translate(boxTrans, firstAtom, lastAtom);
  } // END loop over atom pairs
}

// Image::Ortho()
/** \param Coord Coordinate to image
  * \param bp Box + boundary
  * \param bm Box - boundary
  * \param BoxVec box lengths.
  * \return Vector containing image translation
  */
Vec3 Image::Ortho(Vec3 const& Coord, Vec3 const& bp, Vec3 const& bm, Box const& BoxVec)
{
  Vec3 trans;
  // Determine how far Coord is out of box
  for (int i = 0; i < 3; ++i) {
    trans[i] = 0.0;
    double crd = Coord[i];
    while (crd < bm[i]) {
      crd += BoxVec[i];
      trans[i] += BoxVec[i];
    }
    while (crd > bp[i]) {
      crd -= BoxVec[i];
      trans[i] -= BoxVec[i];
    }
  }
  return trans;
}

// -----------------------------------------------------------------------------
// Image::UnwrapNonortho()
void Image::UnwrapNonortho( Frame& tgtIn, Frame& refIn, PairType const& AtomPairs,
                            Matrix_3x3 const& ucell, Matrix_3x3 const& recip, 
                            bool center, bool useMass ) 
{
  Vec3 vtgt, vref, boxTrans;
  // Loop over atom pairs
  for (PairType::const_iterator atom = AtomPairs.begin();
                                atom != AtomPairs.end(); ++atom)
  {
    int firstAtom = *atom;
    ++atom;
    int lastAtom = *atom;
    if (center) {
      // Use center of coordinates between first and last atoms.
      if (useMass) {
        vtgt = tgtIn.VCenterOfMass(firstAtom, lastAtom);
        vref = refIn.VCenterOfMass(firstAtom, lastAtom); 
      } else {
        vtgt = tgtIn.VGeometricCenter(firstAtom, lastAtom);
        vref = refIn.VGeometricCenter(firstAtom, lastAtom);
      }
    } else {
      // Use position first atom only.
      vtgt = tgtIn.XYZ( firstAtom );
      vref = refIn.XYZ( firstAtom );
    }
    boxTrans.Zero();
    // Calculate original distance from the ref (previous) position. 
    Vec3 vd = vtgt - vref; // dx dy dz
    double minDistanceSquare = vd.Magnitude2();
    // Reciprocal coordinates
    vd = recip * vd ; // recip * dxyz
    double cx = floor(vd[0]);
    double cy = floor(vd[1]);
    double cz = floor(vd[2]);
    // Loop over all possible translations 
    for (int ix = -1; ix < 2; ++ix) {
      for (int iy = -1; iy < 2; ++iy) {
        for (int iz = -1; iz < 2; ++iz) {
          // Calculate the translation.
          Vec3 vcc = ucell.TransposeMult( Vec3( cx+(double)ix, 
                                                cy+(double)iy, 
                                                cz+(double)iz ) ); // ucell^T * ccxyz
          // Calc. the potential new coordinate for tgt
          Vec3 vnew = vtgt - vcc; 
          // Calc. the new distance from the ref (previous) position
          Vec3 vr = vref - vnew; 
          double distanceSquare = vr.Magnitude2();
          // If the orig. distance is greater than the new distance, unwrap. 
          if ( minDistanceSquare > distanceSquare ) {
              minDistanceSquare = distanceSquare;
              boxTrans = vcc;
          }
        }
      }
    }
    // Translate tgt atoms
    boxTrans.Neg();
    tgtIn.Translate( boxTrans, firstAtom, lastAtom );
    // Save new ref positions
    int i3 = firstAtom * 3;
    std::copy( tgtIn.xAddress()+i3, tgtIn.xAddress()+(lastAtom*3), refIn.xAddress()+i3 );
  } // END loop over atom pairs 
}

// Image::UnwrapOrtho()
void Image::UnwrapOrtho( Frame& tgtIn, Frame& refIn, PairType const& AtomPairs,
                         bool center, bool useMass )
{
  Vec3 vtgt, vref, boxTrans;
  Vec3 boxVec = tgtIn.BoxCrd().Lengths();
  // Loop over atom pairs
  for (PairType::const_iterator atom = AtomPairs.begin();
                                atom != AtomPairs.end(); ++atom)
  {
    int firstAtom = *atom;
    ++atom;
    int lastAtom = *atom;
    if (center) {
      // Use center of coordinates between first and last atoms.
      if (useMass) {
        vtgt = tgtIn.VCenterOfMass(firstAtom, lastAtom);
        vref = refIn.VCenterOfMass(firstAtom, lastAtom);
      } else {
        vtgt = tgtIn.VGeometricCenter(firstAtom, lastAtom);
        vref = refIn.VGeometricCenter(firstAtom, lastAtom);
      }
    } else {
      // Use position first atom only.
      vtgt = tgtIn.XYZ( firstAtom );
      vref = refIn.XYZ( firstAtom );
    }
    Vec3 dxyz = vtgt - vref;
    boxTrans[0] = -floor( dxyz[0] / boxVec[0] + 0.5 ) * boxVec[0];
    boxTrans[1] = -floor( dxyz[1] / boxVec[1] + 0.5 ) * boxVec[1];
    boxTrans[2] = -floor( dxyz[2] / boxVec[2] + 0.5 ) * boxVec[2];
    // Translate atoms from first to last
    tgtIn.Translate(boxTrans, firstAtom, lastAtom);
    // Save new ref positions
    int i3 = firstAtom * 3;
    std::copy( tgtIn.xAddress()+i3, tgtIn.xAddress()+(lastAtom*3), refIn.xAddress()+i3 );
  } // END loop over atom pairs
}
