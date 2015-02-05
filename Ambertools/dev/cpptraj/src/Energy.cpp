#include <cmath>
#include "Energy.h"
#include "CpptrajStdio.h"
#include "DistRoutines.h"
#include "TorsionRoutines.h"
#include "Constants.h"

const double Energy_Amber::QFAC = Constants::ELECTOAMBER * Constants::ELECTOAMBER;

// CONSTRUCTOR
Energy_Amber::Energy_Amber() : debug_(0) {}

/** Bond energy */
double Energy_Amber::E_bond(Frame const& fIn, Topology const& tIn, AtomMask const& mask)
{
  time_bond_.Start();
  // Heavy atom bonds
  double Ebond = CalcBondEnergy(fIn, tIn.Bonds(), tIn.BondParm(), mask);
  Ebond += CalcBondEnergy(fIn, tIn.BondsH(), tIn.BondParm(), mask);
  time_bond_.Stop();
  return Ebond;
}

// Energy_Amber::CalcBondEnergy()
double Energy_Amber::CalcBondEnergy(Frame const& fIn, BondArray const& Bonds,
                                    BondParmArray const& BPA, AtomMask const& mask)
{
  double Ebond = 0.0;
  for (BondArray::const_iterator b = Bonds.begin(); b != Bonds.end(); ++b)
  {
    if (mask.AtomInCharMask( b->A1() ) && mask.AtomInCharMask( b->A2() ))
    {
      int bpidx = b->Idx();
      if (bpidx < 0) {
        if (debug_ > 0)
          mprintf("Warning: Bond %i -- %i has no parameters.\n", b->A1()+1, b->A2()+1);
        continue;
      }
      BondParmType const& bp = BPA[bpidx];
      double r2 = DIST2_NoImage( fIn.XYZ(b->A1()), fIn.XYZ(b->A2()) );
      double r = sqrt(r2);
      double rdiff = r - bp.Req();
      double ene = bp.Rk() * (rdiff * rdiff);
      Ebond += ene;
#     ifdef DEBUG_ENERGY
      mprintf("\tBond %4u %4i -- %4i: k= %12.5f  x0= %12.5f  r= %12.5f  E= %12.5e\n",
              b - Bonds.begin(), b->A1()+1, b->A2()+1, bp.Rk(), bp.Req(), r, ene);
#     endif
    }
  }
  return Ebond;
}

// -----------------------------------------------------------------------------
/** Angle energy */
double Energy_Amber::E_angle(Frame const& fIn, Topology const& tIn, AtomMask const& mask)
{
  time_angle_.Start();
  // Heavy atom angles
  double Eang = CalcAngleEnergy(fIn, tIn.Angles(), tIn.AngleParm(), mask);
  Eang += CalcAngleEnergy(fIn, tIn.AnglesH(), tIn.AngleParm(), mask);
  time_angle_.Stop();
  return Eang;
}

// Energy_Amber::CalcAngleEnergy()
double Energy_Amber::CalcAngleEnergy(Frame const& fIn, AngleArray const& Angles,
                                     AngleParmArray const& APA, AtomMask const& mask)
{
  double Eangle = 0.0;
  for (AngleArray::const_iterator a = Angles.begin(); a != Angles.end(); ++a)
  {
    if (mask.AtomInCharMask(a->A1()) &&
        mask.AtomInCharMask(a->A2()) &&
        mask.AtomInCharMask(a->A3()))
    {
      int apidx = a->Idx();
      if (apidx < 0) {
        if (debug_ > 0)
          mprintf("Warning: Angle %i -- %i -- %i has no parameters.\n",
                   a->A1()+1, a->A2()+1, a->A3()+1);
        continue;
      }
      AngleParmType const& ap = APA[apidx];
      double theta = CalcAngle(fIn.XYZ(a->A1()), fIn.XYZ(a->A2()), fIn.XYZ(a->A3()));
      double tdiff = theta - ap.Teq();
      double ene = ap.Tk() * (tdiff * tdiff);
      Eangle += ene;
#     ifdef DEBUG_ENERGY
      mprintf("\tAngle %4u %4i -- %4i -- %4i: k= %12.5f  x0= %12.5f  t= %12.5f  E= %12.5e\n",
              a - Angles.begin(), a->A1()+1, a->A2()+1, a->A3()+1, ap.Tk(), ap.Teq(), theta, ene);
#     endif
    }
  }
  return Eangle;
}

// -----------------------------------------------------------------------------
/** Dihedral energy */
double Energy_Amber::E_torsion(Frame const& fIn, Topology const& tIn, AtomMask const& mask)
{
  time_tors_.Start();
  // Heavy atom dihedrals
  double Edih = CalcTorsionEnergy(fIn, tIn.Dihedrals(), tIn.DihedralParm(), mask);
  Edih += CalcTorsionEnergy(fIn, tIn.DihedralsH(), tIn.DihedralParm(), mask);
  time_tors_.Stop();
  return Edih;
}

// Energy_Amber::CalcTorsionEnergy()
double Energy_Amber::CalcTorsionEnergy(Frame const& fIn, DihedralArray const& Dihedrals,
                                       DihedralParmArray const& DPA, AtomMask const& mask)
{
  double Edih = 0.0;
  for (DihedralArray::const_iterator d = Dihedrals.begin(); d != Dihedrals.end(); d++)
  {
    if (mask.AtomInCharMask(d->A1()) &&
        mask.AtomInCharMask(d->A2()) &&
        mask.AtomInCharMask(d->A3()) &&
        mask.AtomInCharMask(d->A4()))
    {
      int dpidx = d->Idx();
      if (dpidx < 0) {
        if (debug_ > 0)
          mprintf("Warning: Dihedral %i -- %i -- %i -- %i has no parameters.\n",
                   d->A1()+1, d->A2()+1, d->A3()+1, d->A4()+1);
        continue;
      }
      DihedralParmType const& dp = DPA[dpidx];
      double phi = Torsion(fIn.XYZ(d->A1()), fIn.XYZ(d->A2()),
                           fIn.XYZ(d->A3()), fIn.XYZ(d->A4()));
      double ene = dp.Pk() * (1.0 + cos(dp.Pn() * phi - dp.Phase()));
      Edih += ene;
#     ifdef DEBUG_ENERGY
      mprintf("\tDihedral %4u %4i -- %4i -- %4i -- %4i: pk= %12.5f  "
                "pn= %12.5f  phase= %12.5f  p= %12.5f  E= %12.5e\n",
              d - Dihedrals.begin(), d->A1()+1, d->A2()+1, d->A3()+1, d->A4()+1,
              dp.Pk(), dp.Pn(), dp.Phase(), phi, ene);
#     endif
    }
  }
  return Edih;
}

// -----------------------------------------------------------------------------
/** 1-4 nonbond energy */
double Energy_Amber::E_14_Nonbond(Frame const& fIn, Topology const& tIn, AtomMask const& mask,
                                  double& Eq14)
{
  time_14_.Start();
  Eq14 = 0.0;
  // Heavy atom dihedrals
  double Evdw14 = Calc_14_Energy(fIn, tIn.Dihedrals(), tIn.DihedralParm(), tIn, mask, Eq14);
  Evdw14 += Calc_14_Energy(fIn, tIn.DihedralsH(), tIn.DihedralParm(), tIn, mask, Eq14);
  time_14_.Stop();
  return Evdw14;
}

// Energy_Amber::Calc_14_Energy()
double Energy_Amber::Calc_14_Energy(Frame const& fIn, DihedralArray const& Dihedrals,
                                    DihedralParmArray const& DPA, Topology const& tIn,
                                    AtomMask const& mask, double& Eq14)
{
  double Evdw14 = 0.0;
  for (DihedralArray::const_iterator d = Dihedrals.begin(); d != Dihedrals.end(); d++)
  {
    if (d->Type() == DihedralType::NORMAL &&
        mask.AtomInCharMask(d->A1()) && mask.AtomInCharMask(d->A4()))
    {
      int dpidx = d->Idx();
      if (dpidx < 0) {
        if (debug_ > 0)
          mprintf("Warning: 1-4 pair %i -- %i has no parameters.\n", d->A1()+1, d->A4()+1);
        continue;
      }
      DihedralParmType const& dp = DPA[dpidx];
      double rij2 = DIST2_NoImage( fIn.XYZ(d->A1()), fIn.XYZ(d->A4()) );
      double rij = sqrt( rij2 );
      // VDW
      NonbondType const& LJ = tIn.GetLJparam(d->A1(), d->A4());
      double r2    = 1.0 / rij2;
      double r6    = r2 * r2 * r2;
      double r12   = r6 * r6;
      double f12   = LJ.A() * r12;  // A/r^12
      double f6    = LJ.B() * r6;   // B/r^6
      double e_vdw = f12 - f6;      // (A/r^12)-(B/r^6)
      e_vdw /= dp.SCNB();
      Evdw14 += e_vdw;
      // Coulomb
      double qiqj = QFAC * tIn[d->A1()].Charge() * tIn[d->A4()].Charge();
      double e_elec = qiqj / rij;
      e_elec /= dp.SCEE();
      Eq14 += e_elec;
#     ifdef DEBUG_ENERGY
      mprintf("\tEVDW14  %4i -- %4i: A=  %12.5e  B=  %12.5e  r2= %12.5f  E= %12.5e\n",
              d->A1()+1, d->A4()+1, LJ.A(), LJ.B(), rij2, e_vdw);
      mprintf("\tEELEC14 %4i -- %4i: q1= %12.5e  q2= %12.5e  r=  %12.5f  E= %12.5e\n",
              d->A1()+1, d->A4()+1, tIn[d->A1()].Charge(), tIn[d->A4()].Charge(),
              rij, e_elec);
#     endif
    }
  }
  return Evdw14;
}

// -----------------------------------------------------------------------------
// Energy_Amber::E_Nonbond()
double Energy_Amber::E_Nonbond(Frame const& fIn, Topology const& tIn, AtomMask const& mask,
                               double& Eelec)
{
  time_NB_.Start();
  double Evdw = 0.0;
  Eelec = 0.0;
  for (AtomMask::const_iterator maskatom1 = mask.begin();
                                maskatom1 != mask.end();
                                ++maskatom1)
  {
    // Set up coord for this atom
    const double* crd1 = fIn.XYZ( *maskatom1 );
    // Set up exclusion list for this atom
    Atom::excluded_iterator excluded_atom = tIn[*maskatom1].excludedbegin();
    for (AtomMask::const_iterator maskatom2 = maskatom1 + 1;
                                  maskatom2 != mask.end();
                                ++maskatom2)
    {
      // If atom is excluded, just increment to next excluded atom.
      if (excluded_atom != tIn[*maskatom1].excludedend() &&
          *maskatom2 == *excluded_atom)
      {
        ++excluded_atom;
      }
      else
      {
        // TODO: imaged distance
        double rij2 = DIST2_NoImage( crd1, fIn.XYZ( *maskatom2 ) );
        double rij = sqrt( rij2 );
        // VDW
        NonbondType const& LJ = tIn.GetLJparam(*maskatom1, *maskatom2);
        double r2    = 1.0 / rij2;
        double r6    = r2 * r2 * r2;
        double r12   = r6 * r6;
        double f12   = LJ.A() * r12;  // A/r^12
        double f6    = LJ.B() * r6;   // B/r^6
        double e_vdw = f12 - f6;      // (A/r^12)-(B/r^6)
        Evdw += e_vdw;
        // Coulomb
        double qiqj = QFAC * tIn[*maskatom1].Charge() * tIn[*maskatom2].Charge();
        double e_elec = qiqj / rij;
        Eelec += e_elec;
#       ifdef DEBUG_ENERGY
        mprintf("\tEVDW  %4i -- %4i: A=  %12.5e  B=  %12.5e  r2= %12.5f  E= %12.5e\n",
                *maskatom1+1, *maskatom2+1, LJ.A(), LJ.B(), rij2, e_vdw);
        mprintf("\tEELEC %4i -- %4i: q1= %12.5e  q2= %12.5e  r=  %12.5f  E= %12.5e\n",
                *maskatom1, *maskatom2, tIn[*maskatom1].Charge(), tIn[*maskatom2].Charge(),
                rij, e_elec);
#       endif
      }
    }
  }
  time_NB_.Stop();
  return Evdw;
}

void Energy_Amber::PrintTiming() const {
  double total = time_bond_.Total() + time_angle_.Total() +
                 time_tors_.Total() + time_14_.Total() +
                 time_NB_.Total();
  time_bond_.WriteTiming(1, "BOND", total);
  time_angle_.WriteTiming(1, "ANGLE", total);
  time_tors_.WriteTiming(1, "TORSION", total);
  time_14_.WriteTiming(1, "1-4_NONBOND", total);
  time_NB_.WriteTiming(1, "NONBOND", total);
  mprintf("TIME: Total= %.4f\n", total);
}
