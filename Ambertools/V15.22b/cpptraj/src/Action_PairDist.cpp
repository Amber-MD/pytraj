#include <cmath>

#include "Action_PairDist.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h"
#include "DistRoutines.h"



/** Calculate pair distribution function P(r) between two masks.
  * \author Hannes H. Loeffler.
  */

// CONSTRUCTOR
Action_PairDist::Action_PairDist() :
  delta_(0.0),
  maxbin_(0),
  same_mask_(false),
  ub1_(0),
  ub2_(0)
{}

void Action_PairDist::Help()
{
  mprintf("\tout <filename> mask <mask> [mask2 <mask>] [<resolution> delta]\n"
          "  Calculate pair distribution function P(r) between two masks.\n");
}

// Action_PairDist::init()
Action::RetType Action_PairDist::Init(ArgList& actionArgs,
				      TopologyList* PFL, FrameList* FL,
				      DataSetList* DSL, DataFileList* DFL,
				      int debugIn)
{
  InitImaging(true);

  std::string outfileName = actionArgs.GetStringKey("out");

  if (outfileName.empty()) {
    outfileName = "pairdist.dat";
  }

  if (output_.OpenEnsembleWrite(outfileName, DSL->EnsembleNum()) ) {
    mprinterr("Error: PairDist: Could not open output file %s\n",
	      outfileName.c_str());
    return Action::ERR;
  }

  std::string mask1 = actionArgs.GetStringKey("mask");

  if (mask1.empty()) {
    mprinterr("Error: pairdist: No mask1 specified.\n");
    return Action::ERR;
  }

  mask1_.SetMaskString(mask1);

  std::string mask2 = actionArgs.GetStringKey("mask2");

  if (mask2.empty()) {
    same_mask_ = true;
    mask2_.SetMaskString(mask1);
  } else {
    mask2_.SetMaskString(mask2);

    if (mask1_.MaskExpression() != mask2_.MaskExpression() )
      same_mask_ = false;
    else
      same_mask_ = true;
  }

  delta_ = actionArgs.getKeyDouble("delta", 0.01);


  return Action::OK;
}


// Action_PairDist::Setup()
Action::RetType Action_PairDist::Setup(Topology* currentParm,
				       Topology** parmAddress)
{
  if (currentParm->SetupIntegerMask(mask1_) ) return Action::ERR;

  mprintf("\t");
  mask1_.BriefMaskInfo();

  if (mask1_.None()) {
    mprintf("    Error: PairDist::setup: Mask has no atoms.\n");
    return Action::ERR;
  }

  if (currentParm->SetupIntegerMask(mask2_) ) return Action::ERR;

  mask2_.BriefMaskInfo();
  mprintf("\n");

  if (mask2_.None()) {
    mprintf("    Error: PairDist::setup: Mask2 has no atoms.\n");
    return Action::ERR;
  }

  if (mask1_.MaskExpression() != mask2_.MaskExpression() &&
      mask1_.NumAtomsInCommon(mask2_) > 0) {
    mprintf("    Error: PairDist::setup: mask expressions must be either "
	    "exactly the same\n\t(equivalent to mask2 omitted) or masks must "
	    "be non-overlapping.\n");
    return Action::ERR;
  }

  if (same_mask_) {
    ub1_ = mask1_.Nselected() - 1;
    ub2_ = mask1_.Nselected();
  } else {
    ub1_ = mask1_.Nselected();
    ub2_ = mask2_.Nselected();
  }

  SetupImaging(currentParm->BoxType() );

  return Action::OK;
}


// Action_PairDist::action()
Action::RetType Action_PairDist::DoAction(int frameNum,
					  Frame* currentFrame,
					  Frame** frameAddress)
{
  unsigned long bin, j;
  double Dist = 0.0;
  Matrix_3x3 ucell, recip;
  Vec3 a1, a2;
  std::vector<double> tmp;	// per frame histogram


  tmp.resize(histogram_.size() );

  for (unsigned long i = 0; i < ub1_; i++) {
    for (same_mask_ ? j = i + 1 : j = 0; j < ub2_; j++) {
      a1 = currentFrame->XYZ(mask1_[i]);
      a2 = currentFrame->XYZ(mask2_[j]);

      switch (ImageType() ) {
      case NONORTHO:
	currentFrame->BoxCrd().ToRecip(ucell, recip);
	Dist = DIST2_ImageNonOrtho(a1, a2, ucell, recip);
	break;
      case ORTHO:
	Dist = DIST2_ImageOrtho(a1, a2, currentFrame->BoxCrd());
	break;
      case NOIMAGE:
	Dist = DIST2_NoImage(a1, a2);
	break;
      }

      bin = (unsigned long) (sqrt(Dist) / delta_);

      if (bin > maxbin_) {
	maxbin_ = bin;
	tmp.resize(maxbin_ + 1);
	histogram_.resize(maxbin_ + 1);
      }

      tmp[bin]++;
    }
  }

  // FIXME: may be inefficient to call accumulate() on every data point
  // -> pass "array" to accumulate()?
  for (unsigned long i = 0; i < tmp.size(); i++) {
    histogram_[i].accumulate(tmp[i]);
  }

  return Action::OK;
}


// Action_PairDist::print()
void Action_PairDist::Print()
{
  double dist, Pr, sd;


  output_.Printf("# pair-distance distribution P(r)\n"
		 "#distance P(r) stddev\n");

  for (unsigned long i = 0; i < histogram_.size(); i++) {
    Pr = histogram_[i].mean() / delta_;

    if (Pr > 0.0) {
      dist = ((double) i + 0.5) * delta_;
      sd = sqrt(histogram_[i].variance() );
      
      output_.Printf("%10.4f %16.2f %10.2f\n", dist, Pr, sd);
    }
  }
}
