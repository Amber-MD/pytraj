#include <cmath>

#include "Action_OrderParameter.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h"
#include "DistRoutines.h"



/** Calculate pair distribution function P(r) between two masks.
  * \author Hannes H. Loeffler.
  */

const std::string Action_OrderParameter::emptystring = "";
const double Action_OrderParameter::MAXBOND  = 3.0;
const double Action_OrderParameter::MAXBOND2 = 1.5;

// CONSTRUCTOR
Action_OrderParameter::Action_OrderParameter() :
  axis_(DZ),
  delta_(0.1),
  scd_(false),
  maxbin_(0)
{}

void Action_OrderParameter::Help()
{
  mprintf("\tout <filename> [unsat <mask>] [scd] [unsat]\n"
	  "\t[taildist <filename> delta <resolution>\n"
	  "\ttailstart <mask> tailend <mask>]\n"
	  "\t<mask0> ... <maskN>\n");
}

// Action_OrderParameter::init()
Action::RetType Action_OrderParameter::Init(ArgList& actionArgs, TopologyList* PFL, DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  size_t nGroups;
  std::string mask;


  InitImaging(true);

  std::string outfileName = actionArgs.GetStringKey("out");

  if (outfileName.empty()) {
    outfileName = "orderparam.dat";
  }

  if (output_.OpenEnsembleWrite(outfileName, DSL->EnsembleNum()) ) {
    mprinterr("Error: OrderParameter: Could not open output file %s\n",
	      outfileName.c_str());
    return Action::ERR;
  }

  if (actionArgs.hasKey("x") ) {
    axis_ = DX;
  } else if (actionArgs.hasKey("y") ) {
    axis_ = DY;
  } else if (actionArgs.hasKey("z") ) {
    axis_ = DZ;
  }

  std::string taildistName = actionArgs.GetStringKey("taildist");

  if (!taildistName.empty()) {
    if (taildist_.OpenEnsembleWrite(taildistName, DSL->EnsembleNum()) ) {
      mprinterr("Error: OrderParameter: Could not open output file %s\n",
		outfileName.c_str());
      return Action::ERR;
    }

    delta_ = actionArgs.getKeyDouble("delta", 0.01);

    mask = actionArgs.GetStringKey("tailstart");

    if (mask.empty() ) {
      mprinterr("Error: OrderParameter: No tailstart mask specified.\n");
      return Action::ERR;
    }

    tailstart_mask_.SetMaskString(mask);

    mask = actionArgs.GetStringKey("tailend");

    if (mask.empty() ) {
      mprinterr("Error: OrderParameter: No tailend mask specified.\n");
      return Action::ERR;
    }

    tailend_mask_.SetMaskString(mask);
  }

  scd_ = actionArgs.hasKey("scd");

  mask = actionArgs.GetStringKey("unsat");

  if (!mask.empty() ) {
    unsat_mask_.SetMaskString(mask);
  }

  // rest of the command line is the masks for each atom
  // each mask stands for a single atom in a lipid
  while ( (mask = actionArgs.GetMaskNext() ) != emptystring ) {
    masks_.push_back(AtomMask(mask) );
  }

  nGroups = masks_.size();

  if (nGroups < 3) {
    mprinterr("Error: OrderParameter: number of atoms must be at least 3 "
	      "(not %i)\n", nGroups);
    return Action::ERR;
  }

  if (scd_ && (nGroups % 3) ) {
    mprinterr("Error: OrderParameter: scd set but number of masks (%i) "
	      "not a multiple of 3\n", nGroups);
    return Action::ERR;
  }

  orderParameter_.resize(nGroups);

  if (!scd_)
    dbonds_.resize(nGroups);


  return Action::OK;
}


// Action_OrderParameter::Setup()
Action::RetType Action_OrderParameter::Setup(Topology* currentParm,
					     Topology** parmAddress)
{
  int i, nlen1;
  int nlen2 = 0;
  std::vector<AtomMask>::iterator mask;


  SetupImaging(currentParm->BoxType() );


  if (!scd_) {
    if (unsat_mask_.MaskStringSet() && 
	currentParm->SetupIntegerMask(unsat_mask_) )
      return Action::ERR;
  }

  if (taildist_.IsOpen() ) {
    if (tailstart_mask_.MaskStringSet() &&
	currentParm->SetupIntegerMask(tailstart_mask_) )
      return Action::ERR;

    if (tailend_mask_.MaskStringSet() &&
	currentParm->SetupIntegerMask(tailend_mask_) )
      return Action::ERR;

    if (tailstart_mask_.Nselected() != tailend_mask_.Nselected() ) {
      mprinterr("Error: OrderParameter: different number of atoms in\n"
		"tailstart and tailend masks.");
      return Action::ERR;
    }
  }


  for (mask = masks_.begin(), i = 0;
       mask != masks_.end();
       mask++, i++) {

    if (currentParm->SetupIntegerMask(*mask) ) return Action::ERR;

    nlen1 = mask->Nselected();

    if (mask != masks_.begin() && (nlen1 != nlen2) ) {
      mprinterr("Error: OrderParameter: mask groups are not the same length "
		"(%i %i)\n", nlen1, nlen2);
      return Action::ERR;
    }

    nlen2 = nlen1;

    orderParameter_[i].resize(3 * masks_.size() );

    if (!scd_ && !unsat_mask_.None() ) {
      dbonds_[i].resize(nlen1);

      // FIXME: tag masks that are part of a double bond (is this info in
      // topology?)
      for (int j = 0; j < nlen1; j++) {
	for (int k = 0; k < unsat_mask_.Nselected(); k++) {
	  if (masks_[i][j] == unsat_mask_[k]) {
	    dbonds_[i][j] = 1;
	  }
	}
      }
    }
  }

  return Action::OK;
}


// Action_OrderParameter::action()
Action::RetType Action_OrderParameter::DoAction(int frameNum,
						Frame* currentFrame,
						Frame** frameAddress)
{
  int i, j, curr_atom, prev_atom, next_atom;

  unsigned long bin;
  std::vector<double> tmp; 

  double Sx, Sy, Sz, len;
  Vec3 sx, sy, sz, cc1, cc2, c, h1, h2, ca, cb;

  AtomMask::const_iterator C_atom, H1_atom, H2_atom;
  AtomMask::const_iterator start_atom, end_atom;
  std::vector<AtomMask>::const_iterator C_mask, H1_mask, H2_mask;



  if (scd_) {			// C coordinates plus two H coordinates
    // iterate over the all C masks, pointer arithmetics is possible because a
    // vector iterator is a random access iterator
    for (C_mask = masks_.begin(), i = 0;
	 C_mask != masks_.end();
	 C_mask += 3, i += 3) {
      
      H1_mask = C_mask + 1;
      H2_mask = C_mask + 2;

      Sx = Sy = 0.0;

      for (C_atom = C_mask->begin(),
	     H1_atom = H1_mask->begin(),
	     H2_atom = H2_mask->begin();
	     C_atom != C_mask->end() &&
	     H1_atom != H1_mask->end() &&
	     H2_atom != H2_mask->end();
	     C_atom++, H1_atom++, H2_atom++) {

	// C-H1 unit vector
	c = currentFrame->XYZ(*C_atom);
	h1 = currentFrame->XYZ(*H1_atom);

	sx = c - h1;
	len = sqrt(sx.Magnitude2() );

	if (len > MAXBOND2) {
	  mprintf("Warning: OrderParameter: unusual length (> %3.1f A) of "
		  "%f A. Check mask %s %s.\n", MAXBOND2, len,
		  C_mask->MaskString(), H1_mask->MaskString() );
	}

	sx /= len;
	Sx += 0.5 * (3.0 * sx[axis_] * sx[axis_] - 1.0);

	// C-H2 unit vector
	h2 = currentFrame->XYZ(*H2_atom);

	sy = c - h2;
	len = sqrt(sy.Magnitude2() );

	if (len > MAXBOND2) {
	  mprintf("Warning: OrderParameter: unusual length (> %3.1f A) of "
		  "%f A. Check mask %s %s.\n", MAXBOND2, len,
		  C_mask->MaskString(), H2_mask->MaskString() );
	}

	sy /= len;
	Sy += 0.5 * (3.0 * sy[axis_] * sy[axis_] - 1.0);
      }

      orderParameter_[i/3][0].accumulate(Sx / C_mask->Nselected() );
      orderParameter_[i/3][1].accumulate(Sy / C_mask->Nselected() );
    }
  } else {			// C coordinates only
    for (C_mask = masks_.begin() + 1, i = 1;
	 C_mask != masks_.end() - 1;
	 C_mask++, i++) {
      
      Sx = Sy = Sz = 0.0;

      for (j = 0; j < C_mask->Nselected(); j++) {
	prev_atom = (*(C_mask-1))[j];
	curr_atom = (*C_mask)[j];
	next_atom = (*(C_mask+1))[j];

	if (!unsat_mask_.None() && dbonds_[i][j]) {
	  ca = currentFrame->XYZ(curr_atom);
	  cb = currentFrame->XYZ(next_atom);
	} else {
	  ca = currentFrame->XYZ(prev_atom);
	  cb = currentFrame->XYZ(next_atom);
	}

	sz = cb - ca;
	len = sqrt(sz.Magnitude2() );

	if (len > MAXBOND) {
	  mprintf("Warning: OrderParameter: unusual length (> %3.1f A) of "
		  "%f A. Check mask %s.\n", MAXBOND, len, C_mask->MaskString() );
	}

	sz /= len;

	ca = currentFrame->XYZ(curr_atom);
	cb = currentFrame->XYZ(prev_atom);

	cc1 = ca - cb;

	cb = currentFrame->XYZ(next_atom);
	cc2 = ca - cb;

	sx = cc1.Cross(cc2);
	sx.Normalize();

	sy = sz.Cross(sx);
	sy.Normalize();

	// using dotproduct of Sx, Sy and Sz with axis (two components are zero!)
	Sx += 0.5 * (3.0 * sx[axis_] * sx[axis_] - 1.0);
	Sy += 0.5 * (3.0 * sy[axis_] * sy[axis_] - 1.0);
	Sz += 0.5 * (3.0 * sz[axis_] * sz[axis_] - 1.0);
      }

      orderParameter_[i][0].accumulate(Sx / C_mask->Nselected() );
      orderParameter_[i][1].accumulate(Sy / C_mask->Nselected() );
      orderParameter_[i][2].accumulate(Sz / C_mask->Nselected() );
    }
  }

  tmp.resize(tailhist_.size() );

  if (taildist_.IsOpen() ) {
    for (start_atom = tailstart_mask_.begin(), end_atom = tailend_mask_.begin();
	 start_atom != tailstart_mask_.end() && end_atom != tailend_mask_.end();
	 start_atom++, end_atom++) {

      ca = currentFrame->XYZ(*start_atom);
      cb = currentFrame->XYZ(*end_atom);

      c = ca - cb;
      bin = (unsigned long) (sqrt(c.Magnitude2() ) / delta_);
      
      if (bin > maxbin_) {
        maxbin_ = bin;
        tmp.resize(maxbin_ + 1);
        tailhist_.resize(maxbin_ + 1);
      }

      tmp[bin]++;
    }

    for (unsigned long i = 0; i < tmp.size(); i++) {
      tailhist_[i].accumulate(tmp[i]);
    }
  }


  return Action::OK;
}


// Action_OrderParameter::print()
void Action_OrderParameter::Print()
{
  double Sx, Sy, Sz, SCD_1, SCD_2, prob, sd;


  output_.Printf("# order parameters for masks");


  for (std::vector<AtomMask>::iterator mask = masks_.begin();
       mask != masks_.end();
       mask++) {
    output_.Printf(" %s", mask->MaskString() );
  }

  output_.Printf("\n");

  if (scd_) {
    output_.Printf("#Cn %10s %10s %10s %10s\n",
		   "SCD_H1", "sd(SCD_H1)", "SCD_H2", "sd(SCD_H2)");

    for (unsigned int i = 0; i < orderParameter_.size() / 3; i++) {
      Sx = -orderParameter_[i][0].mean();
      Sy = -orderParameter_[i][1].mean();

      output_.Printf("%3u %10.7f %10.7f %10.7f %10.7f\n",
		     i + 1,
		     Sx,
                     sqrt(orderParameter_[i][0].variance()),
		     Sy,
                     sqrt(orderParameter_[i][1].variance()) );
    }
  } else {
    output_.Printf("#Cn %10s %10s %10s %10s %10s %10s %10s %10s\n",
		   "Sx", "sd(Sx)", "Sy", "sd(Sy)", "Sz",  "sd(Sz)", "SCD_z",
		   "SCD_xy");

    for (unsigned int i = 1; i < orderParameter_.size() - 1; i++) {
      Sx = orderParameter_[i][0].mean();
      Sy = orderParameter_[i][1].mean();
      Sz = orderParameter_[i][2].mean();

      SCD_1 = 0.5 * Sz;
      SCD_2 = -(2.0 * Sx + Sy) / 3.0;

      output_.Printf("%3u %10.7f %10.7f %10.7f %10.7f %10.7f %10.7f "
		     "%10.7f %10.7f\n",
		     i + 1,
		     Sx,
                     sqrt(orderParameter_[i][0].variance()),
		     Sy,
                     sqrt(orderParameter_[i][1].variance()),
		     Sz,
                     sqrt(orderParameter_[i][2].variance()),
		     SCD_1,
                     SCD_2);
    }
  }

  if (taildist_.IsOpen() ) {
    taildist_.Printf("# end-to-end distance\n");

    for (unsigned long i = 0; i < tailhist_.size(); i++) {
      prob = tailhist_[i].mean() / delta_;
      sd = sqrt(tailhist_[i].variance() );

      if (prob > 0.0) {
	taildist_.Printf("%10.4f %10.4f %10.7f\n", ((double) i + 0.5) * delta_,
			 prob, sd);
      }
    }
  }
}
