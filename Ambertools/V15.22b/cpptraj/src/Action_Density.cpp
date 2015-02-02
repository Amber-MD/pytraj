#include <cmath>

#include "Action_Density.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h"
#include "DistRoutines.h"


/** Calculate density along a coordinate.
  * \author Hannes H. Loeffler.
  */

const std::string Action_Density::emptystring = "";

// CONSTRUCTOR
Action_Density::Action_Density() :
  axis_(DZ),
  //area_coord_{DZ, DY},		// this is C++11!
  property_(NUMBER),
  delta_(0.0)
{
  area_coord_[0] = DX; area_coord_[1] = DY;
}

void Action_Density::Help()
{
  mprintf("\tout <filename> [delta <resolution>] [x|y|z]\n"
	  "\t[number|mass|charge|electron] <mask1> ... <maskN>\n"
          "  Calculate density along a coordinate.\n");
}

// Action_Density::init()
Action::RetType Action_Density::Init(ArgList& actionArgs,
				     TopologyList* PFL, FrameList* FL,
				     DataSetList* DSL, DataFileList* DFL,
				     int debugIn)
{
  InitImaging(true);

  std::string outfileName = actionArgs.GetStringKey("out");

  if (outfileName.empty()) {
    outfileName = "density.dat";
  }

  if (output_.OpenEnsembleWrite(outfileName, DSL->EnsembleNum()) ) {
    mprinterr("Error: Density: Could not open output file %s\n",
	      outfileName.c_str());
    return Action::ERR;
  }

  if (actionArgs.hasKey("x") ) {
    axis_ = DX;
    area_coord_[0] = DY;
    area_coord_[1] = DZ;
  } else if (actionArgs.hasKey("y") ) {
    axis_ = DY;
    area_coord_[0] = DX;
    area_coord_[1] = DZ;
  } else if (actionArgs.hasKey("z") ) {
    axis_ = DZ;
    area_coord_[0] = DX;
    area_coord_[1] = DY;
  }

  property_ = NUMBER;
  if (actionArgs.hasKey("number") )   property_ = NUMBER;
  if (actionArgs.hasKey("mass") )     property_ = MASS;
  if (actionArgs.hasKey("charge") )   property_ = CHARGE;
  if (actionArgs.hasKey("electron") ) property_ = ELECTRON;

  delta_ = actionArgs.getKeyDouble("delta", 0.01);

  // for compatibility with ptraj, ignored because we rely on the atom code to
  // do the right thing, see Atom.{h,cpp}
 actionArgs.GetStringKey("efile");

  // read the rest of the command line as a series of masks
  std::string maskstr;

  while ( (maskstr = actionArgs.GetMaskNext() ) != emptystring) {
    masks_.push_back( AtomMask(maskstr) );
  }

  minus_histograms_.resize(masks_.size() );
  plus_histograms_.resize(masks_.size() );

  return Action::OK;
}


// Action_Density::Setup()
Action::RetType Action_Density::Setup(Topology* currentParm,
				      Topology** parmAddress)
{
  properties_.resize(0);

  for (std::vector<AtomMask>::iterator mask = masks_.begin();
       mask != masks_.end();
       mask++) {
    if (currentParm->SetupIntegerMask(*mask) ) return Action::ERR;

    std::vector<double> property;

    for (AtomMask::const_iterator idx = mask->begin();
	 idx != mask->end(); idx++) {
      const Atom& atom = (*currentParm)[*idx];

      switch (property_) {
      case NUMBER:
	property.push_back(1.0);
	break;

      case MASS:
	property.push_back(atom.Mass() );
	break;

      case CHARGE:
	property.push_back(atom.Charge() );
	break;

      case ELECTRON:
	property.push_back(atom.AtomicNumber() - atom.Charge() );
	break;
      }
    }

    properties_.push_back(property);

    mprintf("\t");
    mask->BriefMaskInfo();
    mprintf("\n");
  }

  SetupImaging(currentParm->BoxType() );

  return Action::OK;  
}


// Action_Density::action()
Action::RetType Action_Density::DoAction(int frameNum,
					 Frame* currentFrame,
					 Frame** frameAddress)
{
  long slice;
  unsigned long i, j;
  Vec3 coord;
  Box box;


  i = 0;

  for (std::vector<AtomMask>::const_iterator mask = masks_.begin();
       mask != masks_.end();
       mask++) {

    j = 0;

    std::map<long,double> minus_histo, plus_histo;

    for (AtomMask::const_iterator idx = mask->begin();
	 idx != mask->end();
	 idx++) {
      coord = currentFrame->XYZ(*idx);
      slice = (unsigned long) (coord[axis_] / delta_);

      if (coord[axis_] < 0) {
	minus_histo[slice] += properties_[i][j];
      } else {
	plus_histo[slice] += properties_[i][j];
      }

      j++;
    }

    if (minus_histo.size() > 0)
      minus_histograms_[i].accumulate(minus_histo);

    if (plus_histo.size() > 0)
      plus_histograms_[i].accumulate(plus_histo);

    i++;
  }

  box = currentFrame->BoxCrd();
  area_.accumulate(box[area_coord_[0]] * box[area_coord_[1]]);

  return Action::OK;
}


// Action_Density::print()
void Action_Density::Print()
{
  const unsigned int SMALL = 1.0;

  bool first_round, scale_area;
  long minus_minidx = 0, minus_maxidx = 0, plus_minidx = 0, plus_maxidx = 0;
  double density, sd, area;

  std::map<long,double>::iterator first_idx, last_idx;
  statmap curr;



  area = area_.mean();
  sd = sqrt(area_.variance());
  scale_area = (property_ == ELECTRON && area > SMALL);

  mprintf("The average box area in %c/%c is %.2f Angstrom (sd = %.2f).\n",
	  area_coord_[0] + 88, area_coord_[1] + 88, area, sd);

  if (scale_area)
    mprintf("The electron density will be scaled by this area.\n");

  // the search for minimum and maximum indices relies on ordered map
  for (unsigned long i = 0; i < minus_histograms_.size(); i++) {
    first_idx = minus_histograms_[i].mean_begin(); 
    last_idx = minus_histograms_[i].mean_end();

    if (first_idx->first < minus_minidx)
      minus_minidx = first_idx->first;

    if (last_idx != first_idx) {
      last_idx--;
      if (last_idx->first > minus_maxidx)
        minus_maxidx = last_idx->first;
    }
  }

  for (unsigned long i = 0; i < plus_histograms_.size(); i++) {
    first_idx = plus_histograms_[i].mean_begin(); 
    last_idx = plus_histograms_[i].mean_end();

    if (first_idx->first < plus_minidx)
      plus_minidx = first_idx->first;

    if (last_idx != first_idx) {
      last_idx--;
      if (last_idx->first > plus_maxidx)
        plus_maxidx = last_idx->first;
    }
  }

  output_.Printf("# routine version: %s\n#density", ROUTINE_VERSION_STRING);

  for (std::vector<AtomMask>::const_iterator mask = masks_.begin();
       mask != masks_.end();
       mask++) {
    output_.Printf(" %s sd(%s)", mask->MaskString(), mask->MaskString() );
  }

  output_.Printf("\n");

  // make sure we have zero values at beginning and end as this
  // "correctly" integrates the histogram
  minus_minidx--;
  plus_maxidx++;

  for (long i = minus_minidx; i <= minus_maxidx; i++) {
    first_round = true;

    for (unsigned long j = 0; j < minus_histograms_.size(); j++) {
      curr = minus_histograms_[j];

      if (first_round) {
        output_.Printf("%10.4f", -delta_ + ((double) i + 0.5) * delta_);
	first_round = false;
      }
      
      density = curr.mean(i) / delta_;
      sd = sqrt(curr.variance(i) );

      if (scale_area) {
	density /= area;
	sd /= area;
      }

      output_.Printf(" %10.3f %10.5f", density, sd);
    }

    output_.Printf("\n");
  }

  for (long i = plus_minidx; i <= plus_maxidx; i++) {
    first_round = true;

    for (unsigned long j = 0; j < plus_histograms_.size(); j++) {
      curr = plus_histograms_[j];

      if (first_round) {
        output_.Printf("%10.4f", ((double) i + 0.5) * delta_);
	first_round = false;
      }

      density = curr.mean(i) / delta_;
      sd = sqrt(curr.variance(i) );

      if (scale_area) {
	density /= area;
	sd /= area;
      }

      output_.Printf(" %10.3f %10.5f", density, sd);
    }

    output_.Printf("\n");
  }
}

