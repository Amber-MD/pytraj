#include "Action_MRT.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_MRT::Action_MRT() :
  time_(1),
  nStar_(0),
  lowerCutoff2_(0),
  upperCutoff2_(0),
  wSize_(0),
  nOffset_(0),
  idxMaxWin_(0)
{}

void Action_MRT::Help() {

}

// Action_MRT::init()
/** Usage: 
  * mrt out <filename> ([autocorr <filename> [tcorr <time>] [toffset <time>]])
  *    [lower <dist>] [upper <dist>] [time <t>] [tstar <t>] [noimage]
  *    ([solvent <mask> | solute <mask>])
  *    (siteatoms <mask> | onemol <mask> | <sitemask1> ... <sitemaskN>)
  */
Action::RetType Action_MRT::Init(ArgList& actionArgs, TopologyList* PFL, FrameList* FL,
                          DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  useImage_ = !(actionArgs.hasKey("noimage"));
  std::string filename = actionArgs.GetStringKey("out");
  if (filename.empty()) {
    mprinterr("Error: MRT: No output filename specified 'out <filename>'\n");
    return Action::ERR;
  }
  time_ = actionArgs.getKeyDouble("time", 1.0);

  // Time that water can be inside/outside without counting it as
  // having left/entered.
  double tstar = actionArgs.getKeyDouble("tstar", 0.0);
  tstar /= (time_ + 0.5);
  nStar_ = (int)tstar;

  // Lower and upper distance limits
  double lowerCutoff = actionArgs.getKeyDouble("lower", 0.01);
  double upperCutoff = actionArgs.getKeyDouble("upper", 3.50);
  lowerCutoff2_ = lowerCutoff * lowerCutoff;
  upperCutoff2_ = upperCutoff * upperCutoff;

  // If specified, filename for autocorrelation fn
  autoCorr_ = actionArgs.GetStringKey("autocorr");

  // Autocorrelation parameters
  // NOTE: Needed if autoCorr empty?
  double wsize = actionArgs.getKeyDouble("tcorr", 400.0);
  wsize /= (time_ + 0.5);
  wSize_ = (int)wsize;

  double noffset = actionArgs.getKeyDouble("toffset", 10.0);
  noffset /= (time_ + 0.5);
  nOffset_ = (int)noffset;

  if (nOffset_ < 1 || nOffset_ > wSize_) {
    mprinterr("Error: MRT: toffset must be in the range from %8.3f to %8.3f.\n",
               time_, (double)wSize_ * time_);
    return Action::ERR;
  }

  if ( (wSize_ % nOffset_) != 0 ) {
    mprinterr("Error: MRT: tcorr must be multiple of toffset.\n");
    return Action::ERR;
  }

  int nWin = wSize_ / nOffset_;
  idxMaxWin_ = nWin - 1;

  // Get solvent mask expression. 
  // If it is a solute,  mask will be used as a single site.
  solventmask_.SetMaskString( actionArgs.GetStringKey("solvent") );
    
  // Check if MRT to another solute is requested
  solutemask_.SetMaskString( actionArgs.GetStringKey("solute") );

  // Process reference sites. Currently:
  //   1) siteatoms <mask> which expands to mutliple sites per atom in mask
  //   2) onemol <mask> means a single
  //   3) <mask1>...<maskN> multiple sites, use center-of-mass 

  return Action::OK;
}
