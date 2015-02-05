#include <cmath>
#include "Analysis_CrdFluct.h"
#include "CpptrajStdio.h"
#include "Constants.h" // PI
#include "StringRoutines.h" // integerToString

Analysis_CrdFluct::Analysis_CrdFluct() : 
  coords_(0),
  bfactor_(false),
  windowSize_(-1)
{}

void Analysis_CrdFluct::Help() {
  mprintf("\t[crdset <crd set>] [<mask>] [out <filename>] [window <size>] [bfactor]\n"
          "  Calculate atomic positional fluctuations for atoms in <mask>\n"
          "  over windows of specified size.\n"
          "  <crd set> can be created with the 'createcrd' command.\n");
}

// Analysis_CrdFluct::Setup()
Analysis::RetType Analysis_CrdFluct::Setup(ArgList& analyzeArgs, DataSetList* datasetlist,
                            TopologyList* PFLin, DataFileList* DFLin, int debugIn)
{
  bfactor_ = analyzeArgs.hasKey("bfactor");
  // Attempt to get coords dataset from datasetlist
  std::string setname = analyzeArgs.GetStringKey("crdset");
  coords_ = (DataSet_Coords*)datasetlist->FindCoordsSet( setname );
  if (coords_ == 0) {
    mprinterr("Error: crdfluct: Could not locate COORDS set corresponding to %s\n",
              setname.c_str());
    return Analysis::ERR;
  }
  DataFile* outfile = DFLin->AddDataFile( analyzeArgs.GetStringKey("out"), analyzeArgs );
  windowSize_ = analyzeArgs.getKeyInt("window", -1);
  // Get mask
  mask_.SetMaskString( analyzeArgs.GetMaskNext() );

  mprintf("    CRDFLUCT: Atomic fluctuations will be calcd for set %s, mask [%s]\n", 
          coords_->Legend().c_str(), mask_.MaskString());
  if (windowSize_ != -1) mprintf("\tWindow size = %i\n", windowSize_);
  if (outfile != 0) mprintf("\tOutput to %s\n", outfile->DataFilename().base());

  // Set up data sets
  setname = analyzeArgs.GetStringNext();
  if (windowSize_ < 1) {
    // Only one data set for total B-factors
    DataSet* ds = datasetlist->AddSet( DataSet::DOUBLE, setname, "fluct" );
    if (ds == 0) return Analysis::ERR;
    outSets_.push_back( ds );
    if (outfile != 0) outfile->AddSet( ds );
  } else {
    if (coords_->Size() == 0) {
      mprinterr("Error: window size > 0 and COORDS data set %s is empty.\n", 
                 coords_->Legend().c_str());
      mprinterr("Error: Cannot predict how many window data sets will be needed.\n");
      return Analysis::ERR;
    }
    if (setname.empty()) setname = datasetlist->GenerateDefaultName("fluct");
    // Determine how many windows will be needed
    int nwindows = coords_->Size() / windowSize_;
    for (int win = 0; win < nwindows; ++win) {
      int frame = (win + 1) * windowSize_;
      DataSet* ds = datasetlist->AddSetIdx( DataSet::DOUBLE, setname, frame );
      if (ds == 0) return Analysis::ERR;
      ds->SetLegend( "F_" + integerToString( frame ) );
      ds->Dim(Dimension::X).SetLabel("Atom");
      outSets_.push_back( ds );
      if (outfile != 0) outfile->AddSet( ds );
    }
    if ( (coords_->Size() % windowSize_) != 0 ) {
      DataSet* ds = datasetlist->AddSetIdx( DataSet::DOUBLE, setname, coords_->Size() );
      ds->SetLegend("Final");
      outSets_.push_back( ds );
      if (outfile != 0) outfile->AddSet( ds );
    }
    for (SetList::iterator out = outSets_.begin(); out != outSets_.end(); ++out)
      mprintf("\t%s\n", (*out)->Legend().c_str());
  }
  // Setup output file
/*  if (bfactor_)
    outfile->Dim(Dimension::Y).SetLabel("B-factors");
  outfile->Dim(Dimension::X).SetLabel("Atom");*/

  return Analysis::OK;
}

// Analysis_CrdFluct::CalcBfactors()
void Analysis_CrdFluct::CalcBfactors( Frame SumCoords, Frame SumCoords2, double Nsets,
                                      DataSet& outset )
{
  SumCoords.Divide(Nsets);
  SumCoords2.Divide(Nsets);
  //SumCoords2 = SumCoords2 - (SumCoords * SumCoords);
  SumCoords *= SumCoords;
  SumCoords2 -= SumCoords;
  AtomMask::const_iterator maskat = mask_.begin();
  if (bfactor_) {
    // Set up b factor normalization
    // B-factors are (8/3)*PI*PI * <r>**2 hence we do not sqrt the fluctuations
    double bfac = (8.0/3.0)*Constants::PI*Constants::PI;
    for (int i = 0; i < SumCoords2.size(); i+=3) {
      double fluct = (SumCoords2[i] + SumCoords2[i+1] + SumCoords2[i+2]) * bfac;
      outset.Add( *maskat, &fluct );
      ++maskat; 
    }
  } else {
    // Atomic fluctuations
    for (int i = 0; i < SumCoords2.size(); i+=3) {
      double fluct = SumCoords2[i] + SumCoords2[i+1] + SumCoords2[i+2];
      if (fluct > 0)
        outset.Add( *maskat, &fluct );
      ++maskat;
    }
  }
}

// Analysis_CrdFluct::Analyze()
Analysis::RetType Analysis_CrdFluct::Analyze() {
  // Set up mask
  if ( coords_->Top().SetupIntegerMask( mask_ )) return Analysis::ERR;
  mask_.MaskInfo();
  if (mask_.None()) return Analysis::ERR;
  int end = coords_->Size();
  mprintf("\tFluctuation analysis for %i frames (%i atoms each).\n", end, 
          coords_->Top().Natom());
  Frame currentFrame( mask_.Nselected() );
  Frame SumCoords( mask_.Nselected() );
  SumCoords.ZeroCoords();
  Frame SumCoords2( mask_.Nselected() );
  SumCoords2.ZeroCoords();
  int w_count = 0;
  SetList::iterator out = outSets_.begin();
  for (int frame = 0; frame < end; frame++) {
    coords_->GetFrame(frame, currentFrame, mask_);
    SumCoords += currentFrame;
    SumCoords2 += ( currentFrame * currentFrame );
    ++w_count;
    if (w_count == windowSize_) {
      CalcBfactors( SumCoords, SumCoords2, (double)frame, *(*out) );
      ++out;
      w_count = 0;
    }
  }

  if (windowSize_ < 1 || w_count != 0) {
    // For windowSize < 1 this is the only b-factor calc
    CalcBfactors( SumCoords, SumCoords2, (double)end, *(*out) );
    if (w_count != 0) 
      mprintf("Warning: Number of frames (%i) was not evenly divisible by window size.\n",
               end);
  }

  return Analysis::OK;
}
