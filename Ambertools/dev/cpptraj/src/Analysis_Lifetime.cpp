#include "Analysis_Lifetime.h"
#include "CpptrajStdio.h"
#include "ProgressBar.h"
#include "StringRoutines.h" // integerToString

// CONSTRUCTOR
Analysis_Lifetime::Analysis_Lifetime() :
  windowSize_(0),
  fuzzCut_(-1),
  cut_(0.5),
  averageonly_(false),
  cumulative_(false),
  deltaAvg_(false),
  normalizeCurves_(true),
  Compare_(Compare_GreaterThan)
{}

void Analysis_Lifetime::Help() {
  mprintf("\t[out <filename>] <dsetarg0> [ <dsetarg1> ... ]\n"
          "\t[window <windowsize> [name <setname>]] [averageonly]\n"
          "\t[cumulative] [cut <cutoff>] [greater | less] [rawcurve]\n"
          "\t[fuzz <fuzzcut>] [nosort]\n"
          "  Calculate lifetimes for specified data set(s), i.e. time that data is\n"
          "  either greater than or less than <cutoff> (default: > 0.5). If <windowsize>\n"
          "  is given calculate lifetimes over windows of given size.\n");
}

// Analysis_Lifetime::Setup()
Analysis::RetType Analysis_Lifetime::Setup(Array1D const& dsArray, std::string const& fname) {
  if (dsArray.empty()) return Analysis::ERR;
  outfileName_ = fname;
  inputDsets_ = dsArray;
  windowSize_ = -1;
  averageonly_ = false;
  cumulative_ = false;
  deltaAvg_ = false;
  normalizeCurves_ = true;
  cut_ = 0.5;
  fuzzCut_ = -1;
  Compare_ = Compare_GreaterThan;
  return Analysis::OK;
}

inline static int CheckDsetError(DataSet_1D* ds, const char* msg, const char* legend) {
  if (ds == 0) {
    mprinterr("Error: lifetime: Could not allocate %s set for %s\n", msg, legend);
    return 1;
  }
  return 0;
}

// Analysis_Lifetime::Setup()
Analysis::RetType Analysis_Lifetime::Setup(ArgList& analyzeArgs, DataSetList* datasetlist,
                            TopologyList* PFLin, DataFileList* DFLin, int debugIn)
{
  // Get Keywords
  outfileName_ = analyzeArgs.GetStringKey("out");
  std::string setname = analyzeArgs.GetStringKey("name");
  bool sortSets = (!analyzeArgs.hasKey("nosort"));
  windowSize_ = analyzeArgs.getKeyInt("window", -1);
  averageonly_ = analyzeArgs.hasKey("averageonly");
  cumulative_ = analyzeArgs.hasKey("cumulative");
  deltaAvg_ = analyzeArgs.hasKey("delta");
  cut_ = analyzeArgs.getKeyDouble("cut", 0.5);
  fuzzCut_ = analyzeArgs.getKeyInt("fuzz", -1);
  if (fuzzCut_ < 1) fuzzCut_ = -1;
  normalizeCurves_ = !analyzeArgs.hasKey("rawcurve");
  if (analyzeArgs.hasKey("greater"))
    Compare_ = Compare_GreaterThan;
  else if (analyzeArgs.hasKey("less"))
    Compare_ = Compare_LessThan;
  else
    Compare_ = Compare_GreaterThan;
  // Select datasets from remaining args
  if (inputDsets_.AddSetsFromArgs( analyzeArgs.RemainingArgs(), *datasetlist )) {
    mprinterr("Error: lifetime: Could not add data sets.\n");
    return Analysis::ERR;
  }
  // Sort data sets
  if (sortSets) inputDsets_.SortArray1D(); 

  // Create output datasets
  DataFile* outfile = 0;
  DataFile* maxfile = 0;
  DataFile* avgfile = 0;
    if (setname.empty()) 
      setname = datasetlist->GenerateDefaultName( "lifetime" );
  if ( windowSize_ != -1) {
    outfile = DFLin->AddDataFile(outfileName_, analyzeArgs);
    if (!averageonly_ && outfile != 0) {
      maxfile = DFLin->AddDataFile(outfile->DataFilename().DirPrefix() + "max." + 
                                   outfile->DataFilename().Base(), analyzeArgs);
      avgfile = DFLin->AddDataFile(outfile->DataFilename().DirPrefix() + "avg." + 
                                   outfile->DataFilename().Base(), analyzeArgs);
    }
    int didx = 0;
    for (Array1D::const_iterator set = inputDsets_.begin(); set != inputDsets_.end(); ++set)
    {
      DataSet_1D* outSet = (DataSet_1D*)datasetlist->AddSetIdx( DataSet::FLOAT, setname, didx );
      if (CheckDsetError(outSet, "output", (*set)->Legend().c_str())) 
        return Analysis::ERR;
      outSet->SetLegend( (*set)->Legend() );
      outputDsets_.push_back( outSet );
      if (outfile != 0) outfile->AddSet( outSet );
      if (!averageonly_) {
        // MAX
        outSet = (DataSet_1D*)datasetlist->AddSetIdxAspect(DataSet::INTEGER, setname, didx, "max");
        if (CheckDsetError(outSet, "lifetime max", (*set)->Legend().c_str()))
          return Analysis::ERR;
        outSet->SetLegend( (*set)->Legend() );
        maxDsets_.push_back( outSet );
        if (maxfile != 0) maxfile->AddSet( outSet );
        // AVG
        outSet = (DataSet_1D*)datasetlist->AddSetIdxAspect(DataSet::FLOAT, setname, didx, "avg");
        if (CheckDsetError(outSet, "lifetime avg", (*set)->Legend().c_str()))
          return Analysis::ERR;
        outSet->SetLegend( (*set)->Legend() );
        avgDsets_.push_back( outSet );
        if (avgfile != 0) avgfile->AddSet( outSet );
      }
      ++didx;
    }
    // Set step to window size.
    std::string fileArgs = "xstep " + integerToString( windowSize_ ); 
    if (outfile != 0) outfile->ProcessArgs( fileArgs );
    if (maxfile != 0) maxfile->ProcessArgs( fileArgs );
    if (avgfile != 0) avgfile->ProcessArgs( fileArgs );
  }
  // Lifetime curves
  DataFile* crvfile = 0;
  if (!averageonly_) {
    if (!outfileName_.empty()) {
      FileName crvFilename;
      crvFilename.SetFileName( outfileName_ );
      crvfile = DFLin->AddDataFile(crvFilename.DirPrefix() + "crv." + 
                                   crvFilename.Base(), analyzeArgs);
    }
    for (int didx = 0; didx != (int)inputDsets_.size(); didx++)
    {
      DataSet_1D* outSet = (DataSet_1D*)
        datasetlist->AddSetIdxAspect(DataSet::DOUBLE, setname, didx, "curve");
      if (CheckDsetError(outSet, "lifetime curve", inputDsets_[didx]->Legend().c_str()))
        return Analysis::ERR;
      curveSets_.push_back( outSet );
      if (crvfile != 0) crvfile->AddSet( outSet );
    }
  }

  if (!averageonly_)
    mprintf("    LIFETIME: Calculating average lifetime using a cutoff of %f", cut_);
  else
    mprintf("    LIFETIME: Calculating only averages");
  mprintf(" of data in %i sets\n", inputDsets_.size());
  if (!sortSets) mprintf("\tInput data sets will not be sorted.\n");
  if (debugIn > 0)
    for (Array1D::const_iterator set = inputDsets_.begin(); set != inputDsets_.end(); ++set)
      mprintf("\t%s\n", (*set)->Legend().c_str());
  if (Compare_ == Compare_GreaterThan) 
    mprintf("\tValues greater than %f are considered present.\n", cut_);
  else
    mprintf("\tValues less than %f are considered present.\n", cut_);
  if (windowSize_ != -1) {
    mprintf("\tAverage of data over windows will be saved to sets named %s\n",
            setname.c_str());
    mprintf("\tWindow size for averaging: %i\n", windowSize_);
    if (cumulative_)
      mprintf("\tCumulative averages will be saved.\n");
    if (deltaAvg_)
      mprintf("\tChange of average from previous average will be saved.\n");
  }
  if (!outfileName_.empty()) {
    mprintf("\tOutfile: %s", outfileName_.c_str());
    if (!averageonly_ && outfile != 0)
      mprintf(", %s, %s", maxfile->DataFilename().base(), avgfile->DataFilename().base());
    mprintf("\n");
  }
  if (!averageonly_) {
    if (crvfile != 0)
      mprintf("\tLifetime curves output: %s\n", crvfile->DataFilename().base());
    if (normalizeCurves_)
      mprintf("\tLifetime curves will be normalized.\n");
    else
      mprintf("\tLifetime curves will not be normalized.\n");
  }
  if (fuzzCut_ != -1)
    mprintf("\tFuzz value of %i frames will be used.\n", fuzzCut_);
  return Analysis::OK;
}

// TODO: Make most of these class variables. 
inline static void RecordCurrentLifetime(int start, int stop, int length,
                                         int& maximumLifetimeCount,
                                         int& sumLifetimes, int& Nlifetimes,
                                         std::vector<int>& lifetimeCurve)
{
  // Record current lifetime information. 
  if (length > maximumLifetimeCount)
    maximumLifetimeCount = length;
  sumLifetimes += length;
  ++Nlifetimes;
  //mprintf("%i-%i; LC=%i maxLC=%i NL=%i\n", start+1, stop, length, // DEBUG
  //        maximumLifetimeCount, Nlifetimes); // DEBUG
  if (length > (int)lifetimeCurve.size())
    lifetimeCurve.resize( length, 0 );
  for (int lc = 0; lc != length; lc++)
    lifetimeCurve[lc]++;
  // DEBUG: Current lifetime curve
//  mprintf("CRV:");
//  for (std::vector<int>::const_iterator it = lifetimeCurve.begin(); 
//                                        it != lifetimeCurve.end(); ++it)
//    mprintf(" %i", *it);
//  mprintf("\n");
}

// Analysis_Lifetime::Analyze()
Analysis::RetType Analysis_Lifetime::Analyze() {
  float favg;
  int current = 0;
  CpptrajFile standalone_out;
  bool standalone = (!averageonly_ && windowSize_ == -1);
  if (standalone) {
    if (standalone_out.OpenWrite( outfileName_ ))
      return Analysis::ERR;
    standalone_out.Printf("%-10s %10s %10s %10s %10s %s\n","#Set","Nlifetimes",
                          "MaxLT","AvgLT","TotFrames","SetName");
  }
  ProgressBar progress( inputDsets_.size() );
  std::vector<int> lifetimeCurve;
  std::vector<int> localCurve;
  for (unsigned int setIdx = 0; setIdx < inputDsets_.size(); setIdx++) {
    lifetimeCurve.clear();
    localCurve.clear();
    DataSet_1D const& DS = static_cast<DataSet_1D const&>( *inputDsets_[setIdx] );
    if (standalone)
      mprintf("\t\tCalculating lifetimes for set %s\n", DS.Legend().c_str());
    else
      progress.Update( current++ );
    if (DS.Size() < 1) {
      mprintf("Warning: Set %s is empty, skipping.\n", DS.Legend().c_str());
      continue;
    }
    // Loop over all values in set.
    int setSize = (int)DS.Size();
    double sum = 0.0;
    double previous_windowavg = 0.0;
    int windowcount = 0; // Used to trigger averaging
    int Ncount = 0;      // Used in averaging; if !cumulative, == windowcount
    int frame = 0;       // Frame to add data at.
    int maximumLifetimeCount = 0; // Max observed lifetime
    int Nlifetimes = 0;           // # of separate lifetimes observed
    int sumLifetimes = 0;         // sum of lifetimeCount for each lifetime observed
    int potentialLifetimeStart=0;
    int potentialLifetimeStop=0;
    // Are we using fuzz values?
    int startingFuzzCount;
    if (fuzzCut_ < 1)
      startingFuzzCount = -1;
    else
      startingFuzzCount = 0;
    int fuzzCount = startingFuzzCount;
    // Where are we starting
    enum LocationType { OUTSIDE=0, INSIDE, OUTER_FUZZ, INNER_FUZZ };
    //static const char* Lstr[] = {"OUTSIDE", "INSIDE", "OUTER_FUZZ", "INNER_FUZZ"};
    LocationType location;
    if ( Compare_(DS.Dval(0), cut_) )
      location = INSIDE;
    else
      location = OUTSIDE;
    // Loop over all data points
    for (int i = 0; i < setSize; ++i) {
      double dval = DS.Dval(i);
      //mprintf("\t\t\tValue[%i]= %.2f", i,dval);
      if (averageonly_) 
        // Average only
        sum += dval;
      else {
        // Lifetime calculation
        bool frameIsInside = Compare_(dval, cut_);
        bool calculateLifetime = false;
//        mprintf("%8i Location=%s, frameIsInside=%i, fuzzCount=%i\n", // DEBUG
//                i+1, Lstr[location], (int)frameIsInside,  fuzzCount); // DEBUG
        // NOTE: As currently implemented, there must be fuzzCut + 1
        //       consecutive frames for a lifetime to exist.
        switch (location) {
          case OUTSIDE:
            if (frameIsInside) {
              potentialLifetimeStart = i;
              if (fuzzCount == 0) {
                // We have come in from outside but are within fuzz boundary.
                fuzzCount = 1;
                //mprintf("%i: Entered inner fuzz boundary. POTENTIAL LIFETIME START.\n", i+1);
                location = INNER_FUZZ;
              } else {
                // Not doing fuzz calc. We are now inside.
                location = INSIDE;
              }
            }
            break;
          case INSIDE:
            if (!frameIsInside) {
              potentialLifetimeStop = i;
              if (fuzzCount == 0) {
                // We have gone out from inside but are within fuzz boundary.
                fuzzCount = 1;
                //mprintf("%i: Exited to outer fuzz boundary. POTENTIAL LIFETIME STOP.\n", i+1);
                location = OUTER_FUZZ;
              } else {
                // Not doing fuzz calc. We are now outside.
                location = OUTSIDE;
                calculateLifetime = true;
              }
            }
            break;
          case OUTER_FUZZ:
            if (frameIsInside) {
              // We were in outer fuzz but have come back in.
              //mprintf("%i: Back inside from outer fuzz after %i frames. Still a lifetime.\n", i+1, fuzzCount);
              location = INSIDE;
              fuzzCount = startingFuzzCount;
            } else {
              fuzzCount++;
              if (fuzzCount > fuzzCut_) {
                // We have been in outer fuzz too long. Now outside.
                //mprintf("%i: Exited outer fuzz for outside after %i frames.\n", i+1, fuzzCount);
                location = OUTSIDE;
                calculateLifetime = true;
                fuzzCount = startingFuzzCount;
              }
            }
            break;
          case INNER_FUZZ:
            if (!frameIsInside) {
              // We were in inner fuzz but have come back out.
              if (fuzzCount == 0) {
                //mprintf("%i: Exiting inner fuzz for outside. Not a lifetime.\n", i+1);
                location = OUTSIDE;
                fuzzCount = startingFuzzCount;
              } else
                fuzzCount--;
            } else {
              fuzzCount++;
              if (fuzzCount > fuzzCut_) {
                // We have been inside inner fuzz enough.
                //mprintf("%i: Entering inside from inner fuzz after %i frames.\n", i+1, fuzzCount);
                location = INSIDE;
                fuzzCount = startingFuzzCount;
              }
            }
            break;
        } // END switch location
        //mprintf("%8i Start=%i, Stop=%i\n", i+1, potentialLifetimeStart+1, potentialLifetimeStop+1);
        if (calculateLifetime) {
          // Were there enough frames?
          //potentialLifetimeStop++; // Include last frame in lifetime.
          int lifetimeLength = potentialLifetimeStop - potentialLifetimeStart;
          if (lifetimeLength > fuzzCut_) { // TODO: Necessary? Always true if calcLifetime?
            sum += (double)lifetimeLength;
            RecordCurrentLifetime(potentialLifetimeStart, potentialLifetimeStop,
                                  lifetimeLength,
                                  maximumLifetimeCount, sumLifetimes, Nlifetimes,
                                  lifetimeCurve);
          }
        }
      } // END lifetime calc for this frame
      //sum += inputDsets_[setIdx]->Dval(i);
      ++Ncount;
      ++windowcount;
      if (windowcount == windowSize_) {
        // Treat this as the end of an independent run. If lifetime was not
        // calcd this frame, determine if it should be.
        if (!averageonly_) {
          if  (location == INSIDE || location == OUTER_FUZZ) {
            // Were there enough frames?
            potentialLifetimeStop = i + 1; // Include last frame in lifetime.
            int lifetimeLength = potentialLifetimeStop - potentialLifetimeStart;
            if (lifetimeLength > fuzzCut_) { // TODO: Necessary? Always true if calcLifetime?
              sum += (double)lifetimeLength;
              RecordCurrentLifetime(potentialLifetimeStart, potentialLifetimeStop,
                                    lifetimeLength,
                                    maximumLifetimeCount, sumLifetimes, Nlifetimes,
                                    lifetimeCurve);
            }
            // Reset location to prevent potential lifetime trigger next frame
            location = OUTSIDE;
          }
          // If Nlifetimes is 0 then value was never present. 
          if (Nlifetimes == 0)
            favg = 0.0;
          else
            favg = (float)sumLifetimes / (float)Nlifetimes;
          //mprintf("\t\t\t[%i]Max lifetime observed: %i frames\n", frame,maximumLifetimeCount);
          //mprintf("\t\t\t[%i]Avg lifetime: %f frames\n", frame, favg);
          maxDsets_[setIdx]->Add( frame, &maximumLifetimeCount );
          avgDsets_[setIdx]->Add( frame, &favg );
        }
        //mprintf("WINDOW BREAK\n");

        double windowavg = sum / (double)Ncount;
        float fval = (float)(windowavg - previous_windowavg);
        if (deltaAvg_) previous_windowavg = windowavg;
        outputDsets_[setIdx]->Add( frame, &fval );
        //frame += windowcount;
        frame++;
        // Window counter is always reset
        windowcount = 0;
        if (!cumulative_) {
          // Reset average counters
          sum = 0;
          Ncount = 0;
          // Reset lifetime counters
          maximumLifetimeCount = 0;
          Nlifetimes = 0;
          sumLifetimes = 0;
        }
      }
    } // END loop over data points.
    // Print lifetime information if no window
    if ( !averageonly_ && windowSize_ == -1 ) {
      // Update current lifetime total
      if  (location == INSIDE || location == OUTER_FUZZ) {
        // Were there enough frames?
        potentialLifetimeStop = setSize; // Include last frame in lifetime.
        int lifetimeLength = potentialLifetimeStop - potentialLifetimeStart;
        if (lifetimeLength > fuzzCut_) { // TODO: Necessary? Always true if calcLifetime?
          sum += (double)lifetimeLength;
          RecordCurrentLifetime(potentialLifetimeStart, potentialLifetimeStop,
                                lifetimeLength,
                                maximumLifetimeCount, sumLifetimes, Nlifetimes,
                                lifetimeCurve);
        }
      }

      // If Nlifetimes is 0 then value was never present. 
      if (Nlifetimes == 0) 
        favg = 0.0;
      else
        favg = (float)sumLifetimes / (float)Nlifetimes;
      standalone_out.Printf("%10u %10i %10i %10.4f %10.0f %s\n",setIdx,
                            Nlifetimes, maximumLifetimeCount, favg, sum,
                            DS.Legend().c_str());
    }
    // Calculate normalized lifetime curve
    if (!lifetimeCurve.empty() && !curveSets_.empty()) {
      curveSets_[setIdx]->Allocate1D( lifetimeCurve.size() );
      double norm;
      if (normalizeCurves_)
        norm = 1.0 / (double)lifetimeCurve.front();
      else
        norm = 1.0;
      for (unsigned int n = 0; n != lifetimeCurve.size(); n++) {
        double dval = lifetimeCurve[n] * norm;
        curveSets_[setIdx]->Add(n, &dval);
      }
    }
  }
  return Analysis::OK;
}
