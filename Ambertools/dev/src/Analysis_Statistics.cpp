#include <cmath> // sqrt
#include "Analysis_Statistics.h"
#include "DataSet_double.h" // for DISTANCE NOE
#include "CpptrajStdio.h"

// CONSTRUCTOR
Analysis_Statistics::Analysis_Statistics() :
  shift_(0),
  debug_(0),
  NOE_r6_(0),
  NOE_violations_(0),
  NOE_avgViolations_(0),
  NOE_names_(0),
  ignore_negative_violations_(true)
{}

void Analysis_Statistics::Help() {
  mprintf("\t{<name> | all} [shift <value>] [out <filename>] [noeout <filename>]\n"
          "\t [ignorenv]\n"
          "  Calculate various statistical quantities for data in specified data set(s)\n"
          "  based on the data set type (e.g. distance noe, dihedral alpha, etc)\n");
}

// Analysis_Statistics::Setup()
Analysis::RetType Analysis_Statistics::Setup(ArgList& analyzeArgs, DataSetList* DSLin,
                            TopologyList* PFLin, DataFileList* DFLin, int debugIn)
{
  debug_ = debugIn;
  // Get keywords.
  shift_ = analyzeArgs.getKeyDouble("shift", 0);
  filename_ = analyzeArgs.GetStringKey("out");
  ignore_negative_violations_ = analyzeArgs.hasKey("ignorenv");
  DataFile* NOE_out = DFLin->AddDataFile(analyzeArgs.GetStringKey("noeout"), analyzeArgs);
  // Get dataset or all datasets
  bool useAllSets = false;
  if (analyzeArgs.hasKey("all")) {
    useAllSets = true;
    for (DataSetList::const_iterator ds = DSLin->begin(); ds != DSLin->end(); ++ds)
      if ( (*ds)->Ndim() == 1)
        datasets_.push_back( ((DataSet_1D*)*ds) );
  } else {
    // Select datasets from remaining args
    if (datasets_.AddSetsFromArgs( analyzeArgs.RemainingArgs(), *DSLin )) {
      mprinterr("Error: statistics: Could not add data sets\n");
      return Analysis::ERR;
    }
  }
  if (datasets_.empty()) {
    mprinterr("Error: analyze statistics: No 1D datasets to analyze.\n");
    return Analysis::ERR;
  }
  // Count number of NOE data sets
  int numNOEsets = 0;
  for (Array1D::const_iterator set = datasets_.begin(); set != datasets_.end(); ++set)
    if ( (*set)->ScalarMode() == DataSet::M_DISTANCE && (*set)->ScalarType() == DataSet::NOE)
     numNOEsets++;
  if (numNOEsets > 0) {
    std::string dsetName = analyzeArgs.GetStringKey("name");
    if (dsetName.empty())
      dsetName = DSLin->GenerateDefaultName("NOE");
    NOE_r6_ = (DataSet_float*)DSLin->AddSetAspect(DataSet::FLOAT, dsetName, "R6");
    NOE_violations_ = (DataSet_integer*)DSLin->AddSetAspect(DataSet::INTEGER, dsetName, 
                                                            "NViolations");
    NOE_avgViolations_ = (DataSet_float*)DSLin->AddSetAspect(DataSet::FLOAT, dsetName, 
                                                             "AvgViolation");
    NOE_names_ = (DataSet_string*)DSLin->AddSetAspect(DataSet::STRING, dsetName, "NOEnames");
    if (NOE_r6_==0 || NOE_violations_==0 || NOE_avgViolations_==0 || NOE_names_==0) {
      mprinterr("Error: Could not set up NOE data sets.\n");
      return Analysis::ERR;
    }
    NOE_r6_->Dim(0).SetLabel("#NOE");
    if (NOE_out != 0) {
      NOE_out->AddSet( NOE_r6_ );
      NOE_out->AddSet( NOE_violations_ );
      NOE_out->AddSet( NOE_avgViolations_ );
      NOE_out->AddSet( NOE_names_ );
    }
  }
  // INFO
  mprintf("    ANALYZE STATISTICS:");
  if (useAllSets)
    mprintf(" Using all data sets (%zu total).\n", datasets_.size());
  else {
    mprintf(" Using %zu data sets:\n", datasets_.size());
    for (Array1D::const_iterator set = datasets_.begin(); set != datasets_.end(); ++set)
      mprintf("\t%s\n", (*set)->Legend().c_str());
  }
  if (shift_ != 0)
    mprintf("\tShift (about %.2f) is begin applied.\n", shift_);
  if (!filename_.empty())
    mprintf("\tOutput to file %s\n", filename_.c_str());
  if (ignore_negative_violations_)
    mprintf("\tIgnoring negative NOE violations.\n");
  mprintf("# SNB = Values from: Schneider, Neidle, and Berman, \"Conformations of the\n"
          "#       Sugar-Phosphate Backbone in Helical DNA Crystal Structures.\",\n"
          "#       Biopolymers (1997), V.42 (1), pp.113-124.\n");
  return Analysis::OK;
}

// Analysis_Statistics::Analyze()
Analysis::RetType Analysis_Statistics::Analyze() {
  if (outfile_.OpenWrite( filename_ )) return Analysis::ERR;
  for (Array1D::const_iterator ds = datasets_.begin(); ds != datasets_.end(); ++ds)
  {
    mprintf("\t'%s'", (*ds)->Legend().c_str());
    (*ds)->ScalarDescription();
    mprintf("\n");
    DataSet_1D const& data_set = static_cast<DataSet_1D const&>( *(*ds) );
    int Nelements = data_set.Size();
    if (Nelements < 1) {
      mprintf("Warning: analyze statistics: No data in dataset %s, skipping.\n",
              data_set.Legend().c_str());
      continue;
    }

    // Compute average and standard deviation with optional shift.
    bool periodic = data_set.IsTorsionArray();
    double average = 0.0;
    double stddev = 0.0;
    for (int i = 0; i < Nelements; ++i) {
      double value = data_set.Dval( i ) - shift_;
      if (periodic) {
        if (value > 180.0)
          value -= 360.0;
        else if (value < -180.0)
          value += 360.0;
      }
      average += value;
      stddev += (value*value);
    }
    average /= Nelements;
    stddev /= Nelements;
    stddev -= (average * average);
    if (stddev > 0)
      stddev = sqrt( stddev );
    else
      stddev = 0;
    average += shift_;

    // Output average/stddev
    outfile_.Printf("__________________________________________________________________\n\n");
    outfile_.Printf("STATISTICS %6s\n", data_set.Legend().c_str());
    outfile_.Printf("   AVERAGE: %8.4f (%.4f stddev)\n", average, stddev);
    outfile_.Printf("   INITIAL: %8.4f\n   FINAL:   %8.4f\n",
                    data_set.Dval( 0 ), data_set.Dval( Nelements-1 ) );

    // More specific analysis based on MODE
    DataSet::scalarMode mode = data_set.ScalarMode();
    if ( mode == DataSet::M_PUCKER) 
      PuckerAnalysis( data_set, Nelements ); 
    else if ( mode == DataSet::M_TORSION)
      TorsionAnalysis( data_set, Nelements );
    else if ( mode == DataSet::M_DISTANCE)
      DistanceAnalysis( data_set, Nelements );

  } // END loop over DataSets

  return Analysis::OK;
}

const char* Analysis_Statistics::pucker_ss[] = {
  "C3'-endo", "C4'-exo ", "O4'-endo", "C1'-exo ", "C2'-endo", "C3'-exo ",
  "C4'-endo", "O4'-exo ", "C1'-endo", "C2'-exo "
};

// Analysis_Statistics::PuckerAnalysis()
void Analysis_Statistics::PuckerAnalysis( DataSet_1D const& ds, int totalFrames ) {
  int pucker_visits[10];
  int pucker_transitions[10][10];
  double pucker_avg[10];
  double pucker_sd[10];
  int curbin, prevbin;

  for (int j = 0; j < 10; ++j) {
    pucker_visits[j] = 0;
    pucker_avg[j] = 0.0;
    pucker_sd[j] = 0.0;
    for (int k=0;k<10;k++) 
      pucker_transitions[j][k] = 0;
  }

  // Get bin for first frame
  double firstvalue = ds.Dval(0);
  if (firstvalue < 0) firstvalue += 360;
  prevbin = firstvalue / 36;
  // Loop over all frames
  for (int i = 0; i < totalFrames; ++i) {
    double value = ds.Dval( i );
    double dval = value;
    if (dval < 0) dval += 360.0;
    curbin = dval / 36;
    if (curbin < 0 || curbin > 9) {
      mprinterr("Error: stat pucker: frame %i has invalid pucker value.\n", i+1);
    } else {
      ++pucker_visits[curbin];
      pucker_avg[curbin] += value;
      pucker_sd[curbin]  += (value*value);
      if (curbin != prevbin) 
        ++pucker_transitions[prevbin][curbin];
      prevbin = curbin;
    }
  }

  if ( ds.ScalarType() == DataSet::PUCKER)
    outfile_.Printf("\n   This is marked as a nucleic acid sugar pucker phase\n");

  outfile_.Printf("\n            %s %s %s %s %s %s %s %s %s %s\n",
                  pucker_ss[0], pucker_ss[1], pucker_ss[2], pucker_ss[3], pucker_ss[4],
                  pucker_ss[5], pucker_ss[6], pucker_ss[7], pucker_ss[8], pucker_ss[9]);
  outfile_.Printf("           -------------------------------------");
  outfile_.Printf("------------------------------------------------------\n");

  for (int j=0; j < 10; j++) {
    if (pucker_visits[j] > 0) {
      pucker_avg[j] /= pucker_visits[j];
      pucker_sd[j]  /= pucker_visits[j];
      pucker_sd[j] = sqrt(pucker_sd[j] - pucker_avg[j]*pucker_avg[j]);
    }
  }

  // OUTPUT
  outfile_.Printf(" %%occupied |");
  for (int j=0; j < 10; j++) {
    if (pucker_visits[j] > 0) {
      double value = pucker_visits[j]*100.0/totalFrames;
      outfile_.Printf(" %6.1f |", value);
    } else
      outfile_.Printf("        |");
  }
  outfile_.Printf("\n");

  outfile_.Printf(" average   |");
  for (int j=0; j < 10; j++) {
    if (pucker_visits[j] > 0) {
      outfile_.Printf(" %6.1f |", pucker_avg[j]);
    } else
      outfile_.Printf("        |");
  }
  outfile_.Printf("\n");

  outfile_.Printf(" stddev    |");
  for (int j=0; j < 10; j++) {
    if (pucker_visits[j] > 1) {
      outfile_.Printf(" %6.1f |", pucker_sd[j]);
    } else
      outfile_.Printf("        |");
  }
  outfile_.Printf("\n           ----------------------------------------------------------");
  outfile_.Printf("---------------------------------\n");

  if (debug_ > 0) {
  outfile_.Printf("\nTRANSITIONS TABLE: (from/vertical to/horizontal)\n\n");
  outfile_.Printf("           %s %s %s %s %s %s %s %s %s %s\n",
                  pucker_ss[0], pucker_ss[1], pucker_ss[2], pucker_ss[3], pucker_ss[4],
                  pucker_ss[5], pucker_ss[6], pucker_ss[7], pucker_ss[8], pucker_ss[9]);
  outfile_.Printf("           ------------------------------------------");
  outfile_.Printf("-------------------------------------------------\n");
  for (int j=0; j<10; j++) {
    outfile_.Printf("  %s |", pucker_ss[j]);
    for (int k=0; k<10; k++) {
      if (pucker_transitions[j][k] > 0)
        outfile_.Printf(" %6i |", pucker_transitions[j][k]);
      else
        outfile_.Printf("        |");
    }
    outfile_.Printf("\n");
  }
  outfile_.Printf("           ----------------------------------------------------------");
  outfile_.Printf("---------------------------------\n\n");
  }
}

const char* Analysis_Statistics::torsion_ss[] = {
  "g+     ", "a+     ", "t      ", "a-     ", "g-     ", "c      "
};

const double Analysis_Statistics::torsion_offset[] = { 
  0.0, 0.0, 0.0, 0.0, 0.0, 180.0 
};

// Analysis_Statistics::TorsionAnalysis()
void Analysis_Statistics::TorsionAnalysis(DataSet_1D const& ds, int totalFrames) {
  int torsion_visits[6];
  int torsion_transitions[6][6];
  double torsion_avg[6]; 
  double torsion_sd[6];
  int prevbin, curbin;

  for (int j=0;j<6;j++) {
    torsion_visits[j] = 0;
    torsion_avg[j] = 0.0;
    torsion_sd[j] = 0.0;
    for (int k=0;k<6;k++) 
      torsion_transitions[j][k] = 0;
  }
 
  // Get initial bin
  double firstvalue = ds.Dval( 0 );
  if (firstvalue < 0) firstvalue += 360;
  prevbin = (int) (firstvalue - 30.0) / 60;
  // Loop over all frames
  for (int i = 0; i < totalFrames; ++i) {
    double value = ds.Dval( i );
    double dval = value;
    if (dval < 0) dval += 360;
    curbin = (int) (dval - 30.0) / 60;
    if (curbin < 0 || curbin > 5) {
      mprinterr("Error: stat torsion: frame %i has invalid torsion value.\n", i+1);
    } else {
      ++torsion_visits[curbin];
      value += torsion_offset[curbin];
      // Fix for trans averaging
      if (value < -150.0) value += 360.0;
      torsion_avg[curbin] += value;
      torsion_sd[curbin]  += (value*value);
      if (curbin != prevbin) 
        ++torsion_transitions[prevbin][curbin];
      prevbin = curbin;
    }
  }

  // OUTPUT
  outfile_.Printf("\n               %s  %s  %s  %s  %s  %s\n",
          torsion_ss[0], torsion_ss[1], torsion_ss[2],
          torsion_ss[3], torsion_ss[4], torsion_ss[5]);
  outfile_.Printf("           ---------------");
  outfile_.Printf("----------------------------------------\n");

  for (int j=0; j < 6; j++) {
    if (torsion_visits[j] > 0) {
      torsion_avg[j] /= torsion_visits[j];
      torsion_sd[j]  /= torsion_visits[j];
      torsion_sd[j] = sqrt(torsion_sd[j] - torsion_avg[j]*torsion_avg[j]);
      torsion_avg[j] -= torsion_offset[j];
    }
  }

  outfile_.Printf(" %%occupied |");
  for (int j=0; j < 6; j++) {
    if (torsion_visits[j] > 0) {
      double value = torsion_visits[j]*100.0/totalFrames;
      outfile_.Printf(" %6.1f |", value);
    } else
      outfile_.Printf("        |");
  }
  outfile_.Printf("\n");

  outfile_.Printf(" average   |");
  for (int j=0; j < 6; j++) {
    if (torsion_visits[j] > 0) {
      outfile_.Printf(" %6.1f |", torsion_avg[j]);
    } else
      outfile_.Printf("        |");
  }
  outfile_.Printf("\n");

  outfile_.Printf(" stddev    |");
  for (int j=0; j < 6; j++) {
    if (torsion_visits[j] > 1) {
      outfile_.Printf(" %6.1f |", torsion_sd[j]);
    } else
      outfile_.Printf("        |");
  }
  outfile_.Printf("\n           --------------------------");
  outfile_.Printf("-----------------------------\n");

  // Specific torsion types
  switch ( ds.ScalarType() ) {
    case DataSet::ALPHA:
      //              "               g+       a+       t        a-       g-       c 
      outfile_.Printf(" ALPHA       minor             minor            canonical\n");
      outfile_.Printf("\n   O3'-P-O5'-C5', SNB range is 270-300 deg (g-)\n");
      if ( (torsion_visits[0] + torsion_visits[1] + torsion_visits[2] + torsion_visits[5] )
           > (totalFrames * 0.1) )
        outfile_.Printf("   *** > 10%% out of range population detected\n");
      break;

    case DataSet::BETA:
      //              "               g+       a+       t        a-       g-       c 
      outfile_.Printf(" BETA                <-- canonical -->\n");

      outfile_.Printf("\n   P-O5'-C5'-C4', SNB range is 130-200 deg (a+,t)\n");
      if ( (torsion_visits[0] + torsion_visits[3] + torsion_visits[4] + torsion_visits[5] )
           > (totalFrames * 0.05) )
        outfile_.Printf("   *** > 5%% out of range population detected\n");
      break;

    case DataSet::GAMMA:
      //              "               g+       a+       t        a-       g-       c 
      outfile_.Printf(" GAMMA     canonical           minor             minor\n");
      outfile_.Printf("\n   O5'-C5'-C4'-C3', SNB range is 20-80 (g+)\n");
      if (torsion_visits[2] > (totalFrames* 0.1))
        outfile_.Printf("   *** GAMMA trans > 10%% detected!!!\n");
      break;

    case DataSet::DELTA:
      //              "               g+       a+       t        a-       g-       c 
      outfile_.Printf(" DELTA      <------ canonical ------>\n");
      outfile_.Printf("\n   C5'-C4'-C3'-O3', SNB range is 70-180\n");
      outfile_.Printf("   DNA: ~128 with BI (a+), ~144 with BII (a+)\n");
      if ( (torsion_visits[0] + torsion_visits[3] + torsion_visits[4] + torsion_visits[5] )
           > (totalFrames * 0.05) )
        outfile_.Printf("   *** > 5%% out of range population detected\n");
      break;

    case DataSet::EPSILON:
      //              "               g+       a+       t        a-       g-       c 
      outfile_.Printf(" EPSILON                         BI       BII\n");
      outfile_.Printf("\n   C4'-C3'-O3'-P, SNB range is 160-270\n");
      outfile_.Printf("   BI = %6.2f%% (~184), BII = %6.2f%% (~246)\n",
              (torsion_visits[2]*100.0)/totalFrames,
              (torsion_visits[3]*100.0)/totalFrames);
      if ( (torsion_visits[0] + torsion_visits[1] + torsion_visits[4] + torsion_visits[5] )
           > (totalFrames * 0.05) )
        outfile_.Printf("   *** > 5%% out of range population detected\n");
      break;

    case DataSet::ZETA:
      //              "               g+       a+       t        a-       g-       c 
      outfile_.Printf(" ZETA                <----- BII ------------- BI ----->\n");
      outfile_.Printf("\n   C3'-O3'-P-O5', SNB range is 230-300 (BI), 150-210 (BII)\n");
      outfile_.Printf("   BI = %6.2f%% (~265, a-/g-), BII = %6.2f%% (~174, a+/t)\n",
              (torsion_visits[3]+torsion_visits[4])*100.0/totalFrames,
              (torsion_visits[1]+torsion_visits[2])*100.0/totalFrames);
      if ( (torsion_visits[0] + torsion_visits[5] )
           > (totalFrames * 0.05) )
        outfile_.Printf("   *** > 5%% out of range population detected\n");
      break;

    case DataSet::CHI:
      //              "               g+       a+       t        a-       g-       c 
      outfile_.Printf(" CHI                         <-------- anti ------->  <--syn---\n");
      outfile_.Printf("\n   O4'-C1'-NX-CX, SNB range is 200-300\n");
      if ( (torsion_visits[0] + torsion_visits[5] ) > (totalFrames * 0.05) )
        outfile_.Printf(
          "   *** CHI flips; > 5%% out of range populations detected (see table below)\n");
      if ( torsion_visits[1] > (totalFrames * 0.05) )
        outfile_.Printf("   *** Unexpected CHI population in a+ region, > 5%%\n");
      break;

    case DataSet::C2P:
      //              "               g+       a+       t        a-       g-       c 
      outfile_.Printf(" C2' to base      in\n");
      outfile_.Printf("\n   C2'-C1'-NX-CX\n\n");
      break;

    case DataSet::H1P:
      //              "               g+       a+       t        a-       g-       c 
      outfile_.Printf(" H1'       below-plane                           above      in\n");
      outfile_.Printf("\n   H1'-C1'-NX-CX, > 0 H1' below plane (check if sugar in plane)\n\n");
      break;

    default: break;
  } // END switch over dataset scalartype

  if (debug_ > 0) {
    outfile_.Printf("\nTRANSITIONS TABLE: (from/vertical to/horizontal)\n\n");
    outfile_.Printf("              %s  %s  %s  %s  %s  %s\n",
            torsion_ss[0], torsion_ss[1], torsion_ss[2],
            torsion_ss[3], torsion_ss[4], torsion_ss[5]);
    outfile_.Printf("           -----------------------");
    outfile_.Printf("--------------------------------\n");
    for (int j=0; j<6; j++) {
      outfile_.Printf("   %s |", torsion_ss[j]);
      for (int k=0; k<6; k++) {
        if (torsion_transitions[j][k] > 0)
          outfile_.Printf(" %6i |", torsion_transitions[j][k]);
        else
          outfile_.Printf("        |");
      }
      outfile_.Printf("\n");
    }
    outfile_.Printf("           ------------------");
    outfile_.Printf("-------------------------------------\n\n");
  }
}

const char* Analysis_Statistics::distance_ss[] = {
  " < 2.5 ", "2.5-3.5", "3.5-4.5", "4.5-5.5", "5.5-6.5", " > 6.5 "
};

static inline int distbin(double val) {
  int bin = val - 1.5;
  if (bin < 0)
    bin = 0;
  else if (bin > 5)
    bin = 5;
  return bin;
}

// Analysis_Statistics::DistanceAnalysis()
void Analysis_Statistics::DistanceAnalysis( DataSet_1D const& ds, int totalFrames )
{
  int distance_visits[6];
  int distance_transitions[6][6];
  double distance_avg[6];
  double distance_sd[6];
  double average, bound=0.0, boundh=0.0, rexp=0.0;
  int prevbin, curbin, Nb, Nh;

  for (int j=0;j<6;j++) {
    distance_visits[j] = 0;
    distance_avg[j] = 0.0;
    distance_sd[j] = 0.0;
    for (int k=0;k<6;k++) 
      distance_transitions[j][k] = 0;
  }

  // Init for NOE
  bool isNOE = (ds.ScalarType() == DataSet::NOE);
  if (isNOE) {
    outfile_.Printf("   NOE SERIES: S < 2.9, M < 3.5, W < 5.0, blank otherwise.\n    |");
    average = 0;
    Nb = 0;
    Nh = 0;
    DataSet_double const& DsDbl = static_cast<DataSet_double const&>( ds );
    bound = DsDbl.NOE_bound();
    boundh = DsDbl.NOE_boundH();
    rexp = DsDbl.NOE_rexp();
    if (rexp < 0.0) {
      // If lower bound is zero just use boundh, otherwise use avg.
      if (bound > 0.0)
        rexp = (bound + boundh) / 2.0;
      else
        rexp = boundh;
    }
  }

  // Get bin for first value
  prevbin = distbin( ds.Dval(0) );
  // Loop over all frames
  double r6_avg = 0.0;
  int NlowViolations = 0;
  int NhighViolations = 0;
  for (int i = 0; i < totalFrames; ++i) {
    double value = ds.Dval( i );
    curbin = distbin( value );
    ++distance_visits[curbin];
    distance_avg[curbin] += value;
    distance_sd[curbin]  += (value*value);
    if (curbin != prevbin) 
      ++distance_transitions[prevbin][curbin];
    prevbin = curbin;

    // NOE calc
    if (isNOE) {
      r6_avg += pow( value, -6.0 );
      int j = totalFrames / 50.0;
      if (j < 1) j = 1;

      if (value < bound) {
        ++Nb;
        ++NlowViolations;
      }
      if (value < boundh)
        ++Nh;
      else
        ++NhighViolations;

      average += value;
      if (j == 1 || i % j == 1) {
        average /= j;
        if (average < 2.9) {
          outfile_.Printf("S");
        } else if (average < 3.5) {
          outfile_.Printf("M");
        } else if (average < 5.0) {
          outfile_.Printf("W");
        } else {
          outfile_.Printf(" ");
        }
        average = 0.0;
      }
    }
  }

  double d_total = (double)totalFrames;
  // NOE printout
  if (isNOE) {
    outfile_.Printf("|\n");
    if (bound > 0.0) {
      outfile_.Printf("   NOE < %.2f for %.2f%% of the time\n",
              bound, (double)Nb / d_total * 100.0);
    }
    if (boundh > 0.0) {
      outfile_.Printf("   NOE < %.2f for %.2f%% of the time\n",
              boundh, (double)Nh / d_total * 100.0);
    }
    r6_avg /= d_total;
    r6_avg = pow( r6_avg, -1.0/6.0 );
    outfile_.Printf("   NOE <r^-6>^(-1/6)= %.4f\n", r6_avg);
    NOE_r6_->AddElement( (float)r6_avg );
    int total_violations = 0;
    if (bound >= 0.0 && boundh > 0.0) {
      total_violations = NlowViolations+NhighViolations;
      outfile_.Printf("   #Violations: Low= %i High= %i Total= %i\n",
                      NlowViolations, NhighViolations, total_violations);
    }
    NOE_violations_->AddElement( total_violations );
    double avg_violation = 0.0;
    if (rexp > 0.0) {
      avg_violation = r6_avg - rexp;
      if (ignore_negative_violations_ && avg_violation < 0.0)
        avg_violation = 0.0;
      outfile_.Printf("   Rexp= %.4f <Violation>= %.4f\n", rexp, avg_violation);
    }
    NOE_avgViolations_->AddElement( (float)avg_violation );
    NOE_names_->AddElement( "\"" + ds.Legend() + "\"" );
  }

  // OUTPUT
  outfile_.Printf("\n              %s  %s  %s  %s  %s  %s\n",
                  distance_ss[0], distance_ss[1], distance_ss[2],
                  distance_ss[3], distance_ss[4], distance_ss[5]);
  outfile_.Printf("           ---------------");
  outfile_.Printf("----------------------------------------\n");

  for (int j=0; j < 6; j++) {
    if (distance_visits[j] > 0) {
      distance_avg[j] /= distance_visits[j];
      distance_sd[j]  /= distance_visits[j];
      distance_sd[j] = sqrt(distance_sd[j] - distance_avg[j]*distance_avg[j]);
    }
  }

  outfile_.Printf(" %%occupied |");
  for (int j=0; j < 6; j++) {
    if (distance_visits[j] > 0) {
      double value = distance_visits[j]*100.0/d_total;
      outfile_.Printf(" %6.1f |", value);
    } else
      outfile_.Printf("        |");
  }
  outfile_.Printf("\n");

  outfile_.Printf(" average   |");
  for (int j=0; j < 6; j++) {
    if (distance_visits[j] > 0) {
      outfile_.Printf(" %6.3f |", distance_avg[j]);
    } else
      outfile_.Printf("        |");
  }
  outfile_.Printf("\n");

  outfile_.Printf(" stddev    |");
  for (int j=0; j < 6; j++) {
    if (distance_visits[j] > 1) {
      outfile_.Printf(" %6.3f |", distance_sd[j]);
    } else
      outfile_.Printf("        |");
  }
  outfile_.Printf("\n           --------------------------");
  outfile_.Printf("-----------------------------\n");

  if (debug_ > 0) {
    outfile_.Printf("\nTRANSITIONS TABLE: (from/vertical to/horizontal)\n\n");
    outfile_.Printf("            %s  %s  %s  %s  %s  %s\n",
                    distance_ss[0], distance_ss[1], distance_ss[2],
                    distance_ss[3], distance_ss[4], distance_ss[5]);
    outfile_.Printf("           -----------------------");
    outfile_.Printf("--------------------------------\n");
    for (int j=0; j<6; j++) {
      outfile_.Printf("   %s |", distance_ss[j]);
      for (int k=0; k<6; k++) {
        if (distance_transitions[j][k] > 0)
          outfile_.Printf(" %6i |", distance_transitions[j][k]);
        else
          outfile_.Printf("        |");
      }
      outfile_.Printf("\n");
    }
    outfile_.Printf("           ------------------");
    outfile_.Printf("-------------------------------------\n\n");
  }

}
