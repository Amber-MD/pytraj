// Analysis_Clustering
#include "Analysis_Clustering.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // fileExists, integerToString
#include "DataSet_integer.h" // For converting cnumvtime
#include "Trajout.h"
#include "Timer.h"
// Clustering Algorithms
#include "Cluster_HierAgglo.h"
#include "Cluster_DBSCAN.h"
#include "Cluster_Kmeans.h"
#include "Cluster_ReadInfo.h"
#include "Cluster_DPeaks.h"

// CONSTRUCTOR
Analysis_Clustering::Analysis_Clustering() :
  masterDSL_(0),
  coords_(0),
  CList_(0),
  sieve_(1),
  sieveSeed_(-1),
  windowSize_(0),
  drawGraph_(0),
  draw_maxit_(0),
  draw_tol_(0.0),
  cnumvtime_(0),
  clustersVtime_(0),
  cpopvtimefile_(0),
  nofitrms_(false),
  metric_(ClusterList::RMS),
  useMass_(false),
  grace_color_(false),
  norm_pop_(NONE),
  load_pair_(false),
  calc_lifetimes_(false),
  writeRepFrameNum_(false),
  clusterfmt_(TrajectoryFile::UNKNOWN_TRAJ),
  singlerepfmt_(TrajectoryFile::UNKNOWN_TRAJ),
  reptrajfmt_(TrajectoryFile::UNKNOWN_TRAJ),
  debug_(0)
{ } 

// DESTRUCTOR
Analysis_Clustering::~Analysis_Clustering() {
  if (CList_ != 0) delete CList_;
}

void Analysis_Clustering::Help() {
  mprintf("\t[crdset <crd set>]\n");
  mprintf("  Algorithms:\n");
  Cluster_HierAgglo::Help();
  Cluster_DBSCAN::Help();
  Cluster_DPeaks::Help();
  Cluster_Kmeans::Help();
  Cluster_ReadInfo::Help();
  mprintf("  Distance metric options: {rms | srmsd | dme | data}\n"
          "\t{ [[rms | srmsd] [<mask>] [mass] [nofit]] | [dme [<mask>]] |\n"
          "\t   [data <dset0>[,<dset1>,...]] }\n"
          "\t[sieve <#> [random [sieveseed <#>]]] [loadpairdist] [savepairdist] [pairdist <file>]\n"
          "  Output options:\n"
          "\t[out <cnumvtime>] [gracecolor] [summary <summaryfile>] [info <infofile>]\n"
          "\t[summarysplit <splitfile>] [splitframe <comma-separated frame list>]\n"
          "\t[clustersvtime <filename> cvtwindow <window size>]\n"
          "\t[cpopvtime <file> [normpop | normframe]] [lifetime]\n"
          "\t[sil <silhouette file prefix>]\n"
          "  Coordinate output options:\n"
          "\t[ clusterout <trajfileprefix> [clusterfmt <trajformat>] ]\n"
          "\t[ singlerepout <trajfilename> [singlerepfmt <trajformat>] ]\n"
          "\t[ repout <repprefix> [repfmt <repfmt>] [repframe] ]\n"
          "\t[ avgout <avgprefix> [avgfmt <avgfmt>] ]\n"
          "  Experimental options:\n"
          "\t[[drawgraph | drawgraph3d] [draw_tol <tolerance>] [draw_maxit <iterations]]\n"
          "  Cluster structures based on coordinates (RMSD/DME) or given data set(s).\n"
          "  <crd set> can be created with the 'createcrd' command.\n");
}

const char* Analysis_Clustering::PAIRDISTFILE = "CpptrajPairDist";

Analysis::RetType Analysis_Clustering::Setup(ArgList& analyzeArgs, DataSetList* datasetlist,
                            TopologyList* PFLin, DataFileList* DFLin, int debugIn)
{
  debug_ = debugIn;
  // Attempt to get coords dataset from datasetlist
  std::string setname = analyzeArgs.GetStringKey("crdset");
  coords_ = (DataSet_Coords*)datasetlist->FindCoordsSet( setname );
  if (coords_ == 0) {
    mprinterr("Error: clustering: Could not locate COORDS set corresponding to %s\n",
              setname.c_str());
    return Analysis::ERR;
  }
  // Check for DataSet(s) to cluster on, otherwise coords will be used
  cluster_dataset_.clear();
  setname = analyzeArgs.GetStringKey("data");
  metric_ = ClusterList::RMS;
  if (!setname.empty()) {
    ArgList dsnames(setname, ",");
    DataSetList inputDsets;
    for (ArgList::const_iterator name = dsnames.begin(); name != dsnames.end(); ++name) {
      DataSetList tempDSL = datasetlist->GetMultipleSets( *name );
      if (tempDSL.empty()) {
        mprinterr("Error: cluster: %s did not correspond to any data sets.\n");
        return Analysis::ERR;
      }
      inputDsets += tempDSL;
    }
    for (DataSetList::const_iterator ds = inputDsets.begin(); ds != inputDsets.end(); ++ds) {
      // Clustering only allowed on 1D data sets.
      if ( (*ds)->Ndim() != 1 ) {
        mprinterr("Error: Clustering only allowed on 1D data sets, %s is %zuD.\n",
                  (*ds)->Legend().c_str(), (*ds)->Ndim());
        return Analysis::ERR;
      }
      cluster_dataset_.push_back( *ds );
    }
    metric_ = ClusterList::DATA;
  } else {
    int usedme = (int)analyzeArgs.hasKey("dme");
    int userms = (int)analyzeArgs.hasKey("rms");
    int usesrms = (int)analyzeArgs.hasKey("srmsd");
    if (usedme + userms + usesrms > 1) {
      mprinterr("Error: Specify either 'dme', 'rms', or 'srmsd'.\n");
      return Analysis::ERR;
    }
    if      (usedme)  metric_ = ClusterList::DME;
    else if (userms)  metric_ = ClusterList::RMS;
    else if (usesrms) metric_ = ClusterList::SRMSD;
  }
  // Get clustering algorithm
  if (CList_ != 0) delete CList_;
  CList_ = 0;
  if (analyzeArgs.hasKey("hieragglo"))   CList_ = new Cluster_HierAgglo(); 
  else if (analyzeArgs.hasKey("dbscan")) CList_ = new Cluster_DBSCAN();
  else if (analyzeArgs.hasKey("dpeaks")) CList_ = new Cluster_DPeaks();
  else if (analyzeArgs.hasKey("kmeans") ||
           analyzeArgs.hasKey("means" )) CList_ = new Cluster_Kmeans();
  else if (analyzeArgs.hasKey("readinfo") ||
           analyzeArgs.hasKey("readtxt")) CList_ = new Cluster_ReadInfo(); 
  else {
    mprintf("Warning: No clustering algorithm specified; defaulting to 'hieragglo'\n");
    CList_ = new Cluster_HierAgglo();
  }
  if (CList_ == 0) return Analysis::ERR;
  CList_->SetDebug(debug_);
  // Get algorithm-specific keywords
  if (CList_->SetupCluster( analyzeArgs )) return Analysis::ERR; 
  // Get keywords
  useMass_ = analyzeArgs.hasKey("mass");
  sieveSeed_ = analyzeArgs.getKeyInt("sieveseed", -1);
  sieve_ = analyzeArgs.getKeyInt("sieve", 1);
  if (sieve_ < 1) {
    mprinterr("Error: 'sieve <#>' must be >= 1 (%i)\n", sieve_);
    return Analysis::ERR;
  }
  if (analyzeArgs.hasKey("random") && sieve_ > 1)
    sieve_ = -sieve_; // negative # indicates random sieve
  halffile_ = analyzeArgs.GetStringKey("summarysplit");
  if (halffile_.empty()) // For backwards compat.
    halffile_ = analyzeArgs.GetStringKey("summaryhalf");
  if (!halffile_.empty()) {
    ArgList splits( analyzeArgs.GetStringKey("splitframe"), "," );
    if (!splits.empty()) {
      splitFrames_.clear();
      int sf = splits.getNextInteger(-1); // User frame #s start at 1
      while (sf > 0) {
        splitFrames_.push_back( sf );
        sf = splits.getNextInteger(-1);
      }
      if ((int)splitFrames_.size() < splits.Nargs()) {
        mprinterr("Error: Invalid split frame arguments.\n");
        splits.CheckForMoreArgs();
        return Analysis::ERR;
      }
    }
  }
  if (analyzeArgs.hasKey("drawgraph"))
    drawGraph_ = 1;
  else if (analyzeArgs.hasKey("drawgraph3d"))
    drawGraph_ = 2;
  else
    drawGraph_ = 0;
  draw_maxit_ = analyzeArgs.getKeyInt("draw_maxit", 1000);
  draw_tol_ = analyzeArgs.getKeyDouble("draw_tol", 1.0E-5);
  
  DataFile* cnumvtimefile = DFLin->AddDataFile(analyzeArgs.GetStringKey("out"), analyzeArgs);
  DataFile* clustersvtimefile = DFLin->AddDataFile(analyzeArgs.GetStringKey("clustersvtime"),
                                                   analyzeArgs);
  windowSize_ = analyzeArgs.getKeyInt("cvtwindow", 0);
  cpopvtimefile_ = DFLin->AddDataFile(analyzeArgs.GetStringKey("cpopvtime"), analyzeArgs);
  clusterinfo_ = analyzeArgs.GetStringKey("info");
  summaryfile_ = analyzeArgs.GetStringKey("summary");
  nofitrms_ = analyzeArgs.hasKey("nofit");
  grace_color_ = analyzeArgs.hasKey("gracecolor");
  calc_lifetimes_ = analyzeArgs.hasKey("lifetime");
  if (cpopvtimefile_ != 0) {
    if (analyzeArgs.hasKey("normpop"))
      norm_pop_ = CLUSTERPOP;
    else if (analyzeArgs.hasKey("normframe"))
      norm_pop_ = FRAME;
    else
      norm_pop_ = NONE;
  }
  sil_file_ = analyzeArgs.GetStringKey("sil");
  // Options for loading/saving pairwise distance file
  load_pair_ = analyzeArgs.hasKey("loadpairdist");
  bool save_pair = analyzeArgs.hasKey("savepairdist");
  pairdistfile_ = analyzeArgs.GetStringKey("pairdist");
  if ( (load_pair_ || save_pair) && pairdistfile_.empty() )
    pairdistfile_.assign(PAIRDISTFILE);
  else if (!pairdistfile_.empty())
    load_pair_ = true;
  // Output trajectory stuff
  clusterfile_ = analyzeArgs.GetStringKey("clusterout");
  clusterfmt_ = TrajectoryFile::GetFormatFromString( analyzeArgs.GetStringKey("clusterfmt") ); 
  singlerepfile_ = analyzeArgs.GetStringKey("singlerepout");
  singlerepfmt_ = TrajectoryFile::GetFormatFromString( analyzeArgs.GetStringKey("singlerepfmt") );
  reptrajfile_ = analyzeArgs.GetStringKey("repout");
  reptrajfmt_ = TrajectoryFile::GetFormatFromString( analyzeArgs.GetStringKey("repfmt") );
  writeRepFrameNum_ = analyzeArgs.hasKey("repframe");
  avgfile_ = analyzeArgs.GetStringKey("avgout");
  avgfmt_ = TrajectoryFile::GetFormatFromString( analyzeArgs.GetStringKey("avgfmt") );
  // Get the mask string 
  maskexpr_ = analyzeArgs.GetMaskNext();

  // Dataset to store cluster number v time
  cnumvtime_ = datasetlist->AddSet(DataSet::INTEGER, analyzeArgs.GetStringNext(), "Cnum");
  if (cnumvtime_==0) return Analysis::ERR;
  if (cnumvtimefile != 0) cnumvtimefile->AddSet( cnumvtime_ );
  // DataSet for # clusters seen v time
  if (clustersvtimefile != 0) {
    if (windowSize_ < 2) {
      mprinterr("Error: For # clusters seen vs time, cvtwindow must be specified and > 1\n");
      return Analysis::ERR;
    }
    clustersVtime_ = datasetlist->AddSetAspect(DataSet::INTEGER, cnumvtime_->Name(), "NCVT");
    if (clustersVtime_ == 0) return Analysis::ERR;
    clustersvtimefile->AddSet( clustersVtime_ );
  }
  // Save master DSL for Cpopvtime
  masterDSL_ = datasetlist;

  mprintf("    CLUSTER: Using coords dataset %s, clustering using", coords_->Legend().c_str());
  if ( metric_ != ClusterList::DATA ) {
    mprintf(" %s", ClusterList::MetricString( metric_ ));
    if (!maskexpr_.empty())
      mprintf(" (mask [%s])",maskexpr_.c_str());
    else
      mprintf(" (all atoms)");
    if (useMass_)
      mprintf(", mass-weighted");
    if (nofitrms_)
      mprintf(", no fitting");
    else
      mprintf(" best-fit");
  } else {
    if (cluster_dataset_.size() == 1)
      mprintf(" dataset %s", cluster_dataset_[0]->Legend().c_str());
    else
      mprintf(" %u datasets.", cluster_dataset_.size());
  }
  mprintf("\n");
  CList_->ClusteringInfo();
  if (sieve_ > 1)
    mprintf("\tInitial clustering sieve value is %i frames.\n", sieve_);
  else if (sieve_ < -1) {
    mprintf("\tInitial clustering will be randomly sieved (with value %i)", -sieve_);
    if (sieveSeed_ > 0) mprintf(" using random seed %i", sieveSeed_);
    mprintf(".\n");
  }
  if (cnumvtimefile != 0)
    mprintf("\tCluster # vs time will be written to %s\n", cnumvtimefile->DataFilename().base());
  if (clustersvtimefile != 0)
    mprintf("\t# clusters seen vs time will be written to %s\n",
            clustersvtimefile->DataFilename().base());
  if (cpopvtimefile_ != 0) {
    mprintf("\tCluster pop vs time will be written to %s", cpopvtimefile_->DataFilename().base());
    if (norm_pop_==CLUSTERPOP) mprintf(" (normalized by cluster size)");
    else if (norm_pop_==FRAME) mprintf(" (normalized by frame)");
    mprintf("\n");
  }
  if (grace_color_)
    mprintf("\tGrace color instead of cluster number (1-15) will be saved.\n");
  if (calc_lifetimes_)
    mprintf("\tCluster lifetime data sets will be calculated.\n");
  if (load_pair_)
    mprintf("\tPreviously calcd pair distances %s will be used if found.\n",
            pairdistfile_.c_str());
  if (!clusterinfo_.empty())
    mprintf("\tCluster information will be written to %s\n",clusterinfo_.c_str());
  if (!summaryfile_.empty())
    mprintf("\tSummary of cluster results will be written to %s\n",summaryfile_.c_str());
  if (!sil_file_.empty()) {
    mprintf("\tFrame silhouettes will be written to %s.frame.dat, cluster silhouettes\n"
            "\t  will be written to %s.cluster.dat\n", sil_file_.c_str(), sil_file_.c_str());
    if (sieve_ > 1)
      mprintf("\tSilhouette calculation will use sieved frames ONLY.\n");
  }
  if (!halffile_.empty()) {
    mprintf("\tSummary comparing parts of trajectory data for clusters will be written to %s\n",
            halffile_.c_str());
    if (!splitFrames_.empty()) {
      mprintf("\t\tFrames will be split at:");
      for (std::vector<int>::const_iterator f = splitFrames_.begin(); f != splitFrames_.end(); ++f)
        mprintf(" %i", *f);
      mprintf("\n");
    } else
      mprintf("\t\tFrames will be split at the halfway point.\n");
  }
  if (!clusterfile_.empty())
    mprintf("\tCluster trajectories will be written to %s, format %s\n",
            clusterfile_.c_str(), TrajectoryFile::FormatString(clusterfmt_));
  if (!singlerepfile_.empty())
    mprintf("\tCluster representatives will be written to 1 traj (%s), format %s\n",
            singlerepfile_.c_str(), TrajectoryFile::FormatString(singlerepfmt_));
  if (!reptrajfile_.empty()) {
    mprintf("\tCluster representatives will be written to separate trajectories,\n");
    mprintf("\t\tprefix (%s), format %s",reptrajfile_.c_str(), 
            TrajectoryFile::FormatString(reptrajfmt_));
    if (writeRepFrameNum_) mprintf(", with frame #s");
    mprintf("\n");
  }
  if (!avgfile_.empty())
    mprintf("\tAverage structures for clusters will be written to %s, format %s\n",
            avgfile_.c_str(), TrajectoryFile::FormatString(avgfmt_));
  if (drawGraph_ > 0)
    mprintf("\tEXPERIMENTAL: Force-directed graph will be drawn from pairwise distances.\n"
            "\t              Max iterations= %i, min tolerance= %g\n",
                             draw_maxit_, draw_tol_);

  return Analysis::OK;
}

/** This is where the clustering is actually performed. First the distances
  * between each frame are calculated. Then the clustering routine is called.
  */
// TODO: Need to update save to indicate distance type
// NOTE: Should distances be saved only if load_pair?
Analysis::RetType Analysis_Clustering::Analyze() {
  Timer cluster_setup;
  Timer cluster_pairwise;
  Timer cluster_cluster;
  Timer cluster_post;
  Timer cluster_total;
  cluster_total.Start();
  mprintf("\tStarting clustering.\n");
  cluster_setup.Start();
  // Default: USE_FRAMES  - Calculate pair distances from frames.
  //          USE_FILE    - If pairdistfile exists, load pair distances from there.
  // Calculated distances will be saved if not loaded from file.
  ClusterList::DistModeType pairdist_mode = ClusterList::USE_FRAMES; 
  if (load_pair_ && fileExists(pairdistfile_))
    pairdist_mode = ClusterList::USE_FILE;
  // If no dataset specified, use COORDS
  if (cluster_dataset_.empty())
     cluster_dataset_.push_back( (DataSet*)coords_ );
  // Test that cluster data set contains data
  // FIXME make unsigned
  int clusterDataSetSize = (int)cluster_dataset_[0]->Size();
  if (clusterDataSetSize < 1) {
    mprinterr("Error: cluster data set %s does not contain data.\n", 
              cluster_dataset_[0]->Legend().c_str());
    return Analysis::ERR;
  }
  // If more than one data set, make sure they are all the same size.
  for (ClusterDist::DsArray::iterator ds = cluster_dataset_.begin();
                                      ds != cluster_dataset_.end(); ++ds)
  {
    if ((int)(*ds)->Size() != clusterDataSetSize) {
      mprinterr("Error: data set %s size (%i) != first data set %s size (%i)\n",
                (*ds)->Legend().c_str(), (*ds)->Size(), 
                cluster_dataset_[0]->Legend().c_str(), clusterDataSetSize);
      return Analysis::ERR;
    }
  }
  // If no coordinates were specified, disable coordinate output types
  bool has_coords = true;
  if (coords_->Size() < 1) {
    mprintf("Warning: Associated coordinate data set is empty.\n"
            "Warning: Disabling coordinate output.\n");
    has_coords = false;
  }
  cluster_setup.Stop();
  // Calculate distances between frames
  cluster_pairwise.Start();
  if (CList_->CalcFrameDistances( pairdistfile_, cluster_dataset_, pairdist_mode,
                                  metric_, nofitrms_, useMass_, maskexpr_, 
                                  sieve_, sieveSeed_ ))
    return Analysis::ERR;
  cluster_pairwise.Stop();
  // Cluster
  cluster_cluster.Start();
  CList_->Cluster();
  cluster_cluster.Stop();
  cluster_post.Start();
  if (CList_->Nclusters() > 0) {
    // Sort clusters and renumber; also finds centroids for printing
    // representative frames. If sieving, add remaining frames.
    CList_->Renumber( (sieve_ != 1) );
    // DEBUG
    if (debug_ > 0) {
      mprintf("\nFINAL CLUSTERS:\n");
      CList_->PrintClusters();
    }

    // Print ptraj-like cluster info. If no filename is written some info will
    // still be written to STDOUT.
    CList_->PrintClustersToFile(clusterinfo_, clusterDataSetSize);

    // Calculate cluster silhouette
    if (!sil_file_.empty())
      CList_->CalcSilhouette( sil_file_ );

    // Print a summary of clusters
    if (!summaryfile_.empty())
      CList_->Summary(summaryfile_, clusterDataSetSize);

    // Print a summary comparing first half to second half of data for clusters
    if (!halffile_.empty()) {
      // If no split frames were specified, use halfway point.
      if (splitFrames_.empty())
        splitFrames_.push_back( clusterDataSetSize / 2 );
      // Check that none of the split values are invalid.
      std::vector<int> actualSplitFrames;
      for (std::vector<int>::const_iterator f = splitFrames_.begin();
                                            f != splitFrames_.end(); ++f)
        if ( *f < 1 || *f >= clusterDataSetSize )
          mprintf("Warning: split frame %i is out of bounds; ignoring.\n", *f);
        else
          actualSplitFrames.push_back( *f );
      CList_->Summary_Part(halffile_, clusterDataSetSize, actualSplitFrames);
    }

    // Create cluster v time data from clusters.
    CreateCnumvtime( *CList_, clusterDataSetSize );

    // TEST: Draw graph based on point distances
    if (drawGraph_ > 0)
     CList_->DrawGraph( drawGraph_ == 2, cnumvtime_, draw_tol_, draw_maxit_ );

    // Create # clusters seen v time data.
    if (clustersVtime_ != 0)
      NclustersObserved( *CList_, clusterDataSetSize );

    // Create cluster pop v time plots
    if (cpopvtimefile_ != 0)
      CreateCpopvtime( *CList_, clusterDataSetSize );

    // Create cluster lifetime DataSets
    if (calc_lifetimes_)
      ClusterLifetimes( *CList_, clusterDataSetSize );

    // Change cluster num v time to grace-compatible colors if specified.
    if (grace_color_) {
      DataSet_integer& cnum_temp = static_cast<DataSet_integer&>( *cnumvtime_ );
      for (DataSet_integer::iterator ival = cnum_temp.begin();
                                     ival != cnum_temp.end(); ++ival)
      {
        *ival += 1;
        if (*ival > 15) *ival = 15;
      }
    }
    // Coordinate output.
    if (has_coords) {
      // Write clusters to trajectories
      if (!clusterfile_.empty())
        WriteClusterTraj( *CList_ ); 
      // Write all representative frames to a single traj
      if (!singlerepfile_.empty())
        WriteSingleRepTraj( *CList_ );
      // Write all representative frames to separate trajs
      if (!reptrajfile_.empty())
        WriteRepTraj( *CList_ );
      if (!avgfile_.empty())
        WriteAvgStruct( *CList_ );
    }
  } else
    mprintf("\tNo clusters found.\n");
  cluster_post.Stop();
  cluster_total.Stop();
  // Timing data
  mprintf("\tCluster timing data:\n");
  cluster_setup.WriteTiming(1,    "  Cluster Init. :", cluster_total.Total());
  cluster_pairwise.WriteTiming(1, "  Pairwise Calc.:", cluster_total.Total());
  cluster_cluster.WriteTiming(1,  "  Clustering    :", cluster_total.Total());
  cluster_post.WriteTiming(1,     "  Cluster Post. :", cluster_total.Total());
  cluster_total.WriteTiming(1,    "Total:");
  return Analysis::OK;
}

// -----------------------------------------------------------------------------
// Analysis_Clustering::CreateCnumvtime()
/** Put cluster number vs frame into dataset.  */
void Analysis_Clustering::CreateCnumvtime( ClusterList const& CList, int maxFrames ) {
  // FIXME:
  // Cast generic DataSet for cnumvtime back to integer dataset to 
  // access specific integer dataset functions for resizing and []
  // operator. Should this eventually be generic to all atomic DataSets? 
  DataSet_integer& cnum_temp = static_cast<DataSet_integer&>( *cnumvtime_ );
  cnum_temp.Resize( maxFrames );
  // Make all clusters start at -1. This way cluster algorithms that
  // have noise points (i.e. no cluster assigned) will be distinguished.
  std::fill(cnum_temp.begin(), cnum_temp.end(), -1);

  for (ClusterList::cluster_iterator C = CList.begincluster();
                                     C != CList.endcluster(); C++)
  {
    //mprinterr("Cluster %i:\n",CList->CurrentNum());
    int cnum = (*C).Num();
    // Loop over all frames in the cluster
    for (ClusterNode::frame_iterator frame = (*C).beginframe();
                                     frame != (*C).endframe(); frame++)
    {
      //mprinterr("%i,",*frame);
      cnum_temp[ *frame ] = cnum;
    }
    //mprinterr("\n");
    //break;
  }
}

// Analysis_Clustering::CreateCpopvtime()
// NOTE: Should not be called if cpopvtimefile is NULL
void Analysis_Clustering::CreateCpopvtime( ClusterList const& CList, int maxFrames ) {
  std::vector<int> Pop(CList.Nclusters(), 0);
  // Set up output data sets
  std::vector<DataSet*> DSL;
  for (int cnum = 0; cnum < CList.Nclusters(); ++cnum) { 
    DSL.push_back(masterDSL_->AddSetIdxAspect( DataSet::FLOAT, cnumvtime_->Name(), 
                                               cnum, "Pop" ));
    if (DSL.back() == 0) {
      mprinterr("Error: Could not allocate cluster pop v time DataSet.\n");
      return;
    }
    cpopvtimefile_->AddSet( DSL.back() );
  }
  // Set up normalization
  std::vector<double> Norm;
  if (norm_pop_ == CLUSTERPOP) {
    int cnum = 0;
    Norm.resize(CList.Nclusters(), 1.0);
    for (ClusterList::cluster_iterator C = CList.begincluster(); 
                                       C != CList.endcluster(); ++C)
      Norm[cnum++] = (double)((*C).Nframes());
  }
  // Assumes cnumvtime has been calcd and not gracecolor!
  double norm = 1.0;
  DataSet_integer const& cnum_temp = static_cast<DataSet_integer const&>( *cnumvtime_ );
  for (int frame = 0; frame < maxFrames; ++frame) {
    int cluster_num = cnum_temp[frame];
    // Noise points are -1
    if (cluster_num > -1)
      Pop[cluster_num]++;
    for (int cnum = 0; cnum < CList.Nclusters(); ++cnum) {
      // Normalization
      if (norm_pop_ == CLUSTERPOP)
        norm = Norm[cnum];
      else if (norm_pop_ == FRAME)
        norm = (double)(frame + 1);
      //float f = ((double)Pop[cnum] * Norm[cnum]);
      float f = (float)((double)Pop[cnum] / norm);
      DSL[cnum]->Add(frame, &f);
    }
  }
}

// Analysis_Clustering::ClusterLifetimes()
void Analysis_Clustering::ClusterLifetimes( ClusterList const& CList, int maxFrames ) {
  // Set up output data sets. TODO: use ChildDSL
  std::vector<DataSet_integer*> DSL;
  for (int cnum = 0; cnum < CList.Nclusters(); ++cnum) { 
    DSL.push_back((DataSet_integer*)
                  masterDSL_->AddSetIdxAspect( DataSet::INTEGER, cnumvtime_->Name(), 
                                               cnum, "Lifetime" ));
    if (DSL.back() == 0) {
      mprinterr("Error: Could not allocate cluster lifetime DataSet.\n");
      return;
    }
    DSL.back()->Resize( maxFrames );
  }
  // For each frame, assign cluster frame belongs to 1.
  DataSet_integer const& cnum_temp = static_cast<DataSet_integer const&>(*cnumvtime_);
  for (int frame = 0; frame < maxFrames; ++frame) {
    int cluster_num = cnum_temp[frame];
    // Noise points are -1
    if (cluster_num > -1)
      (*DSL[ cluster_num ])[ frame ] = 1;
  }
}

/** Determine how many different clusters are observed within a given time
  * window.
  */
void Analysis_Clustering::NclustersObserved( ClusterList const& CList, int maxFrames ) {
  DataSet_integer const& CVT = static_cast<DataSet_integer const&>( *cnumvtime_ );
  if (CVT.Size() < 1 || CList.Nclusters() < 1) return;
  int dataIdx = 0;
  // True if cluster was observed during window
  std::vector<bool> observed( CList.Nclusters(), false );
  for (int frame = 0; frame < maxFrames; frame++) {
    if (CVT[frame] != -1)
      observed[ CVT[frame] ] = true;
    if ( ((frame+1) % windowSize_) == 0 ) {
      // Count # observed clusters
      int nClustersObserved = 0;
      for (std::vector<bool>::iterator ob = observed.begin(); ob != observed.end(); ++ob)
        if ( *ob ) {
          ++nClustersObserved;
          *ob = false;
        }
      mprintf("DEBUG: WINDOW at frame %i; %i clusters observed\n", frame+1, nClustersObserved);
      clustersVtime_->Add( dataIdx++, &nClustersObserved );
    }
  }
/*
  int currentCluster = CVT[0];
  int nClustersObserved = 1;
  for (int frame = 1; frame < maxFrames; frame++) {
    // Do not count noise as a cluster.
    if (CVT[frame] != currentCluster && CVT[frame] != -1) {
      ++nClustersObserved;
      currentCluster = CVT[frame];
    }
    mprintf("DEBUG: %i %i\n", frame+1, nClustersObserved);
    if ( ((frame+1) % windowSize_) == 0 ) {
      mprintf("DEBUG: WINDOW\n");
      clustersVtime_->Add( dataIdx++, &nClustersObserved );
      nClustersObserved = 1;
    }
  }
*/
  clustersVtime_->SetDim(Dimension::X, Dimension(windowSize_, windowSize_, dataIdx));
} 

// ---------- Cluster Coordinate Output Routines -------------------------------
// Analysis_Clustering::WriteClusterTraj()
/** Write frames in each cluster to a trajectory file.  */
void Analysis_Clustering::WriteClusterTraj( ClusterList const& CList ) {
  Topology* clusterparm = (Topology*)&(coords_->Top()); // TODO: fix cast
  // Loop over all clusters
  for (ClusterList::cluster_iterator C = CList.begincluster();
                                     C != CList.endcluster(); ++C)
  {
    // Create filename based on cluster number.
    int cnum = C->Num();
    std::string cfilename =  clusterfile_ + ".c" + integerToString( cnum );
    // Set up trajectory file 
    Trajout clusterout;
    if (clusterout.InitTrajWrite(cfilename, clusterparm, clusterfmt_)) 
    {
      mprinterr("Error: Could not set up cluster trajectory %s for write.\n",
                cfilename.c_str());
      return;
    }
    // Loop over all frames in cluster
    int set = 0;
    Frame clusterframe = coords_->AllocateFrame();
    for (ClusterNode::frame_iterator fnum = C->beginframe();
                                     fnum != C->endframe(); ++fnum)
    {
      coords_->GetFrame( *fnum, clusterframe );
      clusterout.WriteFrame(set++, clusterparm, clusterframe);
    }
    // Close traj
    clusterout.EndTraj();
  }
}

// Analysis_Clustering::WriteAvgStruct()
void Analysis_Clustering::WriteAvgStruct( ClusterList const& CList ) {
  Topology avgparm = coords_->Top();
  avgparm.SetNframes( 1 );
  // Get extension for representative frame format 
  std::string tmpExt = TrajectoryFile::GetExtensionForType(avgfmt_);
  // Loop over all clusters
  for (ClusterList::cluster_iterator C = CList.begincluster();
                                     C != CList.endcluster(); ++C)
  {
    // Create filename based on cluster number.
    int cnum = C->Num();
    std::string cfilename = avgfile_ + ".c" + integerToString( cnum ) + tmpExt;
    // Set up trajectory file
    Trajout clusterout;
    if (clusterout.InitTrajWrite(cfilename, &avgparm, avgfmt_))
    {
      mprinterr("Error: Could not set up cluster average file %s for write.\n",
                cfilename.c_str());
      return;
    }
    // Get rep frame for rms fitting.
    Frame repframe = coords_->AllocateFrame();
    coords_->GetFrame( C->BestRepFrame(), repframe );
    Vec3 reftrans = repframe.CenterOnOrigin(false);
    // Loop over all frames in cluster
    Frame clusterframe = coords_->AllocateFrame();
    Frame avgframe = clusterframe;
    avgframe.ZeroCoords();
    for (ClusterNode::frame_iterator fnum = C->beginframe();
                                     fnum != C->endframe(); ++fnum)
    {
      coords_->GetFrame( *fnum, clusterframe );
      clusterframe.RMSD_FitToRef( repframe, reftrans );
      avgframe += clusterframe;
    }
    avgframe.Divide( (double)C->Nframes() );
    clusterout.WriteFrame(0, &avgparm, avgframe);
    clusterout.EndTraj();
  }
}
 
// Analysis_Clustering::WriteSingleRepTraj()
/** Write representative frame of each cluster to a trajectory file.  */
void Analysis_Clustering::WriteSingleRepTraj( ClusterList const& CList ) {
  Trajout clusterout;
  // Set up trajectory file. Use parm from COORDS DataSet. 
  Topology *clusterparm = (Topology*)&(coords_->Top()); // TODO: fix cast
  if (clusterout.InitTrajWrite(singlerepfile_, clusterparm, singlerepfmt_)) 
  {
    mprinterr("Error: Could not set up single trajectory for represenatatives %s for write.\n",
                singlerepfile_.c_str());
     return;
  }
  // Set up frame to hold cluster rep coords. 
  Frame clusterframe = coords_->AllocateFrame();
  int framecounter = 0;
  // Write rep frames from all clusters.
  for (ClusterList::cluster_iterator cluster = CList.begincluster(); 
                                     cluster != CList.endcluster(); ++cluster) 
  {
   coords_->GetFrame( cluster->BestRepFrame(), clusterframe );
   clusterout.WriteFrame(framecounter++, clusterparm, clusterframe);
  }
  // Close traj
  clusterout.EndTraj();
}

// Analysis_Clustering::WriteRepTraj()
/** Write representative frame of each cluster to a separate trajectory file,
  * repfile.REPNUM.FMT
  */
void Analysis_Clustering::WriteRepTraj( ClusterList const& CList ) {
  // Get extension for representative frame format 
  std::string tmpExt = TrajectoryFile::GetExtensionForType(reptrajfmt_);
  // Use Topology from COORDS DataSet to set up input frame
  Topology* clusterparm = (Topology*)&(coords_->Top()); // TODO: Fix cast
  Frame clusterframe = coords_->AllocateFrame();
  // Loop over all clusters
  for (ClusterList::cluster_iterator C = CList.begincluster();
                                     C != CList.endcluster(); ++C)
  {
    Trajout clusterout;
    // Get best rep frame # 
    int framenum = C->BestRepFrame();
    // Create filename based on frame #
    std::string cfilename = reptrajfile_ + ".c" + integerToString(C->Num());
    if (writeRepFrameNum_) cfilename += ("." + integerToString(framenum+1));
    cfilename += tmpExt;
    // Set up trajectory file. 
    if (clusterout.InitTrajWrite(cfilename, clusterparm, reptrajfmt_)) 
    {
      mprinterr("Error: Could not set up representative trajectory file %s for write.\n",
                cfilename.c_str());
       return;
    }
    // Write cluster rep frame
    coords_->GetFrame( framenum, clusterframe );
    clusterout.WriteFrame(framenum, clusterparm, clusterframe);
    // Close traj
    clusterout.EndTraj();
  }
}
