#include <cfloat> // DBL_MAX
#include <cmath> // sqrt
#include <vector>
#include <algorithm> // sort
#include "ClusterList.h"
#include "CpptrajStdio.h"
#include "CpptrajFile.h"
#include "Constants.h"
#include "ProgressBar.h"
#ifdef _OPENMP
#  include "omp.h"
#endif
#include "PDBfile.h" // DEBUG

// XMGRACE colors
const char* ClusterList::XMGRACE_COLOR[] = {
  "white", "black", "red", "green", "blue", "yellow", "brown", "grey", "violet",
  "cyan", "magenta", "orange", "indigo", "maroon", "turquoise", "darkgreen"
};

static const char* MetricStringArray[] = {
  "RMSD", "DME", "Symmetry-corrected RMSD", "Data Set(s)"
};

const char* ClusterList::MetricString(DistMetricType dm) {
  return MetricStringArray[dm];
}

// CONSTRUCTOR
ClusterList::ClusterList() : debug_(0), Cdist_(0) {}

// DESTRUCTOR
ClusterList::~ClusterList() {
  if (Cdist_ != 0) delete Cdist_;
}

// ClusterList::SetDebug()
/** Set the debug level */
void ClusterList::SetDebug(int debugIn) {
  debug_ = debugIn;
  if (debug_>0) mprintf("ClusterList debug set to %i\n",debug_);
}

// ClusterList::Renumber()
/** Sort clusters by size and renumber starting from 0, where cluster 0
  * is the largest. Also determine best representative frame and calculate 
  * anything dependent on ClusterDistances since sorting destroys indexing 
  * into ClusterDistances.
  */
void ClusterList::Renumber(bool addSievedFrames) {
  // Before clusters are renumbered, calculate the average distance of 
  // this cluster to every other cluster.
  // Only do this if ClusterDistances has been set.
  if (ClusterDistances_.Nelements() > 0) {
    double numdist = (double) (clusters_.size() - 1);
    for (cluster_it node = clusters_.begin(); node != clusters_.end(); ++node)
    {
      double avgclusterdist = 0.0;
      for (cluster_it node2 = clusters_.begin();
                      node2 != clusters_.end(); node2++)
      {
        if (node == node2) continue;
        //mprintf("DBG:\t\t%i to %i %f\n",(*node).num, (*node2).num, 
        //        ClusterDistances.GetElement( (*node).num, (*node2).num ));
        avgclusterdist += ClusterDistances_.GetCdist( (*node).Num(), (*node2).Num() );
      }
      if (numdist > 0.0)
        avgclusterdist /= numdist;
      //mprintf("DBG:\tCluster %i avg dist = %f\n",(*node).num,avgclusterdist);
      (*node).SetAvgDist( avgclusterdist ); 
    }
  }
  // Update cluster centroids.
  bool centroid_error = false;
  for (cluster_it node = clusters_.begin(); node != clusters_.end(); ++node) {
    node->SortFrameList();
    // Ensure cluster centroid is up-to-date
    node->CalculateCentroid( Cdist_ );
    // Find best representative frame
    if (node->FindBestRepFrame( FrameDistances_ ) == -1) {
      mprinterr("Error: Could not determine represenative frame for cluster %i\n",
                node->Num());
      centroid_error = true;
    }
  }
  // Add back sieved frames
  if (addSievedFrames) {
    if (centroid_error)
      mprinterr("Error: 1 or more centroids not determined. Cannot add sieved frames.\n");
    else {
      mprintf("\tRestoring sieved frames.\n");
      AddSievedFrames();
    }
    // Re-sort cluster frame lists.
    for (cluster_it node = clusters_.begin(); node != clusters_.end(); ++node)
      node->SortFrameList();
  }
  // Sort clusters by population 
  clusters_.sort( );
  // Renumber clusters.
  int newNum = 0;
  for (cluster_it node = clusters_.begin(); node != clusters_.end(); ++node) 
    node->SetNum( newNum++ );
  // TODO: Clear ClusterDistances?
}

// ClusterList::Summary()
/** Print a summary of clusters.  */
void ClusterList::Summary(std::string const& summaryfile, int maxframesIn) {
  CpptrajFile outfile;
  double fmax = (double)maxframesIn;
  if (outfile.OpenWrite(summaryfile)) {
    mprinterr("Error: ClusterList::Summary: Could not set up file.\n");
    return;
  }

  outfile.Printf("%-8s %8s %8s %8s %8s %8s %8s\n","#Cluster","Frames","Frac",
                     "AvgDist","Stdev","Centroid","AvgCDist");
  for (cluster_it node = clusters_.begin();
                  node != clusters_.end(); node++)
  {
    // Since there may be a lot of frames do not calculate SD from the
    // mean (which requires either storing distances or two double loops), 
    // instead use SD = sqrt( (SUM[x^2] - ((SUM[x])^2)/N)/N )
    double internalAvg = 0.0;
    double internalSD = 0.0;
    unsigned int Nelements = 0;
    if ((*node).Nframes() > 1) {
      // Calculate average distance between all frames in this cluster
      ClusterNode::frame_iterator frame2_end = (*node).endframe();
      ClusterNode::frame_iterator frame1_end = frame2_end;
      --frame1_end;
      for (ClusterNode::frame_iterator frm1 = (*node).beginframe();
                                       frm1 != frame1_end; ++frm1)
      {
        // Since this can be called after sieved frames are added back in,
        // need to ensure distances were calcd for these frames.
        if (!FrameDistances_.IgnoringRow(*frm1)) {
          ClusterNode::frame_iterator frm2 = frm1;
          ++frm2;
          for (; frm2 != frame2_end; ++frm2) {
            if (!FrameDistances_.IgnoringRow(*frm2)) {
              double dist = FrameDistances_.GetFdist(*frm1, *frm2);
              internalAvg += dist;
              internalSD += (dist * dist);
              ++Nelements;
            }
          }
        }
      }
      if (Nelements > 0) {
        double norm = 1.0 / ((double)Nelements);
        internalAvg *= norm;
        internalSD *= norm;
        internalSD -= (internalAvg * internalAvg);
        if (internalSD > 0.0)
          internalSD = sqrt( internalSD );
        else
          internalSD = 0.0;
      }
    }
    // OUTPUT
    outfile.Printf("%8i %8i %8.3f %8.3f %8.3f %8i %8.3f\n",
                   node->Num(), node->Nframes(), (double)node->Nframes()/fmax, internalAvg, 
                   internalSD, node->BestRepFrame()+1, node->AvgDist() );
  } // END loop over clusters
  outfile.CloseFile();
}

// ClusterList::Summary_Part
/** Print a summary of clustering for specified portions of the overall traj. 
  */
void ClusterList::Summary_Part(std::string const& summaryfile, int maxframesIn,
                               std::vector<int> const& splitFrames)
{
  const char* nExt[] = {"st", "nd", "rd", "th"};
  if (splitFrames.empty()) return; // Sanity check.
  CpptrajFile outfile;
  double fmax = (double)maxframesIn;
  if (outfile.OpenWrite(summaryfile)) {
    mprinterr("Error: Could not open file '%s'.\n", summaryfile.c_str());
    return;
  }

  // Determine number of frames and traj offset for each part.
  outfile.Printf("# 1st");
  std::vector<double> partMax;
  partMax.reserve( splitFrames.size() + 1 );
  std::vector<int> trajOffset;
  trajOffset.reserve( splitFrames.size() + 1);
  trajOffset.push_back( 0 );
  int lastMax = 0;
  unsigned int eidx = 1;
  for (unsigned int sf = 0; sf < splitFrames.size(); sf++)
  {
    partMax.push_back( (double)(splitFrames[sf] - lastMax) );
    lastMax = splitFrames[sf];
    trajOffset.push_back( lastMax );
    outfile.Printf(" <= %i < %u%s", trajOffset.back(), sf+2, nExt[eidx]);
    if (eidx < 3) ++eidx;
  }
  partMax.push_back( (double)(maxframesIn - lastMax) );
  outfile.Printf("\n# ");
  // Print # of frames in each section
  eidx=0;
  for (std::vector<double>::const_iterator pm = partMax.begin(); pm != partMax.end(); ++pm) {
    if (pm != partMax.begin()) outfile.Printf("  ");
    outfile.Printf("%u%s= %.0f", pm - partMax.begin() + 1, nExt[eidx], *pm);
    if (eidx < 3) ++eidx;
  }
  outfile.Printf("\n");
  // DEBUG
  //mprintf("DEBUG: # Frames (offset):");
  //std::vector<int>::const_iterator of = trajOffset.begin();
  //for (std::vector<double>::const_iterator it = partMax.begin();
  //                                         it != partMax.end(); ++it, ++of)
  //  mprintf(" %.0f (%i)", *it, *of);
  //mprintf("\n");
  // Set up bins
  std::vector<int> numInPart(  splitFrames.size() + 1, 0 );
  std::vector<int> firstFrame( splitFrames.size() + 1, -1);

  // Header
  outfile.Printf("#%-7s %8s %8s %2s %10s", "Cluster", "Total", "Frac", "C#", "Color");
  eidx = 0;
  for (unsigned int pm = 1; pm <= partMax.size(); ++pm) {
    outfile.Printf(" %5s%u%2s", "NumIn", pm, nExt[eidx]);
    if (eidx < 3) ++eidx;
  }
  for (unsigned int pm = 1; pm <= partMax.size(); ++pm)
    outfile.Printf(" %7s%u", "Frac", pm);
  for (unsigned int pm = 1; pm <= partMax.size(); ++pm)
    outfile.Printf(" %7s%u", "First", pm);
  outfile.Printf("\n");
  // LOOP OVER CLUSTERS
  int color = 1; // xmgrace color, 1-15
  for (cluster_it node = clusters_.begin();
                  node != clusters_.end(); node++)
  {
    // Calculate size and fraction of total size of this cluster
    int numframes = (*node).Nframes();
    double frac = (double)numframes / fmax;
    std::fill( numInPart.begin(), numInPart.end(), 0 );
    std::fill( firstFrame.begin(), firstFrame.end(), -1 );
    // DEBUG
    //mprintf("\tCluster %i\n",(*node).num);
    // Count how many frames are in each part. 
    for (ClusterNode::frame_iterator frame1 = (*node).beginframe();
                                     frame1 != (*node).endframe();
                                     frame1++)
    {
      unsigned int bin = splitFrames.size();
      for (unsigned int sf = 0; sf < splitFrames.size(); ++sf) {
        if ( *frame1 < splitFrames[sf] ) {
          bin = sf;
          break;
        }
      }
      if (numInPart[ bin ] == 0)
        firstFrame[ bin ] = *frame1 - trajOffset[ bin ] + 1;
      ++numInPart[ bin ];
    }
    outfile.Printf("%-8i %8i %8.4f %2i %10s", (*node).Num(), numframes, frac,
                   color, XMGRACE_COLOR[color]);
    for (std::vector<int>::const_iterator np = numInPart.begin();
                                          np != numInPart.end(); ++np)
      outfile.Printf(" %8i", *np);
    for (unsigned int pm = 0; pm < partMax.size(); ++pm)
      outfile.Printf(" %8.4f", ((double)numInPart[pm]) / partMax[pm]);
    for (std::vector<int>::const_iterator ff = firstFrame.begin();
                                          ff != firstFrame.end(); ++ff)
      outfile.Printf(" %8i", *ff);
    outfile.Printf("\n");
    if (color<15) ++color;
  }
  outfile.CloseFile();
}

// ClusterList::PrintClustersToFile()
/** Print list of clusters in a style similar to ptraj; each cluster is
  * given a line maxframes characters long, with X for each frame that is
  * in the clusters and . for all other frames. Also print out the
  * representative frame numbers.
  */
void ClusterList::PrintClustersToFile(std::string const& filename, int maxframesIn) {
  CpptrajFile outfile;
  std::string buffer;
  
  if ( outfile.OpenWrite(filename) ) {
    mprinterr("Error: PrintClustersToFile: Could not set up file %s\n",
              filename.c_str());
    return;
  }
  outfile.Printf("#Clustering: %u clusters %i frames\n",
                 clusters_.size(), maxframesIn);
  ComputeDBI( outfile );
  ComputePseudoF( outfile );
  // Call internal info routine.
  ClusterResults( outfile );
  // Do not print trajectory stuff if no filename given (i.e. STDOUT output)
  if (!filename.empty()) {
    for (cluster_it C1_it = clusters_.begin(); 
                    C1_it != clusters_.end(); C1_it++)
    {
      buffer.clear();
      buffer.resize(maxframesIn, '.');
      for (ClusterNode::frame_iterator frame1 = (*C1_it).beginframe();
                                       frame1 != (*C1_it).endframe();
                                       frame1++)
      {
        buffer[ *frame1 ] = 'X';
      }
      buffer += '\n';
      outfile.Write((void*)buffer.c_str(), buffer.size());
    }
  }
  // Print representative frame numbers
  outfile.Printf("#Representative frames:");
  for (cluster_it C = clusters_.begin(); C != clusters_.end(); C++)
    outfile.Printf(" %i", C->BestRepFrame()+1);
  outfile.Printf("\n");
  // Print sieve info if present
  if (FrameDistances_.SieveValue() != 1) {
    if (FrameDistances_.SieveValue() < -1) {
      outfile.Printf("#Sieve value: %i (random)\n#Sieved frames:", -FrameDistances_.SieveValue());
      ClusterSieve::SievedFrames sFrames = FrameDistances_.Sieved();
      for (ClusterSieve::SievedFrames::const_iterator sfrm = sFrames.begin();
                                                      sfrm != sFrames.end(); ++sfrm)
        outfile.Printf(" %i", *sfrm + 1);
      outfile.Printf("\n");
    } else
      outfile.Printf("#Sieve value: %i\n", FrameDistances_.SieveValue());
  }
  outfile.CloseFile();
}

// ClusterList::PrintClusters()
/** Print list of clusters and frame numbers belonging to each cluster.
  */
void ClusterList::PrintClusters() {
  mprintf("CLUSTER: %u clusters, %u frames.\n", clusters_.size(), FrameDistances_.Nframes() );
  for (cluster_it C = clusters_.begin(); C != clusters_.end(); C++) {
    mprintf("\t%8i : ",(*C).Num());
    for (ClusterNode::frame_iterator fnum = (*C).beginframe();
                                     fnum != (*C).endframe(); ++fnum)
      mprintf("%i,",(*fnum)+1);
    mprintf("\n");
  }
}

// -----------------------------------------------------------------------------
// ClusterList::AddCluster()
/** Add a cluster made up of frames specified by the given list of frames to 
  * the cluster list. Cluster # is current cluster list size.
  */
int ClusterList::AddCluster( ClusterDist::Cframes const& framelistIn ) {
  clusters_.push_back( ClusterNode( Cdist_, framelistIn, clusters_.size() ) );
  return 0;
}

// ClusterList::CalcFrameDistances()
int ClusterList::CalcFrameDistances(std::string const& filename, 
                                    ClusterDist::DsArray const& dataSets,
                                    DistModeType mode, DistMetricType metric, bool nofit, 
                                    bool useMass, std::string const& maskexpr,
                                    int sieve, int sieveSeed) 
{
  if (dataSets.empty()) {
    mprinterr("Internal Error: CalcFrameDistances: No DataSets given.\n");
    return 1;
  }
  // Base everything off of the first DataSet
  // TODO: Check all DataSet sizes?
  DataSet* dsIn = dataSets[0];
  // Set up internal cluster disance calculation
  if (metric != DATA) {
    if (!dsIn->IsCoordSet()) {
      mprinterr("Internal Error: Metric is COORDS base but data set is not.\n");
      return 1;
    }
    // Test that the mask expression is valid
    AtomMask testMask( maskexpr );
    Topology const& dsTop = ((DataSet_Coords*)dsIn)->Top();
    if ( dsTop.SetupIntegerMask( testMask ) ) {
      mprinterr("Error: Could not set up mask '%s' for topology %s\n",
                 maskexpr.c_str(), dsTop.c_str());
      return 1;
    }
    testMask.MaskInfo();
    if (testMask.None()) {
      mprinterr("Error: No atoms elected for mask '%s'\n", testMask.MaskString());
      return 1;
    }
    switch (metric) {
      case DME:   Cdist_ = new ClusterDist_DME(dsIn, testMask); break;
      case RMS:   Cdist_ = new ClusterDist_RMS(dsIn, testMask, nofit, useMass); break;
      case SRMSD: Cdist_ = new ClusterDist_SRMSD(dsIn, testMask, nofit, useMass, debug_); break;
      default: return 1; // Sanity check
    }
  } else { // Metric is DATA
    if (dataSets.size() == 1)
      Cdist_ = new ClusterDist_Num(dsIn);
    else // TODO: More than just euclid
      Cdist_ = new ClusterDist_Euclid(dataSets);
  }
  // Attempt to load pairwise distances from file if specified
  if (mode == USE_FILE && !filename.empty()) {
    mprintf("\tLoading pair-wise distances from %s\n", filename.c_str());
    if (FrameDistances_.LoadFile( filename, dsIn->Size() )) {
      mprintf("\tLoading pair-wise distances failed - regenerating from frames.\n");
      mode = USE_FRAMES;
    }
  }
  // Calculate pairwise distances from input DataSet. The ignore array will
  // be set up to ignore sieved frames.
  if (mode == USE_FRAMES) {
    mprintf("\tCalculating pair-wise distances.\n");
    // Set up ClusterMatrix with sieve.
    if (FrameDistances_.SetupWithSieve( dsIn->Size(), sieve, sieveSeed )) {
      mprinterr("Error: Could not setup matrix for pair-wise distances.\n");
      return 1; 
    }
    Cdist_->PairwiseDist(FrameDistances_, FrameDistances_.Sieved() );
  }
  // Save distances if filename specified and not previously loaded.
  if (mode == USE_FRAMES && !filename.empty()) {
    mprintf("\tSaving pair-wise distances to %s\n", filename.c_str());
    FrameDistances_.SaveFile( filename );
  }
  mprintf("\tMemory used by pair-wise matrix: %.4f MB\n",
          (double)FrameDistances_.DataSize() / 1048576);
  // DEBUG - Print Frame distances
  if (debug_ > 1) {
    mprintf("INITIAL FRAME DISTANCES:\n");
    FrameDistances_.PrintElements();
  }
  
  return 0;
}  

// ClusterList::RemoveEmptyClusters()
void ClusterList::RemoveEmptyClusters() {
  cluster_it cnode = clusters_.begin();
  while (cnode != clusters_.end()) {
    if (cnode->Nframes() == 0)
      cnode = clusters_.erase( cnode );
    else
      ++cnode;
  }
}

/** Calculate the distances between each cluster based on centroids. */
void ClusterList::CalcClusterDistances() {
  if (clusters_.empty()) return;
  ClusterDistances_.SetupMatrix( clusters_.size() );
  // Make sure centroid for clusters are up to date
  for (cluster_it C1 = clusters_.begin(); C1 != clusters_.end(); ++C1)
    C1->CalculateCentroid( Cdist_ );
  // Calculate distances between each cluster centroid
  cluster_it Cend = clusters_.end();
  for (cluster_it C1 = clusters_.begin(); C1 != Cend; ++C1) {
    cluster_it C2 = C1;
    ++C2;
    for (; C2 != Cend; ++C2)
      ClusterDistances_.AddElement( Cdist_->CentroidDist( C1->Cent(), C2->Cent() ) );
  }
}

// -----------------------------------------------------------------------------
void ClusterList::AddSievedFramesByCentroid() {
    // NOTE: All cluster centroids must be up to date.
  int frame;
  int nframes = (int)FrameDistances_.Nframes();
  double mindist, dist;
  cluster_it minNode, Cnode;
  ParallelProgress progress( nframes );
# ifdef _OPENMP
  int numthreads, mythread;
  // Need to create a ClusterDist for every thread to ensure memory allocation and avoid clashes
  ClusterDist** cdist_thread;
  // Also need a temp. array to hold which frame goes to which cluster to avoid clashes
  std::vector<cluster_it> frameToCluster( nframes, clusters_.end() );
# pragma omp parallel
  {
    if (omp_get_thread_num()==0)
      numthreads = omp_get_num_threads();
  }
  mprintf("\tParallelizing calculation with %i threads\n", numthreads);
  cdist_thread = new ClusterDist*[ numthreads ];
  for (int i=0; i < numthreads; i++)
    cdist_thread[i] = Cdist_->Copy();
# pragma omp parallel private(mythread, frame, dist, mindist, minNode, Cnode) firstprivate(progress)
{
  mythread = omp_get_thread_num();
  progress.SetThread( mythread );
# pragma omp for schedule(dynamic)
# endif
  for (frame = 0; frame < nframes; ++frame) {
    progress.Update( frame );
    if (FrameDistances_.IgnoringRow(frame)) {
      // Which clusters centroid is closest to this frame?
      mindist = DBL_MAX;
      minNode = clusters_.end();
      for (Cnode = clusters_.begin(); Cnode != clusters_.end(); ++Cnode) {
#       ifdef _OPENMP
        dist = cdist_thread[mythread]->FrameCentroidDist(frame, (*Cnode).Cent());
#       else
        dist = Cdist_->FrameCentroidDist(frame, (*Cnode).Cent());
#       endif
        if (dist < mindist) {
          mindist = dist;
          minNode = Cnode;
        }
      }
      // Add sieved frame to the closest cluster.
#     ifdef _OPENMP
      frameToCluster[frame] = minNode;
#     else
      (*minNode).AddFrameToCluster( frame );
#     endif
    }
  } // END loop over frames
# ifdef _OPENMP
} // END pragma omp parallel
  // Free cdist_thread memory
  for (int i = 0; i < numthreads; i++)
    delete cdist_thread[i];
  delete[] cdist_thread;
  // Now actually add sieved frames to their appropriate clusters
  for (frame = 0; frame < nframes; frame++)
    if (frameToCluster[frame] != clusters_.end())
      (*frameToCluster[frame]).AddFrameToCluster( frame );
# endif
  progress.Finish();
}

// -----------------------------------------------------------------------------
/** The Davies-Bouldin Index (DBI) is a measure of clustering merit; the 
  * smaller the DBI, the better. The DBI is defined as the average, for all 
  * clusters X, of fred, where fred(X) = max, across other clusters Y, of 
  * (Cx + Cy)/dXY ...here Cx is the average distance from points in X to the 
  * centroid, similarly Cy, and dXY is the distance between cluster centroids.
  */
double ClusterList::ComputeDBI(CpptrajFile& outfile) {
  std::vector<double> averageDist;
  averageDist.reserve( clusters_.size() );
  for (cluster_it C1 = clusters_.begin(); C1 != clusters_.end(); ++C1) {
    // Make sure centroid for this cluster is up to date
    C1->CalculateCentroid( Cdist_ );
    // Calculate average distance to centroid for this cluster
    averageDist.push_back( C1->CalcAvgToCentroid( Cdist_ ) );
    if (outfile.IsOpen())
      outfile.Printf("#Cluster %i has average-distance-to-centroid %f\n", 
                     C1->Num(), averageDist.back());
  }
  double DBITotal = 0.0;
  unsigned int nc1 = 0;
  for (cluster_it c1 = clusters_.begin(); c1 != clusters_.end(); ++c1, ++nc1) {
    double MaxFred = 0;
    unsigned int nc2 = 0;
    for (cluster_it c2 = clusters_.begin(); c2 != clusters_.end(); ++c2, ++nc2) {
      if (c1 == c2) continue;
      double Fred = averageDist[nc1] + averageDist[nc2];
      Fred /= Cdist_->CentroidDist( c1->Cent(), c2->Cent() );
      if (Fred > MaxFred)
        MaxFred = Fred;
    }
    DBITotal += MaxFred;
  }
  DBITotal /= (double)clusters_.size();
  if (outfile.IsOpen()) outfile.Printf("#DBI: %f\n", DBITotal);
  return DBITotal;
}

/** The pseudo-F statistic is another measure of clustering goodness. High 
  * values are good. Generally, one selects a cluster-count that gives a peak 
  * in the pseudo-f statistic (or pSF, for short).
  * Formula: A/B, where A = (T - P)/(G-1), and B = P / (n-G). Here n is the 
  * number of points, G is the number of clusters, T is the total distance from
  * the all-data centroid, and P is the sum (for all clusters) of the distances
  * from the cluster centroid.
  */
// NOTE: This calc differs slightly from PTRAJ in that real centroids are used
//       instead of representative structures.
double ClusterList::ComputePseudoF(CpptrajFile& outfile) {
  // Calculation makes no sense with fewer than 2 clusters.
  if (Nclusters() < 2) {
    mprintf("Warning: Fewer than 2 clusters. Not calculating pseudo-F.\n");
    return 0.0;
  }

  // Form a cluster with all points to get a centroid. Use only frames that
  // are in clusters, i.e. ignore noise. Also make sure all cluster centroids
  // are up to date.
  ClusterNode c_all;
  for (cluster_it C1 = clusters_.begin(); C1 != clusters_.end(); ++C1)
  {
    C1->CalculateCentroid( Cdist_ );
    for (ClusterNode::frame_iterator f1 = C1->beginframe(); f1 != C1->endframe(); ++f1)
      c_all.AddFrameToCluster( *f1 );
  }
  // Pseudo-F makes no sense if # clusters == # frames
  if (Nclusters() == c_all.Nframes()) {
    mprintf("Warning: Each frame is in a separate cluster. Not calculating pseudo-F.\n");
    return 0.0;
  }
  c_all.SortFrameList();
  c_all.CalculateCentroid( Cdist_ );

  // Loop over all clusters
  double gss = 0.0; // between-group sum of squares
  double wss = 0.0; // within-group sum of squares
  for (cluster_it C1 = clusters_.begin(); C1 != clusters_.end(); ++C1)
  {
    for (ClusterNode::frame_iterator f1 = C1->beginframe(); f1 != C1->endframe(); ++f1)
    {
      double dist = Cdist_->FrameCentroidDist(*f1, c_all.Cent());
      gss += (dist * dist);
      dist = Cdist_->FrameCentroidDist(*f1, C1->Cent());
      wss += (dist * dist);
    }
  }
  double d_nclusters = (double)Nclusters();
  double d_ntotal = (double)c_all.Nframes();
  double num = (gss - wss) / (d_nclusters - 1.0);
  double den = wss / (d_ntotal - d_nclusters);
  if (den < Constants::SMALL)
    den = Constants::SMALL;
  double pseudof = num / den;
  if (debug_ > 0)
    mprintf("Pseudo-f: Total distance to centroid is %.4f\n"
            "Pseudo-f: Cluster distance to centroid is %.4f\n"
            "Pseudo-f: Numerator %.4f over denominator %.4f gives %.4f\n", 
            gss, wss, num, den, pseudof);
  if (outfile.IsOpen()) outfile.Printf("#pSF: %f\n", pseudof);

  return pseudof;
}

/** The cluster silhouette is a measure of how well each point fits within
  * a cluster. Values of 1 indicate the point is very similar to other points
  * in the cluster, i.e. it is well-clustered. Values of -1 indicate the point
  * is dissimilar and may fit better in a neighboring cluster. Values of 0
  * indicate the point is on a border between two clusters. 
  */
void ClusterList::CalcSilhouette(std::string const& prefix) const {
  mprintf("\tCalculating cluster/frame silhouette.\n");
  CpptrajFile Ffile, Cfile;
  if (Ffile.OpenWrite(prefix + ".frame.dat")) return;
  if (Cfile.OpenWrite(prefix + ".cluster.dat")) return;
  Cfile.Printf("%-8s %10s\n", "#Cluster", "<Si>");
  unsigned int idx = 0;
  for (cluster_iterator Ci = begincluster(); Ci != endcluster(); ++Ci)
  {
    Ffile.Printf("#C%-6i %10s\n", Ci->Num(), "Silhouette");
    double avg_si = 0.0;
    int ci_frames = 0;
    std::vector<double> SiVals;
    for (ClusterNode::frame_iterator f1 = Ci->beginframe(); f1 != Ci->endframe(); ++f1)
    {
      if (FrameDistances_.IgnoringRow(*f1)) continue;
      // Calculate the average dissimilarity of this frame with all other
      // points in this frames cluster.
      double ai = 0.0;
      int self_frames = 0;
      for (ClusterNode::frame_iterator f2 = Ci->beginframe(); f2 != Ci->endframe(); ++f2)
      {
        if (f1 != f2 && !FrameDistances_.IgnoringRow(*f2)) {
          ai += FrameDistances_.GetFdist(*f1, *f2);
          ++self_frames;
        }
      }
      if (self_frames > 0)
        ai /= (double)self_frames;
      //mprintf("\t\tFrame %i cluster %i ai = %g\n", *f1+1, Ci->Num(), ai);
      // Determine lowest average dissimilarity of this frame with all
      // other clusters.
      double min_bi = DBL_MAX;
      for (cluster_iterator Cj = begincluster(); Cj != endcluster(); ++Cj)
      {
        if (Ci != Cj)
        {
          double bi = 0.0;
          int cj_frames = 0;
          // NOTE: ASSUMING NO EMPTY CLUSTERS
          for (ClusterNode::frame_iterator f2 = Cj->beginframe(); f2 != Cj->endframe(); ++f2)
          {
            if (!FrameDistances_.IgnoringRow(*f2)) {
              bi += FrameDistances_.GetFdist(*f1, *f2);
              ++cj_frames;
            }
          }
          bi /= (double)cj_frames;
          //mprintf("\t\tFrame %i to cluster %i bi = %g\n", *f1 + 1, Cj->Num(), bi);
          if (bi < min_bi)
            min_bi = bi;
        }
      }
      double max_ai_bi = std::max( ai, min_bi );
      if (max_ai_bi == 0.0)
        mprinterr("Error: Divide by zero in silhouette calculation for frame %i\n", *f1 + 1);
      else {
        double si = (min_bi - ai) / max_ai_bi;
        SiVals.push_back( si );
        //Ffile.Printf("%8i %10.4f\n", *f1 + 1, si);
        avg_si += si;
        ++ci_frames;
      }
    }
    std::sort( SiVals.begin(), SiVals.end() );
    for (std::vector<double>::const_iterator it = SiVals.begin(); it != SiVals.end(); ++it, ++idx)
      Ffile.Printf("%8i %g\n", idx, *it);
    Ffile.Printf("\n");
    ++idx;
    if (ci_frames > 0)
      avg_si /= (double)ci_frames;
    Cfile.Printf("%8i %g\n", Ci->Num(), avg_si);
  }
}

// -----------------------------------------------------------------------------
void ClusterList::DrawGraph(bool use_z, DataSet* cnumvtime,
                            double min_tol, int max_iteration) const
{
  if (use_z)
    mprintf("\tCreating PDB of graph points based on pairwise distances. B-factor = cluster #.\n");
  else
    mprintf("\tAttempting to draw graph based on pairwise distances.\n");
  unsigned int nframes = FrameDistances_.Nrows();
  std::vector<Vec3> Xarray; // Coords
  std::vector<Vec3> Farray; // Forces
  Xarray.reserve( nframes );
  Farray.assign( nframes, Vec3(0.0) );
  // Initialize coordinates. X and Y only.
  double zcoord = 0.0;
  double theta_deg = 0.0;
  double delta = 360.0 / (double)nframes;
  for (unsigned int n = 0; n != nframes; n++, theta_deg += delta) {
    double theta_rad = Constants::DEGRAD * theta_deg;
    if (use_z)
      zcoord = cos(theta_rad / 2.0);
    Xarray.push_back( Vec3(cos(theta_rad), sin(theta_rad), zcoord) );
  }
  // Write out initial graph
  if (debug_ > 0 && !use_z) {
    CpptrajFile graph0;
    if (graph0.OpenWrite("InitialGraph.dat")) return;
    for (std::vector<Vec3>::const_iterator XV = Xarray.begin();
                                           XV != Xarray.end(); ++XV)
      graph0.Printf("%g %g %u\n", (*XV)[0], (*XV)[1], XV - Xarray.begin() + 1);
    graph0.CloseFile();
  }
  // Degrees of freedom. If Z ever initialized needs to be 3N
  double deg_of_freedom = 2.0 * (double)nframes;
  if (use_z) deg_of_freedom += (double)nframes;
  double fnq = sqrt( deg_of_freedom );
  // Main loop for steepest descent
  const double Rk = 1.0;
  const double dxstm = 1.0E-5;
  const double crits = 1.0E-6;
  double rms = 1.0;
  double dxst = 0.1;
  double last_e = 0.0;
  int iteration = 0;
  mprintf("          \t%8s %12s %12s\n", " ", "ENE", "RMS");
  while (rms > min_tol && iteration < max_iteration) {
    double e_total = 0.0;
    ClusterMatrix::const_iterator Req = FrameDistances_.begin();
    for (unsigned int f1 = 0; f1 != nframes; f1++)
    {
      for (unsigned int f2 = f1 + 1; f2 != nframes; f2++)
      {
        //double Req = FrameDistances_.GetCdist(f1, f2);
        Vec3 V1_2 = Xarray[f1] - Xarray[f2];
        double r2 = V1_2.Magnitude2();
        double s = sqrt(r2);
        double r = 2.0 / s;
        double db = s - *(Req++);
        double df = Rk * db;
        double e = df * db;
        e_total += e;
        df *= r;
        // Apply force
        V1_2 *= df;
        Farray[f1] -= V1_2;
        Farray[f2] += V1_2;
      }
    }
    // Calculate the magnitude of the force vector.
    double sum = 0.0;
    for (std::vector<Vec3>::const_iterator FV = Farray.begin(); FV != Farray.end(); ++FV)
      sum += FV->Magnitude2();
    rms = sqrt( sum ) / fnq;
    // Adjust search step size
    if (dxst < crits) dxst = dxstm;
    dxst = dxst / 2.0;
    if (e_total < last_e) dxst = dxst * 2.4;
    double dxsth = dxst / sqrt( sum );
    last_e = e_total;
    // Update positions and reset force array.
    std::vector<Vec3>::iterator FV = Farray.begin();
    for (std::vector<Vec3>::iterator XV = Xarray.begin();
                                     XV != Xarray.end(); ++XV, ++FV)
    {
      *XV += (*FV * dxsth);
      *FV = 0.0;
    }
    // Write out current E.
    mprintf("Iteration:\t%8i %12.4E %12.4E\n", iteration, e_total, rms);
    iteration++;
  }
  // RMS error 
  ClusterMatrix::const_iterator Req = FrameDistances_.begin();
  double sumdiff2 = 0.0;
  for (unsigned int f1 = 0; f1 != nframes; f1++)
  {
    for (unsigned int f2 = f1 + 1; f2 != nframes; f2++)
    {
      Vec3 V1_2 = Xarray[f1] - Xarray[f2];
      double r1_2 = sqrt( V1_2.Magnitude2() );
      double diff = r1_2 - *Req;
      sumdiff2 += (diff * diff);
      if (debug_ > 0)
        mprintf("\t\t%u to %u: D= %g  Eq= %g  Delta= %g\n",
                f1+1, f2+1, r1_2, *Req, fabs(diff));
      ++Req;
    }
  }
  double rms_err = sqrt( sumdiff2 / (double)FrameDistances_.Nelements() );
  mprintf("\tRMS error of final graph positions: %g\n", rms_err);
  // Write out final graph with cluster numbers.
  std::vector<int> Nums;
  Nums.reserve( nframes );
  if (cnumvtime != 0) {
    ClusterSieve::SievedFrames sievedFrames = FrameDistances_.Sieved();
    DataSet_1D const& CVT = static_cast<DataSet_1D const&>( *cnumvtime );
    for (unsigned int n = 0; n != nframes; n++)
      Nums.push_back( (int)CVT.Dval(sievedFrames[n]) );
  } else
    for (int n = 1; n <= (int)nframes; n++)
      Nums.push_back( n );
  if (!use_z) {
    CpptrajFile graph;
    if (graph.OpenWrite("DrawGraph.dat")) return;
    for (std::vector<Vec3>::const_iterator XV = Xarray.begin();
                                           XV != Xarray.end(); ++XV)
      graph.Printf("%g %g %i \"%u\"\n", (*XV)[0], (*XV)[1], 
                   Nums[XV - Xarray.begin()], XV - Xarray.begin() + 1);
    graph.CloseFile();
  } else {
    // Write out PDB with B-factors
    PDBfile pdbout;
    if (pdbout.OpenWrite("DrawGraph.pdb")) return;
    pdbout.WriteTITLE("Cluster points.");
    for (std::vector<Vec3>::const_iterator XV = Xarray.begin();
                                           XV != Xarray.end(); ++XV)
      pdbout.WriteCoord(PDBfile::HETATM, XV - Xarray.begin() + 1, "HE", "HE", ' ',
                        XV - Xarray.begin() + 1, (*XV)[0], (*XV)[1], (*XV)[2],
                        1.0, Nums[XV - Xarray.begin()], "HE", 0, false);
    pdbout.CloseFile();
  }
}
