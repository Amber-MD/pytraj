#include <cfloat> // DBL_MAX
#include <algorithm> // sort, unique
#include "Cluster_DBSCAN.h"
#include "CpptrajStdio.h"
#include "ProgressBar.h"
#include "StringRoutines.h" // integerToString
#ifdef _OPENMP
#  include "omp.h"
#endif
#include "DataSet_MatrixDbl.h" // DEBUG
#include "DataFile.h" // DEBUG

Cluster_DBSCAN::Cluster_DBSCAN() :
  minPoints_(-1),
  epsilon_(-1.0),
  sieveToCentroid_(true)
{}

void Cluster_DBSCAN::Help() {
  mprintf("\t[dbscan minpoints <n> epsilon <e> [sievetoframe] [kdist <k> [kfile <prefix>]]]\n");
}

int Cluster_DBSCAN::SetupCluster(ArgList& analyzeArgs) {
  kdist_.SetRange(analyzeArgs.GetStringKey("kdist"));
  if (kdist_.Empty()) {
    minPoints_ = analyzeArgs.getKeyInt("minpoints", -1);
    if (minPoints_ < 1) {
      mprinterr("Error: DBSCAN requires minimum # of points to be set and >= 1\n"
                "Error: Use 'minpoints <N>'\n");
      return 1;
    }
    epsilon_ = analyzeArgs.getKeyDouble("epsilon", -1.0);
    if (epsilon_ <= 0.0) {
      mprinterr("Error: DBSCAN requires epsilon to be set and > 0.0\n"
                "Error: Use 'epsilon <e>'\n");
      return 1;
    }
    sieveToCentroid_ = !analyzeArgs.hasKey("sievetoframe");
  } else {
    k_prefix_ = analyzeArgs.GetStringKey("kfile");
    if (!k_prefix_.empty() && k_prefix_.at(k_prefix_.size()-1) != '/')
      k_prefix_ += '/';
  }
  return 0;
}

void Cluster_DBSCAN::ClusteringInfo() {
  mprintf("\tDBSCAN:\n");
  if (!kdist_.Empty()) {
    mprintf("\t\tOnly calculating Kdist graph for K=%s\n", kdist_.RangeArg());
    if (!k_prefix_.empty()) mprintf("\t\tKdist file prefix: %s\n", k_prefix_.c_str());
  } else {
    mprintf("\t\tMinimum pts to form cluster= %i\n", minPoints_);
    mprintf("\t\tCluster distance criterion= %.3f\n", epsilon_);
    if (sieveToCentroid_)
      mprintf("\t\tSieved frames will be added back solely based on their\n"
              "\t\t  closeness to cluster centroids.\n"
              "\t\t  (This option is less accurate but faster.)\n");
    else
      mprintf("\t\tSieved frames will only be added back if they are within\n"
              "\t\t  %.3f of a frame in an existing cluster.\n"
              "\t\t  (This option is more accurate and will identify sieved\n"
              "\t\t  frames as noise but is slower.)\n", epsilon_);
  }
}

// Cluster_DBSCAN::RegionQuery()
void Cluster_DBSCAN::RegionQuery(std::vector<int>& NeighborPts,
                                 std::vector<int> const& FramesToCluster,
                                 int point)
{
  NeighborPts.clear();
  for (std::vector<int>::const_iterator otherpoint = FramesToCluster.begin();
                                        otherpoint != FramesToCluster.end();
                                        ++otherpoint)
  {
    if (point == *otherpoint) continue;
    if ( FrameDistances_.GetFdist(point, *otherpoint) < epsilon_ )
      NeighborPts.push_back( *otherpoint );
  }
}

/** For each point p, calculate function Kdist(p) which is the distance of
  * the Kth nearest point to p.
  */
void Cluster_DBSCAN::ComputeKdist( int Kval, std::vector<int> const& FramesToCluster ) const {
  std::vector<double> dists;
  std::vector<double> Kdist;
  dists.reserve( FramesToCluster.size() ); 
  Kdist.reserve( FramesToCluster.size() );
  std::string outfilename = k_prefix_ + "Kdist." + integerToString(Kval) + ".dat";
  mprintf("\tDBSCAN: Calculating Kdist(%i), output to %s\n", Kval, outfilename.c_str());
  for (std::vector<int>::const_iterator point = FramesToCluster.begin();
                                        point != FramesToCluster.end();
                                        ++point)
  {
    // Store distances from this point
    dists.clear();
    for (std::vector<int>::const_iterator otherpoint = FramesToCluster.begin();
                                          otherpoint != FramesToCluster.end();
                                          ++otherpoint)
      dists.push_back( FrameDistances_.GetFdist(*point, *otherpoint) );
    // Sort distances - first dist should always be 0
    std::sort(dists.begin(), dists.end());
    Kdist.push_back( dists[Kval] );
  }
  std::sort( Kdist.begin(), Kdist.end() );
  CpptrajFile Outfile;
  Outfile.OpenWrite(outfilename);
  Outfile.Printf("%-8s %1i%-11s\n", "#Point", Kval,"-dist");
  // Write out largest to smallest
  unsigned int ik = 0;
  for (std::vector<double>::reverse_iterator k = Kdist.rbegin(); 
                                             k != Kdist.rend(); ++k, ++ik)
    Outfile.Printf("%8u %12.4f\n", ik, *k);
  Outfile.CloseFile();
}

// Cluster_DBSCAN::ComputeKdistMap()
void Cluster_DBSCAN::ComputeKdistMap( Range const& Kvals, 
                                      std::vector<int> const& FramesToCluster ) const
{
  int pt1_idx, pt2_idx, d_idx, point;
  mprintf("\tCalculating Kdist map for %s\n", Kvals.RangeArg());
  double* kdist_array; // Store distance of pt1 to every other point.
  int nframes = (int)FramesToCluster.size();
  // Ensure all Kdist points are within proper range
  Range::const_iterator kval;
  for (kval = Kvals.begin(); kval != Kvals.end(); ++kval)
    if (*kval < 1 || *kval >= nframes) {
      mprinterr("Error: Kdist value %i is out of range (1 <= Kdist < %i)\n",
                 *kval, nframes);
      return;
    }
  int nvals = (int)Kvals.Size();
  double** KMAP; // KMAP[i] has the ith nearest point for each point.
  KMAP = new double*[ nvals ];
  for (int i = 0; i != nvals; i++)
    KMAP[i] = new double[ nframes ];
  ParallelProgress progress( nframes );
# ifdef _OPENMP
# pragma omp parallel private(pt1_idx, pt2_idx, d_idx, kval, point, kdist_array) firstprivate(progress)
  {
  progress.SetThread( omp_get_thread_num() );
#endif
  kdist_array = new double[ nframes ];
# ifdef _OPENMP
# pragma omp for
# endif
  for (pt1_idx = 0; pt1_idx < nframes; pt1_idx++) // X
  {
    progress.Update( pt1_idx );
    point = FramesToCluster[pt1_idx];
    d_idx = 0;
    // Store distances from pt1 to pt2
    for (pt2_idx = 0; pt2_idx != nframes; pt2_idx++)
      kdist_array[d_idx++] = FrameDistances_.GetFdist(point, FramesToCluster[pt2_idx]);
    // Sort distances; will be smallest to largest
    std::sort( kdist_array, kdist_array + nframes );
    // Save the distance of specified nearest neighbors to this point.
    d_idx = 0;
    for (kval = Kvals.begin(); kval != Kvals.end(); ++kval) // Y
      KMAP[d_idx++][pt1_idx] = kdist_array[ *kval ];
  }
  delete[] kdist_array;
# ifdef _OPENMP
  } // END omp parallel
# endif
  progress.Finish();
  // Sort all of the individual kdist plots, smallest to largest.
  for (int i = 0; i != nvals; i++)
    std::sort(KMAP[i], KMAP[i] + nframes);
  // Save in matrix, largest to smallest.
  DataSet_MatrixDbl kmatrix;
  kmatrix.Allocate2D( FramesToCluster.size(), Kvals.Size() );
  for (int y = 0; y != nvals; y++) {
    for (int x = nframes - 1; x != -1; x--)
      kmatrix.AddElement( KMAP[y][x] );
    delete[] KMAP[y];
  }
  delete[] KMAP;
  // Write matrix to file
  DataFile outfile;
  ArgList outargs("usemap");
  outfile.SetupDatafile(k_prefix_ + "Kmatrix.gnu", outargs, debug_);
  outfile.AddSet( (DataSet*)&kmatrix );
  outfile.WriteData();
  // Write out the largest and smallest values for each K.
  // This means for each value of K the point with the furthest Kth-nearest
  // neighbor etc.
  CpptrajFile maxfile;
  if (maxfile.OpenWrite(k_prefix_ + "Kmatrix.max.dat")) return;
  maxfile.Printf("%-12s %12s %12s\n", "#Kval", "MaxD", "MinD");
  d_idx = 0;
  for (kval = Kvals.begin(); kval != Kvals.end(); ++kval, d_idx++)
    maxfile.Printf("%12i %12g %12g\n", *kval, kmatrix.GetElement(0, d_idx),
                   kmatrix.GetElement(nframes-1, d_idx));
  maxfile.CloseFile();
}

// Potential frame statuses.
char Cluster_DBSCAN::UNASSIGNED = 'U';
char Cluster_DBSCAN::NOISE = 'N';
char Cluster_DBSCAN::INCLUSTER = 'C';

/** Ester, Kriegel, Sander, Xu; Proceedings of 2nd International Conference
  * on Knowledge Discovery and Data Mining (KDD-96); pp 226-231.
  */
int Cluster_DBSCAN::Cluster() {
  std::vector<int> NeighborPts;
  std::vector<int> Npts2; // Will hold neighbors of a neighbor
  std::vector<int> FramesToCluster;
  ClusterDist::Cframes cluster_frames;
  // First determine which frames are being clustered.
  // FIXME: Just use sieved array?
  for (int frame = 0; frame < (int)FrameDistances_.Nframes(); ++frame)
    if (!FrameDistances_.IgnoringRow( frame ))
      FramesToCluster.push_back( frame );
  // Calculate Kdist function
  if (!kdist_.Empty()) {
    if (kdist_.Size() == 1)
      ComputeKdist( kdist_.Front(), FramesToCluster );
    else
      ComputeKdistMap( kdist_, FramesToCluster );
    return 0;
  }
  // Set up array to keep track of points that have been visited.
  // Make it the size of FrameDistances so we can index into it. May
  // waste memory during sieving but makes code easier.
  std::vector<bool> Visited( FrameDistances_.Nframes(), false );
  // Set up array to keep track of whether points are noise or in a cluster.
  Status_.assign( FrameDistances_.Nframes(), UNASSIGNED);
  mprintf("\tStarting DBSCAN Clustering:\n");
  ProgressBar cluster_progress(FramesToCluster.size());
  int iteration = 0;
  for (std::vector<int>::iterator point = FramesToCluster.begin();
                                  point != FramesToCluster.end(); ++point)
  {
    if (!Visited[*point]) {
      // Mark this point as visited
      Visited[*point] = true;
      // Determine how many other points are near this point
      RegionQuery( NeighborPts, FramesToCluster, *point );
      if (debug_ > 0) {
        mprintf("\tPoint %i\n", *point + 1);
        mprintf("\t\t%u neighbors:", NeighborPts.size());
      }
      // If # of neighbors less than cutoff, noise; otherwise cluster
      if ((int)NeighborPts.size() < minPoints_) {
        if (debug_ > 0) mprintf(" NOISE\n");
        Status_[*point] = NOISE;
      } else {
        // Expand cluster
        cluster_frames.clear();
        cluster_frames.push_back( *point );
        // NOTE: Use index instead of iterator since NeighborPts may be
        //       modified inside this loop.
        unsigned int endidx = NeighborPts.size();
        for (unsigned int idx = 0; idx < endidx; ++idx) {
          int neighbor_pt = NeighborPts[idx];
          if (!Visited[neighbor_pt]) {
            if (debug_ > 0) mprintf(" %i", neighbor_pt + 1);
            // Mark this neighbor as visited
            Visited[neighbor_pt] = true;
            // Determine how many other points are near this neighbor
            RegionQuery( Npts2, FramesToCluster, neighbor_pt );
            if ((int)Npts2.size() >= minPoints_) {
              // Add other points to current neighbor list
              NeighborPts.insert( NeighborPts.end(), Npts2.begin(), Npts2.end() );
              endidx = NeighborPts.size();
            }
          }
          // If neighbor is not yet part of a cluster, add it to this one.
          if (Status_[neighbor_pt] != INCLUSTER) {
            cluster_frames.push_back( neighbor_pt );
            Status_[neighbor_pt] = INCLUSTER;
          }
        }
        // Remove duplicate frames
        // TODO: Take care of this in Renumber?
        std::sort(cluster_frames.begin(), cluster_frames.end());
        ClusterDist::Cframes::iterator it = std::unique(cluster_frames.begin(), 
                                                        cluster_frames.end());
        cluster_frames.resize( std::distance(cluster_frames.begin(),it) );
        // Add cluster to the list
        AddCluster( cluster_frames );
        if (debug_ > 0) {
          mprintf("\n");
          PrintClusters();
        }
      }
    }
    cluster_progress.Update(iteration++);
  } // END loop over FramesToCluster
  // Calculate the distances between each cluster based on centroids
  CalcClusterDistances();

  return 0;
}

// Cluster_DBSCAN::ClusterResults()
void Cluster_DBSCAN::ClusterResults(CpptrajFile& outfile) const {
  outfile.Printf("#Algorithm: DBSCAN minpoints %i epsilon %g sieveToCentroid %i\n",
                 minPoints_, epsilon_, (int)sieveToCentroid_);
  // List the number of noise points.
  outfile.Printf("#NOISE_FRAMES:");
  unsigned int frame = 1;
  unsigned int numNoise = 0;
  for (std::vector<char>::const_iterator stat = Status_.begin();
                                         stat != Status_.end(); 
                                       ++stat, ++frame)
  {
    if ( *stat == NOISE ) {
      outfile.Printf(" %i", frame);
      ++numNoise;
    }
  }
  outfile.Printf("\n");
  outfile.Printf("#Number_of_noise_frames: %u\n", numNoise); 
}

// Cluster_DBSCAN::AddSievedFrames()
void Cluster_DBSCAN::AddSievedFrames() {
  int n_sieved_noise = 0;
  int Nsieved = 0;
  // NOTE: All cluster centroids must be up to date!
  if (sieveToCentroid_)
    mprintf("\tRestoring sieved frames by closeness to existing centroids.\n");
  else
    mprintf("\tRestoring sieved frames if within %.3f of frame in nearest cluster.\n",
            epsilon_);
  // Vars allocated here in case of OpenMP
  int frame, cidx;
  int nframes = (int)FrameDistances_.Nframes();
  double mindist, dist;
  cluster_it minNode, Cnode;
  bool goodFrame;
  ParallelProgress progress( nframes );
  // Need a temporary array to hold which frame belongs to which cluster. 
  // Otherwise we could be comparoing sieved frames to other sieved frames.
  std::vector<cluster_it> frameToCluster( nframes, clusters_.end() );
# ifdef _OPENMP
  int numthreads, mythread;
  // Need to create a ClusterDist for every thread to ensure memory allocation and avoid clashes
  ClusterDist** cdist_thread;
# pragma omp parallel
  {
    if (omp_get_thread_num()==0)
      numthreads = omp_get_num_threads();
  }
  mprintf("\tParallelizing calculation with %i threads\n", numthreads);
  cdist_thread = new ClusterDist*[ numthreads ];
  for (int i=0; i < numthreads; i++)
    cdist_thread[i] = Cdist_->Copy();
# pragma omp parallel private(mythread, frame, dist, mindist, minNode, Cnode, goodFrame, cidx) firstprivate(progress) reduction(+ : Nsieved, n_sieved_noise)
{
    mythread = omp_get_thread_num();
    progress.SetThread( mythread );
#   pragma omp for schedule(dynamic)
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
      goodFrame = false;
      if ( sieveToCentroid_ || mindist < epsilon_ )
        // Sieving based on centroid only or frame is already within epsilon, accept.
        goodFrame = true;
      else {
        // Check if any frames in the cluster are closer than epsilon to sieved frame.
        for (cidx=0; cidx < (*minNode).Nframes(); cidx++)
        {
#         ifdef _OPENMP
          if ( cdist_thread[mythread]->FrameDist(frame, (*minNode).ClusterFrame(cidx)) < epsilon_ )
#         else
          if ( Cdist_->FrameDist(frame, (*minNode).ClusterFrame(cidx)) < epsilon_ )
#         endif
          {
            goodFrame = true;
            break;
          }
        }
      }
      // Add sieved frame to the closest cluster if closest distance is
      // less than epsilon.
      ++Nsieved;
      if ( goodFrame )
        frameToCluster[frame] = minNode;
      else
        n_sieved_noise++;
    }
  } // END loop over frames
# ifdef _OPENMP
} // END pragma omp parallel
  // Free cdist_thread memory
  for (int i = 0; i < numthreads; i++)
    delete cdist_thread[i];
  delete[] cdist_thread;
# endif
  progress.Finish();
  // Now actually add sieved frames to their appropriate clusters
  for (frame = 0; frame < nframes; frame++)
    if (frameToCluster[frame] != clusters_.end())
      (*frameToCluster[frame]).AddFrameToCluster( frame );
  mprintf("\t%i of %i sieved frames were discarded as noise.\n", 
          n_sieved_noise, Nsieved);
}
