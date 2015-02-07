#include <cmath> // fabs
#include <algorithm> // sort
#include "Cluster_DPeaks.h"
#include "CpptrajStdio.h"
#include "DataSet_Mesh.h"
#include "ProgressBar.h"

Cluster_DPeaks::Cluster_DPeaks() : epsilon_(-1.0), calc_noise_(false) {}

void Cluster_DPeaks::Help() {
  mprintf("\t[dpeaks epsilon <e> [noise] [dvdfile <density_vs_dist_file>]\n"
          "\t  [runavg <runavg_file>] [deltafile <file>]]\n");
}

int Cluster_DPeaks::SetupCluster(ArgList& analyzeArgs) {
  epsilon_ = analyzeArgs.getKeyDouble("epsilon", -1.0);
  if (epsilon_ <= 0.0) {
    mprinterr("Error: DPeaks requires epsilon to be set and > 0.0\n"
              "Error: Use 'epsilon <e>'\n");
    return 1;
  }
  calc_noise_ = analyzeArgs.hasKey("noise");
  dpeaks_ = analyzeArgs.GetStringKey("dvdfile");
  rafile_ = analyzeArgs.GetStringKey("runavg");
  radelta_ = analyzeArgs.GetStringKey("deltafile");
  avg_factor_ = analyzeArgs.getKeyInt("avgfactor", -1);
  if (avg_factor_ != -1 && avg_factor_ < 1) {
    mprinterr("Error: avgfactor must be >= 1.\n");
    return 1;
  }
  return 0;
}

void Cluster_DPeaks::ClusteringInfo() {
  mprintf("\tDPeaks: Cutoff (epsilon) for determining local density is %g\n", epsilon_);
  if (calc_noise_)
    mprintf("\t\tCalculating noise as all points within epsilon of another cluster.\n");
  if (!dpeaks_.empty())
    mprintf("\t\tDensity vs min distance to point with next highest density written to %s\n",
            dpeaks_.c_str());
  if (!rafile_.empty())
    mprintf("\t\tRunning avg of delta vs distance written to %s\n", rafile_.c_str());
  if (avg_factor_ != -1)
    mprintf("\t\tRunning avg window size will be # clustered frames / %i\n", avg_factor_);
  if (!radelta_.empty())
    mprintf("\t\tDelta of distance minus running avg written to %s\n", radelta_.c_str());
}

#ifdef DISABLE
// -----------------------------------------------------------------------------
int Cluster_DPeaks::Cluster() {
  mprintf("\tStarting DPeaks clustering.\n");
  // First determine which frames are being clustered.
  Points_.clear();
  int oidx = 0;
  for (int frame = 0; frame < (int)FrameDistances_.Nframes(); ++frame)
    if (!FrameDistances_.IgnoringRow( frame ))
      Points_.push_back( Cpoint(frame, oidx++) );
  // Sanity check.
  if (Points_.size() < 2) {
    mprinterr("Error: Only 1 frame in initial clustering.\n");
    return 1;
  }

  // Sort distances
  std::vector<float> Distances;
  for (ClusterMatrix::const_iterator mat = FrameDistances_.begin();
                                     mat != FrameDistances_.end(); ++mat)
    Distances.push_back( *mat );
  std::sort( Distances.begin(), Distances.end() );
  unsigned int idx = (unsigned int)((double)Distances.size() * 0.02);
  double bandwidth = (double)Distances[idx];
  mprintf("idx= %u, bandwidth= %g\n", idx, bandwidth);

  // Density via Gaussian kernel
  double maxDist = -1.0;
  for (unsigned int i = 0; i != Points_.size(); i++) {
    for (unsigned int j = i+1; j != Points_.size(); j++) {
      double dist = FrameDistances_.GetFdist(Points_[i].Fnum(), Points_[j].Fnum());
      maxDist = std::max( maxDist, dist );
      dist /= bandwidth;
      double gk = exp(-(dist *dist));
      Points_[i].AddDensity( gk );
      Points_[j].AddDensity( gk );
    }
  }
  mprintf("Max dist= %g\n", maxDist);
  CpptrajFile rhoOut;
  rhoOut.OpenWrite("rho.dat");
  for (unsigned int i = 0; i != Points_.size(); i++)
    rhoOut.Printf("%u %g\n", i+1, Points_[i].Density());
  rhoOut.CloseFile();
  
  // Sort by density, descending
  std::stable_sort( Points_.begin(), Points_.end(), Cpoint::density_sort_descend() );
  CpptrajFile ordrhoOut;
  ordrhoOut.OpenWrite("ordrho.dat");
  for (unsigned int i = 0; i != Points_.size(); i++)
    ordrhoOut.Printf("%u %g %i %i\n", i+1, Points_[i].Density(), Points_[i].Fnum()+1,
                     Points_[i].Oidx()+1);
  ordrhoOut.CloseFile();

  // Determine minimum distances
  int first_idx = Points_[0].Oidx();
  Points_[first_idx].SetDist( -1.0 );
  Points_[first_idx].SetNearestIdx(-1);
  for (unsigned int ii = 1; ii != Points_.size(); ii++) {
    int ord_i = Points_[ii].Oidx();
    Points_[ord_i].SetDist( maxDist );
    for (unsigned int jj = 0; jj != ii; jj++) {
      int ord_j = Points_[jj].Oidx();
      double dist = FrameDistances_.GetFdist(Points_[ord_i].Fnum(), Points_[ord_j].Fnum());
      if (dist < Points_[ord_i].Dist()) {
        Points_[ord_i].SetDist( dist );
        Points_[ord_j].SetNearestIdx( ord_j );
      }
    }
  }
  if (!dpeaks_.empty()) {
    CpptrajFile output;
    if (output.OpenWrite(dpeaks_)) return 1;
    for (Carray::const_iterator point = Points_.begin(); point != Points_.end(); ++point)
      output.Printf("%g %g %i\n", point->Density(), point->Dist(), point->NearestIdx()+1);
    output.CloseFile();
  }
      
  return 1;
}
#endif

// -----------------------------------------------------------------------------
int Cluster_DPeaks::Cluster() {
  mprintf("\tStarting DPeaks clustering.\n");
  Points_.clear();
  // First determine which frames are being clustered.
  for (int frame = 0; frame < (int)FrameDistances_.Nframes(); ++frame)
    if (!FrameDistances_.IgnoringRow( frame ))
      Points_.push_back( Cpoint(frame) );
  // Sanity check.
  if (Points_.size() < 2) {
    mprinterr("Error: Only 1 frame in initial clustering.\n");
    return 1;
  }
  // For each point, determine how many others are within epsilon. Also
  // determine maximum distance between any two points.
  mprintf("\tDetermining local density of each point.\n");
  ProgressBar cluster_progress( Points_.size() );
  double maxDist = -1.0;
  for (Carray::iterator point0 = Points_.begin();
                        point0 != Points_.end(); ++point0)
  {
    cluster_progress.Update(point0 - Points_.begin());
    int density = 0;
    for (Carray::const_iterator point1 = Points_.begin();
                                point1 != Points_.end(); ++point1)
    {
      if (point0 != point1) {
        double dist = FrameDistances_.GetFdist(point0->Fnum(), point1->Fnum());
        maxDist = std::max(maxDist, dist);
        if ( dist < epsilon_ )
          density++;
      }
    }
    point0->SetDensity( density );
  }
  mprintf("DBG: Max dist= %g\n", maxDist);
  // DEBUG: Frame/Density
  CpptrajFile fdout;
  fdout.OpenWrite("fd.dat");
  for (Carray::const_iterator point = Points_.begin(); point != Points_.end(); ++point)
    fdout.Printf("%i %i\n", point->Fnum()+1, point->Density());
  fdout.CloseFile();
  // Sort by density here. Otherwise array indices will be invalid later.
  std::sort( Points_.begin(), Points_.end(), Cpoint::density_sort() );
  // For each point, find the closest point that has higher density. Since 
  // array is now sorted by density the last point has the highest density.
  Points_.back().SetDist( maxDist );
  mprintf("\tFinding closest neighbor point with higher density for each point.\n");
  unsigned int lastidx = Points_.size() - 1;
  cluster_progress.SetupProgress( lastidx );
  for (unsigned int idx0 = 0; idx0 != lastidx; idx0++)
  {
    cluster_progress.Update( idx0 );
    double min_dist = maxDist;
    int nearestIdx = -1; // Index of nearest neighbor with higher density
    Cpoint& point0 = Points_[idx0];
    //mprintf("\nDBG:\tSearching for nearest neighbor to idx %u with higher density than %i.\n",
    //        idx0, point0.Density());
    // Since array is sorted by density we can start at the next point.
    for (unsigned int idx1 = idx0+1; idx1 != Points_.size(); idx1++)
    {
      Cpoint const& point1 = Points_[idx1];
      double dist1_2 = FrameDistances_.GetFdist(point0.Fnum(), point1.Fnum());
      if (point1.Density() > point0.Density())
      {
        if (dist1_2 < min_dist) {
          min_dist = dist1_2;
          nearestIdx = (int)idx1;
          //mprintf("DBG:\t\tNeighbor idx %i is closer (density %i, distance %g)\n",
          //        nearestIdx, point1.Density(), min_dist);
        }
      }
    }
    point0.SetDist( min_dist );
    //mprintf("DBG:\tClosest point to %u with higher density is %i (distance %g)\n",
    //        idx0, nearestIdx, min_dist);
    point0.SetNearestIdx( nearestIdx );
  }
  // Plot density vs distance for each point.
  if (!dpeaks_.empty()) {
    CpptrajFile output;
    if (output.OpenWrite(dpeaks_))
      mprinterr("Error: Could not open density vs distance plot '%s' for write.\n",
                dpeaks_.c_str()); // TODO: Make fatal?
    else {
      output.Printf("%-10s %10s %s %10s %10s\n", "#Density", "Distance",
                    "Frame", "Idx", "Neighbor");
      for (Carray::const_iterator point = Points_.begin();
                                  point != Points_.end(); ++point)
        output.Printf("%-10i %10g \"%i\" %10u %10i\n", point->Density(), point->Dist(),
                      point->Fnum()+1, point-Points_.begin(), point->NearestIdx());
      output.CloseFile();
    }
  }
  // Choose points for which the min distance to point with higher density is
  // anomalously high.
  // Right now all density values are discrete. Try to choose outliers at each
  // value for which there is density.;
 
/*
  // For each point, calculate average distance (X,Y) to points in next and
  // previous density values.
  const double dens_cut = 3.0 * 3.0;
  const double dist_cut = 1.32 * 1.32;
  for (Carray::const_iterator point0 = Points_.begin(); point0 != Points_.end(); ++point0)
  {
    int Npts = 0;
    for (Carray::const_iterator point1 = Points_.begin(); point1 != Points_.end(); ++point1)
    {
      if (point0 != point1) {
        // Only do this for close densities
        double dX = (double)(point0->Density() - point1->Density());
        double dX2 = dX * dX;
        double dY = (point0->Dist() - point1->Dist());
        double dY2 = dY * dY;
        if (dX2 < dens_cut && dY2 < dist_cut) {
          Npts++;
        }
      }
    }
    mprintf("%i %i %i\n", point0->Density(), point0->Fnum()+1, Npts);
  }
*/

/*
  CpptrajFile tempOut;
  tempOut.OpenWrite("temp.dat");
  int currentDensity = -1;
  double distAv = 0.0;
  double distSD = 0.0;
  double sumWts = 0.0;
  int nValues = 0;
  Carray::const_iterator lastPoint = Points_.end() + 1;
  for (Carray::const_iterator point = Points_.begin(); point != lastPoint; ++point)
  {
    if (point == Points_.end() || point->Density() != currentDensity) {
      if (nValues > 0) {
        distAv = distAv / sumWts; //(double)nValues;
        distSD = (distSD / sumWts) - (distAv * distAv);
        if (distSD > 0.0)
          distSD = sqrt(distSD);
        else
          distSD = 0.0;
        //mprintf("Density %i: %i values  Avg= %g  SD= %g  SumWts= %g\n", currentDensity,
        //        nValues, distAv, distSD, sumWts);
        tempOut.Printf("%i %g\n", currentDensity, distAv);
      }
      if (point == Points_.end()) break;
      currentDensity = point->Density();
      distAv = 0.0;
      distSD = 0.0;
      sumWts = 0.0;
      nValues = 0;
    }
    double wt = exp(point->Dist());
    double dval = point->Dist() * wt;
    sumWts += wt;
    distAv += dval;
    distSD += (dval * dval);
    nValues++;
  }
  tempOut.CloseFile(); 
*/

  // BEGIN CALCULATING WEIGHTED DISTANCE AVERAGE
  CpptrajFile tempOut;
  tempOut.OpenWrite("temp.dat");
  DataSet_Mesh weightedAverage;
  Carray::const_iterator cp = Points_.begin();
  // Skip local density of 0.
  //while (cp->Density() == 0 && cp != Points_.end()) ++cp;
  while (cp != Points_.end())
  {
    int densityVal = cp->Density();
    Carray densityArray;
    // Add all points of current density.
    while (cp->Density() == densityVal && cp != Points_.end())
      densityArray.push_back( *(cp++) );
    mprintf("Density value %i has %zu points.\n", densityVal, densityArray.size());
    // Sort array by distance
    std::sort(densityArray.begin(), densityArray.end(), Cpoint::dist_sort());
    // Take the average of the points weighted by their position. 
    double wtDistAv = 0.0;
    double sumWts = 0.0;
    //std::vector<double> weights;
    //weights.reserve( densityArray.size() );
    int maxPt = (int)densityArray.size() - 1;
    for (int ip = 0; ip != (int)densityArray.size(); ++ip) 
    {
      double wt = exp( ip - maxPt );
      //mprintf("\t%10i %10u %10u %10g\n", densityVal, ip, maxPt, wt);
      wtDistAv += (densityArray[ip].Dist() * wt);
      sumWts += wt;
      //weights.push_back( wt );
    }
    wtDistAv /= sumWts;
    // Calculate the weighted sample variance
    //double distSD = 0.0;
    //for (unsigned int ip = 0; ip != densityArray.size(); ++ip) {
    //  double diff = densityArray[ip].Dist() - wtDistAv;
    //  distSD += weights[ip] * (diff * diff);
    //}
    //distSD /= sumWts;
    weightedAverage.AddXY(densityVal, wtDistAv); 
    //tempOut.Printf("%i %g %g %g\n", densityVal, wtDistAv, sqrt(distSD), sumWts);
    tempOut.Printf("%i %g %g\n", densityVal, wtDistAv, sumWts);
/*
    // Find the median.
    double median, Q1, Q3;
    if (densityArray.size() == 1) {
      median = densityArray[0].Dist();
      Q1 = median;
      Q3 = median;
    } else {
      unsigned int q3_beg;
      unsigned int med_idx = densityArray.size() / 2; // Always 0 <= Q1 < med_idx
      if ((densityArray.size() % 2) == 0) {
        median = (densityArray[med_idx].Dist() + densityArray[med_idx-1].Dist()) / 2.0;
        q3_beg = med_idx;
      } else {
        median = densityArray[med_idx].Dist();
        q3_beg = med_idx + 1;
      }
      if (densityArray.size() == 2) {
        Q1 = densityArray[0].Dist();
        Q3 = densityArray[1].Dist();
      } else {
        // Find lower quartile
        unsigned int q1_idx = med_idx / 2;
        if ((med_idx % 2) == 0)
          Q1 = (densityArray[q1_idx].Dist() + densityArray[q1_idx-1].Dist()) / 2.0;
        else
          Q1 = densityArray[q1_idx].Dist();
        // Find upper quartile
        unsigned int q3_size = densityArray.size() - q3_beg;
        unsigned int q3_idx = (q3_size / 2) + q3_beg;
        if ((q3_size %2) == 0)
          Q3 = (densityArray[q3_idx].Dist() + densityArray[q3_idx-1].Dist()) / 2.0;
        else
          Q3 = densityArray[q3_idx].Dist();
      }
    }
    mprintf("\tMedian dist value is %g. Q1= %g   Q3= %g\n", median, Q1, Q3);
*/
  }
  tempOut.CloseFile();
  // END CALCULATING WEIGHTED DISTANCE AVERAGE

/*
  // TEST
  tempOut.OpenWrite("temp2.dat");
  std::vector<double> Hist( Points_.back().Density()+1, 0.0 );
  int gWidth = 3;
  double cval = 3.0;
  double two_c_squared = 2.0 * cval * cval;
  mprintf("DBG: cval= %g, Gaussian denominator is %g\n", cval, two_c_squared);
  for (int wtIdx = 0; wtIdx != (int)weightedAverage.Size(); wtIdx++)
  {
    int bval = weightedAverage.X(wtIdx);
    for (int xval = std::max(bval - gWidth, 0);
             xval != std::min(bval + gWidth + 1, (int)Hist.size()); xval++)
    {
      // a: height (weighted average)
      // b: center (density value)
      // c: width
      // x: density value in histogram 
      //int xval = weightedAverage.X(idx);
      //double bval = weightedAverage.X(wtIdx);
      //double bval = (double)wtIdx;
      double diff = (double)(xval - bval);
      //Hist[xval] += (weightedAverage.Y(wtIdx) * exp( -( (diff * diff) / two_c_squared ) ));
      Hist[xval] = std::max(Hist[xval],
                            weightedAverage.Y(wtIdx) * exp( -( (diff * diff) / two_c_squared ) ));
    }
  }
  for (unsigned int idx = 0; idx != Hist.size(); idx++)
    tempOut.Printf("%u %g\n", idx, Hist[idx]);
  tempOut.CloseFile();
  // END TEST
*/
/*
  // TEST
  // Construct best-fit line segments
  tempOut.OpenWrite("temp2.dat");
  double slope, intercept, correl;
  int segment_length = 3;
  DataSet_Mesh Segment;
  Segment.Allocate1D( segment_length );
  for (int wtIdx = 0; wtIdx != (int)weightedAverage.Size(); wtIdx++)
  {
    Segment.Clear();
    for (int idx = std::max(wtIdx - 1, 0); // TODO: use segment_length
             idx != std::min(wtIdx + 2, (int)weightedAverage.Size()); idx++)
        Segment.AddXY(weightedAverage.X(idx), weightedAverage.Y(idx));
    Segment.LinearRegression(slope, intercept, correl, true);
    for (int idx = std::max(wtIdx - 1, 0); // TODO: use segment_length
             idx != std::min(wtIdx + 2, (int)weightedAverage.Size()); idx++)
    {
      double x = weightedAverage.X(idx);
      double y = slope * x + intercept;
      tempOut.Printf("%g %g %i\n", x, y, weightedAverage.X(wtIdx));
    }
  }
  tempOut.CloseFile(); 
  // END TEST
*/

  // BEGIN WEIGHTED RUNNING AVG/SD OF DISTANCES
  // For each point, determine if it is greater than the average of the
  // weighted average distances of the previous, current, and next densities.
  int width = 2;
  int currentDensity = 0;
  int wtIdx = 0;
  double currentAvg = 0.0;
  double deltaSD = 0.0;
  double deltaAv = 0.0;
  int    Ndelta = 0;
  CpptrajFile raOut;
  if (!rafile_.empty()) raOut.OpenWrite(rafile_);
  CpptrajFile raDelta;
  if (!radelta_.empty()) raDelta.OpenWrite(radelta_);
  std::vector<unsigned int> candidateIdxs;
  std::vector<double> candidateDeltas;
  cp = Points_.begin();
  // Skip over points with zero density
  while (cp != Points_.end() && cp->Density() == 0) ++cp;
  while (weightedAverage.X(wtIdx) != cp->Density() && wtIdx < (int)Points_.size())
    ++wtIdx;
  for (Carray::const_iterator point = cp; point != Points_.end(); ++point)
  {
    if (point->Density() != currentDensity) {
      //currentAvg = weightedAverage.Y(wtIdx);
      // New density value. Determine average.
      currentAvg = 0.0;
     // unsigned int Npt = 0; 
      double currentWt = 0.0;
      for (int idx = std::max(wtIdx - width, 0);
               idx != std::min(wtIdx + width + 1, (int)weightedAverage.Size()); idx++)
      {
        //currentAvg += weightedAverage.Y(idx);
        //Npt++;
        double wt = weightedAverage.Y(idx);
        currentAvg += (weightedAverage.Y(idx) * wt);
        currentWt += wt;
      }
      //currentAvg /= (double)Npt;
      currentAvg /= currentWt;
      //smoothAv += currentAvg;
      //smoothSD += (currentAvg * currentAvg);
      //Nsmooth++;
      currentDensity = point->Density();
      if (raOut.IsOpen())
        raOut.Printf("%i %g %g\n", currentDensity, currentAvg, weightedAverage.Y(wtIdx));
      wtIdx++;
    }
    double delta = (point->Dist() - currentAvg);
    if (delta > 0.0) {
      //delta *= log((double)currentDensity);
      if (raDelta.IsOpen())
        raDelta.Printf("%8i %8.3f %8i %8.3f %8.3f\n",
                       currentDensity, delta, point->Fnum()+1, point->Dist(), currentAvg);
      candidateIdxs.push_back( point - Points_.begin() );
      candidateDeltas.push_back( delta );
      deltaAv += delta;
      deltaSD += (delta * delta);
      Ndelta++;
    }
  }
  raOut.CloseFile();
  deltaAv /= (double)Ndelta;
  deltaSD = (deltaSD / (double)Ndelta) - (deltaAv * deltaAv);
  if (deltaSD > 0.0)
    deltaSD = sqrt(deltaSD);
  else
    deltaSD = 0.0;
  if (raDelta.IsOpen())
    raDelta.Printf("#DeltaAvg= %g  DeltaSD= %g\n", deltaAv, deltaSD);
  raDelta.CloseFile();
  int cnum = 0;
  for (unsigned int i = 0; i != candidateIdxs.size(); i++) {
    if (candidateDeltas[i] > (deltaSD)) {
      Points_[candidateIdxs[i]].SetCluster( cnum++ );
      mprintf("\tPoint %u (frame %i, density %i) selected as candidate for cluster %i\n",
              candidateIdxs[i], Points_[candidateIdxs[i]].Fnum()+1,
              Points_[candidateIdxs[i]].Density(), cnum-1);
    }
  }
  // END WEIGHTED AVG/SD OF DISTANCES

/* 
  // Currently doing this by calculating the running average of density vs 
  // distance, then choosing points with distance > twice the SD of the 
  // running average.
  // NOTE: Store in a mesh data set for now in case we want to spline etc later.
  if (avg_factor_ < 1) avg_factor_ = 10; 
  unsigned int window_size = Points_.size() / (unsigned int)avg_factor_;
  mprintf("\tRunning avg window size is %u\n", window_size);
  // FIXME: Handle case where window_size < frames
  DataSet_Mesh runavg;
  unsigned int ra_size = Points_.size() - window_size + 1;
  runavg.Allocate1D( ra_size );
  double dwindow = (double)window_size;
  double sumx = 0.0;
  double sumy = 0.0;
  for (unsigned int i = 0; i < window_size; i++) {
    sumx += (double)Points_[i].Density();
    sumy += Points_[i].Dist();
  }
  runavg.AddXY( sumx / dwindow, sumy / dwindow );
  for (unsigned int i = 1; i < ra_size; i++) {
    unsigned int nextwin = i + window_size - 1;
    unsigned int prevwin = i - 1;
    sumx = (double)Points_[nextwin].Density() - (double)Points_[prevwin].Density() + sumx;
    sumy =         Points_[nextwin].Dist()    -         Points_[prevwin].Dist()    + sumy;
    runavg.AddXY( sumx / dwindow, sumy / dwindow );
  }
  // Write running average
  if (!rafile_.empty()) {
    CpptrajFile raOut;
    if (raOut.OpenWrite(rafile_))
      mprinterr("Error: Could not open running avg file '%s' for write.\n", rafile_.c_str());
    else {
      for (unsigned int i = 0; i != runavg.Size(); i++)
        raOut.Printf("%g %g\n", runavg.X(i), runavg.Y(i));
      raOut.CloseFile();
    }
  }
  double ra_sd;
  double ra_avg = runavg.Avg( ra_sd );
  // Double stdev to use as cutoff for findning anomalously high peaks.
  ra_sd *= 2.0;
  mprintf("\tAvg of running avg set is %g, SD*2.0 (delta cutoff) is %g\n", ra_avg, ra_sd);
  // For each point in density vs distance plot, determine which running
  // average point is closest. If the difference between the point and the
  // running average point is > 2.0 the SD of the running average data,
  // consider it a 'peak'. 
  CpptrajFile raDelta;
  if (!radelta_.empty())
    raDelta.OpenWrite("radelta.dat");
  if (raDelta.IsOpen())
    raDelta.Printf("%-10s %10s %10s\n", "#Frame", "RnAvgPos", "Delta");
  unsigned int ra_position = 0; // Position in the runavg DataSet
  unsigned int ra_end = runavg.Size() - 1;
  int cnum = 0;
  for (Carray::iterator point = Points_.begin();
                        point != Points_.end(); ++point)
  {
    if (ra_position != ra_end) {
      // Is the next running avgd point closer to this point?
      while (ra_position != ra_end) {
        double dens  = (double)point->Density();
        double diff0 = fabs( dens - runavg.X(ra_position  ) );
        double diff1 = fabs( dens - runavg.X(ra_position+1) );
        if (diff1 < diff0)
          ++ra_position; // Next running avg position is closer for this point.
        else
          break; // This position is closer.
      }
    }
    double delta = point->Dist() - runavg.Y(ra_position);
    if (raDelta.IsOpen())
      raDelta.Printf("%-10i %10u %10g", point->Fnum()+1, ra_position, delta);
    if (delta > ra_sd) {
      if (raDelta.IsOpen())
        raDelta.Printf(" POTENTIAL CLUSTER %i", cnum);
      point->SetCluster(cnum++);
    }
    if (raDelta.IsOpen()) raDelta.Printf("\n");
  }
  raDelta.CloseFile();
*/
  int nclusters = cnum;
  mprintf("\tIdentified %i cluster centers from density vs distance peaks.\n", nclusters);
  // Each remaining point is assigned to the same cluster as its nearest
  // neighbor of higher density. Do this recursively until a cluster
  // center is found.
  cnum = -1;
  for (unsigned int idx = 0; idx != Points_.size(); idx++) {
    if (Points_[idx].Cnum() == -1) {// Point is unassigned.
      AssignClusterNum(idx, cnum);
      //mprintf("Finished recursion for index %i\n\n", idx);
    }
  }
  // Sort by cluster number. NOTE: This invalidates NearestIdx
  std::sort( Points_.begin(), Points_.end(), Cpoint::cnum_sort() );
  // Determine where each cluster starts and stops in Points array
  typedef std::vector<unsigned int> Parray;
  Parray C_start_stop;
  C_start_stop.reserve( nclusters * 2 );
  cnum = -1;
  for (Carray::const_iterator point = Points_.begin(); point != Points_.end(); ++point)
  {
    if (point->Cnum() != cnum) {
      if (!C_start_stop.empty()) C_start_stop.push_back(point - Points_.begin()); // end of cluster
      C_start_stop.push_back(point - Points_.begin()); // beginning of cluster
      cnum = point->Cnum();
    }
  }
  C_start_stop.push_back( Points_.size() ); // end of last cluster
  // Noise calculation.
  if (calc_noise_) {
    mprintf("\tDetermining noise frames from cluster borders.\n");
    // For each cluster find a border region, defined as the set of points
    // assigned to that cluster which are within epsilon of any other
    // cluster.
    // NOTE: Could use a set here to prevent duplicate frames.
    typedef std::vector<Parray> Barray;
    Barray borderIndices( nclusters ); // Hold indices of border points for each cluster.
    for (Parray::const_iterator idx0 = C_start_stop.begin();
                                idx0 != C_start_stop.end(); idx0 += 2)
    {
      int c0 = Points_[*idx0].Cnum();
      //mprintf("Cluster %i\n", c0);
      // Check each frame in this cluster.
      for (unsigned int i0 = *idx0; i0 != *(idx0+1); ++i0)
      {
        Cpoint const& point = Points_[i0];
        // Look at each other cluster
        for (Parray::const_iterator idx1 = idx0 + 2;
                                    idx1 != C_start_stop.end(); idx1 += 2)
        {
          int c1 = Points_[*idx1].Cnum();
          // Check each frame in other cluster
          for (unsigned int i1 = *idx1; i1 != *(idx1+1); i1++)
          {
            Cpoint const& other_point = Points_[i1];
            if (FrameDistances_.GetFdist(point.Fnum(), other_point.Fnum()) < epsilon_) {
              //mprintf("\tBorder frame: %i (to cluster %i frame %i)\n",
              //        point.Fnum() + 1, c1, other_point.Fnum() + 1);
              borderIndices[c0].push_back( i0 );
              borderIndices[c1].push_back( i1 );
            }
          }
        }
      }
    }
    if (debug_ > 0)
      mprintf("Warning: Cluster numbers here may not match final cluster numbers.\n"
              "\tBorder Frames:\n");
    for (Parray::const_iterator idx = C_start_stop.begin();
                                idx != C_start_stop.end(); idx += 2)
    {
      int c0 = Points_[*idx].Cnum();
      if (debug_ > 0)
        mprintf("\tCluster %u: %u frames: %u border frames:", c0, *(idx+1) - *idx,
                borderIndices[c0].size());
      if (borderIndices[c0].empty()) {
        if (debug_ > 0) mprintf(" No border points.\n");
      } else {
        int highestDensity = -1;
        // Find highest density in border region.
        for (Parray::const_iterator bidx = borderIndices[c0].begin();
                                    bidx != borderIndices[c0].end(); ++bidx)
        {
          if (highestDensity == -1)
            highestDensity = Points_[*bidx].Density();
          else
            highestDensity = std::max(highestDensity, Points_[*bidx].Density());
          if (debug_ > 0) mprintf(" %i", Points_[*bidx].Fnum()+1);
        }
        if (debug_ > 0) mprintf(". Highest density in border= %i\n", highestDensity);
        // Mark any point with density <= highest border density as noise.
        for (unsigned int i = *idx; i != *(idx+1); i++)
        {
          Cpoint& point = Points_[i];
          if (point.Density() <= highestDensity) {
            point.SetCluster( -1 );
            if (debug_ > 1)
              mprintf("\t\tMarking frame %i as noise (density %i)\n",
                       point.Fnum()+1, point.Density());
          }
        }
      }
    }
  }
  // Add the clusters.
  for (Parray::const_iterator idx = C_start_stop.begin();
                              idx != C_start_stop.end(); idx += 2)
  {
    ClusterDist::Cframes frames;
    for (unsigned int i = *idx; i != *(idx+1); i++) {
      if (Points_[i].Cnum() != -1)
        frames.push_back( Points_[i].Fnum() );
    }
    if (!frames.empty())
      AddCluster( frames );
  }
  // Calculate the distances between each cluster based on centroids
  CalcClusterDistances();

  return 0;
}

// -----------------------------------------------------------------------------
/** This should never be called for the point with highest density
  * which by definition should be a cluster center.
  */
void Cluster_DPeaks::AssignClusterNum(int idx, int& cnum) {
  // Who is the nearest neighbor with higher density. 
  int neighbor_idx = Points_[idx].NearestIdx();
  //mprintf("Index %i nearest neighbor index %i\n", idx, neighbor_idx);
  // SANITY CHECK
  if (neighbor_idx == -1) {
    mprinterr("Internal Error: In Cluster_DPeaks::AssignClusterNum nearest neighbor is -1.\n");
    return;
  }
  if (Points_[neighbor_idx].Cnum() != -1) {
    // Nearest neighbor has a cluster num assigned.
    cnum = Points_[neighbor_idx].Cnum();
    //mprintf("Neighbor index %i is cluster %i\n", neighbor_idx, cnum);
  } else
    // Ask neighbor to find cluster num.
    AssignClusterNum(neighbor_idx, cnum);
  //mprintf("Index %i cnum %i\n", idx, cnum);
  // At this point cnum should be set. One more sanity check.
  if (cnum == -1) {
    mprinterr("Internal Error: In Cluster_DPeaks::AssignClusterNum could not get"
              " cluster num for index %u.\n", idx);
    return;
  }
  Points_[idx].SetCluster( cnum );
}

void Cluster_DPeaks::ClusterResults(CpptrajFile& outfile) const {
   outfile.Printf("#Algorithm: DPeaks epsilon %g\n", epsilon_);
   if (calc_noise_) {
     outfile.Printf("#NOISE_FRAMES:");
     std::vector<int> noiseFrames;
     for (Carray::const_iterator point = Points_.begin();
                                 point != Points_.end(); ++point)
       if (point->Cnum() == -1) noiseFrames.push_back( point->Fnum()+1 );
    std::sort( noiseFrames.begin(), noiseFrames.end() );
    for (std::vector<int>::const_iterator f = noiseFrames.begin(); f != noiseFrames.end(); ++f)
      outfile.Printf(" %i", *f); 
    outfile.Printf("\n");
    outfile.Printf("#Number_of_noise_frames: %zu\n", noiseFrames.size());
  }
}

void Cluster_DPeaks::AddSievedFrames() {
  mprintf("FIXME: Adding sieved frames not yet supported.\n");
}
