#include "Cluster_ReadInfo.h"
#include "CpptrajStdio.h"
#include "BufferedLine.h"

Cluster_ReadInfo::Cluster_ReadInfo() {}

void Cluster_ReadInfo::Help() {
  mprintf("\t[{readtxt|readinfo} infofile <file>]\n");
}

int Cluster_ReadInfo::SetupCluster(ArgList& analyzeArgs) {
  filename_ = analyzeArgs.GetStringKey("infofile");
  if (filename_.empty()) {
    mprinterr("Error: No cluster info filename given.\n");
    return 1;
  }
  return 0;
}

void Cluster_ReadInfo::ClusteringInfo() {
  mprintf("\tREADINFO: Reading cluster information from previous run, file %s.\n",
          filename_.c_str());
}

static int Err(int code) {
  switch (code) {
    case 0: mprinterr("Error: Could not open info file.\n"); break;
    case 1: mprinterr("Error: Unexpected end of info file.\n"); break;
    case 2: mprinterr("Error: Invalid number of clusters in info file.\n"); break;
    case 3: mprinterr("Error: Invalid number of frames in info file.\n"); break;
  }
  return 1;
}

int Cluster_ReadInfo::Cluster() {
  BufferedLine infile;
  if (infile.OpenFileRead( filename_ )) return Err(0);
  const char* ptr = infile.Line();
  if (ptr == 0) return Err(1);
  ArgList infoLine( ptr, " " );
  int nclusters = infoLine.getKeyInt("#Clustering:", -1);
  if (nclusters == -1) return Err(2);
  int nframes = infoLine.getKeyInt("clusters", -1);
  if (nframes == -1) return Err(3);
  if (nframes != (int)FrameDistances_.Nframes()) {
    mprinterr("Error: # frames in cluster info file (%i) does not match"
              " current # frames (%zu)\n", nframes, FrameDistances_.Nframes());
    return 1;
  }
  // Scan down to clusters
  while (ptr[0] == '#') {
    ptr = infile.Line();
    if (ptr == 0) return Err(1);
    // Save previous clustering info. Includes newline.
    if (ptr[1] == 'A' && ptr[2] == 'l' && ptr[3] == 'g')
      algorithm_.assign( ptr + 12 ); // Right past '#Algorithm: '
  }
  // Read clusters
  ClusterDist::Cframes frames;
  for (int cnum = 0; cnum != nclusters; cnum++) {
    if (ptr == 0) return Err(1);
    frames.clear();
    // TODO: Check for busted lines?
    for (int fidx = 0; fidx != nframes; fidx++) {
      if (ptr[fidx] == 'X')
        frames.push_back( fidx );
    }
    AddCluster( frames );
    mprintf("\tRead cluster %i, %zu frames.\n", cnum, frames.size());
    ptr = infile.Line();
  }
  infile.CloseFile();
  mprintf("\tCalculating the distances between each cluster based on centroids.\n");
  CalcClusterDistances();
  return 0;
}

void Cluster_ReadInfo::ClusterResults(CpptrajFile& outfile) const {
  outfile.Printf("#Algorithm: Read from file '%s'\n", filename_.c_str());
  if (!algorithm_.empty())
    outfile.Printf("#Original algorithm: %s", algorithm_.c_str());
}
