#include <cmath>        // log
#include <iostream>     // ofstream
#include <fstream>      // ofstream
#include "Action_GridFreeEnergy.h"
#include "CpptrajStdio.h" // mprintf, mprinterr
#include "StringRoutines.h" // integerToString
#include "Constants.h" // GASK_KCAL, SMALL

// CONSTRUCTOR
Action_GridFreeEnergy::Action_GridFreeEnergy() :
  maxVoxelOccupancyCount_(600), // NOTE: See header for comments.
  tempInKevin_(293.0),
  grid_(0)
{}

void Action_GridFreeEnergy::Help() {
  mprintf("\t<filename>\n%s\n\t<mask>\n", GridAction::HelpText);
}

// Action_GridFreeEnergy::init()
Action::RetType Action_GridFreeEnergy::Init(ArgList& actionArgs, TopologyList* PFL, DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  // Get output filename
  std::string filename = actionArgs.GetStringNext();
  if (filename.empty()) {
    mprinterr("Error: GridFreeEnergy: no output filename specified.\n");
    return Action::ERR;
  }
  // Get grid options (<nx> <dx> <ny> <dy> <nz> <dz> [box|origin] [negative])
  grid_ = GridInit( "GridFreeEnergy", actionArgs, *DSL );
  if (grid_ == 0) return Action::ERR;

  // Get mask
  std::string maskexpr = actionArgs.GetMaskNext();
  if (maskexpr.empty()) {
    mprinterr("Error: GridFreeEnergy: No mask specified.\n");
    return Action::ERR;
  }
  mask_.SetMaskString(maskexpr);

  // Get extra args
  tempInKevin_ = actionArgs.getKeyDouble("temp", 293.0);

  // Setup output file
  DataFile* outfile = DFL->AddSetToFile(filename, (DataSet*)grid_);
  if (outfile == 0) {
    mprinterr("Error: GridFreeEnergy: Could not set up output file %s\n", filename.c_str());
    return Action::ERR;
  }
  //grid_.PrintXplor( filename_, "", "REMARKS Change in Free energy from bulk solvent with bin normalisation of " + integerToString(currentLargestVoxelOccupancyCount) );

  // Info
  mprintf("    GridFreeEnergy\n");
  GridInfo( *grid_ );
  mprintf("\tGrid will be printed to file %s\n",filename.c_str());
  mprintf("\tMask expression: [%s]\n",mask_.MaskString());
  mprintf("\tTemp is : %f K\n",tempInKevin_);

  // Allocate grid
  //if (GridAllocate()) return 1;

  return Action::OK;
}

// Action_GridFreeEnergy::setup()
Action::RetType Action_GridFreeEnergy::Setup(Topology* currentParm, Topology** parmAddress) {
  // Setup grid, checks box info.
  if (GridSetup( *currentParm )) return Action::ERR;

  // Setup mask
  if (currentParm->SetupIntegerMask( mask_ ))
    return Action::ERR;
  mprintf("\t[%s] %i atoms selected.\n", mask_.MaskString(), mask_.Nselected());
  if (mask_.None()) {
    mprinterr("Error: GridFreeEnergy: No atoms selected for parm %s\n", currentParm->c_str());
    return Action::ERR;
  }

  return Action::OK;
}

// Action_GridFreeEnergy::action()
Action::RetType Action_GridFreeEnergy::DoAction(int frameNum, Frame* currentFrame, 
                                                Frame** frameAddress) 
{
  GridFrame( *currentFrame, mask_, *grid_ );
  return Action::OK;
}

// Action_GridFreeEnergy::print()
void Action_GridFreeEnergy::Print() {
  /* How times does this occupancy count value arise?
   *    i.e. if  
   *                 voxelOccupancyCount[50] = 10 
   * 
   *         then there are 10 voxels with an occupancy count 50
   */
  std::vector<int> voxelOccupancyCount;
  // Largest occupancy count
  int currentLargestVoxelOccupancyCount;
  /* Value of most frequent voxel occupancy count
   *
   * i.e. if 
   *                 voxelOccupancyCount[50] = 10 
   *                 voxelOccupancyCount[51] = 100
   *                 voxelOccupancyCount[52] = 5
   *      then
   *                 mostFrequentVoxelOccupancy would be 51
   */
  int mostFrequentVoxelOccupancy = 0;

  // Zero the voxelOccupancyCount array
  voxelOccupancyCount.assign(maxVoxelOccupancyCount_, 0);

  // Determine the frequency the bin populations
  for (DataSet_GridFlt::iterator gval = grid_->begin(); gval != grid_->end(); ++gval) {
        int bincount = (int) *gval;

        // Increase the array size if needs be
        if (bincount >= (int)voxelOccupancyCount.size()){
          voxelOccupancyCount.resize( bincount+1, 0 );
        }
        // Add one to the counter for the fequency of this bin population
        ++voxelOccupancyCount[ bincount ];
  }



  // Walk the voxelOccupancyCount[] array and determine which
  // occupancy is the most frequent (excluding 0).
  // TODO develop a smart way to work this out....
  currentLargestVoxelOccupancyCount = 0;
  mostFrequentVoxelOccupancy = 0;

  std::ofstream hist;
  hist.open ("hist.dat");
  hist << "#CDBG: i, voxelOccupancyCount[i]" << std::endl;

  // Start from 1, to exclude voxelOccupancyCount[0]
  for (int i = 1; i < (int)voxelOccupancyCount.size(); ++i) {
    hist <<  i << " " << voxelOccupancyCount[i] << std::endl;
    // Determine which voxel has the higest occupancy count
    // i.e. which is the most frequent value.
    if (voxelOccupancyCount[i] > currentLargestVoxelOccupancyCount) {
      currentLargestVoxelOccupancyCount = voxelOccupancyCount[i];
      mostFrequentVoxelOccupancy = i;
    }
  }

  hist.close();
  


  mprintf("CDBG: Most frequent occupancy is %i (%i occurrences)\n", mostFrequentVoxelOccupancy,
          currentLargestVoxelOccupancyCount);
  // The assumption here is that mostFrequentVoxelOccupancy corresponds to the
  // the occupancy count of bulk solvent
  for (DataSet_GridFlt::iterator gval = grid_->begin(); gval != grid_->end(); ++gval) {
    float value = *gval;
    // Avoid log(0) values since this will result in an inf
    if ( (value / mostFrequentVoxelOccupancy) < Constants::SMALL){
      *gval = 0.0;
    } else {
      *gval = (-1.0 * Constants::GASK_KCAL * tempInKevin_ *
                log ((value / mostFrequentVoxelOccupancy)));
    }
  }
}
