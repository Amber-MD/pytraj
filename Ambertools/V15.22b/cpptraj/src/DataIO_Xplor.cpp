#include <cstdio> // sscanf
#include "DataIO_Xplor.h"
#include "CpptrajStdio.h"
#include "DataSet_GridFlt.h"
#include "BufferedLine.h"
#include "ProgressBar.h"

inline int ErrorMsg(const char* msg) {
  mprinterr( msg );
  return 1;
}

// DataIO_Xplor::ReadData()
int DataIO_Xplor::ReadData(std::string const& fname, ArgList& argIn,
                            DataSetList& datasetlist, std::string const& dsname)
{
  // Add grid data set. Default to float for now.
  DataSet* ds = datasetlist.AddSet( DataSet::GRID_FLT, dsname, "GRID" );
  if (ds==0) return 1;
  DataSet_GridFlt& grid = static_cast<DataSet_GridFlt&>( *ds );

  BufferedLine infile;
  if (infile.OpenFileRead(fname)) return 1;
  // First line is ignored
  const char* ptr = infile.Line();
  if (ptr == 0) return ErrorMsg("Error: Unexpected end of file.\n");
  // Next line is number of remark lines?
  int nremarks = 0;
  ptr = infile.Line();
  if (ptr == 0 || sscanf(ptr, "%i", &nremarks) != 1) 
    return ErrorMsg("Error: Could not get # remarks\n");
  mprintf("\t%i remarks\n", nremarks);
  for (int i = 0; i < nremarks; i++) {
    ptr = infile.Line();
    mprintf("\t%s", ptr);
  }
  // Next is num grid points, start grid point, stop grid point for XYZ.
  int GridPts[9];
  ptr = infile.Line();
  if ( sscanf(ptr, "%8i%8i%8i%8i%8i%8i%8i%8i%8i", 
              GridPts  , GridPts+1, GridPts+2,
              GridPts+3, GridPts+4, GridPts+5,
              GridPts+6, GridPts+7, GridPts+8) != 9 )
    return ErrorMsg("Error: Could not read grid dimensions\n");
  // Next is cell x y z alpha beta gamma
  double CellDim[6];
  ptr = infile.Line();
  if ( sscanf(ptr, "%12lf%12lf%12lf%12lf%12lf%12lf",
              CellDim  , CellDim+1, CellDim+2,
              CellDim+3, CellDim+4, CellDim+5) != 6 )
    return ErrorMsg("Error: Could not read cell dimensions.\n");
  // Determine if grid is orthogonal and allocate.
  Box box( CellDim );
  int err = 0;
  if (box.Type() == Box::ORTHO) {
    // Allocate orthogonal grid
    Vec3 spacing( CellDim[0] / (double)GridPts[0],
                  CellDim[1] / (double)GridPts[3],
                  CellDim[2] / (double)GridPts[6] );
    Vec3 origin( (double)GridPts[1] * spacing[0],
                 (double)GridPts[4] * spacing[1],
                 (double)GridPts[7] * spacing[2] );
    err = grid.Allocate_N_O_D(GridPts[0], GridPts[3], GridPts[6], origin, spacing);
  } else {
    // Allocate non-orthogonal grid. Determine where origin is based on ucell
    // and start grid points.
    Matrix_3x3 ucell, recip;
    box.ToRecip(ucell, recip);
    // Turn ucell into delta. Use X axis only.
    Vec3 origin( (ucell[0] / (double)GridPts[0]) * (double)GridPts[1],
                 (ucell[1] / (double)GridPts[0]) * (double)GridPts[1],
                 (ucell[2] / (double)GridPts[0]) * (double)GridPts[1] );
    err = grid.Allocate_N_O_Box(GridPts[0], GridPts[3], GridPts[6], origin, box);
  }
  if (err != 0) return ErrorMsg("Error: Could not allocate grid.\n");
  grid.GridInfo();
  mprintf("\tReading in %zu data elements from XPLOR file.\n", grid.Size());
  // Next is 'ZYX'
  ptr = infile.Line();
  if ( ptr == 0 || (ptr[0] != 'Z' || ptr[1] != 'Y' || ptr[2] != 'X' ))
    return ErrorMsg("Error: Expected 'ZYX'\n"); 
  // Read grid points
  ProgressBar progress( grid.NZ() );
  for (size_t k = 0; k < grid.NZ(); ++k) {
    progress.Update( k );
    ptr = infile.Line(); // Reads starting grid bin, not used. 
    for (size_t j = 0; j < grid.NY(); ++j) {
      size_t i = 0;
      while (i < grid.NX()) {
        ptr = infile.Line();
        if (ptr == 0) {
          mprinterr("Error reading grid value at ijk={%zu %zu %zu}\n", i, j, k);
          return 1;
        }
        size_t nread = (size_t)sscanf(ptr, "%12lf%12lf%12lf%12lf%12lf%12lf",
                                      CellDim  , CellDim+1, CellDim+2,
                                      CellDim+3, CellDim+4, CellDim+5);
        for (size_t n = 0; n < nread; n++)
          grid.SetElement(i++, j, k, (float)CellDim[n]);
      }
    }
  }
  // Set dimensions
  // FIXME: This should be integrated with allocation
  //grid.SetDim(Dimension::X, Dimension(origin[0], spacing[0], GridPts[0], "X"));
  //grid.SetDim(Dimension::Y, Dimension(origin[1], spacing[1], GridPts[1], "Y"));
  //grid.SetDim(Dimension::Z, Dimension(origin[2], spacing[2], GridPts[2], "Z"));

  return 0;
}

// DataIO_Xplor::WriteData3D()
int DataIO_Xplor::WriteData3D(std::string const& fname, DataSetList const& setList)
                              
{
  // Open output file
  CpptrajFile outfile;
  if (outfile.OpenWrite( fname )) {
    mprinterr("Error: Could not open Xplor output file.\n");
    return 1;
  }
  // Warn about writing multiple sets
  if (setList.size() > 1)
    mprintf("Warning: %s: Writing multiple 3D sets in XPLOR format may result in unexpected behavior\n", fname.c_str());
  int err = 0;
  for (DataSetList::const_iterator set = setList.begin(); set != setList.end(); ++set)
    err += WriteSet3D( *(*set), outfile );
  return err;
}

// DataIO_Xplor::WriteSet3D()
int DataIO_Xplor::WriteSet3D( DataSet const& setIn, CpptrajFile& outfile) {
  if (setIn.Ndim() != 3) {
    mprinterr("Internal Error: DataSet %s in DataFile %s has %zu dimensions, expected 3.\n",
              setIn.Legend().c_str(), outfile.Filename().full(), setIn.Ndim());
    return 1;
  }
  DataSet_3D const& set = static_cast<DataSet_3D const&>( setIn );
  // Title
  outfile.Printf("%s\n", title_.c_str());
  // Remarks - Use set legend 
  outfile.Printf("%8i\n%s\n",1,set.Legend().c_str());
  // Header - Num grid points, start grid point, stop grid point
  // NOTE: The commented-out section is how grid was set up in ptraj/cpptraj
  //       before. The new method gives maps that match DX output.
  //outfile.Printf("%8i%8i%8i",   set.NX(), -set.NX()/2 + 1, set.NX()/2 );
  //outfile.Printf("%8i%8i%8i",   set.NY(), -set.NY()/2 + 1, set.NY()/2 );
  //outfile.Printf("%8i%8i%8i\n", set.NZ(), -set.NZ()/2 + 1, set.NZ()/2 );
  // Locate the indices of the absolute origin in order to find starting
  // indices for each axis. FIXME: Is this correct?
  int grid_min_x, grid_min_y, grid_min_z;
  set.BinIndices(0.0, 0.0, 0.0, grid_min_x, grid_min_y, grid_min_z);
  if (grid_min_x != 0) grid_min_x = -grid_min_x;
  if (grid_min_y != 0) grid_min_y = -grid_min_y;
  if (grid_min_z != 0) grid_min_z = -grid_min_z;
  //int grid_min_x = (int)(set.OX()/set.DX());
  outfile.Printf("%8i%8i%8i%8i%8i%8i%8i%8i%8i\n",
    set.NX(), grid_min_x, grid_min_x + set.NX() - 1,
    set.NY(), grid_min_y, grid_min_y + set.NY() - 1,
    set.NZ(), grid_min_z, grid_min_z + set.NZ() - 1);
  // Header - cell x y z alpha beta gamma
  Box box( set.Ucell() );
  outfile.Printf("%12.5f%12.5f%12.5f%12.5f%12.5f%12.5f\n",
                 box[0], box[1], box[2], box[3], box[4], box[5]);
  outfile.Printf("ZYX\n");
  // Print grid bins
  for (size_t k = 0; k < set.NZ(); ++k) {
    outfile.Printf("%8i\n", k);
    for (size_t j = 0; j < set.NY(); ++j) {
      int col = 1;
      for (size_t i = 0; i < set.NX(); ++i) {
        outfile.Printf("%12.5f", set.GetElement(i, j, k));
        if ( (col % 6)==0 )
          outfile.Printf("\n");
        ++col;
      }
      if ( (col-1) % 6 != 0 )
        outfile.Printf("\n");
    }
  }
  outfile.Printf("%8i\n", -9999);
  return 0;
}
