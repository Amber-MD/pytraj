#include <cstdio> // sscanf
#include <cstdlib> // atof
#include "DataIO_OpenDx.h"
#include "CpptrajStdio.h"
#include "DataSet_GridFlt.h"
#include "BufferedLine.h"
#include "ProgressBar.h"

bool DataIO_OpenDx::ID_DataFormat( CpptrajFile& infile ) {
  bool isDX = false;
  if (!infile.OpenFile()) {
    std::string firstLine = infile.GetLine();
    if (!firstLine.empty())
      isDX = (firstLine.compare(0, 28, "object 1 class gridpositions") == 0);
    infile.CloseFile();
  }
  return isDX;
}

// DataIO_OpenDx::ReadData()
int DataIO_OpenDx::ReadData(std::string const& fname, 
                            DataSetList& datasetlist, std::string const& dsname)
{
  // Add grid data set. Default to float for now.
  DataSet* ds = datasetlist.AddSet( DataSet::GRID_FLT, dsname, "GRID" );
  if (ds==0) return 1;
  if (LoadGrid(fname.c_str(), *ds)) {
    // Load failed. Erase grid data set.
    DataSetList::const_iterator last = datasetlist.end();
    --last;
    datasetlist.RemoveSet( last );
    return 1;
  }
  return 0;
}

// DataIO_OpenDx::LoadGrid()
int DataIO_OpenDx::LoadGrid(const char* filename, DataSet& ds)
{
  // TODO: This may need to be changed if new 3D types introduced.
  DataSet_GridFlt& grid = static_cast<DataSet_GridFlt&>( ds );
  // Open file
  BufferedLine infile;
  if (infile.OpenFileRead(filename)) return 1;
  // Skip comments
  std::string line = infile.GetLine();
  while (!line.empty() && line[0] == '#') {
    mprintf("\t%s", line.c_str());
    line = infile.GetLine();
  }
  if (line.empty()) {
    mprinterr("Error: Unexpected EOF in DX file %s\n", filename);
    return 1;
  }
  // object 1 class gridpositions counts nx ny nz
  int nx, ny, nz;
  if (sscanf(line.c_str(), "object 1 class gridpositions counts %d %d %d",
             &nx, &ny, &nz) != 3)
  {
    mprinterr("Error: Reading grid counts from DX file %s\n", filename);
    return 1;
  }
  // origin xmin ymin zmin 
  double oxyz[3];
  line = infile.GetLine();
  if (sscanf(line.c_str(), "origin %lg %lg %lg", oxyz, oxyz+1, oxyz+2) != 3) {
    mprinterr("Error: Reading origin line from DX file %s\n", filename);
    return 1;
  }
  // 3x 'delta hx hy hz'
  double dxyz[3];
  Matrix_3x3 delta(0.0);
  bool isNonortho = false;
  int midx = 0;
  for (int i = 0; i < 3; i++, midx += 3) {
    line = infile.GetLine();
    if (sscanf(line.c_str(), "delta %lg %lg %lg", dxyz, dxyz+1, dxyz+2) != 3) {
      mprinterr("Error: Reading delta line from DX file %s\n", filename);
      return 1;
    }
    // Check that only 1 of the 3 values is non-zero. Otherwise non-ortho.
    if (dxyz[i] != (dxyz[0] + dxyz[1] + dxyz[2]))
      isNonortho = true;
    delta[midx  ] = dxyz[0];
    delta[midx+1] = dxyz[1];
    delta[midx+2] = dxyz[2];
  }
  // object 2 class gridconnections counts nx ny nz
  int nxyz[3];
  line = infile.GetLine();
  if (sscanf(line.c_str(), "object 2 class gridconnections counts %d %d %d",
             nxyz, nxyz+1, nxyz+2) != 3)
  {
    mprinterr("Error: Reading grid connections from DX file %s\n", filename);
    return 1;
  }
  // Sanity check for conflicting grid dimensions
  if (nxyz[0] != nx || nxyz[1] != ny || nxyz[2] != nz) {
    mprinterr("Error: Conflicting grid dimensions in input DX density file %s.\n",
              filename);
    mprinterr("Error: Grid positions: %d %d %d\n", nx, ny, nz);
    mprinterr("Error: Grid connections: %d %d %d\n", nxyz[0], nxyz[1], nxyz[2]);
    return 1;
  }
  // object 3 class array type <type> rank <r> times <i>
  // This line describes whether data will be in binary or ascii format.
  line = infile.GetLine();
  if (line.compare(0, 8, "object 3") != 0) {
    mprinterr("Error: DX file %s; expected 'object 3 ...', got [%s]\n",
              filename, line.c_str());
    return 1;
  }
  if (line.find("binary") != std::string::npos) {
    mprinterr("Error: DX file %s; binary DX files not yet supported.\n", filename);
    return 1;
  }
  // Allocate Grid from dims, origin, and spacing
  int err = 0;
  if (isNonortho) {
    // Create unit cell from delta and bins.
    delta[0] *= (double)nx; delta[1] *= (double)nx; delta[2] *= (double)nx;
    delta[3] *= (double)ny; delta[4] *= (double)ny; delta[5] *= (double)ny;
    delta[6] *= (double)nz; delta[7] *= (double)nz; delta[8] *= (double)nz;
    err = grid.Allocate_N_O_Box(nx,ny,nz, Vec3(oxyz), Box(delta));
  } else
    err = grid.Allocate_N_O_D(nx,ny,nz, Vec3(oxyz), Vec3(delta[0],delta[4],delta[8]));
  if (err != 0) { 
    mprinterr("Error: Could not allocate grid.\n");
    return 1;
  }
  grid.GridInfo();
  // Read in data
  size_t gridsize = grid.Size();
  mprintf("\tReading in %zu data elements from DX file.\n", gridsize); 
  size_t ndata = 0;
  ProgressBar progress( gridsize );
  while (ndata < gridsize) {
    if (infile.Line() == 0) {
      mprinterr("Error: Unexpected EOF hit in %s\n", filename);
      return 1;
    }
    int nTokens = infile.TokenizeLine(" \t");
    for (int j = 0; j < nTokens; j++) {
      if (ndata >= gridsize) {
        mprintf("Warning: Too many grid points found. Only reading %zu grid points.\n", gridsize);
        mprintf("Warning: Check that data region ends with a newline.\n");
        break;
      }
      grid[ndata++] = (float)atof(infile.NextToken());
    }
    progress.Update( ndata );
  }
  // Set dimensions
  // FIXME: This should be integrated with allocation
  //grid.SetDim(Dimension::X, Dimension(oxyz[0], dx, nx, "X"));
  //grid.SetDim(Dimension::Y, Dimension(oxyz[1], dy, ny, "Y"));
  //grid.SetDim(Dimension::Z, Dimension(oxyz[2], dz, nz, "Z"));
  return 0;
}

// DataIO_OpenDx::WriteData3D()
int DataIO_OpenDx::WriteData3D(std::string const& fname, DataSetList const& setList)
{
  // Open output file
  CpptrajFile outfile;
  if (outfile.OpenWrite(fname)) {
    mprinterr("Error: Could not open OpenDX output file.\n");
    return 1;
  }
  // Warn about writing multiple sets
  if (setList.size() > 1)
    mprintf("Warning: %s: Writing multiple 3D sets in OpenDX format may result in unexpected behavior\n", fname.c_str());
  int err = 0;
  for (DataSetList::const_iterator set = setList.begin(); set != setList.end(); ++set)
    err += WriteSet3D( *(*set), outfile );
  return err;
}

// DataIO_OpenDx::WriteSet3D()
int DataIO_OpenDx::WriteSet3D( DataSet const& setIn, CpptrajFile& outfile) {
  if (setIn.Ndim() != 3) {
    mprinterr("Internal Error: DataSet %s in DataFile %s has %zu dimensions, expected 3.\n",
              setIn.Legend().c_str(), outfile.Filename().full(), setIn.Ndim());
    return 1;
  }
  DataSet_3D const& set = static_cast<DataSet_3D const&>( setIn );

  // Print the OpenDX header
  size_t gridsize = set.Size();
  outfile.Printf("object 1 class gridpositions counts %d %d %d\n",
                 set.NX(), set.NY(), set.NZ());
  Vec3 const& oxyz = set.GridOrigin();
  outfile.Printf("origin %lg %lg %lg\n", oxyz[0], oxyz[1], oxyz[2]);
  Matrix_3x3 ucell = set.Ucell();
  double nx = (double)set.NX();
  double ny = (double)set.NY();
  double nz = (double)set.NZ();
  outfile.Printf("delta %lg %lg %lg\n", ucell[0]/nx, ucell[1]/nx, ucell[2]/nx);
  outfile.Printf("delta %lg %lg %lg\n", ucell[3]/ny, ucell[4]/ny, ucell[5]/ny);
  outfile.Printf("delta %lg %lg %lg\n", ucell[6]/nz, ucell[7]/nz, ucell[8]/nz);
  outfile.Printf("object 2 class gridconnections counts %d %d %d\n",
                 set.NX(), set.NY(), set.NZ());
  outfile.Printf(
    "object 3 class array type double rank 0 items %d data follows\n",
    gridsize);

  // Now print out the data. It is already in row-major form (z-axis changes
  // fastest), so no need to do any kind of data adjustment
  for (size_t i = 0UL; i < gridsize - 2UL; i += 3UL)
    outfile.Printf("%g %g %g\n", set[i], set[i+1], set[i+2]);
  // Print out any points we may have missed
  switch (gridsize % 3) {
    case 2: outfile.Printf("%g %g\n", set[gridsize-2], set[gridsize-1]); break;
    case 1: outfile.Printf("%g\n", set[gridsize-1]); break;
  }

  // Print tail
  // TODO: Make this an option
  //if (mode_ == CENTER)
  //  outfile.Printf("\nobject \"density (%s) [A^-3]\" class field\n",
  //                 centerMask_.MaskString());
  //else
    outfile.Printf("\nobject \"density [A^-3]\" class field\n");
  return 0;
}
