#include "GridAction.h"
#include "CpptrajStdio.h"

// GridAction::HelpText
const char* GridAction::HelpText =
  "\t{ data <dsname> | boxref <ref name/tag> <nx> <ny> <nz> |\n"
  "\t  <nx> <dx> <ny> <dy> <nz> <dz> [gridcenter <cx> <cy> <cz>] }\n"
  "\t[box|origin|center <mask>] [negative] [name <gridname>]";

static inline void CheckEven(int& N, char dir) {
  if (N % 2 == 1) {
    ++N;
    mprintf("Warning: number of grid points must be even. Incrementing N%c by 1 to %u\n", dir, N);
  }
}

// GridAction::GridInit()
DataSet_GridFlt* GridAction::GridInit(const char* callingRoutine, ArgList& argIn, DataSetList& DSL) 
{
  DataSet_GridFlt* Grid = 0; 
  bool specifiedCenter = false;
  std::string dsname = argIn.GetStringKey("data");
  std::string refname = argIn.GetStringKey("boxref");
  if (!dsname.empty()) { 
    // Get existing grid dataset
    Grid = (DataSet_GridFlt*)DSL.FindSetOfType( dsname, DataSet::GRID_FLT );
    if (Grid == 0) {
      mprinterr("Error: %s: Could not find grid data set with name %s\n",
                callingRoutine, dsname.c_str());
      return 0;
    }
  } else if (!refname.empty()) {
    // Get grid parameters from reference structure box.
    DataSet_Coords_REF* REF = (DataSet_Coords_REF*)DSL.GetReferenceFrame( refname );
    if (REF == 0) return 0;
    if (REF->Top().ParmBox().Type() == Box::NOBOX) {
      mprinterr("Error: Reference '%s' does not have box information.\n", refname.c_str());
      return 0;
    }
    int nx = argIn.getNextInteger(-1);
    int ny = argIn.getNextInteger(-1);
    int nz = argIn.getNextInteger(-1);
    if (nx < 1 || ny < 1 || nz < 1) {
      mprinterr("Error:  %s: Invalid grid sizes\n", callingRoutine);
      return 0;
    }
    Grid = (DataSet_GridFlt*)DSL.AddSet( DataSet::GRID_FLT, argIn.GetStringKey("name"), "GRID" );
    if (Grid == 0) return 0;
    if (Grid->Allocate_N_O_Box(nx,ny,nz,Vec3(0.0),REF->RefFrame().BoxCrd())) return 0;
  } else {
    // Create new data set.
    // Get nx, dx, ny, dy, nz, dz
    int nx = argIn.getNextInteger(-1);
    double dx = argIn.getNextDouble(-1);
    int ny = argIn.getNextInteger(-1);
    double dy = argIn.getNextDouble(-1);
    int nz = argIn.getNextInteger(-1);
    double dz = argIn.getNextDouble(-1);
    if (nx < 1 || ny < 1 || nz < 1 ||
        dx < 0 || dy < 0 || dz < 0)
    {
      mprinterr("Error: %s: Invalid grid size/spacing.\n", callingRoutine);
      mprinterr("       nx=%i ny=%i nz=%i | dx=%.3f dy=%.3f dz=%.3f\n",
                nx, ny, nz, dx, dy, dz);
      return 0;
    }
    // For backwards compat., enforce even grid spacing.
    CheckEven( nx, 'X' ); CheckEven( ny, 'Y' ); CheckEven( nz, 'Z' );
    Vec3 gridctr(0.0, 0.0, 0.0);
    if (argIn.hasKey("gridcenter")) {
      double cx = argIn.getNextDouble(0.0);
      double cy = argIn.getNextDouble(0.0);
      double cz = argIn.getNextDouble(0.0);
      gridctr.SetVec(cx, cy, cz);
      specifiedCenter = true;
    }
    Grid = (DataSet_GridFlt*)DSL.AddSet( DataSet::GRID_FLT, argIn.GetStringKey("name"), "GRID" );
    if (Grid == 0) return 0;
    // Set up grid from dims, center, and spacing
    // NOTE: # of grid points in each direction with be forced to be even.
    if (Grid->Allocate_N_C_D(nx, ny, nz, gridctr, Vec3(dx,dy,dz))) return 0;
  }
  // Determine offset, Box/origin/center
  mode_ = ORIGIN;
  if (argIn.hasKey("box"))
    mode_ = BOX;
  else if (argIn.hasKey("origin"))
    mode_ = ORIGIN;
  else if (argIn.Contains("center")) {
    std::string maskexpr = argIn.GetStringKey("center");
    if (maskexpr.empty()) {
      mprinterr("Error: 'center' requires <mask>\n");
      return 0;
    }
    centerMask_.SetMaskString( maskexpr );
    mode_ = MASKCENTER;
  }
  if (specifiedCenter) {
    // If center was specified, do not allow an offset
    if (mode_ != ORIGIN)
      mprintf("Warning: Grid offset args (box/center) not allowed with 'gridcenter'.\n"
              "Warning: No offset will be used.\n");
    mode_ = SPECIFIEDCENTER;
  }
  // Negative
  if (argIn.hasKey("negative"))
    increment_ = -1.0;
  else
    increment_ = 1.0;

  return Grid;
}

// GridAction::GridInfo() 
void GridAction::GridInfo(DataSet_GridFlt const& grid) {
  if (mode_ == BOX)
    mprintf("\tOffset for points is box center.\n");
  else if (mode_ == MASKCENTER)
    mprintf("\tOffset for points is center of atoms in mask [%s]\n",
            centerMask_.MaskString());
  if (increment_ > 0)
    mprintf("\tCalculating positive density.\n");
  else
    mprintf("\tCalculating negative density.\n");
  grid.GridInfo();
}

// GridAction::GridSetup()
int GridAction::GridSetup(Topology const& currentParm) {
  // Check box
  if (mode_ == BOX) {
    if (currentParm.BoxType()!=Box::ORTHO) {
      mprintf("Warning: Code to shift to the box center is not yet\n");
      mprintf("Warning: implemented for non-orthorhomibic unit cells.\n");
      mprintf("Warning: Shifting to the origin instead.\n");
      mode_ = ORIGIN;
    }
  } else if (mode_ == MASKCENTER) {
    if ( currentParm.SetupIntegerMask( centerMask_ ) ) return 1;
    centerMask_.MaskInfo();
    if ( centerMask_.None() ) {
      mprinterr("Error: No atoms selected for grid center mask [%s]\n", centerMask_.MaskString());
      return 1;
    }
  }
  return 0;
}
