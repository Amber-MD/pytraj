// DataSetList
// This also includes basic DataSet class and dataType
#include "DataSetList.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // convertToInteger and DigitWidth
#include "ArgList.h"
#include "Range.h"
#include "Constants.h"
// Data types go here
#include "DataSet_double.h"
#include "DataSet_float.h"
#include "DataSet_integer.h"
#include "DataSet_string.h"
#include "DataSet_MatrixDbl.h"
#include "DataSet_MatrixFlt.h"
#include "DataSet_Coords_CRD.h"
#include "DataSet_Vector.h"
#include "DataSet_Modes.h"
#include "DataSet_GridFlt.h"
#include "DataSet_RemLog.h"
#include "DataSet_Mesh.h"
#include "DataSet_Coords_TRJ.h"
#include "DataSet_Coords_REF.h"

// ----- STATIC VARS / ROUTINES ------------------------------------------------
// IMPORTANT: THIS ARRAY MUST CORRESPOND TO DataSet::DataType
const DataSetList::DataToken DataSetList::DataArray[] = {
  { "unknown",     0                           }, // UNKNOWN_DATA
  { "double",        DataSet_double::Alloc     }, // DOUBLE
  { "float",         DataSet_float::Alloc      }, // FLOAT
  { "integer",       DataSet_integer::Alloc    }, // INTEGER
  { "string",        DataSet_string::Alloc     }, // STRING
  { "double matrix", DataSet_MatrixDbl::Alloc  }, // MATRIX_DBL
  { "float matrix",  DataSet_MatrixFlt::Alloc  }, // MATRIX_FLT
  { "coordinates",   DataSet_Coords_CRD::Alloc }, // COORDS
  { "vector",        DataSet_Vector::Alloc     }, // VECTOR
  { "eigenmodes",    DataSet_Modes::Alloc      }, // MODES
  { "float grid",    DataSet_GridFlt::Alloc    }, // GRID_FLT
  { "remlog",        DataSet_RemLog::Alloc     }, // REMLOG
  { "X-Y mesh",      DataSet_Mesh::Alloc       }, // XYMESH
  { "trajectories",  DataSet_Coords_TRJ::Alloc }, // TRAJ
  { "reference",     DataSet_Coords_REF::Alloc }, // REF_FRAME
  { 0, 0 }
};

const char* DataSetList::SetString(DataSet::DataType d) {
  return DataArray[d].Description;
}

// CONSTRUCTOR
DataSetList::DataSetList() :
  debug_(0),
  ensembleNum_(-1),
  hasCopies_(false),
  dataSetsPending_(false) 
{}

// DESTRUCTOR
DataSetList::~DataSetList() {
  Clear();
}

void DataSetList::Clear() {
  if (!hasCopies_)
    for (DataListType::iterator ds = DataList_.begin(); ds != DataList_.end(); ++ds) 
      delete *ds;
  DataList_.clear();
  hasCopies_ = false;
  dataSetsPending_ = false;
} 

DataSetList& DataSetList::operator+=(DataSetList const& rhs) {
  // It is OK if rhs does not have copies, but this should have copies.
  // For now just set hasCopies to true.
  hasCopies_ = true; 
  // Append rhs entries to here
  for (DataListType::const_iterator DS = rhs.begin(); DS != rhs.end(); ++DS)
    DataList_.push_back( *DS );
  return *this;
}

// DataSetList::RemoveSet()
// NOTE: In order to call erase, must use iterator and not const_iterator.
//       Hence, the conversion. The new standard *should* allow const_iterator
//       to be passed to erase(), but this is currently not portable.
/** Erase element pointed to by posIn from the list. */
void DataSetList::RemoveSet( const_iterator posIn ) {
  std::vector<DataSet*>::iterator pos = DataList_.begin() + (posIn - DataList_.begin());
  if (!hasCopies_) delete *pos;
  DataList_.erase( pos ); 
} 

// DataSetList::RemoveSet()
void DataSetList::RemoveSet( DataSet* dsIn ) {
  for (std::vector<DataSet*>::iterator pos = DataList_.begin();
                                       pos != DataList_.end(); ++pos)
  {
    if ( (*pos) == dsIn ) {
      if (!hasCopies_) delete *pos;
      DataList_.erase( pos );
      break;
    }
  }
}

// DataSetList::SetDebug()
void DataSetList::SetDebug(int debugIn) {
  debug_ = debugIn;
  if (debug_>0) 
    mprintf("DataSetList Debug Level set to %i\n",debug_);
}

/** Call Allocate for each 1D DataSet in the list. */
void DataSetList::AllocateSets(long int maxFrames) {
  if (maxFrames < 1L) return;
  for (DataListType::iterator ds = DataList_.begin(); ds != DataList_.end(); ++ds)
  {
    if ((*ds)->Ndim() == 1) // 1D DataSet
      ((DataSet_1D*)(*ds))->Allocate1D( (size_t)maxFrames );
    else if ((*ds)->Ndim() == 4) // COORDS DataSet
      ((DataSet_Coords*)(*ds))->AllocateCoords( (size_t)maxFrames );
  }
}

/* DataSetList::SetPrecisionOfDataSets()
 * Set the width and precision for all datasets in the list.
 */
void DataSetList::SetPrecisionOfDataSets(std::string const& nameIn, int widthIn,
                                         int precisionIn)
{
  if (widthIn < 1)
    mprinterr("Error: Invalid data width (%i)\n", widthIn);
  else {
    DataSetList Sets = GetMultipleSets( nameIn );
    for (DataSetList::const_iterator ds = Sets.begin(); ds != Sets.end(); ++ds) 
      (*ds)->SetPrecision(widthIn, precisionIn);
  }
}

// DataSetList::ParseArgString()
/** Separate argument nameIn specifying DataSet into name, index, and 
  * attribute parts.
  * Possible formats:
  *  - "<name>"         : Plain dataset name.
  *  - "<name>:<index>" : Dataset within larger overall set (e.g. perres:1)
  *  - "<name>[<attr>]" : Dataset with name and given attribute (e.g. rog[max])
  *  - "<name>[<attr>]:<index>" : 
  *       Dataset with name, given attribute, and index (e.g. NA[shear]:1)
  */
std::string DataSetList::ParseArgString(std::string const& nameIn, std::string& idx_arg,
                                        std::string& attr_arg)
{
  std::string dsname( nameIn );
  attr_arg.clear();
  //mprinterr("DBG: ParseArgString called with %s\n", nameIn.c_str());
  // Separate out index arg if present
  size_t idx_pos = dsname.find( ':' );
  if ( idx_pos != std::string::npos ) {
    // Advance to after the ':'
    idx_arg = dsname.substr( idx_pos + 1 );
    //mprinterr("DBG:\t\tIndex Arg [%s]\n", idx_arg.c_str());
    // Drop the index arg
    dsname.resize( idx_pos );
  }

  // Separate out aspect arg if present
  size_t attr_pos0 = dsname.find_first_of( '[' );
  size_t attr_pos1 = dsname.find_last_of( ']' );
  if ( attr_pos0 != std::string::npos && attr_pos1 != std::string::npos ) {
    if ( (attr_pos0 != std::string::npos && attr_pos1 == std::string::npos) ||
         (attr_pos0 == std::string::npos && attr_pos1 != std::string::npos) )
    {
      mprinterr("Error: Malformed attribute ([<attr>]) in dataset name %s\n", nameIn.c_str());
      return 0;
    }
    // Advance to after '[', length is position of ']' minus '[' minus 1 
    attr_arg = dsname.substr( attr_pos0 + 1, attr_pos1 - attr_pos0 - 1 );
    //mprinterr("DBG:\t\tAttr Arg [%s]\n", attr_arg.c_str());
    // Drop the attribute arg
    dsname.resize( attr_pos0 );
  }
  //mprinterr("DBG:\t\tName Arg [%s]\n", dsname.c_str());

  return dsname;
}

// DataSetList::PendingWarning()
void DataSetList::PendingWarning() const {
  if (dataSetsPending_)
    mprintf("Warning: Some Actions currently in Action list need to be run in order to create\n"
            "Warning:   data sets. Try processing currently loaded trajectories with 'run' or\n"
            "Warning:   'go' to generate these data sets.\n");
}

// DataSetList::GetDataSet()
DataSet* DataSetList::GetDataSet( std::string const& nameIn ) const {
  std::string attr_arg;
  std::string idx_arg;
  std::string dsname = ParseArgString( nameIn, idx_arg, attr_arg );
  int idx = -1;
  if (!idx_arg.empty()) idx = convertToInteger(idx_arg); // TODO: Set idx_arg to -1
  DataSet* ds = GetSet( dsname, idx, attr_arg );
  if (ds == 0) {
    mprintf("Warning: Data set '%s' not found.\n", nameIn.c_str());
    PendingWarning();
  }
  return ds;
}

/** The set argument must match EXACTLY, so Data will not return Data:1 */
DataSet* DataSetList::CheckForSet( std::string const& nameIn ) const {
  std::string aspect_arg;
  std::string idx_arg("-1");
  std::string dsname = ParseArgString(nameIn, idx_arg, aspect_arg);
  int idx = convertToInteger(idx_arg);
  return CheckForSet(dsname, idx, aspect_arg);
}

// DataSetList::CheckForSet()
DataSet* DataSetList::CheckForSet(std::string const& dsname, int idx,
                                  std::string const& aspect_arg) const
{
  for (DataListType::const_iterator ds = DataList_.begin(); ds != DataList_.end(); ++ds)
    if ( (*ds)->Name() == dsname )
      if ( (*ds)->Aspect() == aspect_arg )
        if ( (*ds)->Idx() == idx )
          return *ds;
  return 0;
}

// DataSetList::GetMultipleSets()
/** \return a list of all DataSets matching the given argument. */
DataSetList DataSetList::SelectSets( std::string const& nameIn ) const {
  DataSetList dsetOut;
  dsetOut.hasCopies_ = true;
  Range idxrange;

  std::string attr_arg;
  std::string idx_arg;
  std::string dsname = ParseArgString( nameIn, idx_arg, attr_arg );
  //mprinterr("DBG: GetMultipleSets \"%s\": Looking for %s[%s]:%s\n",nameIn.c_str(), dsname.c_str(), attr_arg.c_str(), idx_arg.c_str());
  // If index arg is empty make wildcard (-1)
  if (idx_arg.empty() || idx_arg == "*")
    idxrange.SetRange( -1, 0 ); 
  else
    idxrange.SetRange( idx_arg );
  // If attribute arg not set and name is wildcard, make attribute wildcard.
  if (attr_arg.empty() && dsname == "*")
    attr_arg.assign("*");
  // All start selected
  std::vector<char> SelectedSets(DataList_.size(), 'T');
  // First check name
  std::vector<char>::iterator selected = SelectedSets.begin();
  if ( dsname != "*" ) {
    for (DataListType::const_iterator ds = DataList_.begin(); ds != DataList_.end(); ++ds) {
      if ( (*ds)->Name() != dsname ) *selected = 'F';
      ++selected;
    }
  }
  // Second check aspect
  if ( attr_arg != "*" ) {
    selected = SelectedSets.begin();
    for (DataListType::const_iterator ds = DataList_.begin(); ds != DataList_.end(); ++ds) {
      if ( *selected == 'T' && (*ds)->Aspect() != attr_arg ) *selected = 'F';
      ++selected;
    }
  }
  // Last check index
  if ( !idx_arg.empty() && idx_arg != "*" ) {
    selected = SelectedSets.begin();
    for (DataListType::const_iterator ds = DataList_.begin(); ds != DataList_.end(); ++ds) {
      if ( *selected == 'T' && !idxrange.InRange( (*ds)->Idx() ) ) *selected = 'F';
      ++selected;
    }
  }
  // Add selected DataSets to dsetOut
  selected = SelectedSets.begin();
  for (DataListType::const_iterator ds = DataList_.begin(); ds != DataList_.end(); ++ds)
    if ( *(selected++) == 'T' ) dsetOut.DataList_.push_back( *ds );
  return dsetOut;
}

// DataSetList::GetMultipleSets()
DataSetList DataSetList::GetMultipleSets( std::string const& nameIn ) const {
  DataSetList dsetOut = SelectSets(nameIn);
  if ( dsetOut.empty() ) {
    mprintf("Warning: '%s' selects no data sets.\n", nameIn.c_str());
    PendingWarning();
  }
  return dsetOut;
}

// DataSetList::GetSet()
/** \return Specified Dataset or null if not found.
  */
DataSet* DataSetList::GetSet(std::string const& dsname, int idx, std::string const& aspect) const 
{
  for (DataListType::const_iterator ds = DataList_.begin(); ds != DataList_.end(); ++ds) 
    if ( (*ds)->Matches( dsname, idx, aspect ) ) return *ds;
  return 0;
}

/** Create or append string data set from given array. */
DataSet* DataSetList::AddOrAppendSet(std::string const& nameIn, int idxIn, std::string const& aspectIn,
                                     std::vector<std::string> const& Svals)
{
  if (Svals.empty()) {
    mprinterr("Internal Error: AddOrAppendSet() called with empty array.\n");
    return 0;
  }
  DataSet* ds = CheckForSet(nameIn, idxIn, aspectIn);
  if (ds == 0) {
    ds = AddSetIdxAspect(DataSet::STRING, nameIn, idxIn, aspectIn);
    if (ds == 0) {
      mprinterr("Error: Could not allocate STRING data set %s\n", nameIn.c_str());
      return 0;
    }
    DataSet_string& data = static_cast<DataSet_string&>( *ds );
    data = Svals;
  } else {
    if (ds->Type() != DataSet::STRING) {
      mprinterr("Error: Cannot append string values to set %s type %s\n", ds->Legend().c_str(),
                DataArray[ds->Type()].Description);
      return 0;
    }
    ((DataSet_string*)ds)->Append( Svals );
  }
  return ds;
}

/** Create or append to numerical data set from the given arrays. */
DataSet* DataSetList::AddOrAppendSet(std::string const& nameIn, int idxIn, std::string const& aspectIn,
                                     std::vector<double> const& Xvals, std::vector<double> const& Yvals)
{
  if (Xvals.empty() || Yvals.empty()) {
    mprinterr("Internal Error: AddOrAppendSet() called with empty arrays.\n");
    return 0;
  }
  if (Xvals.size() != Yvals.size()) {
    mprinterr("Internal Error: AddOrAppendSet() called with different size arrays.\n");
    return 0;
  }
  // First determine if X values increase monotonically with a regular step
  DataSet::DataType setType = DataSet::DOUBLE;
  double xstep = 0.0;
  if (Xvals.size() > 1) {
    xstep = Xvals[1] - Xvals[0];
    for (std::vector<double>::const_iterator X = Xvals.begin()+2; X != Xvals.end(); ++X)
      if ((*X - *(X-1)) - xstep > Constants::SMALL) {
        setType = DataSet::XYMESH;
        break;
      }
    //mprintf("DBG: xstep %g format %i\n", xstep, (int)setType);
  }
  // Determine if we are appending or creating a new set.
  DataSet* ds = CheckForSet(nameIn, idxIn, aspectIn);
  if (ds == 0) {
    // Create
    ds = AddSetIdxAspect(setType, nameIn, idxIn, aspectIn);
    if (ds == 0) {
      mprinterr("Error: Could not add set '%s[%s]:%i\n", nameIn.c_str(), aspectIn.c_str(), idxIn);
      return 0;
    }
    if (setType == DataSet::DOUBLE) {
      DataSet_double& data = static_cast<DataSet_double&>( *ds );
      data = Yvals;
      data.SetDim(Dimension::X, Dimension(Xvals.front(), xstep, Xvals.size()));
    } else { // XYMESH
      DataSet_Mesh& data = static_cast<DataSet_Mesh&>( *ds );
      data.SetMeshXY( Xvals, Yvals );
    }
  } else {
    // Append
    if (ds->Type() == DataSet::DOUBLE) {
      if (setType == DataSet::XYMESH) {
        mprinterr("Error: Can not append to double data set when x values are irregular.\n");
        return 0;
      }
      if (xstep != ds->Dim(0).Step()) {
        mprinterr("Error: Can not append to set '%s', X step %g does not match current"
                  " data X step %g\n", ds->Legend().c_str(), ds->Dim(0).Step(), xstep);
        return 0;
      }
      DataSet_double& data = static_cast<DataSet_double&>( *ds );
      data.Append( Yvals );
    } else if (ds->Type() == DataSet::XYMESH) {
      // Doesnt matter if the step is regular or not.
      DataSet_Mesh& data = static_cast<DataSet_Mesh&>( *ds );
      data.Append( Xvals, Yvals );
    } else {
      mprinterr("Error: Can only append to double or mesh data sets.\n");
      return 0;
    }
  }
  return ds;
}

// DataSetList::AddSet()
/** Add a DataSet with given name, or if no name given create a name based on 
  * defaultName and DataSet position.
  */
DataSet* DataSetList::AddSet( DataSet::DataType inType, std::string const& nameIn,
                              const char* defaultName )
{
  if (nameIn.empty()) {
    if (defaultName == 0)
      return AddSetIdxAspect( inType, std::string(), -1, std::string() );
    else
      return AddSetIdxAspect( inType, GenerateDefaultName(defaultName), -1, std::string() ); 
  } else
    return AddSetIdxAspect( inType, nameIn, -1, std::string() );
}

// DataSetList::GenerateDefaultName()
/** Create a name based on the given defaultName and # of DataSets,
  * i.e. defaultName_XXXXX 
  */
std::string DataSetList::GenerateDefaultName(std::string const& defaultName) const {
  // Determine # chars needed to hold text version of set number (min 5).
  size_t extsize = (size_t) DigitWidth( size() );
  if (extsize < 5) extsize = 5;
  if (defaultName.empty())
    return ( "D" + integerToString(size(), extsize) );
  else
    return ( defaultName + "_" + integerToString(size(), extsize) ); 
}

// DataSetList::AddSetIdx()
/** Add DataSet of specified type with given name and index to list. */
DataSet* DataSetList::AddSetIdx(DataSet::DataType inType,
                                std::string const& nameIn, int idxIn)
{
  return AddSetIdxAspect( inType, nameIn, idxIn, std::string() );
}

// DataSetList::AddSetAspect()
/** Add DataSet of specified type with given name and aspect to list. */
DataSet* DataSetList::AddSetAspect(DataSet::DataType inType,
                                   std::string const& nameIn,
                                   std::string const& aspectIn)
{
  return AddSetIdxAspect( inType, nameIn, -1, aspectIn );
}

// DataSetList::AddSetIdxAspect()
DataSet* DataSetList::AddSetIdxAspect(DataSet::DataType inType,
                                      std::string const& nameIn,
                                      int idxIn, std::string const& aspectIn,
                                      std::string const& legendIn)
{
  DataSet* ds = AddSetIdxAspect( inType, nameIn, idxIn, aspectIn );
  if (ds != 0)
    ds->SetLegend( legendIn );
  return ds;
}

// DataSetList::AddSetIdxAspect()
/** Add a DataSet of specified type, set it up and return pointer to it. 
  * \param inType type of DataSet to add.
  * \param nameIn DataSet name.
  * \param idxIn DataSet index, -1 if not specified.
  * \param aspectIn DataSet aspect, empty if not specified.
  * \param MAXin Size to set dataset to; DataSet will be set to maxFrames if < 1.
  * \return pointer to successfully set-up dataset.
  */ 
DataSet* DataSetList::AddSetIdxAspect(DataSet::DataType inType, 
                                     std::string const& nameIn, int idxIn,
                                     std::string const& aspectIn) 
{
  // Do not add to a list with copies
  if (hasCopies_) {
    mprinterr("Internal Error: Adding DataSet %s copy to invalid list.\n", nameIn.c_str());
    return 0;
  }

  // Check if DataSet with same attributes already present.
  DataSet* DS = CheckForSet(nameIn, idxIn, aspectIn);
  if (DS != 0) {
    mprintf("Warning: DataSet '");
    DS->PrintName();
    mprintf("' already present.\n");
    // NOTE: Should return found dataset?
    return 0; 
  }
  TokenPtr token = &(DataArray[inType]);
  if ( token->Alloc == 0) {
    mprinterr("Internal Error: No allocator for DataSet type [%s]\n", token->Description);
    return 0;
  }
  DS = (DataSet*)token->Alloc();
  if (DS==0) {
    mprinterr("Internal Error: DataSet %s memory allocation failed.\n", nameIn.c_str());
    return 0;
  }

  // Set up dataset 
  if ( DS->SetupSet(nameIn, idxIn, aspectIn) ) {
    mprinterr("Error setting up data set %s.\n",nameIn.c_str());
    delete DS;
    return 0;
  }

  DataList_.push_back(DS); 
  //fprintf(stderr,"ADDED dataset %s\n",dsetName);
  return DS;
}

int DataSetList::AddSet( DataSet* dsIn ) {
  if (dsIn == 0 || dsIn->Name().empty()) return 1;
  DataSet* ds = CheckForSet( dsIn->Name(), dsIn->Idx(), dsIn->Aspect() );
  if (ds != 0) {
    mprintf("Warning: DataSet '");
    ds->PrintName();
    mprintf("' already present.\n");
    return 1;
  }
  DataList_.push_back( dsIn );
  return 0;
}

// DataSetList::AddCopyOfSet()
void DataSetList::AddCopyOfSet(DataSet* dsetIn) {
  if (!hasCopies_ && !DataList_.empty()) {
    mprinterr("Internal Error: Adding DataSet (%s) copy to invalid list\n", 
    dsetIn->Legend().c_str());
    return;
  }
  hasCopies_ = true;
  DataList_.push_back( dsetIn );
}

// DataSetList::List()
/** Print information on all data sets in the list, as well as any datafiles
  * that will be written to.
  */
void DataSetList::List() const {
  if (!hasCopies_) { // No copies; this is a Master DSL.
    if (DataList_.empty()) return;
    mprintf("\nDATASETS:\n");
  } else if (DataList_.empty()) {
    mprintf("  No data sets.");
    return;
  }
  if (DataList_.size()==1)
    mprintf("  1 data set:\n");
  else
    mprintf("  %zu data sets:\n", DataList_.size());
  for (unsigned int ds=0; ds<DataList_.size(); ds++) {
    DataSet const& dset = static_cast<DataSet const&>(*DataList_[ds]);
    mprintf("\t");
    dset.PrintName();
    mprintf(" \"%s\"", dset.Legend().c_str());
    mprintf(" (%s", DataArray[dset.Type()].Description);
    dset.ScalarDescription();
    mprintf("), size is %i", dset.Size());
    dset.Info();
    mprintf("\n");
  }
}
#ifdef MPI
// DataSetList::SynchronizeData()
void DataSetList::SynchronizeData() {
  // Sync datasets - does nothing if worldsize is 1
  for (DataListType::iterator ds = DataList_.begin(); ds != DataList_.end(); ++ds) {
    if ( (*ds)->Sync() ) {
      rprintf( "Error syncing dataset %s\n",(*ds)->Legend().c_str());
      //return;
    }
  }
}
#endif
// DataSetList::FindSetOfType()
DataSet* DataSetList::FindSetOfType(std::string const& nameIn, DataSet::DataType typeIn) const
{
  std::string attr_arg;
  std::string idx_arg;
  std::string dsname = ParseArgString( nameIn, idx_arg, attr_arg );
  int idx = -1;
  if (!idx_arg.empty()) idx = convertToInteger(idx_arg);
  for (DataListType::const_iterator ds = DataList_.begin(); ds != DataList_.end(); ++ds) {
    if ( (*ds)->Type() == typeIn ) {
      if ( (*ds)->Matches( dsname, idx, attr_arg ) )
        return (*ds);
    }
  }
  return 0;
}

/** Search for a COORDS DataSet. If no name specified, create a default 
  * COORDS DataSet named _DEFAULTCRD_.
  */
DataSet* DataSetList::FindCoordsSet(std::string const& setname) {
  DataSet* outset = 0;
  if (setname.empty()) {
    // crdset not given, search for the default set
    outset = FindSetOfType("_DEFAULTCRD_", DataSet::COORDS);
    if (outset == 0) {
      // No default set; create one.
      outset = AddSet(DataSet::COORDS, "_DEFAULTCRD_", "CRD");
    }
  } else {
    // crdset specified
    outset = FindSetOfType(setname, DataSet::COORDS);
    // If COORDS not found look for TRAJ
    if (outset == 0)
      outset = FindSetOfType(setname, DataSet::TRAJ);
    // If TRAJ not found, look for REF_FRAME
    if (outset == 0)
      outset = FindSetOfType(setname, DataSet::REF_FRAME);
  }
  return outset;
}

const char* DataSetList::RefArgs = "reference | ref <name> | refindex <#>";

/** Search for a REF_FRAME DataSet by file name/tag. Provided for backwards
  * compatibility with the FrameList::GetFrameByName() routine.
  */
DataSet* DataSetList::GetReferenceFrame(std::string const& refname) const {
  DataSet* ref = 0;
  if (refname[0] == '[') {
    // This is a tag. Handle here since brackets normally reserved for Aspect.
    // FIXME: Will not work if index arg is also provided.
    for (DataListType::const_iterator ds = DataList_.begin(); ds != DataList_.end(); ++ds) {
      if ( (*ds)->Type() == DataSet::REF_FRAME && (*ds)->Name() == refname ) {
        ref = *ds;
        break;
      }
    }
  } else {
    ref = CheckForSet( refname );
    if (ref == 0) {
      // If ref not found, check if base file name was specified instead of
      // full, which is by default how DataSet_Coords_REF chooses name.
      for (DataListType::const_iterator ds = DataList_.begin(); ds != DataList_.end(); ++ds) {
        if ( (*ds)->Type() == DataSet::REF_FRAME) {
          DataSet_Coords_REF const& R = static_cast<DataSet_Coords_REF const&>(*(*ds));
          if (refname == R.FrameName().Base()) {
            ref = *ds;
            break;
          }
        }
      }
    }
    if (ref != 0 && ref->Type() != DataSet::REF_FRAME) {
      mprinterr("Error: Data set '%s' is not a reference frame.\n", refname.c_str());
      ref = 0;
    }
  }
  return ref;
}

/** Search for a REF_FRAME DataSet. Provided for backwards compatibility
  * with the FrameList::GetFrameFromArgs() routine.
  * The keywords in order of precedence are:
  *   - 'ref <name>'  : Get reference frame by full/base filename or tag.
  *   - 'reference'   : First reference frame in list.
  *   - 'refindex <#>': Reference frame at position.
  */
ReferenceFrame DataSetList::GetReferenceFrame(ArgList& argIn) const {
  DataSet* ref = 0;
  // 'ref <name>'
  std::string refname = argIn.GetStringKey("ref");
  if (!refname.empty()) {
    ref = GetReferenceFrame( refname );
    if (ref == 0) {
      mprinterr("Error: Reference '%s' not found.\n", refname.c_str());
      return ReferenceFrame(1);
    } 
  } else {
    int refindex = argIn.getKeyInt("refindex", -1);
    if (argIn.hasKey("reference")) refindex = 0;
    if (refindex > -1) {
      for (DataListType::const_iterator ds = DataList_.begin(); ds != DataList_.end(); ++ds) {
        if ( (*ds)->Type() == DataSet::REF_FRAME) {
          DataSet_Coords_REF const& R = static_cast<DataSet_Coords_REF const&>(*(*ds));
          if (R.RefIndex() == refindex) {
            ref = *ds;
            break;
          }
        }
      }
      if (ref == 0) {
        mprinterr("Error: Reference index %i not found.\n", refindex);
        return ReferenceFrame(1);
      }
    }
  }
  return ReferenceFrame((DataSet_Coords_REF*)ref);
}

// DataSetList::ListReferenceFrames()
void DataSetList::ListReferenceFrames() const {
  // Go through DataSetList, count reference frames, put in temp list.
  std::vector<DataSet_Coords_REF*> refTemp;
  for (DataListType::const_iterator ds = DataList_.begin(); ds != DataList_.end(); ++ds)
    if ( (*ds)->Type() == DataSet::REF_FRAME)
      refTemp.push_back( (DataSet_Coords_REF*)*ds );
  if (!refTemp.empty()) {
    mprintf("\nREFERENCE FRAMES (%zu total):\n", refTemp.size());
    for (std::vector<DataSet_Coords_REF*>::const_iterator ds = refTemp.begin();
                                                          ds != refTemp.end(); ++ds)
      if (!(*ds)->FrameName().empty())
        mprintf("    %i: '%s', frame %i\n", (*ds)->RefIndex(), (*ds)->FrameName().full(),
                (*ds)->Idx());
      else
        mprintf("    (DataSet) '%s'\n", (*ds)->Name().c_str());
  }
}
