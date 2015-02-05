#include <cstdio> // sscanf
#include <algorithm> // sort
#include "DataIO_RemLog.h"
#include "CpptrajStdio.h" 
#include "ProgressBar.h"
#include "StringRoutines.h" // fileExists

// CONSTRUCTOR
DataIO_RemLog::DataIO_RemLog() : n_mremd_replicas_(0), processMREMD_(false),
  searchForLogs_(true)
{
  SetValid( DataSet::REMLOG );
}

// NOTE: Must match ExchgType
const char* DataIO_RemLog::ExchgDescription[] = {
  "Unknown", "Temperature", "Hamiltonian", "MultipleDim", "RXSGLD"
};

/// \return true if char pointer is null.
static inline bool IsNullPtr( const char* ptr ) { 
  if (ptr == 0) {
    mprinterr("Error: Could not read file.\n");
    return true;
  }
  return false;
}

// DataIO_RemLog::ReadRemlogHeader()
int DataIO_RemLog::ReadRemlogHeader(BufferedLine& buffer, ExchgType& type) const {
  int numexchg = -1;
  // Read the first line. Should be '# Replica Exchange log file'
  std::string line = buffer.GetLine();
  if (line.compare(0, 27, "# Replica Exchange log file") != 0) {
    mprinterr("Error: Expected '# Replica Exchange log file', got:\n%s\n", line.c_str());
    return -1;
  }

  // Read past metadata. Save expected number of exchanges.
  while (line[0] == '#') {
    line = buffer.GetLine();
    if (line.empty()) {
      mprinterr("Error: No exchanges in rem log.\n");
      return -1;
    }
    ArgList columns( line );
    // Each line should have at least 2 arguments
    if (columns.Nargs() > 1) {
      if (columns[1] == "exchange")
        break;
      else
      {
        if (debug_ > 0) mprintf("\t%s", line.c_str());
        if (columns[1] == "numexchg")
          numexchg = columns.getNextInteger(-1);
        else if (columns[1] == "Dimension") {
          type = MREMD;
          int ndim = columns.getKeyInt("of", 0);
          if (ndim != (int)GroupDims_.size()) {
            mprinterr("Error: # of dimensions in rem log %i != dimensions in remd.dim file (%u).\n",
                      ndim, GroupDims_.size());
            return -1;
          }
        } else if (columns[1] == "RXSGLD")
          type = RXSGLD;
      }
    }
    if (type == UNKNOWN && columns.hasKey("Rep#,")) {
      if (columns[2] == "Neibr#,") type = HREMD;
      else if (columns[2] == "Velocity") type = TREMD;
    }
  }
  if (numexchg < 1) {
    mprinterr("Error: Invalid number of exchanges (%i) in rem log.\n");
    return -1;
  }
  return numexchg;
}

// DataIO_RemLog::ReadRemdDimFile()
// TODO: Handle cases where groups are not in order.
int DataIO_RemLog::ReadRemdDimFile(std::string const& rd_name) {
  typedef std::map<int,GroupArray> GroupMapType;
  typedef std::pair<GroupMapType::iterator,bool> GroupMapRet;
  typedef std::pair<int,GroupArray> GroupMapElt;
  BufferedLine rd_file;
  if (rd_file.OpenFileRead( rd_name )) {
    mprinterr("Error: Could not read remd dim file '%s'\n", rd_name.c_str());
    return 1;
  }
  const char* separators = " =,()";
  // Read dimension title
  const char* ptr = rd_file.Line();
  if (IsNullPtr( ptr )) return 1;
  // ptr Should end with a newline
  mprintf("\tReplica dimension file '%s' title: %s", rd_name.c_str(), ptr);
  // Read each &multirem section
  GroupDims_.clear();
  DimTypes_.clear();
  ArgList rd_arg;
  while (ptr != 0) {
    rd_arg.SetList( std::string(ptr), separators );
    if ( rd_arg[0] == "&multirem" ) {
      GroupMapType GroupMap;
      std::string desc;
      int n_replicas = 0;
      ExchgType exch_type = UNKNOWN;
      while (ptr != 0) {
        rd_arg.SetList( std::string(ptr), separators );
        if (rd_arg.CommandIs("&end") || rd_arg.CommandIs("/")) break;
        rd_arg.MarkArg(0);
        if ( rd_arg.CommandIs("exch_type") ) {
          if ( rd_arg.hasKey("TEMP") || rd_arg.hasKey("TEMPERATURE") )
            exch_type = TREMD;
          else if ( rd_arg.hasKey("HAMILTONIAN") || rd_arg.hasKey("HREMD") )
            exch_type = HREMD;
          else {
            mprinterr("Error: Unrecognized exch_type: %s\n", rd_arg.ArgLine());
            return 1;
          }
        } else if ( rd_arg.CommandIs("group") ) {
          int group_num = rd_arg.getNextInteger(-1);
          if (group_num < 1) {
            mprinterr("Error: Invalid group number: %i\n", group_num);
            return 1;
          }
          //mprintf("\t\tGroup %i\n", group_num);
          std::vector<int> indices;
          int group_index = rd_arg.getNextInteger(-1);
          while (group_index != -1) {
            indices.push_back( group_index );
            n_replicas++;
            group_index = rd_arg.getNextInteger(-1);
          }
          // Set up partner array for this group
          GroupArray group;
          for (int i = 0; i < (int)indices.size(); i++) {
            int l_idx = i - 1;
            if (l_idx < 0) l_idx = (int)indices.size() - 1;
            int r_idx = i + 1;
            if (r_idx == (int)indices.size()) r_idx = 0;
            group.push_back( GroupReplica(indices[l_idx], indices[i], indices[r_idx]) );
            //mprintf("\t\t\t%i - %i - %i\n", group.back().L_partner(),
            //        group.back().Me(), group.back().R_partner());
          }
          GroupMapRet ret = GroupMap.insert( GroupMapElt(group_num, group) );
          if (ret.second == false) {
            mprinterr("Error: Duplicate group # detected (%i)\n", group_num);
            return 1;
          }
        } else if ( rd_arg.CommandIs("desc") ) {
          desc = rd_arg.GetStringNext(); 
        }
        ptr = rd_file.Line();
      }
      // Place sorted groups for dim into GroupDimType
      GroupDimType Groups;
      for (GroupMapType::const_iterator it = GroupMap.begin(); it != GroupMap.end(); ++it)
        Groups.push_back( it->second );
      mprintf("\tDimension %zu: type '%s', description '%s', groups=%zu, replicas=%i\n", 
              GroupDims_.size() + 1, ExchgDescription[exch_type], desc.c_str(), 
              Groups.size(), n_replicas);
      if (n_mremd_replicas_ == 0)
        n_mremd_replicas_ = n_replicas;
      else if (n_replicas != n_mremd_replicas_) {
        mprinterr("Error: Number of MREMD replicas in dimension (%i) != number of\n"
                  "Error: MREMD replicas in first dimension (%i)\n", n_replicas, n_mremd_replicas_);
        return 1;
      }
      GroupDims_.push_back( Groups );
      DimTypes_.push_back( exch_type );
    }
    ptr = rd_file.Line();
  }
  if (GroupDims_.empty()) {
    mprinterr("Error: No replica dimensions found.\n");
    return 1;
  }
  return 0;
}

/// Get filename up to extension
//TODO May not need to be its own function. Make a general FileName function?
static inline std::string GetPrefix(FileName const& fname) {
  size_t found = fname.Full().rfind( fname.Ext() );
  return fname.Full().substr(0, found);
}

// DataIO_RemLog::SetupTemperatureMap()
/** buffer should be positioned at the first exchange. */
DataSet_RemLog::TmapType 
  DataIO_RemLog::SetupTemperatureMap(BufferedLine& buffer,
                                     std::vector<int>& CrdIdxs) const
{
  DataSet_RemLog::TmapType TemperatureMap;
  std::vector<TlogType> tList;
  TlogType tlog;
  CrdIdxs.clear();
  const char* ptr = buffer.Line();
  while (ptr != 0 && ptr[0] != '#') {
    // For temperature remlog create temperature map. 
    //mprintf("DEBUG: Temp0= %s", ptr+32);
    if ( sscanf(ptr, "%2i%*10f%*10f%*10f%10lf", &tlog.crdidx, &tlog.t0) != 2 ) {
      mprinterr("Error: could not read temperature from T-REMD log.\n"
                "Error: Line: %s", ptr);
      return TemperatureMap;
    }
    tList.push_back( tlog );
    ptr = buffer.Line();
  }
  // Sort temperatures
  std::sort( tList.begin(), tList.end(), TlogType_cmp() );
  // Place sorted temperatures into map starting from replica index 1. Check
  // for duplicate temperatures. Also store the sorted coordinate indices.
  int repidx = 1;
  for (std::vector<TlogType>::const_iterator it = tList.begin();
                                             it != tList.end(); ++it)
  {
    mprintf("\t\tReplica %i => %f (crdidx= %i)\n", repidx, it->t0, it->crdidx); 
    if (it != tList.begin()) {
      if ( it->t0 == (it-1)->t0 ) {
        mprinterr("Error: duplicate temperature %.2f detected in T-REMD remlog\n", it->t0);
        TemperatureMap.clear();
        return TemperatureMap;
      }
    }
    TemperatureMap.insert(std::pair<double,int>(it->t0, repidx++));
    CrdIdxs.push_back( it->crdidx );
  }

  return TemperatureMap;
}

// DataIO_RemLog::CountHamiltonianReps()
/** buffer should be positioned at the first exchange. */
int DataIO_RemLog::CountHamiltonianReps(BufferedLine& buffer) const {
  int n_replicas = 0;
  const char* ptr = buffer.Line();
  while (ptr != 0 && ptr[0] != '#') {
    ptr = buffer.Line();
    ++n_replicas;
  }
  return n_replicas;
}

static inline int MremdNrepsError(int n_replicas, int dim, int groupsize) {
  if (n_replicas != groupsize) {
    mprinterr("Error: Number of replicas in dimension %i (%i) does not match\n"
              "Error: number of replicas in remd.dim file (%u)\n", dim+1, n_replicas,
              groupsize);
    return 1;
  }
  return 0;
}

/** Open all dimension files associated with given replica log. Check that the
  * number of exchanges reported in the files match each other and that each
  * is actually a MREMD log. Will leave all remlogs positioned at the first
  * exchange.
  * \return Expected number of exchanges
  */
int DataIO_RemLog::OpenMremdDims(std::vector<BufferedLine>& buffer, 
                                 Sarray const& dimLogs)
{
  // Sanity check
  if (buffer.size() != dimLogs.size()) {
    mprinterr("Internal Error: File buffer array size %zu != # MREMD logs %zu.\n",
              buffer.size(), dimLogs.size());
    return 1;
  }
  int total_exchanges = -1;
  ExchgType log_type = UNKNOWN;
  // Open remlogs for each dimension as buffered file.
  for (unsigned int dim = 0; dim < GroupDims_.size(); dim++) {
    buffer[dim].CloseFile();
    if (buffer[dim].OpenFileRead( dimLogs[dim]  )) return -1;
    //mprintf("\tOpened %s\n", logname.c_str());
    // Read the remlog header.
    int numexchg = ReadRemlogHeader(buffer[dim], log_type);
    if (numexchg == -1) return -1;
    //mprintf("\t%s should contain %i exchanges\n", dimLogs[dim].c_str(), numexchg);
    if (total_exchanges == -1)
      total_exchanges = numexchg;
    else if (numexchg != total_exchanges) {
      mprinterr("Error: Number of expected exchanges in dimension %i does not match\n"
                "Error: number of expected exchanges in first dimension.\n", dim + 1);
      return -1;
    }
    if ( processMREMD_ ) {
      if (log_type != MREMD) {
        mprinterr("Error: Log type is not MREMD.\n");
        return -1;
      }
    } else {
      // Single dimension. Record log type, ensure subsqeuent opens match.
      if (DimTypes_.empty())
        DimTypes_.push_back( log_type );
      else if ( log_type != DimTypes_.front() ) {
        mprinterr("Error: Log type of %s does not match first log type.\n", 
                  buffer[dim].Filename().full());
        return -1;
      }
    }
  }
  return total_exchanges;
}

/** For non-MREMD runs, setup dimension 1 group. */
void DataIO_RemLog::SetupDim1Group( int group_size ) {
  if (GroupDims_.empty()) GroupDims_.resize( 1 );
  GroupDims_[0].resize( 1 );
  for (int replica = 0; replica < group_size; replica++) {
    int me = replica + 1;
    int l_partner = me - 1;
    if (l_partner < 1) l_partner = group_size;
    int r_partner = me + 1;
    if (r_partner > group_size) r_partner = 1;
    GroupDims_[0][0].push_back( GroupReplica(l_partner, me, r_partner) );
  }
  n_mremd_replicas_ = group_size;
}  

// DataIO_RemLog::ReadHelp()
void DataIO_RemLog::ReadHelp() {
  mprintf("\tnosearch            : Do not automatically search for MREMD dimension logs.\n"
          "\tdimfile <file>      : remd.dim file for processing MREMD logs.\n"
          "\tcrdidx <crd indices>: Use comma-separated list of indices as the initial\n"
          "\t                      coordinate indices.\n"
          "\tMultiple REM logs may be specified.\n");
}

int DataIO_RemLog::processReadArgs(ArgList& argIn) {
  searchForLogs_ = !argIn.hasKey("nosearch");
  // Get dimfile arg
  dimfile_ = argIn.GetStringKey("dimfile");
  crdidx_ = argIn.GetStringKey("crdidx");
  // Unlike other data sets, remlog will find all file names now
  // in case of MREMD. Always add at least one entry to this array
  // for the first file.
  logFilenames_.push_back( std::string("") );
  std::string log_name = argIn.GetStringNext();
  while (!log_name.empty()) {
    if (!fileExists( log_name ))
      mprintf("Warning: '%s' does not exist.\n", log_name.c_str());
    else
      logFilenames_.push_back( log_name );
    log_name = argIn.GetStringNext();
  }
  return 0;
}

// DataIO_RemLog::ReadData()
int DataIO_RemLog::ReadData(std::string const& fname, 
                            DataSetList& datasetlist, std::string const& dsname)
{
  if (!fileExists( fname )) {
    mprinterr("Error: File '%s' does not exist.\n", fname.c_str());
    return 1;
  }
  if (logFilenames_.empty()) // processReadArgs not called
    logFilenames_.push_back( fname );
  else
    logFilenames_[0] = fname;
  if (!dimfile_.empty()) {
    if (ReadRemdDimFile( dimfile_ )) {
      mprinterr("Error: Reading remd.dim file '%s'\n", dimfile_.c_str());
      return 1;
    }
    mprintf("\tExpecting %zu replica dimensions.\n", GroupDims_.size());
  }
  // Get crdidx arg
  ArgList idxArgs( crdidx_, "," );
  mprintf("\tReading from log files:");
  for (Sarray::const_iterator it = logFilenames_.begin(); it != logFilenames_.end(); ++it)
    mprintf(" %s", it->c_str());
  mprintf("\n");
  // ---------------------------------------------
  typedef std::vector<Sarray> LogGroupType;
  LogGroupType logFileGroups;
  if (GroupDims_.empty()) {
    // Not M-REMD; T-REMD or H-REMD. Single dimension.
    processMREMD_ = false;
    GroupDims_.resize( 1 );
    for (Sarray::const_iterator logfile = logFilenames_.begin();
                                logfile != logFilenames_.end(); ++logfile)
      logFileGroups.push_back( Sarray( 1, *logfile ) );
  } else {
    // M-REMD
    processMREMD_ = true;
    // Ensure that each replica log has counterparts for each dimension
    // TODO: Read all headers and check dimensions in log.
    // Two cases; all log files were specified, or only lowest logs were specified.
    if ( searchForLogs_ ) { 
      FileName fname;
      for (Sarray::const_iterator logfile = logFilenames_.begin();
                                  logfile != logFilenames_.end(); ++logfile)
      {
        Sarray dimLogs;
        fname.SetFileName( *logfile );
        // Remove leading '.'
        std::string logExt = fname.Ext();
        if (logExt[0] == '.') logExt.erase(0,1);
        if ( !validInteger(logExt) ) {
          mprinterr("Error: MREMD log %s does not have valid numerical extension.\n", fname.full());
          return 1;
        }
        std::string Prefix = GetPrefix( fname );
        int numericalExt = convertToInteger( logExt );
        if (numericalExt != 1) {
          mprinterr("Error: Must specify MREMD log for dimension 1 (i.e. '%s.1')\n", 
                    Prefix.c_str());
          return 1;
        }
        dimLogs.push_back( *logfile );
        for (int idim = 2; idim <= (int)GroupDims_.size(); idim++) {
          std::string logname = Prefix + "." + integerToString( idim );
          if ( !fileExists(logname) ) {
            mprinterr("Error: MREMD log not found for dimension %i, '%s'\n",
                      idim, logname.c_str());
            return 1;
          }
          dimLogs.push_back( logname );
        }
        logFileGroups.push_back( dimLogs );
      }
    } else {
      // All logs specified. Assume they are given in order.
      Sarray dimLogs;
      Sarray::const_iterator logfile = logFilenames_.begin();
      while (logfile != logFilenames_.end()) {
        dimLogs.clear();
        for (unsigned int dim = 0; dim < GroupDims_.size(); dim++) {
          if (logfile == logFilenames_.end()) {
            mprinterr("Error: Ran out of MREMD logs, run %zu, dimension %u\n",
                      logFileGroups.size() + 1, dim + 1);
            return 1;
          }
          dimLogs.push_back( *(logfile++) );
        }
        logFileGroups.push_back( dimLogs );
      }
    }
  }    
  // Set up temperature maps/coordinate index arrays for each dim/group.
  // Base this on the first set of MREMD replica logs.
  // Open remlogs for each dimension as buffered file.
  std::vector<BufferedLine> buffer( GroupDims_.size() );
  int total_exchanges = OpenMremdDims(buffer, logFileGroups.front());
  if (total_exchanges == -1) return 1;
  // Should now be positioned at the first exchange in each dimension.
  // Set up map/coordinate indices for each group and make sure they match
  // whats in the remd.dim file.
  // Temperature map for dimensions (if needed) 
  std::vector<DataSet_RemLog::TmapType> TemperatureMap( GroupDims_.size() );
  // Coordinate indices for temperature dimensions (if needed)
  std::vector< std::vector<int> > TempCrdIdxs( GroupDims_.size() ); 
  for (int dim = 0; dim < (int)GroupDims_.size(); dim++) {
    int group_size = 0;
    if ( DimTypes_[dim] == TREMD ) {
      TemperatureMap[dim] = SetupTemperatureMap( buffer[dim], TempCrdIdxs[dim] );
      if (TemperatureMap[dim].empty()) return 1;
      group_size = (int)TemperatureMap[dim].size();
      mprintf("\t\tDim %i: %i Temperature reps.\n", dim+1, group_size);
      if (!processMREMD_) SetupDim1Group( group_size );
      for (unsigned int grp = 0; grp < GroupDims_[dim].size(); grp++) {
        if (MremdNrepsError(group_size, dim, GroupDims_[dim][grp].size())) return 1;
      }
    } else if (DimTypes_[dim] == HREMD) {
      group_size = CountHamiltonianReps( buffer[dim] );
      mprintf("\t\tDim %i: %i Hamiltonian reps.\n", dim+1, group_size);
      if (!processMREMD_) SetupDim1Group( group_size );
      for (unsigned int grp = 0; grp < GroupDims_[dim].size(); grp++) {
        if (MremdNrepsError(group_size, dim, GroupDims_[dim][grp].size())) return 1;
      } 
    } else if (DimTypes_[dim] == RXSGLD) {
      group_size = CountHamiltonianReps( buffer[dim] );
      mprintf("\t\tDim %i: %i RXSGLD reps.\n", dim+1, group_size);
      if (!processMREMD_) SetupDim1Group( group_size );
      for (unsigned int grp = 0; grp < GroupDims_[dim].size(); grp++) {
        if (MremdNrepsError(group_size, dim, GroupDims_[dim][grp].size())) return 1;
      }
    } else {
      mprinterr("Error: Unrecognized dimension type.\n");
      return 1;
    }
  } // END loop over replica dimensions
  // DEBUG: Print out dimension layout
  if (debug_ > 0) {
    for (std::vector<GroupDimType>::const_iterator Dim = GroupDims_.begin();
                                                   Dim != GroupDims_.end(); ++Dim)
    {
      mprintf("Dimension %u:\n", Dim - GroupDims_.begin());
      for (GroupDimType::const_iterator Group = Dim->begin();
                                        Group != Dim->end(); ++Group)
      {
        mprintf("\tGroup %u:\n", Group - Dim->begin());
        for (GroupArray::const_iterator Rep = Group->begin();
                                        Rep != Group->end(); ++Rep)
          mprintf("\t\tReplica[%u]= %i\n", Rep - Group->begin(), Rep->Me());
      }
    }
  }
  // Coordinate indices for each replica. Start crdidx = repidx (from 1) for now.
  std::vector<int> CoordinateIndices( n_mremd_replicas_ );
  for (int repidx = 0; repidx < n_mremd_replicas_; repidx++) {
    if (!idxArgs.empty()) {
      // User-specified starting coord indices
      CoordinateIndices[repidx] = idxArgs.getNextInteger(0);
      if (CoordinateIndices[repidx] < 0 || CoordinateIndices[repidx] > n_mremd_replicas_ )
      {
        mprinterr("Error: Given coordinate index out of range or not enough indices given.\n");
        return 1;
      }
    } else if (!processMREMD_ && !TempCrdIdxs.front().empty())
      // Use coordinate indices from 1D T-REMD log
      CoordinateIndices[repidx] = TempCrdIdxs.front()[repidx];
    else
      // Default: starting crdidx = repidx
      CoordinateIndices[repidx] = repidx + 1;
  }
//  if (!idxArgs.empty()) {
    mprintf("\tInitial coordinate indices:");
    for (std::vector<int>::const_iterator c = CoordinateIndices.begin();
                                          c != CoordinateIndices.end(); ++c)
      mprintf(" %i", *c);
    mprintf("\n");
//  }
  // Allocate replica log DataSet
  DataSet* ds = datasetlist.CheckForSet(dsname, -1, "");
  if (ds == 0) {
    // New set
    ds = datasetlist.AddSet( DataSet::REMLOG, dsname, "remlog" );
    if (ds == 0) return 1;
    ((DataSet_RemLog*)ds)->AllocateReplicas(n_mremd_replicas_);
  } else {
    if (ds->Type() != DataSet::REMLOG) {
      mprinterr("Error: Set '%s' is not replica log data.\n", ds->Legend().c_str());
      return 1;
    }
    if ((int)ds->Size() != n_mremd_replicas_) {
      mprinterr("Error: Replica log data '%s' is set up for %zu replicas,"
                " current # replicas is %i\n", ds->Legend().c_str(), ds->Size(),
                n_mremd_replicas_);
      return 1;
    }
  }
  DataSet_RemLog& ensemble = static_cast<DataSet_RemLog&>( *ds );
  // Loop over all remlogs
  for (LogGroupType::const_iterator it = logFileGroups.begin(); it != logFileGroups.end(); ++it)
  { 
    // Open the current remlog, advance to first exchange
    int numexchg = OpenMremdDims(buffer, *it);
    if (numexchg == -1) return 1;
    mprintf("\t%s should contain %i exchanges\n", it->front().c_str(), numexchg);
    // Should now be positioned at 'exchange 1'.
    // Loop over all exchanges.
    ProgressBar progress( numexchg );
    bool fileEOF = false;
    const char* ptr = 0;
    unsigned int current_dim = 0;
    int grp; // Will be set to group number for MREMD or group index otherwise
    for (int exchg = 0; exchg < numexchg; exchg++) {
      progress.Update( exchg );
      // Loop over all groups in the current dimension
      for (unsigned int gidx = 0; gidx < GroupDims_[current_dim].size(); gidx++) {
        if (processMREMD_) {
          if (sscanf(buffer[current_dim].CurrentLine(), "%*s%*s%*i%*s%*s%i", &grp)!=1) {
            mprinterr("Error: Could not get MREMD group number.\n");
            return 1;
          }
          // REMD group indices start from 1
          grp--;
        } else
          grp = gidx;
        // Loop over all replicas in the current group
        //mprintf("--------------------------------------------------------------------------\n");
        for (unsigned int replica = 0; replica < GroupDims_[current_dim][grp].size(); replica++) {
          // Read remlog line.
          ptr = buffer[current_dim].Line();
          if (ptr == 0) {
            mprinterr("Error: reading remlog; unexpected EOF. Dim=%u, Exchg=%i, grp=%u, rep=%u\n",
                      current_dim+1, exchg+1, grp+1, replica+1);
            fileEOF = true;
            // If this is not the first replica remove all partial replicas
            if (replica > 0) ensemble.TrimLastExchange();
            break;
          }
          // ----- T-REMD ----------------------------
          /* Format:
           * '(i2,6f10.2,i8)'
          # Rep#, Velocity Scaling, T, Eptot, Temp0, NewTemp0, Success rate (i,i+1), ResStruct#
            1     -1.00      0.00   -433.24    300.00    300.00      0.00      -1
           * Order during REMD is exchange -> MD, so NewTemp0 is the temp. that gets
           * simulated. TODO: Is that valid?
           */
          if (DimTypes_[current_dim] == TREMD) {
            int tremd_crdidx, current_crdidx; // TODO: Remove tremd_crdidx
            double tremd_scaling, tremd_pe, tremd_temp0, tremd_tempP;
            if ( sscanf(ptr, "%2i%10lf%*10f%10lf%10lf%10lf", &tremd_crdidx, &tremd_scaling,
                        &tremd_pe, &tremd_temp0, &tremd_tempP) != 5 )
            {
              mprinterr("Error reading TREMD line from rem log. Dim=%u, Exchg=%i, grp=%u, rep=%u\n",
                        current_dim+1, exchg+1, grp+1, replica+1);
              mprinterr("Error: Line: %s", ptr);
              return 1;
            }
            // Figure out my position within the group.
            DataSet_RemLog::TmapType::const_iterator tmap = 
              TemperatureMap[current_dim].find( tremd_temp0 );
            if (tmap == TemperatureMap[current_dim].end()) {
              mprinterr("Error: replica temperature %.2f not found in temperature map.\n", 
                        tremd_temp0);
              return 1;
            }
            // What is my actual position? Currently mapped rep nums start from 1
            int tremd_repidx = GroupDims_[current_dim][grp][tmap->second - 1].Me();
            // Who is my partner? ONLY VALID IF EXCHANGE OCCURS
            tmap = TemperatureMap[current_dim].find( tremd_tempP ); // TODO: Make function
            if (tmap == TemperatureMap[current_dim].end()) {
              mprinterr("Error: partner temperature %.2f not found in temperature map.\n",
                        tremd_tempP);
              return 1;
            }
            int tremd_partneridx = GroupDims_[current_dim][grp][tmap->second - 1].Me();
            // Exchange success if velocity scaling is > 0.0
            bool tremd_success = (tremd_scaling > 0.0);
            // If an exchange occured, coordsIdx will be that of partner replica.
            if (tremd_success)
              current_crdidx = CoordinateIndices[tremd_partneridx-1];
            else
              current_crdidx = CoordinateIndices[tremd_repidx-1];
            //mprintf("DEBUG: Exchg %8i Tdim# %2u T=%6.2f group=%2u repidx=%3i partneridx=%3i oldcrdidx=%i newcrdidx=%i\n",
            //        exchg+1, current_dim+1, tremd_temp0, grp+1, tremd_repidx, tremd_partneridx, CoordinateIndices[tremd_repidx-1], current_crdidx);
            // Create replica frame for TREMD
            ensemble.AddRepFrame( tremd_repidx-1,
                                  DataSet_RemLog:: 
                                  ReplicaFrame(tremd_repidx, tremd_partneridx,
                                               current_crdidx,
                                               tremd_success,
                                               tremd_temp0, tremd_pe, 0.0) );
          // ----- H-REMD ----------------------------
          /* Format:
           * '(2i6,5f10.2,4x,a,2x,f10.2)'
       # Rep#, Neibr#, Temp0, PotE(x_1), PotE(x_2), left_fe, right_fe, Success, Success rate (i,i+1)
            1     8    300.00 -25011.03 -24959.58    -27.48      0.00    F        0.00
           */
          } else if (DimTypes_[current_dim] == HREMD) {
            int hremd_grp_repidx, hremd_grp_partneridx, current_crdidx;
            double hremd_temp0, hremd_pe_x1, hremd_pe_x2;
            bool hremd_success;
            if ( sscanf(ptr, "%6i%6i%10lf%10lf%10lf", &hremd_grp_repidx, &hremd_grp_partneridx,
                        &hremd_temp0, &hremd_pe_x1, &hremd_pe_x2) != 5 )
            {
              mprinterr("Error reading HREMD line from rem log. Dim=%u, Exchg=%i, grp=%u, rep=%u\n",
                        current_dim+1, exchg+1, grp+1, replica+1);
              mprinterr("Error: Line: %s", ptr);
              return 1;
            }
            // What is my actual position and who is my actual partner?
            int hremd_repidx = GroupDims_[current_dim][grp][hremd_grp_repidx-1].Me();
            int hremd_partneridx = GroupDims_[current_dim][grp][hremd_grp_partneridx-1].Me();
            //mprintf("DEBUG: Exchg %i Hdim# %u group=%u group_repidx=%i repidx=%i\n",
            //        exchg+1, current_dim+1, grp+1, hremd_grp_repidx, hremd_repidx);
            // Determine if an exchange occurred
            switch ( ptr[66] ) {
              case 'T': hremd_success = true; break;
              case 'F': hremd_success = false; break;
              default: // Should only get here with malformed HREMD log file.
                mprinterr("Error: expected only 'T' or 'F' at character 67, got %c\n", ptr[66]);
                return 1;
            }
            // If an exchange occured, coordsIdx will be that of partner replica.
            if (hremd_success)
              current_crdidx = CoordinateIndices[hremd_partneridx-1];
            else
              current_crdidx = CoordinateIndices[hremd_repidx-1];
            // Create replica frame for HREMD
            ensemble.AddRepFrame( hremd_repidx-1,
                                  DataSet_RemLog::
                                  ReplicaFrame(hremd_repidx, hremd_partneridx,
                                               current_crdidx,
                                               hremd_success,
                                               hremd_temp0, hremd_pe_x1, hremd_pe_x2) );
          // ----- RXSGLD ----------------------------
          /* Format:
           * (i4,i4,2f8.4,2f8.2,e14.6,f8.4) 
           # Rep Stagid Vscale  SGscale Temp Templf Eptot Acceptance(i,i+1)
              1   1  1.0000  1.0000    0.00   21.21 -0.102302E+02  0.0000
           */
          } else if (DimTypes_[current_dim] == RXSGLD) {
            // Consider accept if sgscale is not -1.0.
            int sgld_repidx, sgld_crdidx;
            double sgscale;
            if ( sscanf(ptr, "%4i%4i%*8f%8lf", &sgld_crdidx, &sgld_repidx, &sgscale) != 3 ) {
              mprinterr("Error reading RXSGLD line from rem log. "
                        "Dim=%u, Exchg=%i, grp=%u, rep=%u\n",
                        current_dim+1, exchg+1, grp+1, replica+1);
              mprinterr("Error: Line: %s", ptr);
              return 1;
            }
            bool sgld_success = (sgscale > -1.0);
            // The partner index is not stored in RXSGLD logs.
            bool sgld_up;
            // First exchange for even replicas is up, then down.
            if ((exchg % 2)==0) {
              if ((sgld_repidx % 2)==0) // Exchange up
                sgld_up = true;
              else
                sgld_up = false;
            } else {
              if ((sgld_repidx % 2)==0) // Exchange down
                sgld_up = false;
              else
                sgld_up = true;
            }
            int sgld_partneridx;
            if (sgld_up)
              sgld_partneridx = GroupDims_[current_dim][grp][sgld_repidx-1].R_partner();
            else
              sgld_partneridx = GroupDims_[current_dim][grp][sgld_repidx-1].L_partner();
            // If an exchange occured, coordsIdx will be that of partner replica.
            int current_crdidx;
            if (sgld_success)
              current_crdidx = CoordinateIndices[sgld_partneridx-1];
            else
              current_crdidx = CoordinateIndices[sgld_repidx-1];
            // Create replica frame for SGLD
            ensemble.AddRepFrame( sgld_repidx-1,
                                  DataSet_RemLog::
                                  ReplicaFrame(sgld_repidx, sgld_partneridx,
                                               current_crdidx,
                                               sgld_success,
                                               0.0, 0.0, 0.0) );
          // -----------------------------------------
          } else {
            mprinterr("Error: remlog; unknown type.\n");
            return 1;
          }
          // -----------------------------------------
        } // END loop over replicas in group
        if ( fileEOF ) break; // Error occurred reading replicas, skip rest of groups.
        // Read next group exchange line.
        ptr = buffer[current_dim].Line();
      } // END loop over groups in dimension
      if ( fileEOF ) break; // Error occurred reading replicas, skip rest of exchanges.
      // Update coordinate indices.
      //mprintf("DEBUG: exchange= %i: Updating coordinates\n", exchg + 1);
      for (int repidx = 0; repidx < n_mremd_replicas_; repidx++) {
        //mprintf("DEBUG:\tReplica %i crdidx %i =>", repidx+1, CoordinateIndices[repidx]);
        CoordinateIndices[repidx] = ensemble.LastRepFrame(repidx).CoordsIdx();
        //mprintf(" %i\n", CoordinateIndices[repidx]); // DEBUG
      }
      // Currently each exchange the dimension alternates
      ++current_dim;
      if (current_dim == GroupDims_.size()) current_dim = 0;
    } // END loop over exchanges in remlog
  } // END loop over remlogs
  if (!ensemble.ValidEnsemble()) {
    mprinterr("Error: Ensemble is not valid.\n");
    return 1;
  }
  if (debug_ > 1)
    PrintReplicaStats( ensemble );
  // ---------------------------------------------

  return 0;
}

void DataIO_RemLog::PrintReplicaStats(DataSet_RemLog const& ensemble) {
  mprintf("Replica Stats:\n"
          "%-10s %6s %6s %6s %12s %12s %12s S\n", "#Exchange", "RepIdx", "PrtIdx", "CrdIdx",
          "Temp0", "PE_X1", "PE_X2");
  for (int exchg = 0; exchg < ensemble.NumExchange(); exchg++) {
    for (int replica = 0; replica < (int)ensemble.Size(); replica++) {
      DataSet_RemLog::ReplicaFrame const& frm = ensemble.RepFrame(exchg, replica); 
      mprintf("%10u %6i %6i %6i %12.4f %12.4f %12.4f %1i\n", exchg + 1,
              frm.ReplicaIdx(), frm.PartnerIdx(), frm.CoordsIdx(), frm.Temp0(), 
              frm.PE_X1(), frm.PE_X2(), (int)frm.Success());
    }
  }
}
