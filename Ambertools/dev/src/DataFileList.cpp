// DataFileList
#include "DataFileList.h"
#include "CpptrajStdio.h"
#ifdef MPI
# include "StringRoutines.h" // integerToString
#endif
#ifdef TIMER
# include "Timer.h"
#endif

// CONSTRUCTOR
DataFileList::DataFileList() : 
  debug_(0)
#ifdef MPI
  ,ensembleMode_(-1)
#endif
{}

// DESTRUCTOR
DataFileList::~DataFileList() {
  Clear();
}

// DataFileList::Clear()
void DataFileList::Clear() {
  for (DFarray::iterator it = fileList_.begin(); it != fileList_.end(); it++)
    delete *it;
  fileList_.clear();
}

// DataFileList::RemoveDataFile()
DataFile* DataFileList::RemoveDataFile( DataFile* dfIn ) {
  for (DFarray::iterator it = fileList_.begin(); it != fileList_.end(); ++it) {
    if ( dfIn == *it ) {
      delete *it;
      return (DataFile*)0;
    }
  }
  return dfIn;
}

// DataFileList::RemoveDataSet()
/** Remove given DataSet from any DataFiles in list. */
void DataFileList::RemoveDataSet( DataSet* dsIn ) {
  for (DFarray::iterator df = fileList_.begin(); df != fileList_.end(); ++df)
    (*df)->RemoveSet( dsIn );
}

// DataFileList::SetDebug()
/** Set debug level for DataFileList and all datafiles in it. */
void DataFileList::SetDebug(int debugIn) {
  debug_ = debugIn;
  if (debug_>0)
    mprintf("DataFileList DEBUG LEVEL SET TO %i\n",debug_);
  for (DFarray::iterator df = fileList_.begin(); df != fileList_.end(); ++df)
    (*df)->SetDebug( debug_ );
}

// DataFileList::GetDataFile()
/** Return DataFile specified by given file name if it exists in the list,
  * otherwise return null. Must match full path.
  */
DataFile* DataFileList::GetDataFile(std::string const& nameIn) const {
  if (nameIn.empty()) return 0;
  for (DFarray::const_iterator df = fileList_.begin(); df != fileList_.end(); ++df)
    if (nameIn == (*df)->DataFilename().Full()) return *df;
  return 0;
}

/** Create new DataFile, or return existing DataFile. */
// TODO: Accept const ArgList so arguments are not reset?
DataFile* DataFileList::AddDataFile(std::string const& nameIn, ArgList& argIn) {
  // If no filename, no output desired
  if (nameIn.empty()) return 0;
  std::string name = nameIn;
# ifdef MPI
  if (ensembleMode_ != -1)
    // Ensemble mode, append rank to the output filename.
    name += ("." + integerToString(ensembleMode_));
# endif
  // Check if this filename already in use
  DataFile* Current = GetDataFile(name);
  // If no DataFile associated with name, create new datafile
  if (Current==0) {
    Current = new DataFile();
    if (Current->SetupDatafile(name, argIn, debug_)) {
      mprinterr("Error setting up DataFile %s\n",name.c_str());
      delete Current;
      return 0;
    }
    fileList_.push_back(Current);
  } else {
    // Set debug level
    Current->SetDebug(debug_);
    // Check for keywords that do not match file type
    DataFile::DataFormatType kType = DataFile::GetFormatFromArg( argIn );
    if (kType != DataFile::UNKNOWN_DATA && kType != Current->Type())
      mprintf("Warning: %s is type %s but type %s keyword specified; ignoring keyword.\n",
              Current->DataFilename().full(), Current->FormatString(),
              DataFile::FormatString( kType ));
    // Process Arguments
    if (!argIn.empty())
      Current->ProcessArgs( argIn );
  }
  return Current;
}

// DataFileList::AddDataFile()
DataFile* DataFileList::AddDataFile(std::string const& nameIn) {
  ArgList empty;
  return AddDataFile( nameIn, empty );
}

// DataFileList::AddSetToFile()
/** Add given DataSet to the specified DataFile. If the DataFile does not
  * exist it will be created. Whenever a set is added to a data file
  * reset its writeFile status to true.
  */
DataFile* DataFileList::AddSetToFile(std::string const& nameIn, DataSet* dsetIn) {
  DataFile* DF = AddDataFile( nameIn );
  if (DF == 0) return 0;
  DF->AddSet( dsetIn );
  return DF;
}

// DataFileList::List()
/** Print information on what datasets are going to what datafiles */
void DataFileList::List() const {
  if (fileList_.empty()) {
    //mprintf("NO DATASETS WILL BE OUTPUT\n");
    return;
  }

  mprintf("\nDATAFILES:\n");
  for (DFarray::const_iterator it = fileList_.begin(); it != fileList_.end(); it++) {
    mprintf("  %s (%s): ",(*it)->DataFilename().base(), (*it)->FormatString());
    (*it)->DataSetNames();
    mprintf("\n");
  }
}

// DataFileList::WriteAllDF()
/** Call write for all DataFiles in list for which writeFile is true. Once
  * a file has been written set writeFile to false; it can be reset to
  * true if new DataSets are added to it.
  */
void DataFileList::WriteAllDF() {
  if (fileList_.empty()) return;
# ifdef TIMER
  Timer datafile_time;
  datafile_time.Start();
# endif
  for (DFarray::iterator df = fileList_.begin(); df != fileList_.end(); ++df) {
    if ( (*df)->DFLwrite() ) {
      (*df)->WriteData();
      (*df)->SetDFLwrite( false );
    }
  }
# ifdef TIMER
  datafile_time.Stop();
  mprintf("TIME: Write of all data files took %.4f seconds.\n", datafile_time.Total());
# endif
}

/** Reset writeFile status for all files in list to true. */
void DataFileList::ResetWriteStatus() {
  for (DFarray::iterator df = fileList_.begin(); df != fileList_.end(); ++df)
    (*df)->SetDFLwrite( true );
}

// DataFileList::ProcessDataFileArgs()
/** Process command relating to data files. */
int DataFileList::ProcessDataFileArgs(ArgList& dataArg) {
  // Next string is DataFile name that command will be passed to.
  std::string df_cmd = dataArg.GetStringNext();
  if (df_cmd.empty()) {
    mprintf("Warning: datafile: No filename given.\n");
    return 0;
  }
  // Check for deprecated commands
  if (df_cmd == "create" || df_cmd == "precision") 
    mprintf("Warning: 'datafile %s' is deprecated; use %s instead.\n", 
            df_cmd.c_str(), df_cmd.c_str());
  //mprintf("  [%s]\n",(*dataArg).ArgLine());
  DataFile* df = GetDataFile( df_cmd.c_str() );
  if (df == 0) {
    mprinterr("Error: datafile: File %s not found.\n", df_cmd.c_str());
    return 1;
  }
  // Process command
  int err = df->ProcessArgs( dataArg );
  if (err != 0 || dataArg.CheckForMoreArgs()) return 1;
  return 0;
}
