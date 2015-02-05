// ParmList
#include "TopologyList.h"
#include "CpptrajStdio.h"
#include "ParmFile.h"
#include "StringRoutines.h" // ExpandToFilenames

// CONSTRUCTOR 
TopologyList::TopologyList() : debug_(0) {}

// DESTRUCTOR
TopologyList::~TopologyList() {
  Clear();
}

void TopologyList::SetDebug(int debugIn) {
  debug_ = debugIn;
  if (debug_ > 0)
    mprintf("TopologyList debug level set to %i\n", debug_);
}

void TopologyList::Clear() {
  for (std::vector<Topology*>::iterator top = TopList_.begin();
                                        top != TopList_.end(); top++)
    delete *top;
  TopList_.clear();
}

// TopologyList::GetParm()
/** Return the parm structure with index num. */
Topology* TopologyList::GetParm(int num) const {
  if (num>=(int)TopList_.size() || num<0) return 0;
  return TopList_[num];
}

// TopologyList::GetParmByIndex()
Topology* TopologyList::GetParmByIndex(ArgList& argIn) const {
  int pindex = argIn.getNextInteger(0);
  Topology* parm = GetParm( pindex );
  if ( parm == 0 ) {
    mprinterr("Error: parm index %i not loaded.\n",pindex);
    return 0;
  }
  return parm;
}

const char* TopologyList::ParmArgs = "[parm <parmfile / tag> | parmindex <#>]";

// TopologyList::GetParm()
/** Return the parm structure based on arguments in the given arg list. 
  *   parm <parm name>
  *   parmindex <parm index>
  * \param argIn argument list that contains parm-related keyword
  * \return parm specified by 'parm' or 'parmindex', or the first parm. null on error.
  */
Topology* TopologyList::GetParm(ArgList &argIn) const {
  Topology* ParmOut = 0;
  // Get any parm keywords if present
  int pindex = argIn.getKeyInt("parmindex",0);
  std::string parmfilename = argIn.GetStringKey("parm");
  if (!parmfilename.empty()) {
    for (std::vector<Topology*>::const_iterator tf = TopList_.begin();
                                              tf != TopList_.end(); ++tf)
    {
      if ( (*tf)->Tag()==parmfilename || 
           (*tf)->OriginalFilename().MatchFullOrBase( parmfilename ) )
      {
        ParmOut = *tf;
        break;
      }
    }
  } else {
    ParmOut = GetParm(pindex);
  }
  if (ParmOut==0) {
    mprintf("Warning: Could not get parameter file:\n");
    mprintf("Warning: parmname=%s, pindex=%i\n",parmfilename.c_str(),pindex);
    return 0;
  }
  return ParmOut;
}

// TopologyList::AddParmFile()
/** Add a parameter file to the parm file list. */
int TopologyList::AddParmFile(std::string const& filename) {
  ArgList arg;
  return AddParmFile(filename, arg);
}

// TopologyList::AddParmFile()
/** Add a parameter file to the parm file list with optional tag. */
int TopologyList::AddParmFile(std::string const& filenameIn, ArgList& argIn) 
{
  if (filenameIn.empty()) {
    mprinterr("Error: No topology file name specified.\n");
    return 1;
  }
  StrArray fnames = ExpandToFilenames( filenameIn );
  if (fnames.empty()) {
    mprinterr("Error: '%s' matches no files.\n", filenameIn.c_str());
    return 1;
  }
  std::string ParmTag = argIn.getNextTag();
  int numErr = 0;
  for (StrArray::const_iterator fn = fnames.begin();
                                fn != fnames.end(); ++fn)
  {
    // Check if filename/parmtag already in use
    bool skipFile = false;
    for (std::vector<Topology*>::const_iterator tf = TopList_.begin();
                                                tf != TopList_.end(); ++tf)
    {
      if ( (*tf)->OriginalFilename().Full() == *fn ) {
        mprintf("Warning: Parm '%s' already loaded, skipping.\n",(*fn).c_str());
        skipFile = true;
      }
      if ( !ParmTag.empty() && (*tf)->Tag() == ParmTag ) {
        mprinterr("Error: Parm tag '%s' already in use.\n",ParmTag.c_str());
        numErr++;
        skipFile = true;
      }
    }
    if (skipFile) continue;

    Topology* parm = new Topology();
    ParmFile pfile;
    // NOTE: Arg list will not be modified for multiple parms 
    if (pfile.ReadTopology(*parm, *fn, argIn, debug_)) {
      mprinterr("Error: Could not open topology '%s'\n",(*fn).c_str());
      delete parm;
      numErr++;
      continue;
    }

    if (debug_>0) 
      mprintf("    PARAMETER FILE %zu: %s\n",TopList_.size(),(*fn).c_str());
    // pindex is used for quick identification of the parm file
    parm->SetPindex( TopList_.size() );
    TopList_.push_back(parm);
    parm->SetTag( ParmTag );
    // Only allow the first parm to be tagged.
    ParmTag.clear();
  }
  if (numErr > 0) return 1;
  return 0;
}

// TopologyList::List()
/** Print list of loaded parameter files */
void TopologyList::List() const {
  if (!TopList_.empty()) {
    mprintf("\nPARAMETER FILES:\n");
    for (std::vector<Topology*>::const_iterator top = TopList_.begin();
                                                top != TopList_.end(); top++)
    {
      mprintf(" %i:", (*top)->Pindex());
      (*top)->Brief(0);
      if ((*top)->Nframes() > 0)
        mprintf(", %i frames", (*top)->Nframes());
      mprintf("\n");
    }
  }
}
