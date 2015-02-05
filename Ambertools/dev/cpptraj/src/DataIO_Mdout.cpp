#include <cstdio> // sscanf
#include <cstdlib> // atoi
#include <cstring> // strncmp
#include "DataIO_Mdout.h"
#include "BufferedLine.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // ConvertToDouble

// DataIO_Mdout::ID_DataFormat()
bool DataIO_Mdout::ID_DataFormat(CpptrajFile& infile) {
  if (infile.OpenFile()) return false;
  bool isMdout = false;
  std::string line = infile.GetLine();
  if (line[0] == '\n') {
    line = infile.GetLine();
    if (line.compare(0, 15, "          -----") == 0) {
      line = infile.GetLine();
      if (line.compare(0, 15, "          Amber") == 0)
        isMdout = true;
    }
  }
  infile.CloseFile();
  return isMdout;
}

static inline int EOF_ERROR() {
  mprinterr("Error: Unexpected EOF in MDOUT file.\n");
  return 1;
}

const char* DataIO_Mdout::Enames[] = {
  "Etot",   "EPtot", "GMAX",   "BOND", 
  "ANGLE",  "DIHED",  "VDW",   "EELEC",  "EGB",
  "VDW1-4", "EEL1-4", "RST",   "EAMBER", "Density",
  "RMS",    "EKtot",  "ESURF", "EAMD_BOOST", 0
};

/// \return index of name in Energy[] array, N_FIELDTYPES if not recognized.
DataIO_Mdout::FieldType DataIO_Mdout::getEindex(Sarray const& Name) {
  //mprintf("DEBUG:\tgetEindex(%s,%s)\n", Name[0].c_str(), Name[1].c_str());
  if (Name[0]=="Etot")  return Etot;
  if (Name[0]=="EPtot") return EPtot;
  if (Name[0]=="GMAX") return GMAX; // Not necessary?
  if (Name[0]=="BOND") return BOND;
  if (Name[0]=="ANGLE") return ANGLE;
  if (Name[0]=="DIHED") return DIHED;
  if (Name[0]=="VDWAALS") return VDWAALS;
  if (Name[0]=="EEL" || Name[0]=="EELEC") return EEL;
  if (Name[0]=="EGB") return EGB;
  if ((Name[0]=="1-4" && Name[1]=="VDW") || (Name[0]=="1-4" && Name[1]=="NB")) return VDW14;
  if  (Name[0]=="1-4" && Name[1]=="EEL") return EEL14;
  if (Name[0]=="RESTRAINT") return RESTRAINT;
  if (Name[0]=="EAMBER") return EAMBER;
  if (Name[0]=="Density") return Density;
  if (Name[0]=="RMS") return RMS; // Not necessary?
  if (Name[0]=="EKtot") return EKtot;
  if (Name[0]=="ESURF") return ESURF;
  if (Name[0]=="EAMD_BOOST") return EAMD_BOOST;
  return N_FIELDTYPES;
}

// DataIO_Mdout::ReadData()
int DataIO_Mdout::ReadData(std::string const& fname,
                            DataSetList& datasetlist, std::string const& dsname)
{
  mprintf("\tReading from mdout file: %s\n", fname.c_str());
  // Check if we are going to be appending to sets.
  DataSetList check = datasetlist.SelectSets(dsname + "[*]");
  bool isFirstFile = check.empty();
  // ----- CREATE DATASETS FOR ENERGIES -----
  std::vector<Darray> Esets( N_FIELDTYPES );
  std::vector<double> TimeVals;
  BufferedLine buffer;
  double lastx = 0.0;
  if (!isFirstFile) {
    DataSet_1D const& set = static_cast<DataSet_1D const&>( *(check[0]) );
    lastx = set.Xcrd( set.Size() - 1 );
    mprintf("\tAppending to sets named '%s'; start time is %g\n", dsname.c_str(), lastx);
  }
  unsigned int count = 0; // DataSet index TODO Use TimeVals size?
    if (buffer.OpenFileRead( fname )) return 1;
    // Read first line
    const char* ptr = buffer.Line();
    if (ptr == 0) {
      mprinterr("Error: Nothing in MDOUT file: %s\n", fname.c_str());
      return 1;
    }
    int imin = -1;           // imin for this file
    const char* Trigger = 0; // Trigger must be 8 chars long.
    int frame = 0;           // Frame counter for this file
    double dt = 1.0;         // Timestep for this file (MD)
    double t0 = 0.0;         // Initial time for this file (MD)
    int ntpr = 1;            // Value of ntpr
    int irest = 0;           // Value of irest
    // ----- PARSE THE INPUT SECTION ----- 
    while ( ptr != 0 && strncmp(ptr, "   2.  CONTROL  DATA", 20) != 0 )
      ptr = buffer.Line();
    if (ptr == 0) return EOF_ERROR();
    // Determine whether this is dynamics or minimization, get dt
    ptr = buffer.Line(); // Dashes 
    ptr = buffer.Line(); // Blank 
    ptr = buffer.Line(); // title line
    while ( strncmp(ptr, "   3.  ATOMIC", 13) != 0 ) 
    {
      ArgList mdin_args( ptr, " ,=" ); // Remove commas, equal signs
      // Scan for stuff we want
      //mprintf("DEBUG:\tInput[%i] %s", mdin_args.Nargs(), mdin_args.ArgLine());
      for (int col=0; col < mdin_args.Nargs(); col += 2) {
        int col1 = col + 1;
        if (mdin_args[col] == "imin") {
          imin = convertToInteger( mdin_args[ col1 ] );
          if (debug_ > 0) mprintf("\t\tMDIN: imin is %i\n", imin);
          // Set a trigger for printing. For imin5 this is the word minimization.
          // For imin0 or imin1 this is NSTEP.
          if      (imin==0) Trigger = " NSTEP =";
          else if (imin==1) Trigger = "   NSTEP";
          else if (imin==5) Trigger = "minimiza";
          // Since imin0 and imin1 first trigger has no data, set frame 1 lower.
          if (imin==1 || imin==0) frame = -1;
        } else if (mdin_args[col] == "dt") {
          dt = convertToDouble( mdin_args[ col1 ] );
          if (debug_ > 0) mprintf("\t\tMDIN: dt is %f\n", dt);
        } else if (mdin_args[col] == "t") {
          if (isFirstFile) {
            t0 = convertToDouble( mdin_args[ col1 ] );
            if (debug_ > 0) mprintf("\t\tMDIN: t is %f\n", t0);
          }
        } else if (mdin_args[col] == "ntpr") {
          ntpr = convertToInteger( mdin_args[ col1 ] );
          if (debug_ > 0) mprintf("\t\tMDIN: ntpr is %i\n", ntpr);
        } else if (mdin_args[col] == "irest") {
          irest = convertToInteger( mdin_args[ col1 ] );
          if (debug_ > 0) mprintf("\t\tMDIN: irest is %i\n", irest);
        }
      }
      ptr = buffer.Line();
      if (ptr == 0) return EOF_ERROR();
    }
    if (Trigger == 0) {
      mprinterr("Error: Could not determine whether MDOUT is md, min, or post-process.\n");
      return 1;
    }
    // ----- PARSE THE ATOMIC ... SECTION -----
    while ( ptr != 0 && strncmp(ptr, "   4.  RESULTS", 14) != 0 )
    {
      ptr = buffer.Line();
      // If this is the first mdout file being read and it is a restart,
      // set the initial time value.
      if (isFirstFile && irest == 1) {
        if (strncmp(ptr, " begin time", 11) == 0) {
          sscanf(ptr, " begin time read from input coords = %lf", &lastx);
          if (debug_ > 0) mprintf("\t\tMD restart initial time= %f\n", lastx);
        }
      }
    }
    if (ptr == 0) return EOF_ERROR();
    // ----- PARSE THE RESULTS SECTION -----
    bool finalE = false;
    int nstep;
    int minStep = 0; // For imin=1 only
    if (irest == 0)
      nstep = 0;
    else
      nstep = ntpr;
    double Energy[N_FIELDTYPES];
    std::fill( Energy, Energy+N_FIELDTYPES, 0.0 );
    std::vector<bool> EnergyExists(N_FIELDTYPES, false);
    Sarray Name(2);
    double time = 0.0;
    while (ptr != 0) {
      // Check for end of imin 0 or 1 run; do not record Average and Stdevs
      if ( (imin == 1 && (strncmp(ptr, "                    FINAL", 25) == 0 ||
                          strncmp(ptr, "   5.  TIMINGS",            14) == 0   )) ||
           (imin == 0 && strncmp(ptr, "      A V", 9) == 0))
        finalE = true;
      // Record set for energy post-processing
      if (imin == 5 && strncmp(ptr, "minimizing", 10) == 0)
        nstep = atoi( ptr + 22 );
      // MAIN OUTPUT ROUTINE
      // If the trigger has been reached print output.
      // For imin0 and imin1 the first trigger will have no data.
      // If the end of the file has been reached print then exit.
      if ( strncmp(ptr, Trigger, 8) == 0 || finalE ) {
        if (frame > -1) {
          // Data storage should go here
          for (int i = 0; i < (int)N_FIELDTYPES; i++)
              if (EnergyExists[i]) {
                // NOTE: Since energy terms can appear and vanish over the
                //       course of the mdout file, resize if necessary.
                if (count > Esets[i].size()) Esets[i].resize(count, 0.0);
                Esets[i].push_back( Energy[i] );
              }
          TimeVals.push_back( time );
          count++;
          nstep += ntpr;
        }
        frame++;
        if (finalE) break;
      }
      // Check for NSTEP in minimization or post-processing. Values will be
      // on the next line. NOTE: NSTEP means something different for imin=5.
      if ((imin == 1 || imin == 5) && strncmp(ptr, "   NSTEP", 8) == 0) {
        ptr = buffer.Line(); // Get next line
        //sscanf(ptr, " %6lf    %13lE  %13lE  %13lE", Energy+NSTEP, Energy+EPtot, Energy+RMS, Energy+GMAX);
        sscanf(ptr, " %i %lE %lE %lE", &minStep, Energy+EPtot, Energy+RMS, Energy+GMAX);
        EnergyExists[EPtot] = true;
        EnergyExists[RMS] = true;
        EnergyExists[GMAX] = true;
        ptr = buffer.Line();
      }
      // Tokenize line, scan through until '=' is reached; value after is target.
      int ntokens = buffer.TokenizeLine(" ");
      if (ntokens > 0) {
        int nidx = 0;
        Name[0].clear();
        Name[1].clear();
        for (int tidx = 0; tidx < ntokens; tidx++) {
          const char* tkn = buffer.NextToken();
          if (tkn[0] == '=') {
            FieldType Eindex = getEindex(Name);
            tkn = buffer.NextToken();
            ++tidx;
            if (Eindex != N_FIELDTYPES) {
              Energy[Eindex] = atof( tkn );
              EnergyExists[Eindex] = true;
            }
            nidx = 0;
            Name[0].clear();
            Name[1].clear();
          } else {
            if (nidx > 1) break; // Two tokens, no '=' found. Not an E line.
            Name[nidx++].assign( tkn );
          }
        }
      }
      // Set time
      switch (imin) {
        case 5: time = (double)nstep + lastx; break;
        case 1: time = (double)minStep + lastx; break;
        case 0: time = ((double)nstep * dt) + t0 + lastx; break;
      }
      // Read in next line
      ptr = buffer.Line();
    }
    mprintf("\t%i frames\n", frame);
    lastx = time;
    buffer.CloseFile();
  // ----- SET UP DATA SETS -----
  for (int i = 0; i < (int)N_FIELDTYPES; i++) {
    if (!Esets[i].empty()) {
      DataSet* ds = datasetlist.AddOrAppendSet(dsname, -1, Enames[i], TimeVals, Esets[i]);
      if (ds == 0) return 1;
      ds->SetLegend( dsname + "_" + Enames[i] );
    }
  }
  return 0;
}
