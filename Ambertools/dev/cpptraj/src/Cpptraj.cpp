#include <cstdio> // stdin, fileno
#include <unistd.h> // isatty
#include "Cpptraj.h"
#include "CpptrajStdio.h"
#include "Command.h"
#include "ReadLine.h"
#include "Version.h"
#include "ParmFile.h" // ProcessMask
#include "Timer.h"
#include "StringRoutines.h" // TimeString

void Cpptraj::Usage() {
  mprinterr("\n"
            "Usage: cpptraj [-p <Top0>] [-i <Input0>] [-y <trajin>] [-x <trajout>]\n"
            "               [-c <reference>]\n"
            "               [-h | --help] [-V | --version] [--defines] [-debug <#>]\n"
            "               [--interactive] [--log <logfile>] [-tl]\n"
            "               [-ms <mask>] [-mr <mask>] [--mask <mask>] [--resmask <mask>]\n"
            "       cpptraj <Top> <Input>\n"
            "\t-p <Top0>        : Load <Top0> as a topology file. May be specified more than once.\n"
            "\t-i <Input0>      : Read input from <Input0>. May be specified more than once.\n"
            "\t-y <trajin>      : Read from trajectory file <trajin>; same as input 'trajin <trajin>'.\n"
            "\t-x <trajout>     : Write trajectory file <trajout>; same as input 'trajout <trajout>'.\n"
            "\t-c <reference>   : Read <reference> as reference coordinates; same as input 'reference <reference>'.\n"
            "\t-h | --help      : Print command line help and exit.\n"
            "\t-V | --version   : Print version and exit.\n"
            "\t--defines        : Print compiler defines and exit.\n"
            "\t-debug <#>       : Set global debug level to <#>; same as input 'debug <#>'.\n"
            "\t--interactive    : Force interactive mode.\n"
            "\t--log <logfile>  : Record commands to <logfile> (interactive mode only). Default is 'cpptraj.log'.\n"
            "\t-tl              : Print length of trajectories specified with '-y' to STDOUT.\n"
            "\t-ms <mask>       : Print selected atom numbers to STDOUT.\n"
            "\t-mr <mask>       : Print selected residue numbers to STDOUT.\n"
            "\t--mask <mask>    : Print detailed atom selection to STDOUT.\n"
            "\t--resmask <mask> : Print detailed residue selection to STDOUT.\n\n");
}

void Cpptraj::Intro() {
  mprintf("\nCPPTRAJ: Trajectory Analysis. %s"
# ifdef MPI
          " MPI"
# endif
# ifdef _OPENMP
          " OpenMP"
# endif
          "\n    ___  ___  ___  ___\n     | \\/ | \\/ | \\/ | \n    _|_/\\_|_/\\_|_/\\_|_\n",
          CPPTRAJ_VERSION_STRING);
# ifdef MPI
  mprintf("Running on %i threads\n",CpptrajState::WorldSize());
# endif
  mprintf("\n| Date/time: %s\n", TimeString().c_str());
  double available_mem = AvailableMemory_MB();
  // If < 0 could not be calculated correctly.
  if (available_mem > 0.0)
    mprintf(  "| Available memory: %g MB\n", AvailableMemory_MB());
  mprintf("\n");
}

void Cpptraj::Finalize() {
  mprintf("--------------------------------------------------------------------------------\n"
    "To cite CPPTRAJ use:\n"
    "Daniel R. Roe and Thomas E. Cheatham, III, \"PTRAJ and CPPTRAJ: Software for\n"
    "  Processing and Analysis of Molecular Dynamics Trajectory Data\". J. Chem.\n"
    "  Theory Comput., 2013, 9 (7), pp 3084-3095.\n");
}

int Cpptraj::RunCpptraj(int argc, char** argv) {
  int err = 0;
  Timer total_time;
  total_time.Start();
  Mode cmode = ProcessCmdLineArgs(argc, argv);
  if ( cmode == BATCH ) {
    // If State is not empty, run now. 
    if (!State_.EmptyState())
      err = State_.Run();
  } else if ( cmode == INTERACTIVE ) {
    err = Interactive();
  } else if ( cmode == ERROR ) {
    err = 1;
  }
  total_time.Stop();
  if (cmode != INTERACTIVE)
    mprintf("TIME: Total execution time: %.4f seconds.\n", total_time.Total());
  if (err == 0) Cpptraj::Finalize();
  mprintf("\n");
  return err;
}

/** Process a mask from the command line. */
int Cpptraj::ProcessMask( Sarray const& topFiles, Sarray const& refFiles,
                          std::string const& maskexpr,
                          bool verbose, bool residue ) const
{
  SetWorldSilent(true);
  if (topFiles.empty()) {
    mprinterr("Error: No topology file specified.\n");
    return 1;
  }
  ParmFile pfile;
  Topology parm;
  if (pfile.ReadTopology(parm, topFiles[0], State_.Debug())) return 1;
  if (!refFiles.empty()) {
    DataSet_Coords_REF refFrame;
    if (refFrame.LoadRef( refFiles[0], parm, State_.Debug())) return 1;
    parm.SetReferenceCoords( refFrame.RefFrame() );
  }
  if (!verbose) {
    AtomMask tempMask( maskexpr );
    if (parm.SetupIntegerMask( tempMask )) return 1;
    loudPrintf("Selected=");
    if (residue) {
      int res = -1;
      for (AtomMask::const_iterator atom = tempMask.begin(); 
                                    atom != tempMask.end(); ++atom)
      {
        if (parm[*atom].ResNum() > res) {
          loudPrintf(" %i", parm[*atom].ResNum()+1);
          res = parm[*atom].ResNum();
        }
      }
    } else
      for (AtomMask::const_iterator atom = tempMask.begin();
                                    atom != tempMask.end(); ++atom)
        loudPrintf(" %i", *atom + 1);
    loudPrintf("\n");
  } else {
    if (residue)
      parm.PrintResidueInfo( maskexpr );
    else
      parm.PrintAtomInfo( maskexpr );
  }
  return 0;
}

/** Read command line args. */
Cpptraj::Mode Cpptraj::ProcessCmdLineArgs(int argc, char** argv) {
  bool hasInput = false;
  bool interactive = false;
  Sarray inputFiles;
  Sarray topFiles;
  Sarray trajinFiles;
  Sarray trajoutFiles;
  Sarray refFiles;
  for (int i = 1; i < argc; i++) {
    std::string arg(argv[i]);
    if ( arg == "--help" || arg == "-h" ) {
      // --help, -help: Print usage and exit
      SetWorldSilent(true);
      Usage();
      return QUIT;
    }
    if ( arg == "-V" || arg == "--version" ) {
      // -V, --version: Print version number and exit
      SetWorldSilent( true );
      loudPrintf("CPPTRAJ: Version %s\n", CPPTRAJ_VERSION_STRING);
      return QUIT;
    }
    if ( arg == "--internal-version" ) {
      // --internal-version: Print internal version number and quit.
      SetWorldSilent( true );
      loudPrintf("CPPTRAJ: Internal version # %s\n", CPPTRAJ_INTERNAL_VERSION);
      return QUIT;
    }
    if ( arg == "--defines" ) {
      // --defines: Print information on compiler defines used and exit
      SetWorldSilent( true );
      loudPrintf("Compiled with:");
#     ifdef DEBUG
      loudPrintf(" -DDEBUG");
#     endif
#     ifdef HASBZ2
      loudPrintf(" -DHASBZ2");
#     endif
#     ifdef HASGZ
      loudPrintf(" -DHASGZ");
#     endif
#     ifdef BINTRAJ
      loudPrintf(" -DBINTRAJ");
#     endif
#     ifdef MPI
      loudPrintf(" -DMPI");
#     endif
#     ifdef _OPENMP
      loudPrintf(" -D_OPENMP");
#     endif
#     ifdef NO_MATHLIB
      loudPrintf(" -DNO_MATHLIB");
#     endif
#     ifdef TIMER
      loudPrintf(" -DTIMER");
#     endif
      loudPrintf("\n");
      return QUIT;
    }
    if (arg == "-tl") {
      // -tl: Trajectory length
      if (topFiles.empty()) {
        mprinterr("Error: No topology file specified.\n");
        return ERROR;
      }
      SetWorldSilent( true );
      if (State_.TrajLength( topFiles[0], trajinFiles )) return ERROR;
      return QUIT;
    }
    if ( arg == "--interactive" )
      interactive = true;
    else if ( arg == "-debug" && i+1 != argc) {
      // -debug: Set overall debug level
      ArgList dbgarg( argv[++i] );
      State_.SetListDebug( dbgarg );
    } else if ( arg == "--log" && i+1 != argc)
      // --log: Set up log file for interactive mode
      logfilename_ = argv[++i];
    else if ( arg == "-p" && i+1 != argc) {
      // -p: Topology file
      topFiles.push_back( argv[++i] );
    } else if ( arg == "-y" && i+1 != argc) {
      // -y: Trajectory file in
      trajinFiles.push_back( argv[++i] );
    } else if ( arg == "-x" && i+1 != argc) {
      // -x: Trajectory file out
      trajoutFiles.push_back( argv[++i] );
    } else if ( arg == "-c" && i+1 != argc) {
      // -c: Reference file
      refFiles.push_back( argv[++i] );
    } else if (arg == "-i" && i+1 != argc) {
      // -i: Input file(s)
      inputFiles.push_back( argv[++i] );
    } else if (arg == "-ms" && i+1 != argc) {
      // -ms: Parse mask string, print selected atom #s
      if (ProcessMask( topFiles, refFiles, std::string(argv[++i]), false, false )) return ERROR;
      return QUIT;
    } else if (arg == "-mr" && i+1 != argc) {
      // -mr: Parse mask string, print selected res #s
      if (ProcessMask( topFiles, refFiles, std::string(argv[++i]), false, true )) return ERROR;
      return QUIT;
    } else if (arg == "--mask" && i+1 != argc) {
      // --mask: Parse mask string, print selected atom details
      if (ProcessMask( topFiles, refFiles, std::string(argv[++i]), true, false )) return ERROR;
      return QUIT;
    } else if (arg == "--resmask" && i+1 != argc) {
      // --resmask: Parse mask string, print selected residue details
      if (ProcessMask( topFiles, refFiles, std::string(argv[++i]), true, true )) return ERROR;
      return QUIT;
    } else if ( i == 1 ) {
      // For backwards compatibility with PTRAJ; Position 1 = TOP file
      topFiles.push_back( argv[i] );
    } else if ( i == 2 ) {
      // For backwards compatibility with PTRAJ; Position 2 = INPUT file
      inputFiles.push_back( argv[i] );
    } else {
      // Unrecognized
      mprintf("  Unrecognized input on command line: %i: %s\n", i,argv[i]);
      Usage();
      return ERROR;
    }
  }
  Cpptraj::Intro();
  // Add all topology files specified on command line.
  for (Sarray::const_iterator topFilename = topFiles.begin();
                              topFilename != topFiles.end();
                              ++topFilename)
    if (State_.PFL()->AddParmFile( *topFilename )) return ERROR;
  // Add all reference trajectories specified on command line.
  for (Sarray::const_iterator refName = refFiles.begin();
                              refName != refFiles.end();
                              ++refName)
    if (State_.AddReference( *refName )) return ERROR;
  // Add all input trajectories specified on command line.
  for (Sarray::const_iterator trajinName = trajinFiles.begin();
                              trajinName != trajinFiles.end();
                              ++trajinName)
    if (State_.AddTrajin( *trajinName )) return ERROR;
  // Add all output trajectories specified on command line.
  if (!trajoutFiles.empty()) {
    hasInput = true; // This allows direct traj conversion with no other input 
    for (Sarray::const_iterator trajoutName = trajoutFiles.begin();
                                trajoutName != trajoutFiles.end();
                                ++trajoutName)
      if (State_.AddTrajout( *trajoutName )) return ERROR;
  }
  // Process all input files specified on command line.
  if ( !inputFiles.empty() ) {
    hasInput = true;
    for (Sarray::const_iterator inputFilename = inputFiles.begin();
                                inputFilename != inputFiles.end();
                                ++inputFilename)
    {
      Command::RetType c_err = Command::ProcessInput( State_, *inputFilename );
      if (c_err == Command::C_ERR && State_.ExitOnError()) return ERROR;
      if (c_err == Command::C_QUIT) return QUIT;
    }
  }
  // Determine whether to enter interactive mode
  if (!hasInput || interactive) {
    // Test if input is really from a console
    if ( isatty(fileno(stdin)) )
      return INTERACTIVE;
    else {
      // "" means read from STDIN
      Command::RetType c_err = Command::ProcessInput( State_, "" ); 
      if (c_err == Command::C_ERR && State_.ExitOnError()) return ERROR;
      if (c_err == Command::C_QUIT) return QUIT;
    }
  }
  return BATCH;
}

// Cpptraj::Interactive()
int Cpptraj::Interactive() {
  ReadLine inputLine;
  // By default when interactive do not exit on errors
  State_.SetNoExitOnError();
  // Open log file. If no name has been set, use default.
  CpptrajFile logfile_;
  if (logfilename_.empty())
    logfilename_.assign("cpptraj.log");
  if (fileExists(logfilename_)) {
    // Load previous history.
    if (logfile_.OpenRead(logfilename_)==0) {
      mprintf("\tLoading previous history from log '%s'\n", logfile_.Filename().full());
      std::string previousLine = logfile_.GetLine();
      while (!previousLine.empty()) {
        if (previousLine[0] != '#') {
          // Remove any newline chars.
          std::size_t found = previousLine.find_first_of("\r\n");
          if (found != std::string::npos)
            previousLine[found] = '\0';
          inputLine.AddHistory( previousLine.c_str() );
        }
        previousLine = logfile_.GetLine();
      }
      logfile_.CloseFile();
    }
  }
  logfile_.OpenAppend(logfilename_);
  if (logfile_.IsOpen())
    logfile_.Printf("# %s\n", TimeString().c_str());
  Command::RetType readLoop = Command::C_OK;
  while ( readLoop != Command::C_QUIT ) {
    if (inputLine.GetInput()) {
      // EOF (Ctrl-D) specified. If state is not empty, ask before exiting.
      if (!State_.EmptyState()) {
        if (inputLine.YesNoPrompt("EOF (Ctrl-D) specified but there are actions/"
                                  "analyses/trajectories queued. Really quit? [y/n]> "))
          break;
      } else
        break;
    }
    if (!inputLine.empty()) {
      readLoop = Command::Dispatch( State_, *inputLine );
      if (logfile_.IsOpen() && readLoop != Command::C_ERR) {
        logfile_.Printf("%s\n", inputLine.c_str());
        logfile_.Flush();
      }
    }
    // If state is not empty, ask before exiting.
    if (readLoop == Command::C_QUIT && !State_.EmptyState()) {
      if (inputLine.YesNoPrompt("There are actions/analyses/trajectories queued. "
                                "Really quit? [y/n]> "))
        break;
      else
        readLoop = Command::C_OK;
    }
  }
  logfile_.CloseFile();
  if (readLoop == Command::C_ERR) return 1;
  return 0;
}
