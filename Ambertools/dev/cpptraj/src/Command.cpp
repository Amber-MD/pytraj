#include <cstdio> // for ProcessInput
#include <cstdlib> // system
#include "Command.h"
#include "CpptrajStdio.h"
#include "DistRoutines.h" // GenerateAmberRst
#include "TorsionRoutines.h" // GenerateAmberRst
#include "Constants.h" // GenerateAmberRst
#include "DataSet_Coords_TRJ.h" // LoadTraj
#include "DataSet_double.h" // DataSetCmd
#include "ParmFile.h" // ReadOptions, WriteOptions
#include "StringRoutines.h" // tildeExpansion
#include "Timer.h"
#include "RPNcalc.h" // Calc
// INC_ACTION==================== ALL ACTION CLASSES GO HERE ===================
#include "Action_Distance.h"
#include "Action_Rmsd.h"
#include "Action_Dihedral.h"
#include "Action_Angle.h"
#include "Action_AtomMap.h"
#include "Action_Strip.h"
#include "Action_DSSP.h"
#include "Action_Center.h"
#include "Action_Hbond.h"
#include "Action_Image.h"
#include "Action_Surf.h"
#include "Action_Radgyr.h"
#include "Action_Mask.h"
#include "Action_Closest.h"
#include "Action_NAstruct.h"
#include "Action_Pucker.h"
#include "Action_Outtraj.h"
#include "Action_Average.h"
#include "Action_Radial.h"
#include "Action_DistRmsd.h"
#include "Action_Jcoupling.h"
#include "Action_Pairwise.h"
#include "Action_Molsurf.h"
#include "Action_CheckStructure.h"
#include "Action_DihedralScan.h"
#include "Action_Rotdif.h"
#include "Action_RunningAvg.h"
#include "Action_AtomicFluct.h"
#include "Action_Watershell.h"
#include "Action_Contacts.h"
#include "Action_Vector.h"
#include "Action_Principal.h"
#include "Action_Matrix.h"
#include "Action_LIE.h"
#include "Action_Grid.h"
#include "Action_GridFreeEnergy.h"
#include "Action_Dipole.h"
#include "Action_Projection.h"
#include "Action_ClusterDihedral.h"
#include "Action_Unwrap.h"
#include "Action_Diffusion.h"
#include "Action_DNAionTracker.h"
#include "Action_Scale.h"
#include "Action_RandomizeIons.h"
#include "Action_AutoImage.h"
#include "Action_STFC_Diffusion.h"
#include "Action_AtomicCorr.h"
#include "Action_Bounds.h"
#include "Action_Rotate.h"
#include "Action_Translate.h"
#include "Action_Box.h"
#include "Action_CreateCrd.h"
#include "Action_MultiDihedral.h"
#include "Action_MakeStructure.h"
#include "Action_SymmetricRmsd.h"
#include "Action_Volmap.h"
#include "Action_Spam.h"
#include "Action_Temperature.h"
#include "Action_Gist.h"
#include "Action_CreateReservoir.h"
#include "Action_Density.h"
#include "Action_PairDist.h"
#include "Action_OrderParameter.h"
#include "Action_FixAtomOrder.h"
#include "Action_NMRrst.h"
#include "Action_FilterByData.h"
#include "Action_LESsplit.h"
#include "Action_NativeContacts.h"
#include "Action_VelocityAutoCorr.h"
#include "Action_SetVelocity.h"
#include "Action_MultiVector.h"
#include "Action_MinImage.h"
#include "Action_ReplicateCell.h"
#include "Action_AreaPerMol.h"
#include "Action_Energy.h"
#include "Action_CheckChirality.h"
#include "Action_Channel.h" // EXPERIMENTAL

// INC_ANALYSIS================= ALL ANALYSIS CLASSES GO HERE ==================
#include "Analysis_Hist.h"
#include "Analysis_Corr.h"
#include "Analysis_Matrix.h"
#include "Analysis_Timecorr.h"
#include "Analysis_IRED.h"
#include "Analysis_Modes.h"
#include "Analysis_CrankShaft.h"
#include "Analysis_Statistics.h"
#include "Analysis_CrossCorr.h"
#include "Analysis_AutoCorr.h"
#include "Analysis_Lifetime.h"
#include "Analysis_FFT.h"
#include "Analysis_CrdFluct.h"
#include "Analysis_RmsAvgCorr.h"
#include "Analysis_Rms2d.h"
#include "Analysis_Clustering.h"
#include "Analysis_RunningAvg.h"
#include "Analysis_MeltCurve.h"
#include "Analysis_Overlap.h"
#include "Analysis_AmdBias.h"
#include "Analysis_RemLog.h"
#include "Analysis_Integrate.h"
#include "Analysis_Spline.h"
#include "Analysis_Average.h"
#include "Analysis_KDE.h"
#include "Analysis_MultiHist.h"
#include "Analysis_Divergence.h"
#include "Analysis_VectorMath.h"
#include "Analysis_Regression.h"
#include "Analysis_LowestCurve.h"
#include "Analysis_CurveFit.h"
#include "Analysis_PhiPsi.h"
// ---- Command Functions ------------------------------------------------------
/// Warn about deprecated commands.
void Command::WarnDeprecated(TokenPtr token)
{
  mprinterr("Error: '%s' is deprecated.\n", token->Cmd);
  if (token->Help != 0)
    token->Help();
}

/** Search Commands list for a specific type of command. */
Command::TokenPtr Command::SearchTokenType(CommandType dtype,
                                           ArgList const& argIn)
{
  for (TokenPtr token = Commands; token->Type != NONE; ++token)
  {
    if (token->Type == DEPRECATED && argIn.CommandIs( token->Cmd )) {
      WarnDeprecated( token );
      return 0;
    }
    if (dtype != token->Type) continue;
    if (argIn.CommandIs( token->Cmd )) return token;
  }
  mprinterr("'%s': Command not found.\n", argIn.Command());
  return 0;
}

/// Strings that correspond to CommandType
const char* Command::CommandTitle[] = { 0, "Topology", "Trajectory", "Action",
  "Analysis", "General", "Deprecated" };

/** List all commands of the given type, or all commands if type
  * is NONE.
  */
void Command::ListCommands(CommandType dtype) {
  std::string Line;
  CommandType lastType = NONE;
  for (TokenPtr token = Commands; token->Type != DEPRECATED; ++token)
  {
    CommandType currentType = token->Type;
    if (dtype != NONE && dtype != currentType) continue;
    // Command group type title
    if (currentType != lastType) {
      if (!Line.empty()) {
        mprintf("%s\n", Line.c_str());
        Line.clear();
      }
      mprintf("%s Commands:\n", CommandTitle[currentType]);
      lastType = currentType;
    }
    if (Line.empty()) Line.assign("        ");
    std::string Command(token->Cmd);
    Command.append(" ");
    if ( Line.size() + Command.size() > 80 ) {
      mprintf("%s\n", Line.c_str());
      Line.assign("        ");
    }
    Line.append(Command);
  }
  if (!Line.empty())
    mprintf("%s\n", Line.c_str());
}

/** Search the Commands list for given command.
  * \return the token if found, 0 if not.
  */
Command::TokenPtr Command::SearchToken(ArgList& argIn) {
  // SPECIAL CASE: For backwards compat. remove analyze prefix
  if (argIn.CommandIs("analyze")) {
    argIn.RemoveFirstArg();
    argIn.MarkArg(0); // Mark new first arg as command
    return (SearchTokenType(ANALYSIS, argIn));
  }
  // Search for command.
  for (TokenPtr token = Commands; token->Type != NONE; ++token)
    if (argIn.CommandIs( token->Cmd )) {
      if (token->Type == DEPRECATED) {
        WarnDeprecated( token );
        return 0; 
      } else
        return token;
    }
  //mprinterr("'%s': Command not found.\n", argIn.Command());
  return 0;
}

/** Search for the given command and execute it. */
Command::RetType Command::Dispatch(CpptrajState& State,
                                   std::string const& commandIn)
{
  ArgList cmdArg( commandIn );
  cmdArg.MarkArg(0); // Always mark the first arg as the command 
  TokenPtr cmdToken = SearchToken( cmdArg );
  Command::RetType ret_val = Command::C_OK;
  if (cmdToken == 0) {
    // Try to evaluate the expression.
    RPNcalc calc;
    calc.SetDebug( State.Debug() );
    if (calc.ProcessExpression( commandIn ))
      ret_val = Command::C_ERR;
    else {
      if (calc.Evaluate(*State.DSL()))
        ret_val = Command::C_ERR;
    }
    if (ret_val == Command::C_ERR)
      mprinterr("'%s': Invalid command or expression.\n", commandIn.c_str());
  } else
    ret_val = cmdToken->Fxn( State, cmdArg, cmdToken->Alloc );
  return ret_val;
}

/// Used by ProcessInput to determine when line ends.
static inline bool EndChar(char ptr) {
  if (ptr=='\n' || ptr=='\r' || ptr=='\0' || ptr==EOF) return true;
  return false;
}

/** Read commands from an input file, or from STDIN if given filename
  * is empty. '#' indicates the beginning of a comment, backslash at the 
  * end of a line indicates continuation (otherwise indicates 'literal').
  * \return 0 if successfully read, 1 on error.
  */
Command::RetType Command::ProcessInput(CpptrajState& State, 
                                       std::string const& inputFilename)
{
  FILE* infile; // TODO: CpptrajFile
  if (inputFilename.empty()) {
    mprintf("INPUT: Reading Input from STDIN\n");
    infile = stdin;
  } else {
    std::string fname = tildeExpansion(inputFilename);
    if (fname.empty()) return C_ERR;
    mprintf("INPUT: Reading Input from file %s\n", fname.c_str());
    if ( (infile=fopen(fname.c_str(),"r"))==0 ) {
      rprinterr("Error: Could not open input file %s\n", fname.c_str());
      return C_ERR;
    }
  }
  // Read in each line of input. Newline or null terminates. \ continues line.
  std::string inputLine;
  unsigned int idx = 0;
  char lastchar = '0';
  char ptr = 0;
  int nInputErrors = 0;
  RetType cmode = C_OK;
  while ( ptr != EOF ) {
    ptr = (char)fgetc(infile);
    // Skip leading whitespace
    if (idx == 0 && isspace(ptr)) {
      while ( (ptr = (char)fgetc(infile))!=EOF )
        if ( !isspace(ptr) ) break;
    }
    // If '#' is encountered, skip the rest of the line
    if (ptr=='#')
      while (!EndChar(ptr)) ptr=(char)fgetc(infile);
    // newline, null, or EOF terminates the line
    if (EndChar(ptr)) {
      // If no chars in string continue
      if (inputLine.empty()) continue;
      // Print the input line that will be sent to dispatch
      mprintf("  [%s]\n",inputLine.c_str());
      // Call Dispatch to convert input to arglist and process.
      cmode = Command::Dispatch(State, inputLine);
      if (cmode == C_ERR) {
        nInputErrors++;
        if (State.ExitOnError()) break;
      } else if (cmode == C_QUIT)
        break;
      // Reset Input line
      inputLine.clear();
      idx = 0;
      continue;
    }
    // Any consecutive whitespace is skipped
    if (idx > 0) lastchar = inputLine[idx-1];
    if (isspace(ptr) && isspace(lastchar)) continue;
    // Backslash followed by newline continues to next line. Otherwise backslash
    // followed by next char will be inserted. 
    if (ptr=='\\') {
      ptr = (char)fgetc(infile);
      if ( ptr == EOF ) break;
      if (ptr == '\n' || ptr == '\r') continue;
      inputLine += "\\";
      inputLine += ptr;
      idx += 2;
      continue;
    }
    // Add character to input line
    inputLine += ptr;
    ++idx;
  }
  if (!inputFilename.empty())
    fclose(infile);
  if (nInputErrors > 0) {
    mprinterr("\t%i errors encountered reading input.\n", nInputErrors);
    return C_ERR;
  }
  return cmode;
}

// ====================== CPPTRAJ COMMANDS HELP ================================
static void Help_Help() {
  mprintf("\t{[<cmd>] | General | Action | Analysis | Topology | Trajectory}\n"
          "  With no arguments list all known commands, otherwise display help for\n"
          "  command <cmd>. If General/Action/Analysis/Topology/Trajectory specified\n"
          "  list commands only in that category.\n");
}

static void Help_System() { mprintf("  Call command from system.\n"); }

static void Help_NoProgress() {
  mprintf("  Do not print progress while reading in trajectories.\n");
}

static void Help_NoExitOnError() {
  mprintf("  Do not exit when errors are encountered. This is the default\n"
          "  in interactive mode.\n");
}

static void Help_Run() {
  mprintf("  Process all trajectories currently in input trajectory list.\n"
          "  All actions in action list will be run on each frame.\n"
          "  If not processing ensemble input, all analyses in analysis\n"
          "  list will be run after trajectory processing.\n");
}

static void Help_Quit() { mprintf("  Exit CPPTRAJ\n"); }

static void Help_List() {
  mprintf("\t[<type>] (<type> =%s)\n"
          "  List currently loaded objects of the specified type. If no type is given\n"
          "  then list all loaded objects.\n", CpptrajState::PrintListKeys().c_str());
}

static void Help_Debug() {
  mprintf("\t[<type>] <#> (<type> =%s)\n", CpptrajState::PrintListKeys().c_str());
  mprintf("  Set debug level for new objects of the specified type. If no type is given\n"
          "  then set debug level for all new objects. Does not affect current objects.\n");
}

static void Help_Clear() {
  mprintf("\t[ {all | <type>} ] (<type> =%s)\n", CpptrajState::PrintListKeys().c_str());
  mprintf("  Clear currently loaded objects of the specified type. If 'all' is specified\n"
          "  then clear all loaded objects.\n");
}

static void Help_RemoveData() {
  mprintf("\t[<arg>]\n"
          "  Remove data sets(s) corresponding to <arg> from data set list.\n");
}

static void Help_ActiveRef() {
  mprintf("\t<#>\n"
          "  Set the reference structure to be used for coordinate-based mask parsing.\n"
          "  <#> starts from 0 (first loaded reference).\n");
}

static void Help_Create_DataFile() {
  mprintf("\t<filename> <dataset0> [<dataset1> ...]\n"
          "  Add a file with specified data sets to the data file list. Does not\n"
          "  immediately write the data.\n");
  DataFile::WriteHelp();
  DataFile::WriteOptions();
}

static void Help_DataFile() {
  mprintf("\t<data filename> <datafile cmd>\n"
          "  Pass <datafile cmd> to specified data file currently in data file list.\n");
  DataFile::WriteHelp();
  DataFile::WriteOptions();
}

static void Help_DataSetCmd() {
  mprintf("\t{ legend <legend> <set> | \n"
          "\t  [mode <mode>] [type <type>] <set arg1> [<set arg 2> ...] }\n"
          "\tOptions for 'type noe':\n"
          "\t  %s\n"
          "  Either set the legend for a single data set, or change the mode/type for\n"
          "  one or more data sets.\n", Action_Distance::NOE_Help);
}

static void Help_ReadData() {
  mprintf("\t<filename> [name <dsname>] [as <fmt>] [<format options>]\n"
          "  Read data from <filename> into data sets.\n");
  DataFile::ReadOptions();
}

static void Help_ReadInput() {
  mprintf("\t<filename>\n"
          "  Read commands from input file <filename>\n");
}

static void Help_Write_DataFile() {
  mprintf("\t[<filename> <dataset0> [<dataset1> ...]]\n");
  DataFile::WriteHelp();
  mprintf("  With no arguments, write all files currently in the data file list.\n"
          "  Otherwise, write specified data sets to <filename> immediately.\n");
  DataFile::WriteOptions();
}

static void Help_Precision() {
  mprintf("\t{<filename> | <dataset arg>} [<width>] [<precision>]\n"
          "  Set precision for all datasets in datafile <filename> or dataset(s)\n"
          "  specified by <dataset arg> to <width>.<precision>. If width/precision\n"
          "  is not specified then default to 12.4\n");
}

static void Help_Select() {
  mprintf("\t[<parmindex>] <mask>\n"
          "  Show atom numbers selected by <mask> for parm <parmindex>\n"
          "  (default first parm)\n");
}

static void Help_SelectDS() {
  mprintf("\t<dataset selection>\n"
          "  Show results of data set selection. Data set selection format is:\n"
          "\t<name>[<aspect]:<idx range>\n"
          "  Where '<name>' is the data set name, '[<aspect>]' is the data set aspect,\n"
          "  and <idx range> is a numerical range specifying data set indices (i.e. 2-5,7 etc).\n"
          "  The aspect and index portions may be optional. An asterisk '*' may be used as\n"
          "  a wildcard. E.g. 'selectds R2', 'selectds RoG[Max]', 'selectds PR[res]:2-12'\n");
}

static void Help_Trajin() {
  mprintf("\t<filename> {[<start>] [<stop> | last] [offset]} | lastframe\n"
          "\t           %s\n", TopologyList::ParmArgs);
  mprintf("\t           [ remdtraj [remdtrajtemp <T> | remdtrajidx <#>]\n"
          "\t             [trajnames <rep1>,<rep2>,...,<repN> ] ]\n"
          "  Load trajectory specified by <filename> to the input trajectory list.\n");
  TrajectoryFile::ReadOptions();
}

static void Help_Ensemble() {
  mprintf("\t<file0> {[<start>] [<stop> | last] [offset]} | lastframe\n"
          "\t        %s\n", TopologyList::ParmArgs);
  mprintf("\t        [trajnames <file1>,<file2>,...,<fileN>\n"
          "\t        [remlog <remlogfile> [nstlim <nstlim> ntwx <ntwx>]]\n"
          "  Load an ensemble of trajectories starting with <file0> that will be\n"
          "  processed together as an ensemble.\n");
}

static void Help_Trajout() {
  mprintf("\t<filename> [<fileformat>] [append] [nobox]\n"
          "\t           %s [onlyframes <range>] [title <title>]\n", TopologyList::ParmArgs);
  mprintf("\t           %s\n", ActionFrameCounter::HelpText);
  mprintf("\t           [ <Format Options> ]\n"
          "  Write frames after all actions have been processed to output trajectory\n"
          "  specified by <filename>.\n");
  TrajectoryFile::WriteOptions();
}

static void Help_Reference() {
  mprintf("\t<filename> [<frame#>] [<mask>] [TAG] [lastframe]\n"
          "  Load trajectory <filename> as a reference frame.\n");
}

static void Help_Parm() {
  mprintf("\t<filename> [<tag>] [nobondsearch | bondsearch [<offset>]]\n"
          "  Add <filename> to the topology list.\n");
  ParmFile::ReadOptions();
}
static void Help_ParmInfo() {
  mprintf("\t[<parmindex>] [<mask>]\n"
          "  Print information on topology <parmindex> (0 by default).\n");
}

static void Help_ParmWrite() {
  mprintf("\tout <filename> [<parmindex>] [<fmt>] [nochamber]\n"
          "  Write topology <parmindex> to <filename>.\n");
  ParmFile::WriteOptions();
}

static void Help_ParmStrip() {
  mprintf("\t<mask> [<parmindex>]\n"
          "  Strip atoms in mask from topology <parmindex>.\n");
}

static void Help_ParmBox() {
  mprintf("\t[<parmindex>] [x <xval>] [y <yval>] [z <zval>]\n"
          "\t              [alpha <a>] [beta <b>] [gamma <g>] [nobox]\n"
          "  Set the specified topology box info to what is specified. If nobox is\n"
          "  specified, remove box info.\n");
}

static void Help_Solvent() {
  mprintf("\t[<parmindex>] { <mask> | none }\n"
          "  Set solvent for the specified topology (default 0) based on <mask>.\n"
          "  If 'none' specified, remove all solvent information.\n");
}

static void Help_AtomInfo() {
  mprintf("\t[<parmindex>] [<mask>]\n"
          "  Print information on atoms in <mask> for topology <parmindex> (0 by default).\n");
}

static void Help_BondInfo() {
  mprintf("\t[<parmindex>] [<mask>]\n"
          "  Print bond information of atoms in <mask> for topology <parmindex> (0 by default).\n");
}

static void Help_AngleInfo() {
  mprintf("\t[<parmindex>] [<mask>]\n"
      "  Print angle information of atoms in <mask> for topology <parmindex> (0 by default).\n");
}

static void Help_DihedralInfo() {
  mprintf("\t[<parmindex>] [<mask>]\n"
      "  Print dihedral information of atoms in <mask> for topology <parmindex> (0 by default).\n");
}

static void Help_ChargeInfo() {
  mprintf("\t[<parmindex>] <mask>\n"
          "  Print the total charge of atoms in <mask> for topology <parmindex> (0 by default).\n");
}

static void Help_MassInfo() {
  mprintf("\t[<parmindex>] <mask>\n"
          "  Print the total mass of atoms in <mask> for topology <parmindex> (0 by default).\n");
}

static void Help_ResInfo() {
  mprintf("\t[<parmindex>] [<mask>] [short]\n"
          "  Print information for residues in <mask> for topology <parmindex> (0 by default).\n");
}
static void Help_MolInfo() {
  mprintf("\t[<parmindex>] [<mask>]\n"
          "  Print information for molecules in <mask> for topology <parmindex> (0 by default).\n");
}

static void Help_LoadCrd() {
  mprintf("\t<filename> %s [<trajin args>] [<name>]\n", TopologyList::ParmArgs);
  mprintf("  Load trajectory <filename> as a COORDS data set named <name> (default <filename>).\n");
}

static void Help_LoadTraj() {
  mprintf("\tname <setname> [<filename>]\n"
          "  Create/add to TRAJ data set named <setname>. If no <filename> given, convert\n"
          "  currently loaded input trajectories to TRAJ data set; otherwise add <filename>\n"
          "  to TRAJ data set <setname>\n");
}

static void Help_CrdAction() {
  mprintf("\t<crd set> <actioncmd> [<action args>] [crdframes <start>,<stop>,<offset>]\n"
          "  Perform action <actioncmd> on COORDS data set <crd set>.\n");
}

static void Help_CrdOut() {
  mprintf("\t<crd set> <filename> [<trajout args>] [crdframes <start>,<stop>,<offset>]\n"
          "  Write COORDS data set <crd set> to trajectory file <filename>\n");
}

static void Help_RunAnalysis() {
  mprintf("\t[<analysis> [<analysis args>]]\n"
          "  If specified alone, run all analyses in the analysis list.\n"
          "  Otherwise run the specified analysis immediately.\n");
}

// ---------- Information on Deprecated commands -------------------------------
static void Deprecate_MinDist() {
  mprinterr("  Use the 'nativecontacts' action instead.\n");
}

static void Deprecate_Hbond() {
  mprinterr("  Hydrogen bond acceptors and donors are defined within the 'hbond' action.\n");
}

static void Deprecate_TopSearch() {
  mprinterr("  Bonds and/or molecules are automatically searched for if needed.\n");
}

static void Deprecate_ParmBondInfo() {
  mprinterr("  Use bonds, bondinfo, or printbonds instead.\n");
}

static void Deprecate_ParmResInfo() {
  mprinterr("  Use resinfo instead.\n");
}

static void Deprecate_ParmMolInfo() {
  mprinterr("  Use molinfo instead.\n");
}

static void Deprecate_AvgCoord() {
  mprinterr("  Use 'vector center' (optionally with keyword 'magnitude') instead.\n");
}

// ---------- GENERAL COMMANDS -------------------------------------------------
/// Set active reference for distance-based masks etc.
Command::RetType ActiveRef(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  ReferenceFrame REF = State.DSL()->GetReferenceFrame( argIn );
  if (!REF.error() && REF.empty()) {
    ArgList singleIntArg(integerToString(argIn.getNextInteger(0)));
    REF = State.DSL()->GetReferenceFrame( singleIntArg );
  }
  if (REF.error() || REF.empty()) return Command::C_ERR;
  State.SetActiveReference( REF.RefPtr() );
  return Command::C_OK; 
}

/// Clear data in specified lists
Command::RetType ClearList(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  return (Command::RetType)State.ClearList( argIn );
}

Command::RetType RemoveData(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  return (Command::RetType)State.RemoveDataSet( argIn );
}

/// Set debug value for specified list(s)
Command::RetType SetListDebug(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  return (Command::RetType)State.SetListDebug( argIn );
}

/// List all members of specified list(s)
Command::RetType ListAll(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  return (Command::RetType)State.ListAll( argIn );
}

static void Help_SilenceActions() { mprintf("Silence Actions Init/Setup output.\n"); }
/// Silence Actions Init/Setup output.
Command::RetType SilenceActions(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{ State.SetActionSilence( true ); return Command::C_OK; }

/// Perform action on given COORDS dataset
Command::RetType CrdAction(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  std::string setname = argIn.GetStringNext();
  if (setname.empty()) {
    mprinterr("Error: %s: Specify COORDS dataset name.\n", argIn.Command());
    return Command::C_ERR;
  }
  DataSet_Coords* CRD = (DataSet_Coords*)State.DSL()->FindCoordsSet( setname );
  if (CRD == 0) {
    mprinterr("Error: %s: No COORDS set with name %s found.\n", argIn.Command(), setname.c_str());
    return Command::C_ERR;
  }
  mprintf("\tUsing set '%s'\n", CRD->Legend().c_str());
  Timer total_time;
  total_time.Start();
  // Start, stop, offset
  int start, stop, offset;
  ArgList crdarg( argIn.GetStringKey("crdframes"), "," );
  if (Trajin::CheckFrameArgs(crdarg, CRD->Size(), start, stop, offset)) return Command::C_ERR;
  if (State.Debug() > 0) mprintf("\tDBG: Frames %i to %i, offset %i\n", start+1, stop, offset);
  ArgList actionargs = argIn.RemainingArgs();
  actionargs.MarkArg(0);
  Command::TokenPtr tkn = Command::SearchTokenType( Command::ACTION, actionargs);
  if ( tkn == 0 ) return Command::C_ERR;
  Action* act = (Action*)tkn->Alloc();
  if (act == 0) return Command::C_ERR;
  if ( act->Init( actionargs, State.PFL(), State.DSL(), State.DFL(), State.Debug() ) != Action::OK ) {
    delete act;
    return Command::C_ERR;
  }
  actionargs.CheckForMoreArgs();
  // Set up frame and parm for COORDS.
  Topology* originalParm = new Topology();
  *originalParm = CRD->Top();
  Frame* originalFrame = new Frame( CRD->Top().Atoms() );
  // Set up for this topology
  Topology* currentParm = originalParm;
  if ( act->Setup( currentParm, &currentParm ) == Action::ERR ) {
    delete act;
    return Command::C_ERR;
  }
  // Loop over all frames in COORDS.
  ProgressBar progress( stop - start );
  int set = 0;
  for (int frame = start; frame < stop; frame += offset) {
    progress.Update( set );
    CRD->GetFrame( frame, *originalFrame );
    Frame* currentFrame = originalFrame;
    if (act->DoAction( set, currentFrame, &currentFrame ) == Action::ERR) {
      mprinterr("Error: crdaction: Frame %i, set %i\n", frame + 1, set + 1);
      break;
    }
    // Check if frame was modified. If so, update COORDS.
    // TODO: Have actions indicate whether they will modify coords
    //if ( currentFrame != originalFrame ) 
      CRD->SetCRD( frame, *currentFrame );
    set++;
  }
  // Check if parm was modified. If so, update COORDS.
  if ( currentParm != originalParm ) {
    mprintf("Info: crdaction: Parm for %s was modified by action %s\n",
            CRD->Legend().c_str(), actionargs.Command());
    CRD->SetTopology( *currentParm );
  }
  act->Print();
  State.MasterDataFileWrite();
  delete originalFrame;
  delete originalParm;
  delete act;
  total_time.Stop();
  mprintf("TIME: Total action execution time: %.4f seconds.\n", total_time.Total());
  return Command::C_OK;
}

/// Write out COORDS dataset
Command::RetType CrdOut(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  std::string setname = argIn.GetStringNext();
  if (setname.empty()) {
    mprinterr("Error: crdout: Specify COORDS dataset name.\n");
    return Command::C_ERR;
  }
  DataSet_Coords* CRD = (DataSet_Coords*)State.DSL()->FindCoordsSet( setname );
  if (CRD == 0) {
    mprinterr("Error: crdout: No COORDS set with name %s found.\n", setname.c_str());
    return Command::C_ERR;
  }
  mprintf("\tUsing set '%s'\n", CRD->Legend().c_str());
  setname = argIn.GetStringNext();
  // Start, stop, offset
  int start, stop, offset;
  ArgList crdarg( argIn.GetStringKey("crdframes"), "," );
  if (Trajin::CheckFrameArgs(crdarg, CRD->Size(), start, stop, offset)) return Command::C_ERR;
  if (State.Debug() > 0) mprintf("\tDBG: Frames %i to %i, offset %i\n", start+1, stop, offset);
  Trajout outtraj;
  Topology* currentParm = (Topology*)&(CRD->Top()); // TODO: Fix cast
  if (outtraj.InitTrajWrite( setname, argIn, currentParm, TrajectoryFile::UNKNOWN_TRAJ)) {
    mprinterr("Error: crdout: Could not set up output trajectory.\n");
    return Command::C_ERR;
  }
  outtraj.PrintInfo( 1 );
  Frame currentFrame = CRD->AllocateFrame(); 
  ProgressBar progress( stop );
  for (int frame = start; frame < stop; frame += offset) {
    progress.Update( frame );
    CRD->GetFrame( frame, currentFrame );
    if ( outtraj.WriteFrame( frame, currentParm, currentFrame ) ) {
      mprinterr("Error writing %s to output trajectory, frame %i.\n",
                CRD->Legend().c_str(), frame + 1);
      break;
    }
  }
  return Command::C_OK;
}

/// Load single trajectory as DataSet_Coords
Command::RetType LoadCrd(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  // Get parm
  Topology* parm = State.PFL()->GetParm( argIn );
  if (parm == 0) {
    mprinterr("Error: loadcrd: No parm files loaded.\n");
    return Command::C_ERR;
  }
  // Load trajectory
  Trajin_Single trajin;
  trajin.SetDebug( State.Debug() );
  if (trajin.SetupTrajRead(argIn.GetStringNext(), argIn, parm)) {
    mprinterr("Error: loadcrd: Could not set up input trajectory.\n");
    return Command::C_ERR;
  }
  // Create input frame
  Frame frameIn;
  frameIn.SetupFrameV(parm->Atoms(), trajin.TrajCoordInfo());
  // Get output set name; use base file name as set name if none specified. 
  // NOTE: Default name should NEVER get used.
  std::string setname = argIn.GetStringNext();
  if (setname.empty())
    setname = trajin.TrajFilename().Base();
  // Check if set already present
  DataSet_Coords* coords = (DataSet_Coords*)State.DSL()->FindSetOfType(setname, DataSet::COORDS);
  if (coords == 0) {
    // Create Set 
    coords = (DataSet_Coords*)State.DSL()->AddSet(DataSet::COORDS, setname, "__DCRD__");
    if (coords == 0) {
      mprinterr("Error: loadcrd: Could not set up COORDS data set.\n");
      return Command::C_ERR;
    }
    coords->SetTopology( *parm );
    mprintf("\tLoading trajectory '%s' as '%s'\n", trajin.TrajFilename().full(), setname.c_str());
  } else {
    // Check that topology matches. For now just check # atoms.
    if (parm->Natom() != coords->Top().Natom()) {
      mprinterr("Error: Trajectory '%s' # atoms %i does not match COORDS data set '%s' (%i)\n",
                trajin.TrajFilename().full(), parm->Natom(),
                coords->Legend().c_str(), coords->Top().Natom());
      return Command::C_ERR;
    }
    mprintf("\tAppending trajectory '%s' to COORDS data set '%s'\n", 
            trajin.TrajFilename().full(), coords->Legend().c_str());
  }
  // Read trajectory
  trajin.BeginTraj(true);
  trajin.PrintInfoLine();
  while (trajin.GetNextFrame( frameIn ))
    coords->AddFrame( frameIn );
  trajin.EndTraj();
  return Command::C_OK;
}

/// Convert input traj list to TRAJ data set
Command::RetType LoadTraj(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  // Get Keywords
  std::string setname = argIn.GetStringKey("name");
  if (setname.empty()) {
    mprinterr("Error: Must provide data set name ('name <setname>')\n");
    return Command::C_ERR;
  }
  DataSet_Coords_TRJ* trj = (DataSet_Coords_TRJ*)
                            State.DSL()->FindSetOfType(setname, DataSet::TRAJ);
  if (trj == 0)
    trj = (DataSet_Coords_TRJ*)
          State.DSL()->AddSet(DataSet::TRAJ, setname, "__DTRJ__");
  if (trj == 0) {
    mprinterr("Error: Could not set up TRAJ data set.\n");
    return Command::C_ERR;
  }
  std::string trajname = argIn.GetStringNext();
  if (trajname.empty()) {
    // Add all existing input trajectories
    if (State.InputTrajList().empty()) {
      mprinterr("Error: No input trajectories loaded.\n");
      return Command::C_ERR;
    }
    if (State.InputTrajList().Mode() != TrajinList::NORMAL) {
      mprinterr("Error: Cannot convert ensemble input trajectories to data.\n");
      return Command::C_ERR;
    }
    mprintf("\tSaving currently loaded input trajectories as data set with name '%s'\n",
            setname.c_str());
    for (TrajinList::const_iterator Trajin = State.InputTrajList().begin();
                                    Trajin != State.InputTrajList().end(); ++Trajin)
      if (trj->AddInputTraj( *Trajin )) return Command::C_ERR;
    // TODO: Clear input trajectories from trajinList?
  } else {
    // Add the named trajectory
    if (trj->AddSingleTrajin( trajname, argIn, State.PFL()->GetParm(argIn) ))
      return Command::C_ERR;
  }
  return Command::C_OK;
}
// -----------------------------------------------------------------------------
static void Help_CombineCoords() {
  mprintf("\t<crd1> <crd2> ... [parmname <topname>] [crdname <crdname>]\n"
          "  Combined two COORDS data sets.\n");
}

static inline void CombinedCoords_AddBondArray(Topology* top, BondArray const& barray,
                                               int atomOffset)
{
  for (BondArray::const_iterator bond = barray.begin(); bond != barray.end(); ++bond)
  {
    //mprintf("DBG:\t\tBonding %i and %i\n", bond->A1() + atomOffset + 1, bond->A2() + atomOffset + 1);
    top->AddBond( bond->A1() + atomOffset, bond->A2() + atomOffset );
  }
}

/// Combine two COORDS DataSets
Command::RetType CombineCoords(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  std::string parmname = argIn.GetStringKey("parmname");
  std::string crdname  = argIn.GetStringKey("crdname");
  // Get COORDS DataSets.
  std::vector<DataSet_Coords*> CRD;
  std::string setname = argIn.GetStringNext();
  while (!setname.empty()) {
    DataSet_Coords* ds = (DataSet_Coords*)State.DSL()->FindCoordsSet( setname );
    if (ds == 0) {
      mprinterr("Error: %s: No COORDS set with name %s found.\n", argIn.Command(), setname.c_str());
      return Command::C_ERR;
    }
    CRD.push_back( ds );
    setname = argIn.GetStringNext();
  }
  if (CRD.size() < 2) {
    mprinterr("Error: %s: Must specify at least 2 COORDS data sets\n", argIn.Command());
    return Command::C_ERR;
  }
  Topology* CombinedTop = new Topology();
  if (CombinedTop == 0) return Command::C_ERR;
  if (parmname.empty())
    parmname = CRD[0]->Top().ParmName() + "_" + CRD[1]->Top().ParmName();
  CombinedTop->SetParmName( parmname, FileName() );
  // TODO: Check Parm box info.
  size_t minSize = CRD[0]->Size();
  for (unsigned int setnum = 0; setnum != CRD.size(); ++setnum) {
    if (CRD[setnum]->Size() < minSize)
      minSize = CRD[setnum]->Size();
    CombinedTop->AppendTop( CRD[setnum]->Top() );
  }
  CombinedTop->Brief("Combined parm:");
  State.PFL()->AddParm( CombinedTop );
  // Combine coordinates
  if (crdname.empty())
    crdname = CRD[0]->Legend() + "_" + CRD[1]->Legend();
  mprintf("\tCombining %zu frames from each set into %s\n", minSize, crdname.c_str());
  DataSet_Coords* CombinedCrd = (DataSet_Coords*)State.DSL()->AddSet(DataSet::COORDS, crdname, "CRD");
  if (CombinedCrd == 0) {
    mprinterr("Error: Could not create COORDS data set.\n");
    return Command::C_ERR;
  }
  CombinedCrd->SetTopology( *CombinedTop );
  // FIXME: Only copying coords for now
  Frame CombinedFrame( CombinedTop->Natom() * 3 );
  std::vector<Frame> InputFrames;
  for (unsigned int setnum = 0; setnum != CRD.size(); ++setnum)
    InputFrames.push_back( CRD[setnum]->AllocateFrame() );
  for (size_t nf = 0; nf != minSize; ++nf) {
    CombinedFrame.ClearAtoms();
    for (unsigned int setnum = 0; setnum != CRD.size(); ++setnum)
    {
      CRD[setnum]->GetFrame( nf, InputFrames[setnum] );
      for (int atnum = 0; atnum < CRD[setnum]->Top().Natom(); atnum++)
        CombinedFrame.AddXYZ( InputFrames[setnum].XYZ(atnum) );
    }
    CombinedCrd->AddFrame( CombinedFrame );
  }
/* FIXME: This code is fast but only works for DataSet_Coords_CRD
  Frame::CRDtype CombinedFrame( CombinedTop->Natom() * 3 );
  for (size_t nf = 0; nf != minSize; ++nf) {
    size_t offset = 0;
    for (unsigned int setnum = 0; setnum != CRD.size(); ++setnum)
    {
      size_t crd_offset = (size_t)CRD[setnum]->Top().Natom() * 3;
      std::copy( CRD[setnum]->CRD(nf).begin(), CRD[setnum]->CRD(nf).begin() + crd_offset,
                 CombinedFrame.begin() + offset );
      offset += crd_offset;
    }
    CombinedCrd->AddCRD( CombinedFrame );
  }
*/
  return Command::C_OK;
}

// -----------------------------------------------------------------------------
static void Help_GenerateAmberRst() {
  mprintf("\t<mask1> <mask2> [<mask3>] [<mask4>]\n"
          "\tr1 <r1> r2 <r2> r3 <r3> r4 <r4> rk2 <rk2> rk3 <rk3>\n"
          "\t{%s}\n"
          "\t[{%s} [offset <off>] [width <width>]\n"
          "\t[out <outfile>] [overwrite]\n"
          "  Generate Amber-format restraint from 2 or more mask expressions.\n",
          TopologyList::ParmArgs, DataSetList::RefArgs);
}

/// Generate amber restraints from given masks.
Command::RetType GenerateAmberRst(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  // Get parm
  Topology* parm = State.PFL()->GetParm( argIn );
  if (parm == 0) {
    mprinterr("Error: No parm files loaded.\n");
    return Command::C_ERR;
  }
  // Get optional reference coords
  ReferenceFrame RefCrd = State.DSL()->GetReferenceFrame(argIn);
  // Get arguments
  bool overwrite = argIn.hasKey("overwrite");
  double r1 = argIn.getKeyDouble("r1", 0.0);
  double r2 = argIn.getKeyDouble("r2", 0.0);
  double r3 = argIn.getKeyDouble("r3", 0.0);
  double r4 = argIn.getKeyDouble("r4", 0.0);
  double rk2 = argIn.getKeyDouble("rk2", 0.0);
  double rk3 = argIn.getKeyDouble("rk3", 0.0);
  // crddist will be !RefCrd.empty()
  double offset = argIn.getKeyDouble("offset", 0.0);
  double width = argIn.getKeyDouble("width", 0.5);
  std::string outName = argIn.GetStringKey("out");
  CpptrajFile outfile;
  int err = 0;
  if (overwrite)
    err = outfile.OpenWrite( outName );
  else
    err = outfile.OpenAppend( outName );
  if (err != 0) {
    mprinterr("Error: Could not open output file.\n");
    return Command::C_ERR;
  }
  // TODO: else if (strcmp(pch,"nobb")==0) nobb=1;
  // Assume everything else is a mask
  std::vector<AtomMask> rstMasks;
  std::string maskExpr = argIn.GetMaskNext();
  while (!maskExpr.empty()) {
    if ( rstMasks.size() >= 4 )
      mprintf("Warning: 4 masks already defined. Skipping '%s'\n", maskExpr.c_str());
    else {
      AtomMask tmpmask( maskExpr );
      int maskerr = 0;
      if (!RefCrd.empty())
        maskerr = parm->SetupIntegerMask( tmpmask, RefCrd.Coord() );
      else
        maskerr = parm->SetupIntegerMask( tmpmask );
      if ( maskerr != 0 ) {
        mprinterr("Error: Could not set up mask '%s'\n", tmpmask.MaskString());
        return Command::C_ERR;
      }
      tmpmask.MaskInfo();
      if ( tmpmask.None() ) {
        mprinterr("Error: '%s' corresponds to no atoms.\n", tmpmask.MaskString());
        return Command::C_ERR;
      }
      rstMasks.push_back( tmpmask );
    }
    maskExpr = argIn.GetMaskNext();
  }
  if (State.Debug() > 0) mprintf("\tDefined %zu rst masks.\n", rstMasks.size());
  if (rstMasks.size() < 2) {
    mprinterr("Error: Must specify at least 2 masks for restraint.\n");
    return Command::C_ERR;
  }
  // TODO: Remove backbone atoms for 'nobb'?
  // If a reference frame was specified and distance restraint, use center of
  // mass distance/angle/torsion between masks as r2.
  if ( !RefCrd.empty() ) {
    if ( RefCrd.Parm().Pindex() != parm->Pindex() )
      mprintf("Warning: Reference topology does not match specified topology.\n");
    Vec3 a1 = RefCrd.Coord().VCenterOfMass( rstMasks[0] );
    Vec3 a2 = RefCrd.Coord().VCenterOfMass( rstMasks[1] );
    if (rstMasks.size() == 2)
      r2 = DIST_NoImage( a1, a2 );
    else if (rstMasks.size() == 3) {
      Vec3 a3 = RefCrd.Coord().VCenterOfMass( rstMasks[2] );
      r2 = CalcAngle(a1.Dptr(), a2.Dptr(), a3.Dptr()) * Constants::RADDEG;
    } else if (rstMasks.size() == 4) {
      Vec3 a3 = RefCrd.Coord().VCenterOfMass( rstMasks[2] );
      Vec3 a4 = RefCrd.Coord().VCenterOfMass( rstMasks[3] );
      r2 = Torsion(a1.Dptr(), a2.Dptr(), a3.Dptr(), a4.Dptr()) * Constants::RADDEG;
    }
    r2 += offset;
    r3 = r2;
    r1 = r2 - width;
    r4 = r3 + width;
    mprintf("\tCoM value from ref will be used, r1=%f, r2=%f, r3=%f, r4=%f\n", r1,r2,r3,r4);
  } 
  // Print restraint header 
  outfile.Printf(" &rst iat=");
  for (std::vector<AtomMask>::const_iterator M = rstMasks.begin();
                                             M != rstMasks.end(); ++M)
  {
    if ((*M).Nselected() == 1)
      outfile.Printf("%i,", (*M)[0] + 1);
    else
      outfile.Printf("-1,");
  }
  outfile.Printf("0\n");
  // Print Restraint boundaries and constants
  outfile.Printf("   r1=%f, r2=%f, r3=%f, r4=%f, rk2=%f, rk3=%f,\n",
                 r1, r2, r3, r4, rk2, rk3);
  // Print out atom groups if necessary
  unsigned int group = 1;
  for (std::vector<AtomMask>::const_iterator M = rstMasks.begin();
                                             M != rstMasks.end(); ++M, group++)
  {
    if ((*M).Nselected() > 1) {
      outfile.Printf("   ");
      unsigned int j = 1;
      for (AtomMask::const_iterator atom = (*M).begin();
                                    atom != (*M).end(); ++atom, j++)
        outfile.Printf("IGR%u(%u)=%i,", group, j, (*atom) + 1);
      outfile.Printf("\n");
    }
  }
  // Finish restraint
  outfile.Printf("   nstep1=0, nstep2=0,\n &end\n");
  outfile.CloseFile();
  return Command::C_OK;
}

// -----------------------------------------------------------------------------
/// Add DataSets specified by arguments to given DataFile.
// NOTE: Used byt Create_DataFile and Write_DataFile
// TODO: Put in DataFile?
static int AddSetsToDataFile(DataFile& df, ArgList const& dsetArgs, DataSetList& DSL)
{
  int err = 0;
  for (ArgList::const_iterator dsa = dsetArgs.begin(); dsa != dsetArgs.end(); ++dsa) {
    DataSetList Sets = DSL.GetMultipleSets( *dsa );
    if (Sets.empty())
      mprintf("Warning: %s does not correspond to any data sets.\n", (*dsa).c_str());
    for (DataSetList::const_iterator set = Sets.begin(); set != Sets.end(); ++set) {
      mprintf(" %s", (*set)->Legend().c_str());
      if ( df.AddSet(*set) ) {
        mprinterr("Error: Could not add data set %s to file.\n", (*set)->Legend().c_str());
        ++err;
      }
    }
  }
  mprintf("\n");
  return err;
}

/// Add a new DataFile to DFL with specified DataSets, to be written later.
Command::RetType Create_DataFile(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  // Next string is datafile that command pertains to.
  std::string name1 = argIn.GetStringNext();
  if (name1.empty()) {
    mprinterr("Error: No filename given.\n");
    return Command::C_ERR;
  }
  DataFile* df = State.DFL()->AddDataFile(name1, argIn);
  if (df == 0) return Command::C_ERR;
  return (Command::RetType)( AddSetsToDataFile(*df, argIn.RemainingArgs(), *(State.DSL())) );
}

/// Write DataFile with specified DataSets immediately, or force write of all DataFiles in State
Command::RetType Write_DataFile(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  // Next string is datafile that command pertains to.
  std::string name1 = argIn.GetStringNext();
  if (name1.empty()) {
    State.DFL()->ResetWriteStatus();
    State.MasterDataFileWrite();
    return Command::C_OK;
  }
  DataFile* df = new DataFile();
  if (df == 0) return Command::C_ERR;
  if (df->SetupDatafile( name1, argIn, State.Debug() )) {
    delete df;
    return Command::C_ERR;
  }
  mprintf("\tWriting sets to %s, format '%s'\n", df->DataFilename().full(), df->FormatString());
  int err = AddSetsToDataFile(*df, argIn.RemainingArgs(), *(State.DSL()));
  if (err == 0) df->WriteData();
  delete df;
  return (Command::RetType)err;
}

/// Process DataFile-specific command
Command::RetType DataFileCmd(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  return (Command::RetType)( State.DFL()->ProcessDataFileArgs( argIn ) );
}

/// Process DataSet-specific command
Command::RetType DataSetCmd(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  if (argIn.Contains("legend")) { // Set legend for one data set
    std::string legend = argIn.GetStringKey("legend");
    DataSet* ds = State.DSL()->GetDataSet( argIn.GetStringNext() );
    if (ds == 0) return Command::C_ERR;
    mprintf("\tChanging legend '%s' to '%s'\n", ds->Legend().c_str(), legend.c_str());
    ds->SetLegend( legend );
    return Command::C_OK;
  }
  // Change mode/type for one or more sets.
  std::string modeKey = argIn.GetStringKey("mode");
  std::string typeKey = argIn.GetStringKey("type");
  if (modeKey.empty() && typeKey.empty()) {
    mprinterr("Error: No valid keywords specified.\n");
    return Command::C_ERR;
  }
  // First determine mode if specified.
  DataSet::scalarMode dmode = DataSet::UNKNOWN_MODE;
  if (!modeKey.empty()) {
    dmode = DataSet::ModeFromKeyword( modeKey );
    if (dmode == DataSet::UNKNOWN_MODE) {
      mprinterr("Error: Invalid mode keyword '%s'\n", modeKey.c_str());
      return Command::C_ERR;
    }
  }
  // Next determine type if specified.
  DataSet::scalarType dtype = DataSet::UNDEFINED;
  if (!typeKey.empty()) {
    dtype = DataSet::TypeFromKeyword( typeKey, dmode );
    if (dtype == DataSet::UNDEFINED) {
      mprinterr("Error: Invalid type keyword '%s'\n", typeKey.c_str());
      return Command::C_ERR;
    }
  }
  // Additional options for type 'noe'
  double noe_lbound=0.0, noe_ubound=0.0, noe_rexp = -1.0;
  if (dtype == DataSet::NOE) {
    if (Action_Distance::NOE_Args(argIn, noe_lbound, noe_ubound, noe_rexp))
      return Command::C_ERR;
  }
  if (dmode != DataSet::UNKNOWN_MODE)
    mprintf("\tDataSet mode = %s\n", DataSet::Smodes[dmode]);
  if (dtype != DataSet::UNDEFINED)
    mprintf("\tDataSet type = %s\n", DataSet::Stypes[dtype]);
  // Loop over all DataSet arguments 
  std::string ds_arg = argIn.GetStringNext();
  while (!ds_arg.empty()) {
    DataSetList dsl = State.DSL()->GetMultipleSets( ds_arg );
    for (DataSetList::const_iterator ds = dsl.begin(); ds != dsl.end(); ++ds)
    {
      if ( (*ds)->Ndim() != 1 )
        mprintf("Warning:\t\t'%s': Can only set mode/type for 1D data sets.\n",
                (*ds)->Legend().c_str());
      else {
        if ( dtype == DataSet::NOE ) {
          if ( (*ds)->Type() != DataSet::DOUBLE )
            mprintf("Warning:\t\t'%s': Can only set NOE parameters for 'double' data sets.\n",
                (*ds)->Legend().c_str());
          else
            ((DataSet_double*)(*ds))->SetNOE(noe_lbound, noe_ubound, noe_rexp);
        }
        mprintf("\t\t'%s'\n", (*ds)->Legend().c_str());
        (*ds)->SetScalar( dmode, dtype );
      }
    }
    ds_arg = argIn.GetStringNext();
  }
  return Command::C_OK;
}

// -----------------------------------------------------------------------------
static void Help_DataFilter() {
  mprintf("\t<dataset arg> min <min> max <max> [out <file> [name <setname>]]\n"
          "  Create a data set (optionally named <setname>) containing 1 for\n"
          "  data within given <min> and <max> criteria for each specified\n"
          "  data set. There must be at least one <min> and <max> argument,\n"
          "  and can be as many as there are specified data sets.\n");
}

/// Use the filter command on DataSets outside trajectory processing.
Command::RetType DataFilter(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  Action_FilterByData filterAction;
  if (filterAction.Init(argIn, State.PFL(), State.DSL(), State.DFL(), State.Debug()) != Action::OK)
    return Command::C_ERR;
  size_t nframes = filterAction.DetermineFrames();
  if (nframes < 1) {
    mprinterr("Error: No data to filter. All sets must contain some data.\n");
    return Command::C_ERR;
  }
  ProgressBar progress( nframes );
  for (size_t frame = 0; frame != nframes; frame++) {
    progress.Update( frame );
    filterAction.DoAction(frame, (Frame*)0, (Frame**)0); // Filter does not need frame.
  }
  // Trigger master datafile write just in case
  State.MasterDataFileWrite();
  return Command::C_OK;
}
// -----------------------------------------------------------------------------

/// Read data from file into master DataSetList
Command::RetType ReadData(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  DataFile dataIn;
  dataIn.SetDebug( State.DFL()->Debug() );
  std::string filenameIn = argIn.GetStringNext();
  StrArray fnames = ExpandToFilenames( filenameIn );
  if (fnames.empty()) {
    mprinterr("Error: '%s' matches no files.\n", filenameIn.c_str());
    return Command::C_ERR;
  }
  int err = 0;
  for (StrArray::const_iterator fn = fnames.begin(); fn != fnames.end(); ++fn) {
    if (dataIn.ReadDataIn( *fn, argIn, *State.DSL() )!=0) {
      mprinterr("Error: Could not read data file '%s'.\n", fn->c_str());
      err++;
    }
  }
  if (err > 0) return Command::C_ERR;
  return Command::C_OK;
}

/// Exit
Command::RetType Quit(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  return Command::C_QUIT;
}

/// Run a system command
Command::RetType SystemCmd(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  system( argIn.ArgLine() );
  return Command::C_OK;
}

/// Find help for command/topic
Command::RetType Help(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  ArgList arg = argIn;
  arg.RemoveFirstArg();
  if (arg.empty())
    // NONE in this context means list all commands
    Command::ListCommands(Command::NONE);
  else if (arg.CommandIs("General"))
    Command::ListCommands(Command::GENERAL);
  else if (arg.CommandIs("Topology"))
    Command::ListCommands(Command::PARM);
  else if (arg.CommandIs("Action"))
    Command::ListCommands(Command::ACTION);
  else if (arg.CommandIs("Analysis"))
    Command::ListCommands(Command::ANALYSIS);
  else if (arg.CommandIs("Trajectory"))
    Command::ListCommands(Command::TRAJ);
  else {
    Command::TokenPtr dispatchToken = Command::SearchToken( arg );
    if (dispatchToken == 0 || dispatchToken->Help == 0)
      mprinterr("No help found for %s\n", arg.Command());
    else
      dispatchToken->Help();
  }
  return Command::C_OK;
}

/// Run the current State
Command::RetType RunState(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  return (Command::RetType)State.Run();
}

/// Read input from a file.
Command::RetType ReadInput(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  // Next arg should be a filename. Not allowed to be blank in command.
  std::string inputFilename = argIn.GetStringNext();
  if (inputFilename.empty()) {
    mprinterr("Error: No input filename given.\n");
    return Command::C_ERR;
  }
  return Command::ProcessInput(State, inputFilename);
}

/// Tell CpptrajState to ignore errors if possible
Command::RetType NoExitOnError(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  State.SetNoExitOnError();
  mprintf("\tAttempting to ignore errors if possible.\n");
  return Command::C_OK;
}

/// Tell CpptrajState not to use a progress bar during Run.
Command::RetType NoProgress(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  State.SetNoProgress();
  mprintf("\tProgress bar will not be used during Run.\n");
  return Command::C_OK;
}

///  Set precision for specific set or all sets in specified DataFile
Command::RetType Precision(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  // Next string is DataSet(s)/DataFile that command pertains to.
  std::string name1 = argIn.GetStringNext();
  if (name1.empty()) {
    mprinterr("Error: No filename/setname given.\n");
    return Command::C_ERR;
  }
  // This will break if dataset name starts with a digit...
  int width = argIn.getNextInteger(12);
  if (width < 1) {
    mprintf("Error: Cannot set width < 1 (%i).\n", width);
    return Command::C_ERR;
  }
  int precision = argIn.getNextInteger(4);
  if (precision < 0) precision = 0;
  DataFile* df = State.DFL()->GetDataFile(name1);
  if (df != 0) {
    mprintf("\tSetting precision for all sets in %s to %i.%i\n", df->DataFilename().base(),
            width, precision);
    df->SetDataFilePrecision(width, precision);
  } else {
    State.DSL()->SetPrecisionOfDataSets( name1, width, precision );
  }
  return Command::C_OK;
}

/// Run specified analysis or all analyses in State.
Command::RetType RunAnalysis(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  // If only 1 arg (the command) run all analyses in list
  if (argIn.Nargs() == 1) {
    int eval = State.RunAnalyses();
    State.MasterDataFileWrite();
    if (eval == 0)
      return Command::C_OK;
    else
      return Command::C_ERR;
  }
  // Run specified analysis
  // FIXME: Use RemoveFirstArg
  ArgList analyzeargs = argIn.RemainingArgs();
  analyzeargs.MarkArg(0);
  Command::TokenPtr tkn = Command::SearchTokenType( Command::ANALYSIS, analyzeargs );
  if ( tkn == 0 ) return Command::C_ERR;
  Analysis* ana = (Analysis*)tkn->Alloc();
  if (ana == 0) return Command::C_ERR;
  Timer total_time;
  total_time.Start();
  Command::RetType err = Command::C_ERR;
  if ( ana->Setup( analyzeargs, State.DSL(), State.PFL(), State.DFL(), State.Debug() ) == 
                   Analysis::OK )
  {
    analyzeargs.CheckForMoreArgs();
    if (ana->Analyze() != Analysis::ERR) {
      err = Command::C_OK;
      State.MasterDataFileWrite();
    }
  }
  delete ana;
  total_time.Stop();
  mprintf("TIME: Total analysis execution time: %.4f seconds.\n", total_time.Total());
  return err;
}

/// Show results of mask expression
Command::RetType SelectAtoms(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  AtomMask tempMask( argIn.GetMaskNext() );
  Topology* parm = State.PFL()->GetParmByIndex( argIn );
  if (parm == 0) return Command::C_ERR;
  if (parm->SetupIntegerMask( tempMask )) return Command::C_ERR;
  mprintf("Selected %i atoms.\n", tempMask.Nselected());
  if (!argIn.hasKey("total"))
    tempMask.PrintMaskAtoms("Selected");
  return Command::C_OK;
}

/// Show results of DataSet expression
Command::RetType SelectDataSets(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  std::string dsarg = argIn.GetStringNext();
  DataSetList dsets = State.DSL()->GetMultipleSets( dsarg );
  if (!dsets.empty()) {
    mprintf("SelectDS: Arg '%s':", dsarg.c_str());
    dsets.List();
  }
  return Command::C_OK;
}

// -----------------------------------------------------------------------------
static void Help_PrintData() {
  mprintf("\t<data set>\n"
          "  Print data from data set to screen.\n");
}

Command::RetType PrintData(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  std::string ds_arg = argIn.GetStringNext();
  if (ds_arg.empty()) {
    mprinterr("Error: No data set arg specified.\n");
    return Command::C_ERR;
  }
  DataSet* ds = State.DSL()->GetDataSet( ds_arg );
  if (ds == 0) return Command::C_ERR;

  DataFile ToStdout;
  ToStdout.SetupStdout(argIn, State.Debug());
  ToStdout.AddSet( ds );
  ToStdout.WriteData();
  return Command::C_OK;
}

// -----------------------------------------------------------------------------
static void Help_Calc() {
  mprintf("\t<expression>\n"
          "  Evaluate the given mathematical expression.\n");
}

/// Parse a mathematical expression.
Command::RetType Calc(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  RPNcalc calc;
  calc.SetDebug( State.Debug() );
  // Do NOT include command in expression.
  if (calc.ProcessExpression( argIn.ArgString().substr(argIn[0].size()) ))
    return Command::C_ERR;
  if (calc.Evaluate(*State.DSL())) return Command::C_ERR;
  return Command::C_OK;
}
  
// ---------- TRAJECTORY COMMANDS ----------------------------------------------
/// Add output trajectory to State
Command::RetType Trajout(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  return (Command::RetType)State.AddTrajout( argIn );
}

/// Add input trajectory to State
Command::RetType Trajin(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  return (Command::RetType)State.AddTrajin( argIn, false );
}

/// Add ensemble of input trajectories to State.
Command::RetType Ensemble(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  return (Command::RetType)State.AddTrajin( argIn, true );
}

/// Add reference trajectory to State
Command::RetType Reference(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  return (Command::RetType)State.AddReference( argIn.GetStringNext(), argIn );
}

// ---------- TOPOLOGY COMMANDS ------------------------------------------------
/// Load topology to State
Command::RetType LoadParm(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  return (Command::RetType)State.PFL()->AddParmFile(argIn.GetStringNext(), argIn);
}

/// Print info for specified parm or atoms in specified parm.
Command::RetType ParmInfo(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  Topology* parm = State.PFL()->GetParmByIndex( argIn );
  if (parm == 0) return Command::C_ERR;
  parm->Summary();
  return Command::C_OK;
}

/// Print info for atoms in mask.
Command::RetType AtomInfo(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  Topology* parm = State.PFL()->GetParmByIndex( argIn );
  if (parm == 0) return Command::C_ERR;
  parm->PrintAtomInfo( argIn.GetMaskNext() );
  return Command::C_OK;
}

/// Print bond info for atoms in mask.
Command::RetType BondInfo(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  Topology* parm = State.PFL()->GetParmByIndex( argIn );
  if (parm == 0) return Command::C_ERR;
  parm->PrintBondInfo( argIn.GetMaskNext() );
  return Command::C_OK;
}

/// Print angle info for atoms in mask.
Command::RetType AngleInfo(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  Topology* parm = State.PFL()->GetParmByIndex( argIn );
  if (parm == 0) return Command::C_ERR;
  parm->PrintAngleInfo( argIn.GetMaskNext() );
  return Command::C_OK;
}

/// Print dihedral info for atoms in mask.
Command::RetType DihedralInfo(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  Topology* parm = State.PFL()->GetParmByIndex( argIn );
  if (parm == 0) return Command::C_ERR;
  parm->PrintDihedralInfo( argIn.GetMaskNext() );
  return Command::C_OK;
}

/// Print residue info for atoms in mask.
Command::RetType ResInfo(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  Topology* parm = State.PFL()->GetParmByIndex( argIn );
  if (parm == 0) return Command::C_ERR;
  bool printShort = argIn.hasKey("short");
  if (printShort)
    parm->PrintShortResInfo( argIn.GetMaskNext(), argIn.getKeyInt("maxwidth",50) );
  else
    parm->PrintResidueInfo( argIn.GetMaskNext() );
  return Command::C_OK;
}

/// Print molecule info for atoms in mask.
Command::RetType MolInfo(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  Topology* parm = State.PFL()->GetParmByIndex( argIn );
  if (parm == 0) return Command::C_ERR;
  parm->PrintMoleculeInfo( argIn.GetMaskNext() );
  return Command::C_OK;
}

/// Print the total charge of atoms in mask
Command::RetType ChargeInfo(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  Topology* parm = State.PFL()->GetParmByIndex( argIn );
  if (parm == 0) return Command::C_ERR;
  if (parm->PrintChargeMassInfo( argIn.GetMaskNext(), 0 )) return Command::C_ERR;
  return Command::C_OK;
}

/// Print the total charge of atoms in mask
Command::RetType MassInfo(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  Topology* parm = State.PFL()->GetParmByIndex( argIn );
  if (parm == 0) return Command::C_ERR;
  if (parm->PrintChargeMassInfo( argIn.GetMaskNext(), 1 )) return Command::C_ERR;
  return Command::C_OK;
}

/// Modify specified parm box info
Command::RetType ParmBox(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  Box pbox;
  bool nobox = false;
  if ( argIn.hasKey("nobox") )
    nobox = true;
  else {
    pbox.SetX( argIn.getKeyDouble("x",0) );
    pbox.SetY( argIn.getKeyDouble("y",0) );
    pbox.SetZ( argIn.getKeyDouble("z",0) );
    pbox.SetAlpha( argIn.getKeyDouble("alpha",0) );
    pbox.SetBeta(  argIn.getKeyDouble("beta",0)  );
    pbox.SetGamma( argIn.getKeyDouble("gamma",0) );
  }
  Topology* parm = State.PFL()->GetParmByIndex( argIn );
  if (parm == 0) return Command::C_ERR;
  if (nobox)
    mprintf("\tRemoving box information from parm %i:%s\n", parm->Pindex(), parm->c_str());
  else
    // Fill in missing parm box information from specified parm
    pbox.SetMissingInfo( parm->ParmBox() );
  parm->SetParmBox( pbox );
  return Command::C_OK;
}

/// Strip atoms from specified parm
Command::RetType ParmStrip(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  Topology* parm = State.PFL()->GetParmByIndex( argIn );
  if (parm == 0) return Command::C_ERR;
  // Check if this topology has already been used to set up an input
  // trajectory, as this will break the traj read.
  for (TrajinList::const_iterator tIn = State.InputTrajList().begin();
                                  tIn != State.InputTrajList().end(); ++tIn)
    if ( (*tIn)->TrajParm() == parm ) {
      mprinterr("Error: Topology '%s' has already been used to set up trajectory '%s'.\n"
                "Error:   To strip this topology use the 'strip' action.\n",
                parm->c_str(), (*tIn)->TrajFilename().full());
      return Command::C_ERR;
    }
  AtomMask tempMask( argIn.GetMaskNext() );
  // Since want to keep atoms outside mask, invert selection
  tempMask.InvertMask();
  if (parm->SetupIntegerMask( tempMask )) return Command::C_ERR;
  mprintf("\tStripping atoms in mask [%s] (%i) from %s\n",tempMask.MaskString(),
           parm->Natom() - tempMask.Nselected(), parm->c_str());
  Topology* tempParm = parm->modifyStateByMask(tempMask);
  if (tempParm==0) {
    mprinterr("Error: %s: Could not strip parm.\n", argIn.Command());
    return Command::C_ERR;
  } else {
    // Replace parm with stripped version
    *parm = *tempParm;
    parm->Brief("Stripped parm:");
    delete tempParm;
  }
  return Command::C_OK;
}

/// Write parm to Amber Topology file.
Command::RetType ParmWrite(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  std::string outfilename = argIn.GetStringKey("out");
  if (outfilename.empty()) {
    mprinterr("Error: No output filename specified (use 'out <filename>').\n");
    return Command::C_ERR;
  }
  int err = 0;
  ParmFile pfile;
  // Check if a COORDS data set was specified.
  std::string crdset = argIn.GetStringKey("crdset");
  if (crdset.empty()) {
    Topology* parm = State.PFL()->GetParmByIndex( argIn );
    if (parm == 0) return Command::C_ERR;
    err = pfile.WriteTopology( *parm, outfilename, argIn, ParmFile::UNKNOWN_PARM, State.Debug() );
  } else {
    DataSet_Coords* ds = (DataSet_Coords*)State.DSL()->FindCoordsSet(crdset);
    if (ds == 0) return Command::C_ERR;
    mprintf("\tUsing topology from data set '%s'\n", ds->Legend().c_str());
    err = pfile.WriteTopology(ds->Top(), outfilename, argIn, ParmFile::UNKNOWN_PARM, State.Debug());
  }
  if (err != 0)
    return Command::C_ERR;
  return Command::C_OK;
}

/// Modify parm solvent information
Command::RetType ParmSolvent(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  std::string maskexpr;
  if (!argIn.hasKey("none")) {
    maskexpr = argIn.GetMaskNext();
    if ( maskexpr.empty() ) {
      mprinterr("Error: solvent: No mask specified.\n");
      return Command::C_ERR;
    }
  }
  // Get parm index
  Topology* parm = State.PFL()->GetParmByIndex( argIn );
  if (parm == 0) return Command::C_ERR;
  parm->SetSolvent( maskexpr );
  return Command::C_OK;
}

/// Scale dihedral force constants in specfied parm by factor.
Command::RetType ScaleDihedralK(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  Topology* parm = State.PFL()->GetParm( argIn );
  if (parm == 0) return Command::C_ERR;
  double scale_factor = argIn.getNextDouble(1.0);
  mprintf("\tScaling dihedral force constants in %s by %f\n", parm->c_str(), scale_factor);
  parm->ScaleDihedralK( scale_factor );
  return Command::C_OK;
}

// ---------- DISPATCHABLE COMMANDS --------------------------------------------
/// Add an action to the State ActionList
Command::RetType AddAction(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  return ( (Command::RetType)State.AddAction( Alloc, argIn ) );
}

/// Add an action to the State AnalysisList
Command::RetType AddAnalysis(CpptrajState& State, ArgList& argIn, Command::AllocType Alloc)
{
  return ( (Command::RetType)State.AddAnalysis( Alloc, argIn ) );
}

// ================ LIST OF ALL COMMANDS =======================================
/** Ideally keep this array first sorted by type (1st field), then 
  * alphabetically by command string (2nd field).
  */
const Command::Token Command::Commands[] = {
  // GENERAL COMMANDS
  { GENERAL, "activeref",     0, Help_ActiveRef,       ActiveRef       },
  { GENERAL, "calc",          0, Help_Calc,            Calc            },
  { GENERAL, "clear",         0, Help_Clear,           ClearList       },
  { GENERAL, "combinecrd",    0, Help_CombineCoords,   CombineCoords   },
  { GENERAL, "crdaction",     0, Help_CrdAction,       CrdAction       },
  { GENERAL, "crdout",        0, Help_CrdOut,          CrdOut          },
  { GENERAL, "create",        0, Help_Create_DataFile, Create_DataFile },
  { GENERAL, "datafile",      0, Help_DataFile,        DataFileCmd     },
  { GENERAL, "datafilter",    0, Help_DataFilter,      DataFilter      },
  { GENERAL, "dataset",       0, Help_DataSetCmd,      DataSetCmd      },
  { GENERAL, "debug",         0, Help_Debug,           SetListDebug    },
  { GENERAL, "exit" ,         0, Help_Quit,            Quit            },
  { GENERAL, "gnuplot",       0, Help_System,          SystemCmd       },
  { GENERAL, "go",            0, Help_Run,             RunState        },
  { GENERAL, "head",          0, Help_System,          SystemCmd       },
  { GENERAL, "help",          0, Help_Help,            Help            },
  { GENERAL, "less",          0, Help_System,          SystemCmd       },
  { GENERAL, "list",          0, Help_List,            ListAll         },
  { GENERAL, "loadcrd",       0, Help_LoadCrd,         LoadCrd         },
  { GENERAL, "loadtraj",      0, Help_LoadTraj,        LoadTraj        },
  { GENERAL, "ls",            0, Help_System,          SystemCmd       },
  { GENERAL, "noexitonerror", 0, Help_NoExitOnError,   NoExitOnError   },
  { GENERAL, "noprogress",    0, Help_NoProgress,      NoProgress      },
  { GENERAL, "precision",     0, Help_Precision,       Precision       },
  { GENERAL, "printdata",     0, Help_PrintData,       PrintData       },
  { GENERAL, "prnlev",        0, Help_Debug,           SetListDebug    },
  { GENERAL, "pwd",           0, Help_System,          SystemCmd       },
  { GENERAL, "quit" ,         0, Help_Quit,            Quit            },
  { GENERAL, "readdata",      0, Help_ReadData,        ReadData        },
  { GENERAL, "readinput",     0, Help_ReadInput,       ReadInput       },
  { GENERAL, "removedata",    0, Help_RemoveData,      RemoveData      },
  { GENERAL, "rst"   ,        0, Help_GenerateAmberRst,GenerateAmberRst},
  { GENERAL, "run"   ,        0, Help_Run,             RunState        },
  { GENERAL, "runanalysis",   0, Help_RunAnalysis,     RunAnalysis     },
  { GENERAL, "select",        0, Help_Select,          SelectAtoms     },
  { GENERAL, "selectds",      0, Help_SelectDS,        SelectDataSets  },
  { GENERAL, "silenceactions",0, Help_SilenceActions,  SilenceActions  },
  { GENERAL, "write",         0, Help_Write_DataFile,  Write_DataFile  },
  { GENERAL, "writedata",     0, Help_Write_DataFile,  Write_DataFile  },
  { GENERAL, "xmgrace",       0, Help_System,          SystemCmd       },
  // TRAJECTORY COMMANDS
  { TRAJ,    "ensemble",      0, Help_Ensemble,        Ensemble        },
  { TRAJ,    "reference",     0, Help_Reference,       Reference       },
  { TRAJ,    "trajin",        0, Help_Trajin,          Trajin          },
  { TRAJ,    "trajout",       0, Help_Trajout,         Trajout         },
  // TOPOLOGY COMMANDS
  { PARM,    "angleinfo",     0, Help_AngleInfo,       AngleInfo       },
  { PARM,    "angles",        0, Help_AngleInfo,       AngleInfo       },
  { PARM,    "atominfo",      0, Help_AtomInfo,        AtomInfo        },
  { PARM,    "atoms",         0, Help_AtomInfo,        AtomInfo        },
  { PARM,    "bondinfo",      0, Help_BondInfo,        BondInfo        },
  { PARM,    "bonds",         0, Help_BondInfo,        BondInfo        },
  { PARM,    "charge",        0, Help_ChargeInfo,      ChargeInfo      },
  { PARM,    "dihedralinfo",  0, Help_DihedralInfo,    DihedralInfo    },
  { PARM,    "dihedrals",     0, Help_DihedralInfo,    DihedralInfo    },
  { PARM,    "mass",          0, Help_MassInfo,        MassInfo        },
  { PARM,    "molinfo",       0, Help_MolInfo,         MolInfo         },
  { PARM,    "parm",          0, Help_Parm,            LoadParm        },
  { PARM,    "parmbox",       0, Help_ParmBox,         ParmBox         },
  { PARM,    "parminfo",      0, Help_ParmInfo,        ParmInfo        },
  { PARM,    "parmstrip",     0, Help_ParmStrip,       ParmStrip       },
  { PARM,    "parmwrite",     0, Help_ParmWrite,       ParmWrite       },
  { PARM,    "printangles",   0, Help_AngleInfo,       AngleInfo       },
  { PARM,    "printatoms",    0, Help_AtomInfo,        AtomInfo        },
  { PARM,    "printbonds",    0, Help_BondInfo,        BondInfo        },
  { PARM,    "printdihedrals",0, Help_DihedralInfo,    DihedralInfo    },
  { PARM,    "resinfo",       0, Help_ResInfo,         ResInfo         },
  { PARM,    "scaledihedralk",0, 0,                    ScaleDihedralK  },
  { PARM,    "solvent",       0, Help_Solvent,         ParmSolvent     },
  // INC_ACTION: ACTION COMMANDS
  { ACTION, "angle", Action_Angle::Alloc, Action_Angle::Help, AddAction },
  { ACTION, "areapermol", Action_AreaPerMol::Alloc, Action_AreaPerMol::Help, AddAction },
  { ACTION, "atomiccorr", Action_AtomicCorr::Alloc, Action_AtomicCorr::Help, AddAction },
  { ACTION, "atomicfluct", Action_AtomicFluct::Alloc, Action_AtomicFluct::Help, AddAction },
  { ACTION, "atommap", Action_AtomMap::Alloc, Action_AtomMap::Help, AddAction },
  { ACTION, "autoimage", Action_AutoImage::Alloc, Action_AutoImage::Help, AddAction },
  { ACTION, "average", Action_Average::Alloc, Action_Average::Help, AddAction },
  { ACTION, "bounds", Action_Bounds::Alloc, Action_Bounds::Help, AddAction },
  { ACTION, "box", Action_Box::Alloc, Action_Box::Help, AddAction },
  { ACTION, "center", Action_Center::Alloc, Action_Center::Help, AddAction },
  { ACTION, "channel", Action_Channel::Alloc, Action_Channel::Help, AddAction },
  { ACTION, "check", Action_CheckStructure::Alloc, Action_CheckStructure::Help, AddAction },
  { ACTION, "checkchirality", Action_CheckChirality::Alloc, Action_CheckChirality::Help, AddAction },
  { ACTION, "checkoverlap", Action_CheckStructure::Alloc, Action_CheckStructure::Help, AddAction },
  { ACTION, "checkstructure", Action_CheckStructure::Alloc, Action_CheckStructure::Help, AddAction },
  { ACTION, "closest", Action_Closest::Alloc, Action_Closest::Help, AddAction },
  { ACTION, "closestwaters", Action_Closest::Alloc, Action_Closest::Help, AddAction },
  { ACTION, "clusterdihedral", Action_ClusterDihedral::Alloc, Action_ClusterDihedral::Help, AddAction },
  { ACTION, "contacts", Action_Contacts::Alloc, Action_Contacts::Help, AddAction },
  { ACTION, "createcrd", Action_CreateCrd::Alloc, Action_CreateCrd::Help, AddAction },
  { ACTION, "createreservoir", Action_CreateReservoir::Alloc, Action_CreateReservoir::Help, AddAction },
  { ACTION, "density", Action_Density::Alloc, Action_Density::Help, AddAction },
  { ACTION, "diffusion", Action_Diffusion::Alloc, Action_Diffusion::Help, AddAction },
  { ACTION, "dihedral", Action_Dihedral::Alloc, Action_Dihedral::Help, AddAction },
  { ACTION, "dihedralscan", Action_DihedralScan::Alloc, Action_DihedralScan::Help, AddAction },
  { ACTION, "dipole", Action_Dipole::Alloc, Action_Dipole::Help, AddAction },
  { ACTION, "distance", Action_Distance::Alloc, Action_Distance::Help, AddAction },
//  { ACTION, "dnaiontracker", Action_DNAionTracker::Alloc, Action_DNAionTracker::Help, AddAction },
  { ACTION, "drms", Action_DistRmsd::Alloc, Action_DistRmsd::Help, AddAction },
  { ACTION, "drmsd", Action_DistRmsd::Alloc, Action_DistRmsd::Help, AddAction },
  { ACTION, "dssp", Action_DSSP::Alloc, Action_DSSP::Help, AddAction },
  { ACTION, "energy", Action_Energy::Alloc, Action_Energy::Help, AddAction },
  { ACTION, "filter", Action_FilterByData::Alloc, Action_FilterByData::Help, AddAction },
  { ACTION, "fixatomorder", Action_FixAtomOrder::Alloc, Action_FixAtomOrder::Help, AddAction },
  { ACTION, "gist", Action_Gist::Alloc, Action_Gist::Help, AddAction },
//  { ACTION, "gfe", Action_GridFreeEnergy::Alloc, Action_GridFreeEnergy::Help, AddAction },
  { ACTION, "grid", Action_Grid::Alloc, Action_Grid::Help, AddAction },
  { ACTION, "hbond", Action_Hbond::Alloc, Action_Hbond::Help, AddAction },
  { ACTION, "image", Action_Image::Alloc, Action_Image::Help, AddAction },
  { ACTION, "jcoupling", Action_Jcoupling::Alloc, Action_Jcoupling::Help, AddAction },
  { ACTION, "lessplit", Action_LESsplit::Alloc, Action_LESsplit::Help, AddAction },
  { ACTION, "lie", Action_LIE::Alloc, Action_LIE::Help, AddAction },
  { ACTION, "lipidorder", Action_OrderParameter::Alloc, Action_OrderParameter::Help, AddAction },
  { ACTION, "makestructure", Action_MakeStructure::Alloc, Action_MakeStructure::Help, AddAction },
  { ACTION, "mask", Action_Mask::Alloc, Action_Mask::Help, AddAction },
  { ACTION, "matrix", Action_Matrix::Alloc, Action_Matrix::Help, AddAction },
  { ACTION, "minimage", Action_MinImage::Alloc, Action_MinImage::Help, AddAction },
  { ACTION, "molsurf", Action_Molsurf::Alloc, Action_Molsurf::Help, AddAction },
  { ACTION, "multidihedral", Action_MultiDihedral::Alloc, Action_MultiDihedral::Help, AddAction },
  { ACTION, "multivector", Action_MultiVector::Alloc, Action_MultiVector::Help, AddAction },
  { ACTION, "nastruct", Action_NAstruct::Alloc, Action_NAstruct::Help, AddAction },
  { ACTION, "nativecontacts", Action_NativeContacts::Alloc, Action_NativeContacts::Help, AddAction },
  { ACTION, "nmrrst", Action_NMRrst::Alloc, Action_NMRrst::Help, AddAction },
  { ACTION, "outtraj", Action_Outtraj::Alloc, Action_Outtraj::Help, AddAction },
  { ACTION, "pairdist", Action_PairDist::Alloc, Action_PairDist::Help, AddAction },
  { ACTION, "pairwise", Action_Pairwise::Alloc, Action_Pairwise::Help, AddAction },
  { ACTION, "principal", Action_Principal::Alloc, Action_Principal::Help, AddAction },
  { ACTION, "projection", Action_Projection::Alloc, Action_Projection::Help, AddAction },
  { ACTION, "pucker", Action_Pucker::Alloc, Action_Pucker::Help, AddAction },
  { ACTION, "radgyr", Action_Radgyr::Alloc, Action_Radgyr::Help, AddAction },
  { ACTION, "radial", Action_Radial::Alloc, Action_Radial::Help, AddAction },
  { ACTION, "randomizeions", Action_RandomizeIons::Alloc, Action_RandomizeIons::Help, AddAction },
  { ACTION, "replicatecell", Action_ReplicateCell::Alloc, Action_ReplicateCell::Help, AddAction },
  { ACTION, "rms", Action_Rmsd::Alloc, Action_Rmsd::Help, AddAction },
  { ACTION, "rmsd", Action_Rmsd::Alloc, Action_Rmsd::Help, AddAction },
  { ACTION, "rog", Action_Radgyr::Alloc, Action_Radgyr::Help, AddAction },
  { ACTION, "rotate", Action_Rotate::Alloc, Action_Rotate::Help, AddAction },
  { ACTION, "rotdif", Action_Rotdif::Alloc, Action_Rotdif::Help, AddAction },
  { ACTION, "runavg", Action_RunningAvg::Alloc, Action_RunningAvg::Help, AddAction },
  { ACTION, "runningaverage", Action_RunningAvg::Alloc, Action_RunningAvg::Help, AddAction },
  { ACTION, "scale", Action_Scale::Alloc, Action_Scale::Help, AddAction },
  { ACTION, "secstruct", Action_DSSP::Alloc, Action_DSSP::Help, AddAction },
  { ACTION, "setvelocity", Action_SetVelocity::Alloc, Action_SetVelocity::Help, AddAction },
  { ACTION, "spam", Action_Spam::Alloc, Action_Spam::Help, AddAction },
  { ACTION, "stfcdiffusion", Action_STFC_Diffusion::Alloc, Action_STFC_Diffusion::Help, AddAction },
  { ACTION, "strip", Action_Strip::Alloc, Action_Strip::Help, AddAction },
  { ACTION, "surf", Action_Surf::Alloc, Action_Surf::Help, AddAction },
  { ACTION, "symmrmsd", Action_SymmetricRmsd::Alloc, Action_SymmetricRmsd::Help, AddAction },
  { ACTION, "temperature", Action_Temperature::Alloc, Action_Temperature::Help, AddAction },
  { ACTION, "trans", Action_Translate::Alloc, Action_Translate::Help, AddAction },
  { ACTION, "translate", Action_Translate::Alloc, Action_Translate::Help, AddAction },
  { ACTION, "unstrip", Action_Unstrip::Alloc, Action_Unstrip::Help, AddAction },
  { ACTION, "unwrap", Action_Unwrap::Alloc, Action_Unwrap::Help, AddAction },
  { ACTION, "vector", Action_Vector::Alloc, Action_Vector::Help, AddAction },
  { ACTION, "velocityautocorr", Action_VelocityAutoCorr::Alloc, Action_VelocityAutoCorr::Help, AddAction },
  { ACTION, "volmap", Action_Volmap::Alloc, Action_Volmap::Help, AddAction},
  { ACTION, "watershell", Action_Watershell::Alloc, Action_Watershell::Help, AddAction },
  // INC_ANALYSIS: ANALYSIS COMMANDS
  { ANALYSIS, "2drms", Analysis_Rms2d::Alloc, Analysis_Rms2d::Help, AddAnalysis },
  { ANALYSIS, "amdbias", Analysis_AmdBias::Alloc, Analysis_AmdBias::Help, AddAnalysis },
  { ANALYSIS, "autocorr", Analysis_AutoCorr::Alloc, Analysis_AutoCorr::Help, AddAnalysis },
  { ANALYSIS, "avg", Analysis_Average::Alloc, Analysis_Average::Help, AddAnalysis },
  { ANALYSIS, "cluster", Analysis_Clustering::Alloc, Analysis_Clustering::Help, AddAnalysis },
  { ANALYSIS, "corr", Analysis_Corr::Alloc, Analysis_Corr::Help, AddAnalysis },
  { ANALYSIS, "correlationcoe", Analysis_Corr::Alloc, Analysis_Corr::Help, AddAnalysis },
  { ANALYSIS, "crank", Analysis_CrankShaft::Alloc, Analysis_CrankShaft::Help, AddAnalysis },
  { ANALYSIS, "crankshaft", Analysis_CrankShaft::Alloc, Analysis_CrankShaft::Help, AddAnalysis },
  { ANALYSIS, "crdfluct", Analysis_CrdFluct::Alloc, Analysis_CrdFluct::Help, AddAnalysis },
  { ANALYSIS, "crosscorr", Analysis_CrossCorr::Alloc, Analysis_CrossCorr::Help, AddAnalysis },
  { ANALYSIS, "curvefit", Analysis_CurveFit::Alloc, Analysis_CurveFit::Help, AddAnalysis },
  { ANALYSIS, "diagmatrix", Analysis_Matrix::Alloc, Analysis_Matrix::Help, AddAnalysis },
  { ANALYSIS, "divergence", Analysis_Divergence::Alloc, Analysis_Divergence::Help, AddAnalysis },
  { ANALYSIS, "fft", Analysis_FFT::Alloc, Analysis_FFT::Help, AddAnalysis },
  { ANALYSIS, "hist", Analysis_Hist::Alloc, Analysis_Hist::Help, AddAnalysis },
  { ANALYSIS, "histogram", Analysis_Hist::Alloc, Analysis_Hist::Help, AddAnalysis },
  { ANALYSIS, "integrate", Analysis_Integrate::Alloc, Analysis_Integrate::Help, AddAnalysis },
  { ANALYSIS, "ired", Analysis_IRED::Alloc, Analysis_IRED::Help, AddAnalysis },
  { ANALYSIS, "kde", Analysis_KDE::Alloc, Analysis_KDE::Help, AddAnalysis },
  { ANALYSIS, "lifetime", Analysis_Lifetime::Alloc, Analysis_Lifetime::Help, AddAnalysis },
  { ANALYSIS, "lowestcurve", Analysis_LowestCurve::Alloc, Analysis_LowestCurve::Help, AddAnalysis },
  { ANALYSIS, "matrix", Analysis_Matrix::Alloc, Analysis_Matrix::Help, AddAnalysis },
  { ANALYSIS, "meltcurve", Analysis_MeltCurve::Alloc, Analysis_MeltCurve::Help, AddAnalysis },
  { ANALYSIS, "modes", Analysis_Modes::Alloc, Analysis_Modes::Help, AddAnalysis },
  { ANALYSIS, "multihist", Analysis_MultiHist::Alloc, Analysis_MultiHist::Help, AddAnalysis },
  { ANALYSIS, "overlap", Analysis_Overlap::Alloc, Analysis_Overlap::Help, AddAnalysis },
  { ANALYSIS, "phipsi", Analysis_PhiPsi::Alloc, Analysis_PhiPsi::Help, AddAnalysis },
  { ANALYSIS, "regress", Analysis_Regression::Alloc, Analysis_Regression::Help, AddAnalysis },
  { ANALYSIS, "remlog", Analysis_RemLog::Alloc, Analysis_RemLog::Help, AddAnalysis },
  { ANALYSIS, "rms2d", Analysis_Rms2d::Alloc, Analysis_Rms2d::Help, AddAnalysis },
  { ANALYSIS, "rmsavgcorr", Analysis_RmsAvgCorr::Alloc, Analysis_RmsAvgCorr::Help, AddAnalysis },
  { ANALYSIS, "runningavg", Analysis_RunningAvg::Alloc, Analysis_RunningAvg::Help, AddAnalysis },
  { ANALYSIS, "spline", Analysis_Spline::Alloc, Analysis_Spline::Help, AddAnalysis },
  { ANALYSIS, "stat", Analysis_Statistics::Alloc, Analysis_Statistics::Help, AddAnalysis },
  { ANALYSIS, "statistics", Analysis_Statistics::Alloc, Analysis_Statistics::Help, AddAnalysis },
  { ANALYSIS, "timecorr", Analysis_Timecorr::Alloc, Analysis_Timecorr::Help, AddAnalysis },
  { ANALYSIS, "vectormath", Analysis_VectorMath::Alloc, Analysis_VectorMath::Help, AddAnalysis },
  // DEPRECATED COMMANDS
  { DEPRECATED, "acceptor",     0, Deprecate_Hbond,        0 },
  { DEPRECATED, "avgcoord",     0, Deprecate_AvgCoord,     0 },
  { DEPRECATED, "bondsearch",   0, Deprecate_TopSearch,    0 },
  { DEPRECATED, "donor",        0, Deprecate_Hbond,        0 },
  { DEPRECATED, "maxdist",      0, Deprecate_MinDist,      0 },
  { DEPRECATED, "mindist",      0, Deprecate_MinDist,      0 },
  { DEPRECATED, "molsearch",    0, Deprecate_TopSearch,    0 },
  { DEPRECATED, "nobondsearch", 0, Deprecate_TopSearch,    0 },
  { DEPRECATED, "nomolsearch",  0, Deprecate_TopSearch,    0 },
  { DEPRECATED, "parmbondinfo", 0, Deprecate_ParmBondInfo, 0 },
  { DEPRECATED, "parmmolinfo",  0, Deprecate_ParmMolInfo,  0 },
  { DEPRECATED, "parmresinfo",  0, Deprecate_ParmResInfo,  0 },
  { NONE      , 0,              0, 0,                      0 }
};
