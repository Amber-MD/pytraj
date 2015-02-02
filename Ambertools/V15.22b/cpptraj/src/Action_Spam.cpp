// Action_Spam
#include <cmath> // sqrt
#include <cstdio> // sscanf, sprintf
#include "Action_Spam.h"
#include "Box.h"
#include "CpptrajFile.h"
#include "CpptrajStdio.h"
#include "Constants.h" // ELECTOAMBER
#include "DataSet_integer.h"
#include "DistRoutines.h"
#include "StringRoutines.h" // fileexists

// CONSTRUCTOR
Action_Spam::Action_Spam() :
  ensembleNum_(-1),
  bulk_(0.0),
  purewater_(false),
  reorder_(false),
  cut2_(144.0),
  onecut2_(1.0 / 144.0),
  doublecut_(24.0),
  site_size_(2.5),
  sphere_(false),
  Nframes_(0),
  overflow_(false)
{ }

void Action_Spam::Help() {
  mprintf("\t<filename> [solv <solvname>] [reorder] [name <name>] [bulk <value>]\n"
          "\t[purewater] [cut <cut>] [info <infofile>] [summary <summary>]\n"
          "\t[site_size <size>] [sphere] [out <datafile>]\n\n"
          "    <filename> : File with the peak locations present (XYZ- format)\n"
          "    <solvname> : Name of the solvent residues\n"
          "    <cut>      : Non-bonded cutoff for energy evaluation\n"
          "    <value>    : SPAM free energy of the bulk solvent\n"
          "    <infofile> : File with stats about which sites are occupied when.\n"
          "    <size>     : Size of the water site around each density peak.\n"
          "    [sphere]   : Treat each site like a sphere.\n"
          "    [purewater]: The system is pure water---used to parametrize the bulk values.\n"
          "    [reorder]  : The solvent should be re-ordered so the same solvent molecule\n"
          "                 is always in the same site.\n"
          "    <summary>  : File with the summary of all SPAM results. If not specified,\n"
          "                 no SPAM energies will be calculated.\n"
          "    <datafile> : Data file with all SPAM energies for each snapshot.\n");
}

// Action_Spam::init()
Action::RetType Action_Spam::Init(ArgList& actionArgs, TopologyList* PFL,
                      FrameList *FL, DataSetList *DSL, DataFileList *DFL,
                      int debugIn)
{
  ensembleNum_ = DSL->EnsembleNum();
  // Always use imaged distances
  InitImaging(true);
  // This is needed everywhere in this function scope
  std::string filename;

  // See if we're doing pure water. If so, we don't need a peak file
  purewater_ = actionArgs.hasKey("purewater");

  if (purewater_) {
    // We still need the cutoff
    double cut = actionArgs.getKeyDouble("cut", 12.0);
    cut2_ = cut * cut;
    doublecut_ = 2 * cut;
    onecut2_ = 1 / cut2_;
    // See if we write to a data file
    datafile_ = actionArgs.GetStringKey("out");
    // Generate the data set name, and hold onto the master data set list
    std::string ds_name = actionArgs.GetStringKey("name");
    if (ds_name.empty())
      ds_name = myDSL_.GenerateDefaultName("SPAM");
    // We only have one data set averaging over every water. Add it here
    myDSL_.AddSet(DataSet::DOUBLE, ds_name, NULL);
    solvname_ = actionArgs.GetStringKey("solv");
    if (solvname_.empty())
      solvname_ = std::string("WAT");
  }else {
    // Get the file name with the peaks defined in it
    filename = actionArgs.GetStringNext();

    if (filename.empty() || !fileExists(filename)) {
      mprinterr("Spam: Error: Peak file [%s] does not exist!\n", filename.c_str());
      return Action::ERR;
    }

    // Get the remaining optional arguments
    solvname_ = actionArgs.GetStringKey("solv");
    if (solvname_.empty())
      solvname_ = std::string("WAT");
    reorder_ = actionArgs.hasKey("reorder");
    bulk_ = actionArgs.getKeyDouble("bulk", 0.0);
    double cut = actionArgs.getKeyDouble("cut", 12.0);
    cut2_ = cut * cut;
    doublecut_ = 2 * cut;
    onecut2_ = 1 / cut2_;
    infoname_ = actionArgs.GetStringKey("info");
    if (infoname_.empty())
      infoname_ = std::string("spam.info");
    // The default maskstr is the Oxygen atom of the solvent
    summaryfile_ = actionArgs.GetStringKey("summary");
    // Divide site size by 2 to make it half the edge length (or radius)
    site_size_ = actionArgs.getKeyDouble("site_size", 2.5) / 2.0;
    sphere_ = actionArgs.hasKey("sphere");
    // If it's a sphere, square the radius to compare with
    if (sphere_)
      site_size_ *= site_size_;
    datafile_ = actionArgs.GetStringKey("out");
    std::string ds_name = actionArgs.GetStringKey("name");
    if (ds_name.empty())
      ds_name = myDSL_.GenerateDefaultName("SPAM");

    // Parse through the peaks file and extract the peaks
    CpptrajFile peakfile;
    if (peakfile.OpenRead(filename)) {
      mprinterr("SPAM: Error: Could not open %s for reading!\n", filename.c_str());
      return Action::ERR;
    }
    std::string line = peakfile.GetLine();
    int npeaks = 0;
    while (!line.empty()) {
      if (sscanf(line.c_str(), "%d", &npeaks) != 1) {
        line = peakfile.GetLine();
        continue;
      }
      line = peakfile.GetLine();
      break;
    }
    while (!line.empty()) {
      double x, y, z, dens;
      if (sscanf(line.c_str(), "C %lg %lg %lg %lg", &x, &y, &z, &dens) != 4) {
        line = peakfile.GetLine();
        continue;
      }
      line = peakfile.GetLine();
      peaks_.push_back(Vec3(x, y, z));
    }
    peakfile.CloseFile();
    // Check that our initial number of peaks matches our parsed peaks. Warn
    // otherwise
    if (npeaks != (int)peaks_.size())
      mprinterr("SPAM: Warning: %s claims to have %d peaks, but really has %d!\n",
                filename.c_str(), npeaks, peaks_.size());
    // Now add all of the data sets
    for (int i = 0; i < (int)peaks_.size(); i++) {
      myDSL_.AddSetAspect(DataSet::DOUBLE, ds_name,
                                      integerToString(i+1).c_str());
      // Add a new list of integers to keep track of omitted frames
      std::vector<int> vec;
      peakFrameData_.push_back(vec);
    }
  }

  // Print info now
  if (purewater_) {
    mprintf("SPAM: Calculating bulk value for pure solvent\n");
    if (!datafile_.empty())
      mprintf("SPAM: Printing solvent energies to %s\n", datafile_.c_str());
    mprintf("SPAM: Using a %.2f Angstrom non-bonded cutoff with shifted EEL.\n",
            sqrt(cut2_));
    if (reorder_)
      mprintf("SPAM: Warning: Re-ordering makes no sense for pure solvent.\n");
    if (!summaryfile_.empty())
      mprintf("SPAM: Printing solvent SPAM summary to %s\n", summaryfile_.c_str());
  }else {
    mprintf("SPAM: Solvent [%s] density peaks taken from %s.\n",
            solvname_.c_str(), filename.c_str());
    mprintf("SPAM: %d density peaks will be analyzed from %s.\n",
            peaks_.size(), filename.c_str());
    mprintf("SPAM: Occupation information printed to %s.\n", infoname_.c_str());
    mprintf("SPAM: Sites are ");
    if (sphere_)
      mprintf("spheres with diameter %.3lf\n", site_size_);
    else
      mprintf("boxes with edge length %.3lf\n", site_size_);
    if (reorder_) {
      mprintf("SPAM: Re-ordering trajectory so each site always has ");
      mprintf("the same water molecule.\n");
    }
    if (summaryfile_.empty() && datafile_.empty()) {
      if (!reorder_) {
        mprinterr("SPAM: Error: Not re-ordering trajectory or calculating energies. ");
        mprinterr("Nothing to do!\n");
        return Action::ERR;
      }
      mprintf("SPAM: Not calculating any SPAM energies\n");
    }else {
      mprintf("SPAM: Using a non-bonded cutoff of %.2lf Ang. with a EEL shifting function.\n",
              sqrt(cut2_));
      mprintf("SPAM: Bulk solvent SPAM energy taken as %.3lf kcal/mol\n", bulk_);
    }
  }
  mprintf("#Citation: Cui, G.; Swails, J.M.; Manas, E.S.; \"SPAM: A Simple Approach\n"
          "#          for Profiling Bound Water Molecules\"\n"
          "#          J. Chem. Theory Comput., 2013, 9 (12), pp 5539–5549.\n");

  return Action::OK;
}

// Action_Spam::setup()
Action::RetType Action_Spam::Setup(Topology* currentParm, Topology** parmAddress) {

  // We need box info
  if (currentParm->BoxType() == Box::NOBOX) {
    mprinterr("Error: SPAM: Must have explicit solvent with periodic boundaries!");
    return Action::ERR;
  }

  // See if our box dimensions are too small for our cutoff...
  if (currentParm->ParmBox().BoxX() < doublecut_ ||
      currentParm->ParmBox().BoxY() < doublecut_ ||
      currentParm->ParmBox().BoxZ() < doublecut_) {
    mprinterr("Error: SPAM: The box appears to be too small for your cutoff!\n");
    return Action::ERR;
  }
  // Set up the solvent_residues_ vector
  int resnum = 0;
  for (Topology::res_iterator res = currentParm->ResStart();
       res != currentParm->ResEnd(); res++) {
    if (res->Name().Truncated() == solvname_) {
      solvent_residues_.push_back(*res);
      // Tabulate COM
      double mass = 0.0;
      for (int i = res->FirstAtom(); i < res->LastAtom(); i++)
        mass += (*currentParm)[i].Mass();
    }
    resnum++;
  }

  // DEBUG
  mprintf("SPAM: Found %d solvent residues [%s]\n", solvent_residues_.size(),
          solvname_.c_str());

  // Set up the charge array and check that we have enough info
  if (SetupParms(currentParm)) return Action::ERR;

  // Back up the parm
  // NOTE: This is a full copy - use reference instead?
  CurrentParm_ = *currentParm;

  return Action::OK;
}

// Action_Spam::SetupParms
/** Sets the temporary charge array and makes sure that we have the necessary
  * parameters in our topology to calculate nonbonded energy terms
  */
int Action_Spam::SetupParms(Topology* ParmIn) {
  // Store the charges
  atom_charge_.clear();
  atom_charge_.reserve( ParmIn->Natom() );
  for (Topology::atom_iterator atom = ParmIn->begin();
       atom != ParmIn->end(); ++atom)
    atom_charge_.push_back( atom->Charge() * Constants::ELECTOAMBER );
  if (!ParmIn->Nonbond().HasNonbond()) {
    mprinterr("Error: SPAM: Parm does not have LJ information.\n");
    return 1;
  }
  return 0;
}

double Action_Spam::Calculate_Energy(Frame *frameIn, Residue const& res) {

  // The first atom of the solvent residue we want the energy from
  double result = 0;
  /* Now loop through all atoms in the residue and loop through the pairlist to
   * get the energies
   */
  for (int i = res.FirstAtom(); i < res.LastAtom(); i++) {
    Vec3 atm1 = Vec3(frameIn->XYZ(i));
    for (int j = 0; j < CurrentParm_.Natom(); j++) {
      if (j >= res.FirstAtom() && j < res.LastAtom()) continue;
      Vec3 atm2 = Vec3(frameIn->XYZ(j));
      double dist2;
      // Get imaged distance
      Matrix_3x3 ucell, recip;
      switch( ImageType() ) {
        case NONORTHO:
          frameIn->BoxCrd().ToRecip(ucell, recip);
          dist2 = DIST2_ImageNonOrtho(atm1, atm2, ucell, recip);
          break;
        case ORTHO:
          dist2 = DIST2_ImageOrtho(atm1, atm2, frameIn->BoxCrd());
          break;
        default:
          dist2 = DIST2_NoImage(atm1, atm2);
      }
      if (dist2 < cut2_) {
        double qiqj = atom_charge_[i] * atom_charge_[j];
        NonbondType const& LJ = CurrentParm_.GetLJparam(i, j);
        double r2 = 1 / dist2;
        double r6 = r2 * r2 * r2;
                  // Shifted electrostatics: qiqj/r * (1-r/rcut)^2 + VDW
        double shift = (1 - dist2 * onecut2_);
        result += qiqj / sqrt(dist2) * shift * shift + LJ.A() * r6 * r6 - LJ.B() * r6;
      }
    }
  }
  return result;
}

// Action_Spam::action()
Action::RetType Action_Spam::DoAction(int frameNum, Frame* currentFrame,
                                     Frame ** frameAddress)
{

  Nframes_++;

  // Check that our box is still big enough...
  overflow_ = overflow_ || currentFrame->BoxCrd().BoxX() < doublecut_ ||
                           currentFrame->BoxCrd().BoxY() < doublecut_ ||
                           currentFrame->BoxCrd().BoxZ() < doublecut_;
  if (purewater_)
    return DoPureWater(frameNum, currentFrame);
  else
    return DoSPAM(frameNum, currentFrame);

}

// Action_Spam::DoPureWater
/** Carries out SPAM analysis for pure water to parametrize bulk */
Action::RetType Action_Spam::DoPureWater(int frameNum, Frame* currentFrame) {
  /* This is relatively simple... We have to loop through every water molecule
   * for every frame, calculate the energy of that solvent molecule, and add
   * that to our one data set. Therefore we will have NFRAMES * NWATER data
   * points
   */
  int wat = 0;
  int basenum = frameNum * solvent_residues_.size();
  for (std::vector<Residue>::const_iterator res = solvent_residues_.begin();
        res != solvent_residues_.end(); res++) {
    double ene = Calculate_Energy(currentFrame, *res);
    myDSL_[0]->Add(basenum + wat, &ene);
    wat++;
  }
  return Action::OK;
}

// Action_Spam::DoSPAM
/** Carries out SPAM analysis on a typical system */
Action::RetType Action_Spam::DoSPAM(int frameNum, Frame* currentFrame) {
  // Set up a function pointer to see if we are inside
  bool (*inside)(Vec3, Vec3, double);
  if (sphere_)
    inside = &inside_sphere;
  else
    inside = &inside_box;

  /* A list of all solvent residues and the sites that they are reserved for. An
   * unreserved solvent residue has an index -1. At the end, we will go through
   * and re-order the frame if requested.
   */
  std::vector<int> reservations(solvent_residues_.size(), -1);
  // Tabulate all of the COMs
  std::vector<Vec3> comlist;
  for (std::vector<Residue>::const_iterator res = solvent_residues_.begin();
       res != solvent_residues_.end(); res++) {
    comlist.push_back(currentFrame->VCenterOfMass(res->FirstAtom(), res->LastAtom()));
  }
  // Loop through each peak and then scan through every residue, and assign a
  // solvent residue to each peak
  int pknum = 0;
  for (std::vector<Vec3>::const_iterator pk = peaks_.begin();
       pk != peaks_.end(); pk++) {
    int resnum = 0;
    for (std::vector<Residue>::const_iterator res = solvent_residues_.begin();
         res != solvent_residues_.end(); res++) {
      // If we're inside, make sure this residue is not already `claimed'. If it
      // is, assign it to the closer peak center
      if (inside(*pk, comlist[resnum], site_size_)) {
        if (reservations[resnum] > 0) {
          Vec3 diff1 = comlist[resnum] - *pk;
          Vec3 diff2 = comlist[resnum] - peaks_[ reservations[resnum] ];
          // If we are closer, update. Otherwise do nothing
          if (diff1.Magnitude2() < diff2.Magnitude2())
            reservations[resnum] = pknum;
        }else
          reservations[resnum] = pknum;
      }
      resnum++;
    }
    pknum++;
  }

  /* Now we have a vector of reservations. We want to make sure that each site
   * is occupied once and only once. If a site is unoccupied, add frameNum to
   * this peak's data set in peakFrameData_. If a site is double-occupied, add
   * -frameNum to this peak's data set in peakFrameData_.
   */
  std::vector<bool> occupied(peaks_.size(), false);
  std::vector<bool> doubled(peaks_.size(), false); // to avoid double-additions
  for (std::vector<int>::const_iterator it = reservations.begin();
       it != reservations.end(); it++) {
    if (*it > -1) {
      if (!occupied[*it])
        occupied[*it] = true;
      else if (!doubled[*it]) {
        peakFrameData_[*it].push_back(-frameNum); // double-occupied, frameNum will be ignored
        doubled[*it] = true;
      }
    }
  }
  // Now loop through and add all non-occupied sites
  for (unsigned int i = 0; i < peaks_.size(); i++)
    if (!occupied[i]) 
      peakFrameData_[i].push_back(frameNum);
  // Now adjust the occupied vectors to only contain 'true' for sites we need to
  // analyze (i.e., make all doubled points 'unoccupied')
  for (unsigned int i = 0; i < peaks_.size(); i++)
    if (doubled[i])
      occupied[i] = false;

  // If we have to calculate energies, do that here
  if (!summaryfile_.empty() || !datafile_.empty()) {
    /* Loop through every peak, then loop through the water molecules to find
     * which one is in that site, and calculate the LJ and EEL energies for that
     * water molecule within a given cutoff.
     */
    for (int peak = 0; peak < (int)peaks_.size(); peak++) {
      // Skip unoccupied peaks
      if (!occupied[peak]) {
        double val = 0;
        myDSL_[peak]->Add(frameNum, &val);
        continue;
      }
      for (unsigned int i = 0; i < reservations.size(); i++) {
        if (reservations[i] != peak) continue;
        /* Now we have our residue number. Create a pairlist for each solvent
         * molecule that can be used for each atom in that residue. Should
         * provide some time savings.
         */
        double ene = Calculate_Energy(currentFrame, solvent_residues_[i]);
        myDSL_[peak]->Add(frameNum, &ene);
        break;
      }
    }
  }

  // If we have to re-order trajectories, do that here
  if (reorder_) {
    /* Loop over every occupied site and swap the atoms so the same solvent
     * residue is always in the same site
     */
    for (int i = 0; i < (int)peaks_.size(); i++) {
      // Skip unoccupied sites
      if (!occupied[i]) continue;
      for (unsigned int j = 0; j < solvent_residues_.size(); j++) {
        // This is the solvent residue in our site
        if (reservations[j] == i) {
          for (int k = 0; k < solvent_residues_[j].NumAtoms(); k++)
            currentFrame->SwapAtoms(solvent_residues_[i].FirstAtom()+k,
                                    solvent_residues_[j].FirstAtom()+k);
          // Since we swapped solvent_residues_ of 2 solvent atoms, we also have
          // to swap reservations[i] and reservations[j]...
          int tmp = reservations[j];
          reservations[j] = reservations[i]; reservations[i] = tmp;
        }
      }
    }
  }

  return Action::OK;

}

void Action_Spam::Print() {
  // Print the spam info file if we didn't do pure water
  if (!purewater_) {
    CpptrajFile info;
    if (info.OpenEnsembleWrite(infoname_, ensembleNum_)) {
      mprinterr("Error: SPAM: Could not open %s for writing.\n",
                infoname_.c_str());
      return;
    }

    // Warn about any overflows
    if (overflow_)
      mprinterr("Warning: SPAM: Some frames had a box too small for the cutoff.\n");

    // Print information about each missing peak
    info.Printf("# There are %d density peaks and %d frames\n\n",
                (int)peaks_.size(), Nframes_);
    // Loop over every Data set
    for (unsigned int i = 0; i < peakFrameData_.size(); i++) {
      // Skip peaks with 0 unoccupied sites
      if (peakFrameData_[i].size() == 0) continue;
      // Find out how many double-occupied frames there are
      int ndouble = 0;
      for (unsigned int j = 0; j < peakFrameData_[i].size(); j++)
        if (peakFrameData_[i][j] < 0) ndouble++;
      info.Printf("# Peak %u has %d omitted frames (%d double-occupied)\n",
                  i, peakFrameData_[i].size(), ndouble);
      for (unsigned int j = 0; j < peakFrameData_[i].size(); j++) {
        if (j > 0 && j % 10 == 0) info.Printf("\n");
        info.Printf("%7d ", peakFrameData_[i][j]);
      }
      info.Printf("\n\n");
    }

    // Now close the info file
    info.CloseFile();
  }

  // Print the summary file with the calculated SPAM energies
  if (!summaryfile_.empty()) {
    // Not enabled yet -- just print out the data files.
    mprinterr("Warning: SPAM: SPAM calculation not yet enabled.\n");
    if (datafile_.empty()) datafile_ = summaryfile_;
  }
  // Now print the energy data
  if (!datafile_.empty()) {
    // Now write the data file with all of the SPAM energies
    DataFile dfl;
    ArgList dummy;
    dfl.SetupDatafile(datafile_, dummy, 0);
    for (int i = 0; i < (int)myDSL_.size(); i++) {
      dfl.AddSet(myDSL_[i]);
    }
    dfl.WriteData();
  }
}

bool inside_box(Vec3 gp, Vec3 pt, double edge) {
  return (gp[0] + edge > pt[0] && gp[0] - edge < pt[0] &&
          gp[1] + edge > pt[1] && gp[1] - edge < pt[1] &&
          gp[2] + edge > pt[2] && gp[2] - edge < pt[2]);
}

bool inside_sphere(Vec3 gp, Vec3 pt, double rad2) {
  return ( (gp[0]-pt[0])*(gp[0]-pt[0]) + (gp[1]-pt[1])*(gp[1]-pt[1]) +
           (gp[2]-pt[2])*(gp[2]-pt[2]) < rad2 );
}
