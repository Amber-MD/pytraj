#include "Action_CreateReservoir.h"
#include "CpptrajStdio.h"

// CONSTRUCTOR
Action_CreateReservoir::Action_CreateReservoir() :
  original_trajparm_(0),
  ene_(0),
  bin_(0),
# ifdef BINTRAJ
  reservoirT_(0.0),
# endif
  iseed_(0),
  trajIsOpen_(false),
  nframes_(0)
{}

void Action_CreateReservoir::Help() {
  mprintf("\t<filename> ene <energy data set> [bin <cluster bin data set>]\n"
          "\ttemp0 <temp0> iseed <iseed> [velocity]\n"
          "\t[parm <parmfile> | parmindex <#>] [title <title>]\n"
          "  Create structure reservoir for use with reservoir REMD simulations using\n"
          "  energies in <energy data set>, temperature <temp0> and random seed <iseed>\n"
          "  Include velocities if [velocity] is specified. If <cluster bin data set> is\n"
          "  specified from e.g. a previous 'clusterdihedral' command, the reservoir can\n"
          "  be used for non-Boltzmann reservoir REMD (rremd==3).\n");
}

// Action_CreateReservoir::Init()
Action::RetType Action_CreateReservoir::Init(ArgList& actionArgs, TopologyList* PFL, DataSetList* DSL, DataFileList* DFL, int debugIn)
{
# ifdef BINTRAJ
  // Get keywords
  filename_ = actionArgs.GetStringNext();
  if (filename_.empty()) {
    mprinterr("Error: createreservoir: No filename specified.\n");
    return Action::ERR;
  }
  reservoirT_ = actionArgs.getKeyDouble("temp0", -1.0);
  if (reservoirT_ < 0.0) {
    mprinterr("Error: Reservoir temperature must be specified and cannot be < 0.0\n");
    return Action::ERR;
  }
  iseed_ = actionArgs.getKeyInt("iseed", 0);
  if (iseed_ < 1) {
    mprinterr("Error: Reservoir random seed must be specified and > 0\n");
    return Action::ERR;
  }
  reservoir_.SetVelocity( actionArgs.hasKey("velocity") );
  // Get parm for reservoir traj
  original_trajparm_ = PFL->GetParm( actionArgs );
  if (original_trajparm_ == 0) {
    mprinterr("Error: createreservoir: no topology.\n");
    return Action::ERR;
  }
  // Get energy data set
  std::string eneDsname = actionArgs.GetStringKey("ene");
  DataSet* dstmp = DSL->GetDataSet( eneDsname );
  if (dstmp == 0) {
    mprinterr("Error: could not get energy data set %s\n", eneDsname.c_str());
    return Action::ERR;
  }
  if (dstmp->Type() != DataSet::FLOAT && dstmp->Type() != DataSet::DOUBLE) {
    mprinterr("Error: energy data set %s must be type FLOAT or DOUBLE.\n", eneDsname.c_str());
    return Action::ERR;
  }
  if (dstmp->Ndim() != 1) {
    mprinterr("Error: energy data set is not 1D (%u)\n", dstmp->Ndim());
    return Action::ERR;
  }
  ene_ = static_cast<DataSet_1D*>( dstmp );
  // Get bin data set
  std::string binDSname = actionArgs.GetStringKey("bin");
  if (!binDSname.empty()) {
    dstmp = DSL->GetDataSet( binDSname );
    if (dstmp == 0) {
      mprinterr("Error: could not get bin data set %s\n", binDSname.c_str());
      return Action::ERR;
    } else if (dstmp->Ndim() != 1) {
      mprinterr("Error: bin data set must be one dimensional.\n");
      return Action::ERR;
    }
    bin_ = static_cast<DataSet_1D*>( dstmp );
  }
  trajIsOpen_ = false;
  nframes_ = 0;
  // Setup output reservoir file
  reservoir_.SetDebug( debugIn );
  // Set title
  std::string title = actionArgs.GetStringKey("title");
  if (title.empty())
    title.assign("Cpptraj generated structure reservoir");
  reservoir_.SetTitle( title );
  // Process additional netcdf traj args
  //reservoir_.processWriteArgs( actionArgs );

  mprintf("    CREATERESERVOIR: %s, energy data %s", filename_.c_str(),
          ene_->Legend().c_str());
  if (bin_ != 0)
    mprintf(", bin data %s", bin_->Legend().c_str());
  mprintf("\n\tReservoir temperature= %.2f, random seed= %i\n", reservoirT_, iseed_);
  if (reservoir_.HasV())
    mprintf("\tVelocities will be written to reservoir.\n");
  mprintf("\tTopology: %s\n", original_trajparm_->c_str());
  return Action::OK;
# else
  mprinterr("Error: reservoir requires NetCDF. Reconfigure with NetCDF enabled.\n");
  return Action::ERR;
# endif
}

// Action_CreateReservoir::Setup()
Action::RetType Action_CreateReservoir::Setup(Topology* currentParm, Topology** parmAddress)
{
# ifdef BINTRAJ
  // Check that input parm matches current parm
  if (original_trajparm_->Pindex() != currentParm->Pindex()) {
    mprintf("Info: createreservoir was set up for topology %s\n", original_trajparm_->c_str());
    mprintf("Info: skipping topology %s\n", currentParm->c_str());
    return Action::ERR;
  }
  if (!trajIsOpen_) {
    mprintf("\tCreating reservoir file %s\n", filename_.c_str());
    // Use parm to set up box info for the reservoir. FIXME: Necessary?
    reservoir_.SetBox( currentParm->ParmBox() );
    // Set up write and open - no append.
    if (reservoir_.setupTrajout( filename_, currentParm, currentParm->Nframes(), false))
      return Action::ERR;
    // Add reservoir vars to netcdf traj
    if (reservoir_.createReservoir((bin_!=0), reservoirT_, iseed_)) {
      mprinterr("Error: Could not add reservoir variables to netcdf trajectory.\n");
      return Action::ERR;
    }

    trajIsOpen_ = true;
    nframes_ = 0;
  }
  return Action::OK;
# else
  return Action::ERR;
# endif 
}

// Action_CreateReservoir::DoAction()
Action::RetType Action_CreateReservoir::DoAction(int frameNum, Frame* currentFrame, 
                                                 Frame** frameAddress) 
{
# ifdef BINTRAJ
  int bin = -1;
  if (bin_ != 0) bin = (int)bin_->Dval(frameNum);
  if (reservoir_.writeReservoir(nframes_++, *currentFrame, ene_->Dval(frameNum), bin))
    return Action::ERR;
  return Action::OK;
# else
  return Action::ERR;
# endif
}

// Action_CreateReservoir::Print()
void Action_CreateReservoir::Print() {
# ifdef BINTRAJ
  mprintf("\tReservoir %s: %u frames.\n", filename_.c_str(), nframes_);
  reservoir_.closeTraj();
  trajIsOpen_=false;
# endif
}
