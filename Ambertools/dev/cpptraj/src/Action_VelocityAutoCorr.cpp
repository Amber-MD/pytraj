#include "Action_VelocityAutoCorr.h"
#include "CpptrajStdio.h"
#include "ProgressBar.h"
#include "DataSet_Mesh.h"
#include "DataSet_double.h"
#include "Constants.h"
#include "Corr.h"
#ifdef _OPENMP
#  include "omp.h"
#endif

// CONSTRUCTOR
Action_VelocityAutoCorr::Action_VelocityAutoCorr() :
  useVelInfo_(false), useFFT_(false), normalize_(false), VAC_(0), tstep_(0.0),
  maxLag_(0) {}

void Action_VelocityAutoCorr::Help() {
  mprintf("\t[<set name>] [<mask>] [usevelocity] [out <filename>]\n"
          "\t[maxlag <time>] [tstep <timestep>] [direct] [norm]\n"
          "  Calculate velocity auto-correlation function for atoms in <mask>\n");
}

// Action_VelocityAutoCorr::Init()
Action::RetType Action_VelocityAutoCorr::Init(ArgList& actionArgs, TopologyList* PFL, DataSetList* DSL, DataFileList* DFL, int debugIn)
{
  useVelInfo_ = actionArgs.hasKey("usevelocity");
  mask_.SetMaskString( actionArgs.GetMaskNext() );
  DataFile* outfile =  DFL->AddDataFile( actionArgs.GetStringKey("out"), actionArgs );
  maxLag_ = actionArgs.getKeyInt("maxlag", -1);
  tstep_ = actionArgs.getKeyDouble("tstep", 1.0);
  useFFT_ = !actionArgs.hasKey("direct");
  normalize_ = actionArgs.hasKey("norm");
  // Set up output data set
  VAC_ = DSL->AddSet(DataSet::DOUBLE, actionArgs.GetStringNext(), "VAC");
  if (VAC_ == 0) return Action::ERR;
  VAC_->Dim(0).SetStep( tstep_ );
  if (outfile != 0) outfile->AddSet( VAC_ ); 

  mprintf("    VELOCITYAUTOCORR:\n"
          "\tCalculate velocity auto-correlation function for atoms in mask '%s'\n",
          mask_.MaskString());
  if (useVelInfo_)
    mprintf("\tUsing velocity information present in frames.\n");
  else
    mprintf("\tCalculating velocities between consecutive frames.\n");
  if (outfile != 0)
    mprintf("\tOutput data set '%s' to '%s'\n", VAC_->Legend().c_str(), 
            outfile->DataFilename().full());
  if (maxLag_ < 1)
    mprintf("\tMaximum lag will be half total # of frames");
  else
    mprintf("\tMaximum lag is %i frames", maxLag_);
  mprintf(", time step is %f ps\n", tstep_);
  if (useFFT_)
    mprintf("\tUsing FFT to calculate autocorrelation function.\n");
  else
    mprintf("\tUsing direct method to calculate autocorrelation function.\n");
  if (normalize_)
    mprintf("\tNormalizing autocorrelation function to 1.0\n");
  return Action::OK;
}

// Action_VelocityAutoCorr::Setup()
/** For this to be valid the same # of atoms should be selected each time. */
Action::RetType Action_VelocityAutoCorr::Setup(Topology* currentParm,
                                               Topology** parmAddress)
{
  if (currentParm->SetupIntegerMask( mask_ )) return Action::ERR;
  mask_.MaskInfo();
  if (mask_.None()) {
    mprintf("Warning: No atoms selected by mask.\n");
    return Action::ERR;
  }
  // If using velocity info, check that it is present.
  if (useVelInfo_ && !currentParm->ParmCoordInfo().HasVel()) {
    mprinterr("Error: 'usevelocity' specified but no velocity info assocated with %s\n",
              currentParm->c_str());
    return Action::ERR;
  }
  // If we have already started recording velocities, check that the current 
  // number of selected atoms has remained the same.
  if (!Vel_.empty()) {
    if ((int)Vel_.size() != mask_.Nselected()) {
      mprinterr("Error: # of selected atoms %i has changed (previously %zu)\n",
                mask_.Nselected(), Vel_.size());
      return Action::ERR;
    }
  } else
    Vel_.resize( mask_.Nselected() );
  return Action::OK;
}

// Action_VelocityAutoCorr::DoAction()
Action::RetType Action_VelocityAutoCorr::DoAction(int frameNum, 
                                                  Frame* currentFrame,
                                                  Frame** frameAddress)
{
  if (!useVelInfo_) {
    // Calculate pseudo-velocities between frames.
    if (!previousFrame_.empty()) {
      // This is the first frame which we can calculate pseudo-velocity.
      VelArray::iterator vel = Vel_.begin();
      for (AtomMask::const_iterator atom = mask_.begin();
                                    atom != mask_.end(); 
                                  ++atom, ++vel)
        vel->AddVxyz( (Vec3(currentFrame->XYZ(*atom)) - Vec3(previousFrame_.XYZ(*atom))) / tstep_ );
    }
    previousFrame_ = *currentFrame;
  } else {
    // Use velocity information in the frame.
    // FIXME: Eventually dont assume V is in Amber units.
    VelArray::iterator vel = Vel_.begin();
    for (AtomMask::const_iterator atom = mask_.begin();
                                  atom != mask_.end(); 
                                ++atom, ++vel)
      vel->AddVxyz( Vec3(currentFrame->VXYZ( *atom )) * Constants::AMBERTIME_TO_PS );
  }
  return Action::OK;
}

// Action_VelocityAutoCorr::Print()
void Action_VelocityAutoCorr::Print() {
  if (Vel_.empty()) return;
  mprintf("    VELOCITYAUTOCORR:\n");
  mprintf("\t%zu vectors have been saved, total length of each = %zu\n",
          Vel_.size(), Vel_[0].Size());
  int maxlag;
  if (maxLag_ <= 0) {
    maxlag = (int)Vel_[0].Size() / 2;
    mprintf("\tSetting maximum lag to 1/2 total time (%i)\n", maxlag);
  } else if (maxLag_ > (int)Vel_[0].Size()) {
    maxlag = (int)Vel_[0].Size();
    mprintf("\tSpecified maximum lag > total length, setting to %i\n", maxlag);
  } else
    maxlag = maxLag_;
  // DEBUG
  //for (VelArray::iterator vel = Vel_.begin(); vel != Vel_.end(); ++vel) {
  //  mprintf("Vector %u:\n", vel - Vel_.begin());
  //  for (DataSet_Vector::iterator vec = vel->begin(); vec != vel->end(); ++vec)
  //    mprintf("\t%u {%f %f %f}\n", vec - vel->begin(), (*vec)[0], (*vec)[1], (*vec)[2]);
  //}
  // Allocate space for output correlation function values.
  DataSet_double& Ct = static_cast<DataSet_double&>( *VAC_ );
  Ct.Resize( maxlag );
  if (!useFFT_) {
    // DIRECT METHOD 
    ParallelProgress progress( maxlag );
    int t;
    unsigned int dtmax, dt;
#   ifdef _OPENMP
#   pragma omp parallel private(t, dtmax, dt) firstprivate(progress)
    {
      progress.SetThread(omp_get_thread_num());
#     pragma omp for schedule(dynamic)
#   endif
      for (t = 0; t < maxlag; ++t)
      {
        progress.Update( t );
        dtmax = Vel_[0].Size() - t;
        for (dt = 0; dt < dtmax; ++dt)
        {
          for (VelArray::const_iterator vel = Vel_.begin(); vel != Vel_.end(); ++vel)
            Ct[t] += (*vel)[dt] * (*vel)[dt + t];
        }
        Ct[t] /= (double)(dtmax * Vel_.size());
        //mprintf("\tCt[%i]= %f\n", t, Ct[t]); // DEBUG
      }
#   ifdef _OPENMP
    } // END pragma omp parallel
#   endif
    progress.Finish();
  } else {
    // FFT METHOD
    // Since FFT is cyclic, unroll vectors into a 1D array; in the resulting
    // transformed array after FFT, every 3rd value will be the correlation
    // via dot products that we want (once it is normalized).
    unsigned int total_length = Vel_[0].Size() * 3;
    CorrF_FFT pubfft( total_length );
    ComplexArray data1 = pubfft.Array();
    //mprintf("Complex Array Size is %i (%i actual)\n", data1.size(), data1.size()*2);
    ProgressBar progress( Vel_.size() );
    unsigned int nvel = 0;
    for (VelArray::iterator vel = Vel_.begin(); vel != Vel_.end(); ++vel, ++nvel)
    {
      //mprintf("Vector %u\n", vel - Vel_.begin()); // DEBUG
      progress.Update( nvel );
      // Place vector from each frame into 1D array
      unsigned int nd = 0; // Will be used to index complex data
      for (DataSet_Vector::const_iterator vec = vel->begin(); vec != vel->end(); ++vec, nd+=6)
      {
        //mprintf("\tFrame %u assigned to complex array %u\n", vec - vel->begin(), nd); // DEBUG
        data1[nd  ] = (*vec)[0]; data1[nd+1] = 0.0;
        data1[nd+2] = (*vec)[1]; data1[nd+3] = 0.0;
        data1[nd+4] = (*vec)[2]; data1[nd+5] = 0.0;
        //mprintf("\t  Complex[%u]= %f  Complex[%u]= %f  Complex[%u]= %f\n", // DEBUG
        //        nd, data1[nd], nd+2, data1[nd+2], nd+4, data1[nd+4]);
      }
      data1.PadWithZero( total_length );
      //mprintf("Before FFT:\n"); // DEBUG
      //for (unsigned int cd = 0; cd < (unsigned int)data1.size()*2; cd += 2)
      //  mprintf("\t%u: %f + %fi\n", cd/2, data1[cd], data1[cd+1]);
      pubfft.AutoCorr( data1 );
      //mprintf("Results of FFT:\n"); // DEBUG
      //for (unsigned int cd = 0; cd < (unsigned int)data1.size()*2; cd += 2)
      //  mprintf("\t%u: %f + %fi\n", cd/2, data1[cd], data1[cd+1]);
      // Increment nd by 3 here so it can be used for normalization.
      nd = 0;
      for (int t = 0; t < maxlag; t++, nd += 3) {
        //mprintf("\tdata1[%u] = %f", nd*2, data1[nd*2]); // DEBUG
        //mprintf("  norm= %u", total_length - nd); // DEBUG
        //Ct[t] += data1[nd*2] * 3.0 / (double)(total_length - nd);
        Ct[t] += data1[nd*2];
        //mprintf("  Ct[%i]= %f\n", t, Ct[t]); // DEBUG
      }
    }
    // Normalization
    for (int t = 0, nd = 0; t < maxlag; t++, nd += 3)
      Ct[t] *= ( 3.0 / (double)((total_length - nd) * Vel_.size()) );
      //Ct[t] /= (double)Vel_.size();
  }
  // Integration to get diffusion coefficient.
  mprintf("\tIntegrating data set %s, step is %f\n", VAC_->Legend().c_str(), VAC_->Dim(0).Step());
  DataSet_Mesh mesh;
  mesh.SetMeshXY( static_cast<DataSet_1D const&>(*VAC_) );
  double total = mesh.Integrate_Trapezoid();
  const double ANG2_PS_TO_CM2_S = 10.0 / 6.0;
  mprintf("\t3D= %g Å^2/ps, %g x10^-5 cm^2/s\n", total, total * ANG2_PS_TO_CM2_S);
  mprintf("\t D= %g Å^2/ps, %g x10^-5 cm^2/s\n", total / 3.0, total * ANG2_PS_TO_CM2_S / 3.0);
  mprintf("\t6D= %g Å^2/ps, %g x10^-5 cm^2/s\n", total * 2.0, total * ANG2_PS_TO_CM2_S * 2.0);
  if (normalize_) {
    // Normalize VAC fn to 1.0
    mprintf("\tNormalizing VAC function to 1.0, C[0]= %g\n", Ct[0]);
    double norm = 1.0 / Ct[0];
    for (int t = 0; t < maxlag; ++t)
      Ct[t] *= norm;
  }
}
