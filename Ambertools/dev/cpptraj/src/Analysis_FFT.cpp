#include <cmath> // sqrt
#include "Analysis_FFT.h"
#include "CpptrajStdio.h"
#include "PubFFT.h"

// CONSTRUCTOR
Analysis_FFT::Analysis_FFT() : dt_(0.0) {}

void Analysis_FFT::Help() {
  mprintf("\t<dset0> [<dset1> ...] [out <outfile>] [name <outsetname>] [dt <samp_int>]\n"
          "  Perform forward fast-Fourier transformation of data set(s)\n");
}

// Analysis_FFT::Setup()
Analysis::RetType Analysis_FFT::Setup(ArgList& analyzeArgs, DataSetList* datasetlist,
                            TopologyList* PFLin, DataFileList* DFLin, int debugIn)
{
  std::string setname = analyzeArgs.GetStringKey("name");
  DataFile* outfile = DFLin->AddDataFile(analyzeArgs.GetStringKey("out"), analyzeArgs);
  dt_ = analyzeArgs.getKeyDouble("dt",1.0);
  // Select datasets from remaining args
  if (input_dsets_.AddSetsFromArgs( analyzeArgs.RemainingArgs(), *datasetlist )) {
    mprinterr("Error: Could not add data sets.\n");
    return Analysis::ERR;
  }
  if (input_dsets_.empty()) {
    mprinterr("Error: No input data sets.\n");
    return Analysis::ERR;
  }
  // If setname is empty generate a default name
  if (setname.empty())
    setname = datasetlist->GenerateDefaultName( "FFT" );
  // Setup output datasets.
  int idx = 0;
  if ( input_dsets_.size() == 1 )
    idx = -1; // Only one input set, no need to refer to it by index
  for ( Array1D::const_iterator DS = input_dsets_.begin(); 
                                DS != input_dsets_.end(); ++DS) 
  {
    DataSet* dsout = datasetlist->AddSetIdx( DataSet::DOUBLE, setname, idx++ );
    if (dsout==0) return Analysis::ERR;
    dsout->SetLegend( (*DS)->Legend() );
    output_dsets_.push_back( (DataSet_1D*)dsout );
    if (outfile != 0) outfile->AddSet( dsout );
  }

  mprintf("    FFT: Calculating FFT for %u data sets.\n", input_dsets_.size());
  mprintf("\tTime step: %f\n", dt_);
  if ( !setname.empty() )
    mprintf("\tSet name: %s\n", setname.c_str() );
  if ( outfile != 0 )
    mprintf("\tOutfile name: %s\n", outfile->DataFilename().base());

  return Analysis::OK;
}

/** Calculate FFT for input DataSets. FFT magnitude is reported. Magnitude
  * is normalized by N / 2. Only data up to the Nyquist frequency is used. 
  */
//TODO: Deal with vectors
Analysis::RetType Analysis_FFT::Analyze() {
  // Ensure input data sets have the same number of points.
  size_t maxsize_ = 0;
  std::vector<bool> skipSet(input_dsets_.size(), true);
  std::vector<bool>::iterator skip = skipSet.begin();
  for (Array1D::const_iterator DS = input_dsets_.begin();
                               DS != input_dsets_.end(); 
                             ++DS, ++skip)
  {
    // Check for empty set
    if ( (*DS)->Empty() ) {
      mprintf("Warning: Set %s is empty, skipping.\n", (*DS)->Legend().c_str() );
      continue;
    }
    if ( maxsize_ == 0 )
      maxsize_ = (*DS)->Size();
    else if ( (*DS)->Size() != maxsize_ ) {
      mprintf("Warning: Set %s does not have same size (%u) as initial set (%u). Skipping.\n",
              (*DS)->Legend().c_str(), (*DS)->Size(), maxsize_ );
      continue;
    }
    *skip = false;
  }
  // Setup FFT
  PubFFT pubfft;
  pubfft.SetupFFTforN( maxsize_ );
  //mprintf("DEBUG: FFT size is %i\n",pubfft.size());
  // Set up complex number array
  ComplexArray data1( pubfft.size() );

  double sr = 1.0 / dt_;              // 1 / sampling interval, sampling rate (freq)
  double fnyquist = sr / 2.0;         // Nyquist frequency
  double total_time = dt_ * (double)maxsize_; // Total time (fundamental period)
  double f0 = 1.0 / total_time;       // Fundamental frequency (first harmonic)
  Dimension Xdim(0.0, f0, maxsize_, "Freq.");
  mprintf("\tReporting FFT magnitude, normalized by N/2.\n"
          "\tOnly data up to the Nyquist frequency will be used.\n");
  mprintf("\tSampling rate= %f ps^-1, Nyquist freq.= %f ps^-1\n", sr, fnyquist);
  mprintf("\tPoints= %zu, Fundamental period= %f ps, fundamental freq.= %f ps^-1\n", 
          maxsize_, total_time, f0);
  double norm = (double)maxsize_ / 2;

  skip = skipSet.begin(); 
  Array1D::const_iterator dsout = output_dsets_.begin();
  for (Array1D::const_iterator DS = input_dsets_.begin(); 
                               DS != input_dsets_.end();
                             ++DS, ++dsout, ++skip)
  {
    if (*skip) continue;
    mprintf("\t\tCalculating FFT for set %s\n", (*DS)->Legend().c_str());
    // Reset data1 so it is padded with zeros
    data1.PadWithZero(0);
    // Place data from DS in real spots in data1
    int datasize =  (*DS)->Size();
    for (int i = 0; i < datasize; ++i)
      data1[i*2] = ((DataSet_1D*)(*DS))->Dval(i);
    // DEBUG
    //for (int i = 0; i < pubfft.size()*2; i+=2)
    //  mprintf("\t\t\t%i FFTinR=%f  FFTinI=%f\n",i/2,data1[i],data1[i+1]);
    // Perform FFT
    pubfft.Forward( data1 );
    // Place real data from FFT in output Data up to the Nyquist frequency
    int i2 = 0;
    for (int i1 = 0; i1 < datasize; ++i1) {
     double freq = i1 * f0;
     if (freq > fnyquist) break;
     double magnitude = sqrt(data1[i2]*data1[i2] + data1[i2+1]*data1[i2+1]);
     magnitude /= norm;
     //mprintf("\t\t\tReal=%f  Img=%f  Mag=%f\n",data1[i2],data1[i2+1],magnitude);
     (*dsout)->Add( i1, &magnitude );
     i2 += 2;
    }
    (*dsout)->SetDim(Dimension::X, Xdim);
  }
  return Analysis::OK;
}
