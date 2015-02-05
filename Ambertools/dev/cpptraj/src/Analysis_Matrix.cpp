#include "Analysis_Matrix.h"
#include "CpptrajStdio.h"
#include "DataSet_MatrixDbl.h"

// CONSTRUCTOR
Analysis_Matrix::Analysis_Matrix() :
  matrix_(0),
  modes_(0),
  thermo_temp_(298.15),
  nevec_(0),
  thermopt_(false),
  reduce_(false),
  eigenvaluesOnly_(false),
  nmwizopt_(false),
  nmwizvecs_(0)
{}

void Analysis_Matrix::Help() {
  mprintf("\t<name> [out <filename>] [thermo [outthermo <filename>] [temp <T>]]\n"
          "\t[vecs <#>] [name <modesname>] [reduce]\n"
          "\t[ nmwiz [nmwizvecs <n>] [nmwizfile <file>] %s\n"
          "\t  nmwizmask <mask> ]\n"
          "  Diagonalize given symmetric matrix to obtain eigenvectors\n"
          "  and eigenvalues.\n", TopologyList::ParmArgs);
}

// Analysis_Matrix::Setup()
Analysis::RetType Analysis_Matrix::Setup(ArgList& analyzeArgs, DataSetList* DSLin,
                            TopologyList* PFLin, DataFileList* DFLin, int debugIn)
{
#ifdef NO_MATHLIB
  mprinterr("Error: Compiled without LAPACK routines.\n");
  return Analysis::ERR;
#else
  // Get matrix name
  std::string mname = analyzeArgs.GetStringNext();
  if (mname.empty()) {
    mprinterr("Error: Missing matrix name (first argument).\n");
    return Analysis::ERR;
  }
  // Find matrix in DataSetList.
  matrix_ = (DataSet_2D*)DSLin->FindSetOfType( mname, DataSet::MATRIX_DBL );
  if (matrix_ == 0)
    matrix_ = (DataSet_2D*)DSLin->FindSetOfType( mname, DataSet::MATRIX_FLT );
  if (matrix_ == 0) {
    mprinterr("Error: Could not find matrix named %s\n",mname.c_str());
    return Analysis::ERR;
  }
  // Check that matrix is symmetric (half-matrix incl. diagonal).
  if (matrix_->Kind() != DataSet_2D::HALF) {
    mprinterr("Error: Only works for symmetric matrices (i.e. no mask2)\n");
    return Analysis::ERR;
  }
  
  // nmwiz flag
  nmwizopt_ = analyzeArgs.hasKey("nmwiz");
  if (nmwizopt_) { 
    nmwizvecs_ = analyzeArgs.getKeyInt("nmwizvecs", 20);
    if (nmwizvecs_ < 1) {
      mprinterr("Error: nmwizvecs must be >= 1\n");
      return Analysis::ERR;
    }
    nmwizfile_ = analyzeArgs.GetStringKey("nmwizfile");
    Topology* parmIn = PFLin ->GetParm( analyzeArgs);
    if (parmIn == 0) {
      mprinterr("Error: nmwiz: No topology specified.\n");
      return Analysis::ERR;
    }
    AtomMask nmwizMask( analyzeArgs.GetStringKey("nmwizmask") );
    if (parmIn->SetupIntegerMask( nmwizMask )) return Analysis::ERR;
    nmwizMask.MaskInfo();
    Topology* nparm = parmIn->partialModifyStateByMask( nmwizMask );
    if (nparm == 0) return Analysis::ERR;
    nmwizParm_ = *nparm;
    delete nparm;
    nmwizParm_.Brief("nmwiz topology");
  }
  
  // Filenames
  DataFile* outfile = DFLin->AddDataFile( analyzeArgs.GetStringKey("out"), analyzeArgs);
  // Thermo flag
  thermopt_ = analyzeArgs.hasKey("thermo");
  outthermo_ = analyzeArgs.GetStringKey("outthermo");
  thermo_temp_ = analyzeArgs.getKeyDouble("temp", 298.15);
  if (thermopt_ && matrix_->Type()!=DataSet_2D::MWCOVAR) {
    mprinterr("Error: Parameter 'thermo' only works for mass-weighted covariance matrix ('mwcovar').\n");
    return Analysis::ERR;
  }
  // Number of eigenvectors; allow "0" only in case of 'thermo'
  nevec_ = analyzeArgs.getKeyInt("vecs",0);
  if (thermopt_) {
    if (nevec_ < 0) nevec_ = 0;
  } else if (nevec_ <= 0) {
    mprintf("Warning: # of eigenvectors specified is < 1 (%i) and 'thermo' not specified.\n",
            nevec_);
    mprintf("Warning: Specify # eigenvectors with 'vecs <#>'. Setting to 1.\n");
    nevec_ = 1;
  }
  // Reduce flag
  reduce_ = analyzeArgs.hasKey("reduce");
  // Set up DataSet_Modes
  std::string modesname = analyzeArgs.GetStringKey("name");
  modes_ = (DataSet_Modes*)DSLin->AddSet( DataSet::MODES, modesname, "Modes" );
  if (modes_==0) return Analysis::ERR;
  // Output string for writing modes file.
  modes_->SetType( matrix_->Type() );
  if (outfile != 0) outfile->AddSet( modes_ );

  // Print Status
  mprintf("    DIAGMATRIX: Diagonalizing matrix %s",matrix_->Legend().c_str());
  if (outfile != 0)
    mprintf(" and writing modes to %s", outfile->DataFilename().full());
  mprintf("\n\tCalculating %i eigenvectors", nevec_);
  if (thermopt_) {
    mprintf(" and thermodynamic data at %.2f K, output to", thermo_temp_);
    if (!outthermo_.empty())
      mprintf(" %s", outthermo_.c_str());
    else
      mprintf(" STDOUT");
  }
  mprintf("\n");
  if (nmwizopt_) {
    mprintf("\tWriting %i modes to NMWiz file", nmwizvecs_);
    if (!nmwizfile_.empty())
      mprintf(" %s\n", nmwizfile_.c_str());
    else
      mprintf(" STDOUT\n");
  }
  if (nevec_>0 && reduce_)
    mprintf("\tEigenvectors will be reduced\n");
  if (!modesname.empty())
    mprintf("\tStoring modes with name: %s\n",modesname.c_str());
  return Analysis::OK;
#endif
}

// Analysis_Matrix::Analyze()
Analysis::RetType Analysis_Matrix::Analyze() {
  modes_->SetAvgCoords( *matrix_ );
  // Calculate eigenvalues / eigenvectors
  if (modes_->CalcEigen( *matrix_, nevec_ )) return Analysis::ERR;
  if (matrix_->Type() == DataSet_2D::MWCOVAR) {
    DataSet_MatrixDbl* Dmatrix = static_cast<DataSet_MatrixDbl*>( matrix_ );
    if ( Dmatrix->Mass().empty() ) {
      mprinterr("Error: MWCOVAR Matrix %s does not have mass info.\n", matrix_->Legend().c_str());
      return Analysis::ERR;
    }
    // Convert eigenvalues to cm^-1
    if (modes_->EigvalToFreq(thermo_temp_)) return Analysis::ERR;
    // Mass-wt eigenvectors // TODO Do not pass in Mass again, done above in SetAvgCoords
    if (modes_->MassWtEigvect( Dmatrix->Mass() )) return Analysis::ERR;
    // Calc thermo-chemistry if specified
    if (thermopt_) {
      CpptrajFile outfile;
      outfile.OpenWrite(outthermo_);
      modes_->Thermo( outfile, 1, thermo_temp_, 1.0 );
      outfile.CloseFile();
    }
  }
  // Print nmwiz file if specified
  if (nmwizopt_) NMWizOutput(); 

  if (reduce_) {
    if (modes_->ReduceVectors()) return Analysis::ERR;
  }
  //modes_->PrintModes(); // DEBUG

  return Analysis::OK;
}

int Analysis_Matrix::NMWizOutput() const {
  CpptrajFile outfile; 
  if (outfile.OpenWrite(nmwizfile_)) return 1;
  // Check # vecs
  int nvecs;
  if (nmwizvecs_ <= modes_->Nmodes())
    nvecs = nmwizvecs_;
  else {
    mprintf("Warning: nmwizvecs > # eigenvectors, only writing %i vecs.\n",
            modes_->Nmodes());
    nvecs = modes_->Nmodes();
  }
  // Check # atoms
  if (nmwizParm_.Natom() * 3 != modes_->VectorSize()) {
    mprinterr("Error: nmwiz topology size %i does not match eigenvector size %i.\n",
              nmwizParm_.Natom() * 3, modes_->VectorSize());
    return 1;
  }
  
  outfile.Printf("nmwiz_load %s\n", nmwizfile_.c_str());

  outfile.Printf("name default_name\n");  //TODO: get from optionally provided pdb file

  outfile.Printf("atomnames ");
  for (Topology::atom_iterator atom = nmwizParm_.begin(); atom != nmwizParm_.end(); ++atom)
        outfile.Printf("%s ", atom->c_str());
  outfile.Printf("\n");

  outfile.Printf("resnames ");
  for (Topology::atom_iterator atom = nmwizParm_.begin(); atom != nmwizParm_.end(); ++atom)
    outfile.Printf("%s ", nmwizParm_.Res(atom->ResNum()).c_str());
  outfile.Printf("\n");

  outfile.Printf("resids ");
  for (Topology::atom_iterator atom = nmwizParm_.begin(); atom != nmwizParm_.end(); ++atom)
    outfile.Printf("%d ", atom->ResNum()+1);
  outfile.Printf("\n");

  outfile.Printf("chainids \n");    //TODO: get from optionally provided pdb file

  outfile.Printf("bfactors \n");    //TODO: get from optionally provided pdb file

  outfile.Printf("coordinates ");
  for (int i = 0; i < modes_->NavgCrd(); ++i)
          outfile.Printf("%8.3f ", modes_->AvgCrd()[i]);
  outfile.Printf("\n");

  for (int vec = 0; vec < nvecs; ++vec){
    outfile.Printf("mode %i %12.10f ", vec+1, 1/modes_->Eigenvalue(vec));
    const double* Vec = modes_->Eigenvector(vec);
    for (int i = 0 ; i < modes_->VectorSize(); ++i)
      outfile.Printf("%12.5f ", Vec[i]);
    outfile.Printf("\n");
  }
  return 0;
}
