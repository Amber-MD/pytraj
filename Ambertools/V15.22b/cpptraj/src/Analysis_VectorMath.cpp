#include "Analysis_VectorMath.h"
#include "CpptrajStdio.h"
#include "Constants.h"
#include "DataSet_double.h"

/// Strings corresponding to modes, used in output.
const char* Analysis_VectorMath::ModeString[] = {
  "Dot product", "Angle from dot product", "Cross product" };

// CONSTRUCTOR
Analysis_VectorMath::Analysis_VectorMath() :
  mode_(DOTPRODUCT), vinfo1_(0), vinfo2_(0), DataOut_(0), norm_(false) {}

void Analysis_VectorMath::Help() {
  mprintf("\tvec1 <vecname1> vec2 <vecname2> [out <filename>] [norm] [name <setname>]\n"
          "\t[ dotproduct | dotangle | crossproduct ]\n"
          "  Calculate dot product, angle from dot product (degrees), or cross product\n"
          "  for specified vectors.\n");
}

// Analysis_VectorMath::Setup()
Analysis::RetType Analysis_VectorMath::Setup(ArgList& analyzeArgs, DataSetList* DSLin,
                            TopologyList* PFLin, DataFileList* DFLin, int debugIn)
{
  // Get Vectors
  vinfo1_ = (DataSet_Vector*)DSLin->FindSetOfType( analyzeArgs.GetStringKey("vec1"),
                                                   DataSet::VECTOR );
  vinfo2_ = (DataSet_Vector*)DSLin->FindSetOfType( analyzeArgs.GetStringKey("vec2"),
                                                   DataSet::VECTOR );
  if (vinfo1_ == 0 ) {
    mprinterr("Error: 'vec1' not found.\n");
    return Analysis::ERR;
  }
  if (vinfo2_ == 0) {
    mprinterr("Error: 'vec2' not found.\n");
    return Analysis::ERR;
  }
  std::string setname = analyzeArgs.GetStringKey("name");
  norm_ = analyzeArgs.hasKey("norm");
  // Check for dotproduct/crossproduct keywords
  DataOut_ = 0;
  if (analyzeArgs.hasKey("dotproduct")) {
    mode_ = DOTPRODUCT;
    if ((DataOut_ = DSLin->AddSet(DataSet::DOUBLE, setname, "Dot")) == 0) return Analysis::ERR;
  } else if (analyzeArgs.hasKey("dotangle")) {
    mode_ = DOTANGLE;
    norm_ = true; // Vecs must be normalized for angle calc to work
    if ((DataOut_ = DSLin->AddSet(DataSet::DOUBLE, setname, "Angle")) == 0) return Analysis::ERR;
  } else if (analyzeArgs.hasKey("crossproduct")) {
    mode_ = CROSSPRODUCT;
    if ((DataOut_ = DSLin->AddSet(DataSet::VECTOR, setname, "Cross")) == 0) return Analysis::ERR;
  } else
    mode_ = DOTPRODUCT;
  // Set up output file in DataFileList if necessary
  DataFile* outfile = DFLin->AddDataFile( analyzeArgs.GetStringKey("out"), analyzeArgs );
  if (outfile != 0) outfile->AddSet( DataOut_ );

  // Print Status
  mprintf("    VECTORMATH: Calculating %s of vectors %s and %s\n", 
            ModeString[mode_], vinfo1_->Legend().c_str(), vinfo2_->Legend().c_str());
  if (norm_) mprintf("\tVectors will be normalized.\n");
  if (outfile != 0)
    mprintf("\tResults are written to %s\n", outfile->DataFilename().full());

  return Analysis::OK;
}

// Analysis_VectorMath::DotProduct()
int Analysis_VectorMath::DotProduct()
{
  DataSet_double& Out = static_cast<DataSet_double&>( *DataOut_ );
  DataSet_Vector& V1 = static_cast<DataSet_Vector&>( *vinfo1_ );
  DataSet_Vector& V2 = static_cast<DataSet_Vector&>( *vinfo2_ );
  Out.Resize( V1.Size() );
  for (unsigned int v = 0; v < V1.Size(); ++v) {
    if (norm_) {
      V1[v].Normalize();
      V2[v].Normalize();
    }
    if (mode_ == DOTPRODUCT)
      Out[v] = V1[v] * V2[v];
    else // DOTANGLE
      Out[v] = V1[v].Angle( V2[v] ) * Constants::RADDEG;
  }
  return 0;
}

// Analysis_VectorMath::CrossProduct()
int Analysis_VectorMath::CrossProduct()
{
  DataSet_Vector& Out = static_cast<DataSet_Vector&>( *DataOut_ );
  DataSet_Vector& V1 = static_cast<DataSet_Vector&>( *vinfo1_ );
  DataSet_Vector& V2 = static_cast<DataSet_Vector&>( *vinfo2_ );
  Out.ReserveVecs( V1.Size() );
  for (unsigned int v = 0; v < V1.Size(); ++v) {
    if (norm_) {
      V1[v].Normalize();
      V2[v].Normalize();
    }
    Out.AddVxyz( V1[v].Cross( V2[v] ) );
  }
  return 0;
}

// Analysis_VectorMath::Analyze()
Analysis::RetType Analysis_VectorMath::Analyze() {
  // Ensure vectors have the same # of frames
  if (vinfo1_->Size() != vinfo2_->Size()) {
    mprinterr("Error: # Frames in vec %s (%i) != # Frames in vec %s (%i)\n",
              vinfo1_->Legend().c_str(), vinfo1_->Size(),
              vinfo2_->Legend().c_str(), vinfo2_->Size());
    return Analysis::ERR;
  }
  int err = 0;
  if (mode_ == CROSSPRODUCT)
    err = CrossProduct();
  else // DOTPRODUCT || DOTANGLE
    err = DotProduct();
  if (err != 0) return Analysis::ERR;
  return Analysis::OK;
}
