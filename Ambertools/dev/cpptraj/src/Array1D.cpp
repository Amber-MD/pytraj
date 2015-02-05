#include <algorithm> // sort
#include "Array1D.h"
#include "CpptrajStdio.h"

// COPY CONSTRUCTOR
Array1D::Array1D(const Array1D& rhs) : array_(rhs.array_) {}

// ASSIGNMENT
Array1D& Array1D::operator=(const Array1D& rhs) {
  if (this == &rhs) return *this;
  array_ = rhs.array_;
  return *this;
}

// CONSTRUCTOR
Array1D::Array1D(DataSetList const& SetList) {
  AddDataSets( SetList );
  if (array_.empty())
    mprinterr("Internal Error: No 1D data sets present.");
}

// Array1D::push_back()
int Array1D::push_back( DataSet_1D* const& val ) {
  // Save blank pointers, no error.
  //FIXME: This is only done for DataIO_Std reads
  if (val == 0)
    array_.push_back( val );
  else if (val->Ndim() == 1)
    array_.push_back( val );
  else
    return 1;
  return 0;
}

void Array1D::SortArray1D() {
  std::sort( array_.begin(), array_.end(), DataSet::DS_PtrCmp() );
}

// Array1D::AddDataSets()
int Array1D::AddDataSets(DataSetList const& SetList) {
  for (DataSetList::const_iterator ds = SetList.begin(); ds != SetList.end(); ++ds)
    if ( push_back( (DataSet_1D*)*ds ) ) {
      mprinterr("Error: Only 1D data sets allowed.");
      array_.clear();
      return 1;
    }
  return 0;
}

// Array1D::AddTorsionSets()
int Array1D::AddTorsionSets(DataSetList const& SetList) {
  for (DataSetList::const_iterator ds = SetList.begin(); ds != SetList.end(); ++ds) {
    // Ensure data sets are 1D and periodic
    if ( (*ds)->Ndim() == 1 ) {
      DataSet_1D* ds1 = (DataSet_1D*)(*ds);
      if ( ds1->IsTorsionArray() )
        array_.push_back( ds1 );
      else
        mprintf("Warning: Set '%s' is not periodic, skipping.\n", (*ds)->Legend().c_str());
    } else
      mprintf("Warning: Set '%s' is not 1D, skipping.\n", (*ds)->Legend().c_str());
  }
  return 0;
}

// Array1D::AddSetsFromArgs()
int Array1D::AddSetsFromArgs(ArgList const& dsetArgs, DataSetList const& DSLin) {
  DataSetList input_dsl;
  for (ArgList::const_iterator dsa = dsetArgs.begin(); dsa != dsetArgs.end(); ++dsa)
    input_dsl += DSLin.GetMultipleSets( *dsa );
  if (input_dsl.empty()) {
    mprinterr("Error: No data sets selected.\n");
    return 1;
  }
  // Add to main list
  array_.clear();
  if (AddDataSets( input_dsl ))
    return 1;
  return 0;
}

// Array1D::DetermineMax() 
size_t Array1D::DetermineMax() const {
  size_t maxFrames = 0L;
  for (std::vector<DataSet_1D*>::const_iterator set = array_.begin(); set != array_.end(); ++set)
    if ( (*set)->Size() > maxFrames )
      maxFrames = (*set)->Size();
  return maxFrames;
}

// Array1D::CheckXDimension()
int Array1D::CheckXDimension() const {
  int err = 0;
  Dimension const& Xdim = static_cast<Dimension const&>(array_[0]->Dim(0));
  for (std::vector<DataSet_1D*>::const_iterator set = array_.begin(); set != array_.end(); ++set)
  {
    if ((*set)->Dim(0) != Xdim) {
      mprinterr("Error: X Dimension of %s != %s\n", (*set)->Legend().c_str(),
                array_[0]->Legend().c_str());
      mprinterr("Error:  %s: Min=%f Step=%f\n", (*set)->Legend().c_str(),
                (*set)->Dim(0).Min(), (*set)->Dim(0).Step());
      mprinterr("Error:  %s: Min=%f Step=%f\n", array_[0]->Legend().c_str(),
                Xdim.Min(), Xdim.Step());
      ++err;
    }
  }
  return err;
}
