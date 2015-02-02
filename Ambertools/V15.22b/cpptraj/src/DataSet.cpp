// DataSet
#include "DataSet.h"
#include "CpptrajStdio.h"
#include "StringRoutines.h" // SetStringFormatString etc

// CONSTRUCTOR
DataSet::DataSet() :
  data_format_(0),
  idx_(-1),
  dType_(UNKNOWN_DATA),
  colwidth_(0),
  width_(0),
  precision_(0),
  leftAlign_(false),
  scalarmode_(UNKNOWN_MODE),
  scalartype_(UNDEFINED)
{ }

/// CONSTRUCTOR - Take type, width, precision, and dimension
DataSet::DataSet(DataType typeIn, int widthIn, int precisionIn, int dimIn) :
  data_format_(0),
  idx_(-1),
  dType_(typeIn),
  dim_(dimIn),
  colwidth_(widthIn),
  width_(widthIn),
  precision_(precisionIn),
  leftAlign_(false),
  scalarmode_(UNKNOWN_MODE),
  scalartype_(UNDEFINED)
{
  SetDataSetFormat(leftAlign_);
}

// COPY CONSTRUCTOR
DataSet::DataSet(const DataSet& rhs) :
  data_format_(0),
  name_(rhs.name_),
  idx_(rhs.idx_),
  aspect_(rhs.aspect_),
  legend_(rhs.legend_),
  dType_(rhs.dType_),
  dim_(rhs.dim_),
  colwidth_(rhs.colwidth_),
  width_(rhs.width_),
  precision_(rhs.precision_),
  leftAlign_(rhs.leftAlign_),
  format_(rhs.format_),
  scalarmode_(rhs.scalarmode_),
  scalartype_(rhs.scalartype_)
{
  if (!format_.empty())
    data_format_ = format_.c_str();
}

// ASSIGNMENT
DataSet& DataSet::operator=(const DataSet& rhs) {
  if (this == &rhs) return *this;
  name_ = rhs.name_;
  idx_ = rhs.idx_;
  aspect_ = rhs.aspect_;
  legend_ = rhs.legend_;
  dType_ = rhs.dType_;
  dim_ = rhs.dim_;
  colwidth_ = rhs.colwidth_;
  width_ = rhs.width_;
  precision_ = rhs.precision_;
  leftAlign_ = rhs.leftAlign_;
  format_ = rhs.format_;
  if (!format_.empty()) 
    data_format_ = format_.c_str();
  scalarmode_ = rhs.scalarmode_;
  scalartype_ = rhs.scalartype_;
  return *this;
}

// DataSet::SetWidth()
/** Set only DataSet width */
void DataSet::SetWidth(int widthIn) {
  width_ = widthIn;
  SetDataSetFormat( leftAlign_ );
}

// DataSet::SetPrecision()
/** Set DataSet width and precision; recalc. output format string.
  */
void DataSet::SetPrecision(int widthIn, int precisionIn) {
  width_ = widthIn;
  precision_ = precisionIn;
  SetDataSetFormat( leftAlign_ );
}

// DataSet::SetupSet()
/** Set up DataSet name and optionally index and/or aspect. Also create
  * default legend if not already set.
  * \param nameIn The DataSet name; should be unique and not empty.
  * \param idxIn DataSet index; if no index should be -1.
  * \param aspectIn DataSet aspect; if no aspect should be empty.
  */
int DataSet::SetupSet(std::string const& nameIn, int idxIn, std::string const& aspectIn)
{
  // Dataset name
  if (nameIn.empty()) {
    mprinterr("Internal Error: DataSet has no name.\n");
    return 1;
  }
  name_ = nameIn;
  // Set index and aspect if given
  if (idxIn != -1) idx_ = idxIn;
  if (!aspectIn.empty()) aspect_ = aspectIn;
  // If no legend set yet create a default one. Possible formats are:
  //  - Name[Aspect]
  //  - Name:idx
  //  - Aspect:Idx
  //  - Name
  if (legend_.empty()) {
    if (!aspect_.empty() && idx_ == -1)
      legend_ = name_ + "[" + aspect_ + "]";
    else if (aspect_.empty() && idx_ != -1)
      legend_ = name_ + ":" + integerToString( idx_ );
    else if (!aspect_.empty() && idx_ != -1)
      legend_ = aspect_ + ":" + integerToString( idx_ );
    else
      legend_ = name_;
  }
  return 0;
}

// DataSet::SetDataSetFormat()
/** Sets the output format strings for DataSet data and name.
  * \param leftAlign if true the data and header will be left-aligned,
  *        otherwise they will be preceded by a space.
  * \return 0 on success, 1 on error.
  */
// TODO: Each data set has their own 
int DataSet::SetDataSetFormat(bool leftAlignIn) {
  leftAlign_ = leftAlignIn;
  // Set data format string.
  // NOTE: According to C++ std 4.7/4 (int)true == 1
  colwidth_ = width_ + (int)(!leftAlign_);
  switch (dType_) {
    case MODES :
    case REMLOG:
    case MATRIX_DBL:
    case XYMESH:
    case TRAJ   :
    case DOUBLE : format_ = SetDoubleFormatString(width_, precision_, 0); break;
    case MATRIX_FLT:
    case GRID_FLT  :
    case COORDS : 
    case FLOAT  : format_ = SetDoubleFormatString(width_, precision_, 1); break;
    case INTEGER: format_ = SetIntegerFormatString(width_); break;
    case STRING : format_ = SetStringFormatString(width_, leftAlign_); break;
    case VECTOR:
      format_ = SetDoubleFormatString(width_, precision_, 0); 
      colwidth_ = (width_ + 1) * 6; // Vx Vy Vz Ox Oy Oz
      break;
    default:
      mprinterr("Error: No format string defined for this data type (%s).\n", 
                Legend().c_str());
      return 1;
  }
  // If we are not left-aligning prepend a space to the format string.
  if (!leftAlign_) format_ = " " + format_;
  // Assign format to a constant ptr to avoid continuous calls to c_str
  data_format_ = format_.c_str();
  return 0;
}

// DataSet::Matches()
bool DataSet::Matches( std::string const& dsname, int idxnum, std::string const& aspect ) const
{
  /*mprintf("DEBUG: Input: %s[%s]:%i  This Set: %s[%s]:%i\n",
          dsname.c_str(), aspect.c_str(), idxnum, 
          name_.c_str(), aspect_.c_str(), idx_);*/
  if ( dsname != name_ && dsname != "*") return false;
  // Currently match any index if not specified.
  if (idxnum != -1 && idxnum != idx_) return false;
  // If aspect specified make sure it matches. 
  if (!aspect.empty() && (aspect != aspect_ && aspect != "*")) return false;
  // If no aspect specified but dataset has aspect do not match.
  if (aspect.empty() && !aspect_.empty()) return false;
  //mprintf("\tMATCH\n");
  return true;
}

const char* DataSet::Smodes[] = {"distance","angle","torsion","pucker","rms",0};
const char* DataSet::Stypes[] = {"alpha","beta","gamma",
  "delta","epsilon","zeta","pucker","chi","h1p","c2p",
  "phi","psi","pchi","omega","noe",0};
const DataSet::scalarMode DataSet::TypeModes[] = {M_TORSION,M_TORSION,M_TORSION,
  M_TORSION,M_TORSION,M_TORSION,M_PUCKER,M_TORSION,M_TORSION,M_TORSION,
  M_TORSION,M_TORSION,M_TORSION,M_TORSION,M_DISTANCE,UNKNOWN_MODE}; 

// DataSet::ScalarDescription()
void DataSet::ScalarDescription() const {
  if (scalarmode_ != UNKNOWN_MODE) {
    mprintf(", %s", Smodes[scalarmode_]);
    if (scalartype_ != UNDEFINED)
      mprintf("(%s)", Stypes[scalartype_]);
  }
}

// DataSet::ModeFromKeyword()
DataSet::scalarMode DataSet::ModeFromKeyword(std::string const& key) {
  for (int i = 0; i != (int)UNKNOWN_MODE; i++)
    if (key.compare( Smodes[i] ) == 0) return (scalarMode)i;
  return UNKNOWN_MODE;
}

DataSet::scalarType DataSet::TypeFromKeyword(std::string const& key, scalarMode const& mIn)
{
  scalarMode dm = mIn;
  return TypeFromKeyword(key, dm);
}

// DataSet::TypeFromKeyword()
DataSet::scalarType DataSet::TypeFromKeyword(std::string const& key, scalarMode& modeIn) {
  for (int i = 0; i != (int)UNDEFINED; i++)
    if (key.compare( Stypes[i] ) == 0) {
      if (modeIn != UNKNOWN_MODE) {
        // Is type valid for given mode?
        if (modeIn != TypeModes[i]) {
          mprinterr("Error: Type '%s' not valid for mode '%s'\n",
                    Stypes[i], Smodes[TypeModes[i]]);
          return UNDEFINED;
        }
      } else
        modeIn = TypeModes[i];
      return (scalarType)i;
    }
  return UNDEFINED;
}
