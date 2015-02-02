#include <cstdio> // sscanf
#include "DataIO_Grace.h"
#include "CpptrajStdio.h"
#include "BufferedLine.h"
#include "Array1D.h"

// TODO: Set dimension labels
// Dont assume anything about set ordering
int DataIO_Grace::ReadData(std::string const& fname, ArgList& argIn,
                           DataSetList& datasetlist, std::string const& dsname)
{
  ArgList dataline;
  int setnum = 0;
  int frame = 0;
  DataSet_1D* dset = 0;
  Array1D Dsets;
  std::vector<std::string> labels;
  double dval, xval;
  const char* linebuffer;
  std::vector<double> Xvals;
  
  // Allocate and set up read buffer
  BufferedLine buffer;
  if (buffer.OpenFileRead( fname )) return 1;

  // Read chunks from file
  while ( (linebuffer = buffer.Line()) != 0 ) {
    if (linebuffer[0] == '@') {
      // Command: create command line without the @
      dataline.SetList(linebuffer+1, " \t");
      if ( !dataline.CommandIs("legend") && dataline.Contains("legend") ) {
        // Legend keyword that does not come first. Dont store blanks.
        std::string lbl = dataline.GetStringKey("legend");
        if (!lbl.empty())
          labels.push_back( lbl );
      } else if (dataline.CommandIs("target")) {
        if (dset != 0) {
          // Set was previously allocated. Figure out X dimension.
          dset->SetDim(Dimension::X, DataIO::DetermineXdim(Xvals));
          Xvals.clear();
        }
        // Indicates dataset will be read soon. Allocate new set.
        dset = (DataSet_1D*)datasetlist.AddSetIdx( DataSet::DOUBLE, dsname, setnum++);
        if (dset == 0) {
          mprinterr("Error: %s: Could not allocate data set.\n", buffer.Filename().full());
          return 1;
        }
        Dsets.push_back( dset );
        frame = 0;
      }
    } else if (linebuffer[0] != '#' && linebuffer[0] != '&') { // Skip comments and dataset end
      // Data
      if (dset==0) {
        mprinterr("Error: %s: Malformed grace file. Expected 'target' before data.\n", 
                  buffer.Filename().full());
        return 1;
      }
      sscanf(linebuffer,"%lf %lf", &xval, &dval);
      dset->Add( frame++, &dval );
      Xvals.push_back( xval );
    }
  } // END loop over file
  buffer.CloseFile();
  // Figure out X dimension for last set read
  if (dset != 0)
    dset->SetDim(Dimension::X, DataIO::DetermineXdim(Xvals));

  // Set DataSet legends if specified
  if (!labels.empty()) {
    if (Dsets.size() == labels.size()) {
      mprintf("\tLabels:\n");
      Array1D::const_iterator set = Dsets.begin();
      for (std::vector<std::string>::iterator lbl = labels.begin();
                                              lbl != labels.end(); ++lbl)
      {
        mprintf("\t\t%s\n", (*lbl).c_str());
        (*set)->SetLegend( *lbl );
        ++set;
      }
    } else {
      mprintf("Warning: Number of labels (%zu) does not match number of sets (%zu)\n",
              labels.size(), Dsets.size());
    }
  }
    
  return 0;
}
// -----------------------------------------------------------------------------
void DataIO_Grace::WriteHelp() {
  mprintf("\tinvert: Flip X/Y axes.\n");
}

// DataIO_Grace::processWriteArgs()
int DataIO_Grace::processWriteArgs(ArgList &argIn) {
  isInverted_ = argIn.hasKey("invert");
  return 0;
}

// DataIO_Grace::WriteData()
int DataIO_Grace::WriteData(std::string const& fname, DataSetList const& SetList)
{
  int err = 0;
  // Open output file.
  CpptrajFile file;
  if (file.OpenWrite( fname )) return 1;
  if (isInverted_)
    err = WriteDataInverted(file, SetList);
  else
    err = WriteDataNormal(file, SetList);
  file.CloseFile();
  return err;
}

// DataIO_Grace::WriteDataNormal()
int DataIO_Grace::WriteDataNormal(CpptrajFile& file, DataSetList const& SetList) {
  // Hold all 1D data sets.
  Array1D Sets( SetList );
  if (Sets.empty()) return 1;
  // Determine size of largest DataSet.
  //size_t maxFrames = Sets.DetermineMax();
  // Grace header. Use first data set for labels
  // TODO: DataFile should pass in axis information 
/*  file.Printf(
    "@with g0\n@  xaxis label \"%s\"\n@  yaxis label \"%s\"\n@  legend 0.2, 0.995\n@  legend char size 0.60\n",
    Dim[0].Label().c_str(), Dim[1].Label().c_str()
  );*/
  file.Printf("@with g0\n@  xaxis label \"%s\"\n@  yaxis label \"%s\"\n"
              "@  legend 0.2, 0.995\n@  legend char size 0.60\n",
              Sets[0]->Dim(0).Label().c_str(), "");
  // Loop over DataSets
  unsigned int setnum = 0;
  for (Array1D::const_iterator set = Sets.begin(); set != Sets.end(); ++set) {
    size_t maxFrames = (*set)->Size();
    // Set information
    file.Printf("@  s%u legend \"%s\"\n@target G0.S%u\n@type xy\n",
                   setnum, (*set)->Legend().c_str(), setnum );
    // Setup set X coord format.
    std::string x_col_format = SetupCoordFormat( maxFrames, (*set)->Dim(0), 8, 3 );
    // Write Data for set
    for (size_t frame = 0L; frame < maxFrames; frame++) {
      file.Printf(x_col_format.c_str(), (*set)->Xcrd(frame));
      (*set)->WriteBuffer(file, frame);
      file.Printf("\n");
    }
    ++setnum;
  }
  return 0;
}

// DataIO_Grace::WriteDataInverted()
int DataIO_Grace::WriteDataInverted(CpptrajFile& file, DataSetList const& SetList)
{
  // Hold all 1D data sets.
  Array1D Sets( SetList );
  if (Sets.empty()) return 1;
  // Determine size of largest DataSet.
  size_t maxFrames = Sets.DetermineMax();
  // Grace header. Use first DataSet for axis labels.
  // TODO: DataFile should pass in axis info.
/*  file.Printf(
    "@with g0\n@  xaxis label \"%s\"\n@  yaxis label \"%s\"\n@  legend 0.2, 0.995\n@  legend char size 0.60\n",
    Dim[1].Label().c_str(), Dim[0].Label().c_str() 
  );*/
  file.Printf("@with g0\n@  xaxis label \"%s\"\n@  yaxis label \"%s\"\n"
              "@  legend 0.2, 0.995\n@  legend char size 0.60\n",
              "", Sets[0]->Dim(0).Label().c_str());
  // Setup set X coord format. 
  Dimension Xdim;
  Xdim.SetStep( 1 );
  std::string x_col_format = SetupCoordFormat( Sets.size(), Xdim, 8, 3 );
  // Loop over frames
  for (size_t frame = 0L; frame < maxFrames; frame++) {
    // Set information
    file.Printf("@target G0.S%u\n@type xy\n",frame);
    // Loop over all Set Data for this frame
    unsigned int setnum = 0;
    for (Array1D::const_iterator set = Sets.begin(); set != Sets.end(); ++set, ++setnum) {
      file.Printf(x_col_format.c_str(), Xdim.Coord( setnum ));
      (*set)->WriteBuffer(file, frame);
      file.Printf(" \"%s\"\n", (*set)->Legend().c_str());
    }
  }
  return 0;
}
