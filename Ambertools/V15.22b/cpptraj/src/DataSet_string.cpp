// DataSet_string
#include "DataSet_string.h"
#include "MpiRoutines.h"

// DataSet_string::Allocate()
/** Reserve space in the Data and Frames arrays. */
int DataSet_string::Allocate1D( size_t sizeIn ) {
  Data_.reserve( sizeIn );
  return 0;
}

// DataSet_string::Add()
/** Insert data vIn at frame. */
void DataSet_string::Add(size_t frame, const void* vIn) {
  if (frame > Data_.size())
    Data_.resize( frame, "NoData" );
  std::string Temp( (const char*)vIn );
  // Check string width.
  if ( (int)Temp.size() > Width() )
    SetPrecision(Temp.size(), 0);
  // Always insert at the end
  // NOTE: No check for duplicate frame values.
  Data_.push_back( Temp );
}

// DataSet_string::WriteBuffer()
/** Write data at frame to CharBuffer. If no data for frame write 0.0.
  */
void DataSet_string::WriteBuffer(CpptrajFile &cbuffer, size_t frame) const {
  if (frame >= Data_.size())
    cbuffer.Printf(data_format_, "NoData");
  else
    cbuffer.Printf(data_format_, Data_[frame].c_str());
}

// DataSet_string::Sync()
/** First, non-master threads convert their vectors into C-arrays.
  * These arrays are then sent to the master, where they are put 
  * into the master arrays. It is assumed that master (rank 0) has 
  * first chunk of data, rank 1 has next and so on.
  */
int DataSet_string::Sync() {
  unsigned int dataSize;
  unsigned int masterStringSize = 0;
  unsigned int stringSize;
  char* values = 0;

  if (worldsize==1) return 0;

  for ( int rank = 1; rank < worldsize; ++rank) {
    if ( worldrank == rank ) {
      // ----- RANK -------
      // Get size of data on rank.
      dataSize = Data_.size();
      // Send rank size to master
      parallel_sendMaster(&dataSize, 1, rank, PARA_INT);
      // If size is 0 on rank, skip this rank.
      if (dataSize == 0) continue;
      // Get sum size of each string on rank (incl. null char).
      stringSize = 0;
      for ( std::vector<std::string>::iterator str_it = Data_.begin(); 
                                          str_it != Data_.end(); ++str_it)
        stringSize += ( (*str_it).size() + 1 ); // +1 for null char.
      // Send sum string size to master
      parallel_sendMaster(&stringSize, 1, rank, PARA_INT);
      // Allocate space on rank
      values = new char[ stringSize ];
      // Copy each string (incl. null char) to the char array
      char* ptr = values;
      for ( std::vector<std::string>::iterator str_it = Data_.begin(); 
                                          str_it != Data_.end(); ++str_it) 
      {
        size_t length = (*str_it).copy( ptr, (*str_it).size() + 1 );
        ptr += length;
      }
      // Send arrays to master
      //parallel_sendMaster(frames, dataSize, rank, PARA_INT);
      parallel_sendMaster(values, stringSize, rank, PARA_CHAR);
      // Free arrays on rank
      delete[] values;
    } else if (worldrank == 0) {
      // ----- MASTER -----
      // Master receives size from rank
      parallel_sendMaster(&dataSize, 1, rank, PARA_INT);
      // If size was 0 on rank, skip rank.
      if (dataSize == 0) continue;
      // Master receives sum string size from rank
      parallel_sendMaster(&stringSize, 1, rank, PARA_INT);
      // Reallocate if necessary
      //if (dataSize > masterSize) {
      //  if ( frames != 0 ) delete[] frames;
      //  frames = new int[ dataSize ];
      //  masterSize = dataSize;
      //}
      if (stringSize > masterStringSize) {
        if ( values != 0 ) delete[] values;
        values = new char[ stringSize ];
        masterStringSize = stringSize;
      }
      // Master receives arrays
      //parallel_sendMaster(frames, dataSize, rank, PARA_INT);
      parallel_sendMaster(values, stringSize, rank, PARA_CHAR);
      // Insert frames and values to master arrays
      char* ptr = values;
      for (unsigned int i = 0; i < dataSize; ++i) {
        Data_.push_back( ptr );
        ptr += ( Data_.back().size() + 1 );
      }
    }
  } // End loop over ranks > 0

  // Free master array
  if (worldrank == 0 && values != 0 ) delete[] values;

  return 0;
}
