// DataSet_integer
#include "DataSet_integer.h"
#include "MpiRoutines.h"

// DataSet_integer::Allocate()
/** Reserve space in the Data and Frames arrays. */
int DataSet_integer::Allocate1D( size_t sizeIn ) {
  Data_.reserve( sizeIn );
  return 0;
}

// DataSet_integer::Add()
/** Insert data vIn at frame. */
void DataSet_integer::Add(size_t frame, const void* vIn) {
  if (frame > Data_.size())
    Data_.resize( frame, 0 );
  // Always insert at the end
  // NOTE: No check for duplicate frame values.
  Data_.push_back( *((int*)vIn) );
}

// DataSet_integer::WriteBuffer()
/** Write data at frame to CharBuffer. If no data for frame write 0.0.
  */
void DataSet_integer::WriteBuffer(CpptrajFile &cbuffer, size_t frame) const {
  if (frame >= Data_.size())
    cbuffer.Printf(data_format_, 0);
  else
    cbuffer.Printf(data_format_, Data_[frame]);
}

// DataSet_integer::Sync()
/** First, non-master threads convert their vectors into C-arrays.
  * These arrays are then sent to the master, where they are put 
  * into the master arrays. It is assumed that master (rank 0) has 
  * first chunk of data, rank 1 has next and so on.
  */
int DataSet_integer::Sync() {
  unsigned int dataSize;
  unsigned int masterSize = 0;
  int* values = 0;

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
      // Allocate space for temp array on rank, put Data_ into values.
      values = new int[ dataSize ];
      std::copy(Data_.begin(), Data_.end(), values);
      //frames = new int[ dataSize ];
      // Send arrays to master
      //parallel_sendMaster(frames, dataSize, rank, PARA_INT);
      parallel_sendMaster(values, dataSize, rank, PARA_INT);
      // Free arrays on rank
      delete[] values;
    } else if (worldrank == 0) {
      // ----- MASTER -----
      // Master receives size from rank
      parallel_sendMaster(&dataSize, 1, rank, PARA_INT);
      // If size was 0 on rank, skip rank.
      if (dataSize == 0) continue;
      // Reallocate temp array on master if necessary
      if (dataSize > masterSize) {
        if ( values != 0 ) delete[] values;
        values = new int[ dataSize ];
        masterSize = dataSize;
      }
      // Master receives arrays
      //parallel_sendMaster(frames, dataSize, rank, PARA_INT);
      parallel_sendMaster(values, dataSize, rank, PARA_INT);
      // Insert frames and values to master arrays
      for (unsigned int i = 0; i < dataSize; ++i) {
        //Frames_.push_back( frames[i] );
        Data_.push_back( values[i] );
      }
    }
  } // End loop over ranks > 0

  // Free master array
  if (worldrank == 0 && values != 0 ) delete[] values;

  return 0;
}
