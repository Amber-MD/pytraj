#ifndef INC_CLUSTERMATRIX_H
#define INC_CLUSTERMATRIX_H
#include <string>
#include <vector>
#include "Matrix.h"
#include "ClusterSieve.h"
// NOTE: This only needs to inherit directly if going on DataSetList.
class ClusterMatrix {
  public:
    ClusterMatrix() {}
    /// Set up matrix with sieve value.
    int SetupWithSieve(size_t, size_t, int);
    ClusterMatrix(const ClusterMatrix&);
    ClusterMatrix& operator=(const ClusterMatrix&);

    int SaveFile(std::string const&) const;
    int LoadFile(std::string const&, int);
    int SetupMatrix(size_t);
    /// Indicate given row/col should be ignored.
    void Ignore(int row)            { ignore_[row] = true;   }
    /// \return true if given row/col has been ignored.
    bool IgnoringRow(int row) const { return ignore_[row];   }
    /// \return Number of frames (original nrows)
    size_t Nframes()          const { return ignore_.size(); }
    /// \return Actual number of rows in matrix.
    size_t Nrows()            const { return Mat_.Nrows();   }
    /// \return An array containing sieved frame numbers.
    ClusterSieve::SievedFrames Sieved() const { return sievedFrames_.Frames(); }
    /// \return Sieve value
    int SieveValue()                    const { return sievedFrames_.Sieve();  }
    /// Set the row and column of the smallest element.
    double FindMin(int&, int&) const;
    void PrintElements() const;
    /// \return an element indexed by sievedFrames.
    inline double GetFdist(int, int) const;
    /// \return an element.
    inline double GetCdist(int c, int r) const { return Mat_.element(c,r); }
    inline void SetElement(int, int, double);
    size_t Nelements()        const { return Mat_.size();               }
    int AddElement(double d)        { return Mat_.addElement((float)d); }
    size_t DataSize() const;
    typedef Matrix<float>::iterator const_iterator;
    const_iterator begin() const { return Mat_.begin(); }
    const_iterator end()   const { return Mat_.end();   }
  private:
    static const unsigned char Magic_[];
    /// For reading/writing 8 byte unsigned integers
    typedef unsigned long long int uint_8;
    /// For reading/writing 8 byte signed integers
    typedef long long int sint_8;
    /// If true, ignore the row/col when printing/searching etc
    std::vector<bool> ignore_;
    Matrix<float> Mat_;
    ClusterSieve sievedFrames_;
};
// Inline functions
double ClusterMatrix::GetFdist(int col, int row) const {
  // row and col are based on original; convert to reduced
  // FIXME: This assumes GetElement will never be called for a sieved frame.
  return Mat_.element(sievedFrames_.FrameToIdx(col), 
                      sievedFrames_.FrameToIdx(row));
}

void ClusterMatrix::SetElement(int col, int row, double val) {
  Mat_.setElement(col, row, val);
}
#endif
