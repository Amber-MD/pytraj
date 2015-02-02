#ifndef INC_HUNGARIAN_H
#define INC_HUNGARIAN_H
#include <vector>
#include "Matrix.h"
/// Use Hungarian algorithm to perform minimum-matching on a matrix.
class Hungarian {
  public:
    Hungarian() : nrows_(0), ncols_(0) {}
    /// Initialize NxN full matrix for Hungarian algorithm.
    int Initialize(size_t);
    /// Add an element to matrix for Hungarian algorithm.
    void AddElement(double d) { matrix_.addElement( d ); }
    /// \return Array containing Map[col] = row
    std::vector<int> Optimize();
  private:
    int AssignRowsToColumns();
    void CoverZeroElements();
    void UpdateMatrix();
#   ifdef DEBUG_HUNGARIAN
    void PrintLines(const char*);
    void PrintMatrix(const char*);
#   endif
    Matrix<double> matrix_;            ///< Working matrix. 
    std::vector<bool> lineThroughRow_; ///< True if specified row is marked.
    std::vector<bool> lineThroughCol_; ///< True if specified col is marked.
    std::vector<int> assignRowToCol_;  ///< map[col] = row
    std::vector<int> assignColToRow_;  ///< map[row] = col
    // TODO: Make these size_t
    int nrows_;                        ///< # of rows in matrix.
    int ncols_;                        ///< # of cols in matrix.
};
#endif
