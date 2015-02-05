#ifndef INC_DIMENSION_H
#define INC_DIMENSION_H
#include <string>
/// Holds information about a coordinate dimension.
class Dimension {
  public:
    enum DimIdxType { X = 0, Y, Z };
    Dimension();
    /// CONSTRUCTOR - Min, step, bins
    Dimension(double, double, int);
    /// CONSTRUCTOR - Min, step, bins, label
    Dimension(double, double, int, std::string const&);
    Dimension(const Dimension&);
    Dimension& operator=(const Dimension&);
 
    void SetLabel(std::string const& l) { label_ = l;  }
    void SetMin(double m)               { min_ = m; minIsSet_ = true; }
    void SetMax(double m)               { max_ = m; maxIsSet_ = true; }
    void SetStep(double s)              { step_ = s;   }
    void SetBins(int b)                 { bins_ = b;   }

    std::string const& Label() const { return label_;    }
    double Min()               const { return min_;      }
    double Max()               const { return max_;      }
    double Step()              const { return step_;     }
    int Bins()                 const { return bins_;     }
    bool MinIsSet()            const { return minIsSet_; }
    bool MaxIsSet()            const { return maxIsSet_; }
    double Coord(size_t i)     const { return ((step_ * (double)i) + min_); }
    inline bool operator!=(const Dimension&) const;
    /// Attempt to set up bins or step.
    int CalcBinsOrStep();
    void PrintDim() const;
  private:
    std::string label_;
    double min_;
    double max_;
    double step_;
    int bins_;
    bool minIsSet_;     ///< True if SetMin has been called.
    bool maxIsSet_;     ///< True if SetMax has been called.
};
// ----- INLINE FUNCTIONS ------------------------------------------------------
// FIXME: Check everything??
bool Dimension::operator!=(const Dimension& rhs) const {
  if (min_ != rhs.min_ || step_ != rhs.step_) return true;
  return false;
}
#endif
