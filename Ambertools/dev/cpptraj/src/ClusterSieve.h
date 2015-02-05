#ifndef INC_CLUSTERSIEVE_H
#define INC_CLUSTERSIEVE_H
#include <vector>
#include <cstddef>
/// Used to map actual frame numbers to ClusterMatrix internal indices.
class ClusterSieve {
  public:
    enum SieveType { NONE=0, REGULAR, RANDOM };
    typedef std::vector<int> SievedFrames;
    ClusterSieve();
    /// Setup no sieve, regular sieve, or random sieve.
    int SetSieve(int, size_t, int);
    /// Setup sieve from previously obtained ignore array.
    int SetSieve(int, std::vector<bool> const&);
    /// \return an array of sieved frame numbers.
    SievedFrames Frames() const;
    /// \return size of data in bytes
    size_t DataSize() const;
    /// \return an array index corresponding to a sieved frame.
    inline int FrameToIdx(int frame) const { return frameToIdx_[frame]; }
    inline size_t MaxFrames()        const { return frameToIdx_.size(); }
    inline int Sieve()               const { return sieve_;             }
    inline SieveType Type()          const { return type_;              }
  private:
    inline void DetermineTypeFromSieve(int);
    SieveType type_;
    int sieve_; ///< Sieve value; > 1 is regular, < -1 is random.
    std::vector<int> frameToIdx_;
};
#endif
