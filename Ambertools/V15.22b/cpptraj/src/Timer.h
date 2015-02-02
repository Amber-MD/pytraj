#ifndef INC_TIMER_H
#define INC_TIMER_H
/// Class used to get timing information.
/** Under the hood the function gettimeofday is used by default, which has
  * microsecond resolution. If TIMER is defined, clock_gettime is used instead
  * which has nanosecond resolution and is superior to gettimeofday in many 
  * ways, but requires linking to another library for some versions of GLIB
  * an so is less portable.
  */
class Timer {
  public:
    Timer();
    void Start() { GetWallTime(start_sec_, start_ns_); }
    void Stop();
    double Total() const { return total_; }
    void WriteTiming(int, const char*, double) const;
    void WriteTiming(int i, const char* h) const {
      return WriteTiming(i, h, 0.0);
    }
  private:
    void GetWallTime(int&, int&);
    int start_sec_;
    int start_ns_;
    double total_;
};
#endif
