#ifndef INC_PROGRESSBAR_H
#define INC_PROGRESSBAR_H
/// Used to print progress to screen
class ProgressBar {
  public:
    ProgressBar();
    ProgressBar(int);
    void SetupProgress(int);
    void Update(int);
  private:
    int unknown_;
    int max_;
    float C_over_max_;
    float targetPercent_;
    bool unknownframes_;
};
/// Used to track progress in parallel
class ParallelProgress {
  public:
    ParallelProgress()      : C_over_max_(1.0), tgt_(0.0), thread_(0) {}
    ParallelProgress(int m) : C_over_max_(100.0/(float)m), tgt_(0.0), thread_(0) {}
    ParallelProgress(const ParallelProgress&);
    void SetThread(int t) { thread_ = t;                       }
    void Update(int it)   { if (thread_==0) printProgress(it); }
    void Finish();
  private:
    void printProgress(int);
    float C_over_max_;
    float tgt_;
    int thread_;
}; 
#endif
