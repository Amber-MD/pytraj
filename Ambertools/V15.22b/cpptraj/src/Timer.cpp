#ifdef TIMER
# include <time.h>
#else
# include <cstddef>
# include <sys/time.h>
#endif
#ifdef _MSC_VER
   // tw struct timeval is in there on windows
#  include <Winsock2.h>
#endif
#include "Timer.h"
#include "CpptrajStdio.h"

Timer::Timer() : start_sec_(0), start_ns_(0), total_(0.0) {}

#ifdef _MSC_VER 
/* tw 
 * replacement for gettimeofday based on  
 * http://www.cpp-programming.net/c-tidbits/gettimeofday-function-for-windows/
 */
#if defined(_MSC_VER) || defined(_MSC_EXTENSIONS)
#define DELTA_EPOCH_IN_MICROSECS  11644473600000000Ui64
#else
#define DELTA_EPOCH_IN_MICROSECS  11644473600000000ULL
#endif

struct timezone
{
        int  tz_minuteswest; /* minutes W of Greenwich */
        int  tz_dsttime;     /* type of dst correction */
};

int gettimeofday(struct timeval *tv, struct timezone *tz)
{
        FILETIME ft;
        unsigned __int64 tmpres = 0;
        static int tzflag;

        if (NULL != tv)
        {
                GetSystemTimeAsFileTime(&ft);

                tmpres |= ft.dwHighDateTime;
                tmpres <<= 32;
                tmpres |= ft.dwLowDateTime;

                /*converting file time to unix epoch*/
                tmpres /= 10;  /*convert into microseconds*/
                tmpres -= DELTA_EPOCH_IN_MICROSECS;
                tv->tv_sec = (long)(tmpres / 1000000UL);
                tv->tv_usec = (long)(tmpres % 1000000UL);
        }
        return 0;
}
/* end time.h from tw */
#endif

void Timer::GetWallTime(int& sec, int& ns) {
# ifdef TIMER
  struct timespec wall_time;
  clock_gettime(CLOCK_REALTIME, &wall_time);
  sec = wall_time.tv_sec;
  ns = wall_time.tv_nsec;
# else
  struct timeval wall_time;
  gettimeofday(&wall_time, NULL);
  sec = wall_time.tv_sec;
  ns = wall_time.tv_usec;
# endif
};

void Timer::Stop() {
  int stop_sec, stop_ns;
  GetWallTime(stop_sec,  stop_ns);
  double seconds = (double)( stop_sec - start_sec_ );
  double nano = (double)( stop_ns - start_ns_ );
# ifdef TIMER
  total_ += ( seconds + (nano / 1000000000) );
# else
  total_ += ( seconds + (nano / 1000000) );
#endif
}

void Timer::WriteTiming(int indents, const char* header, double FracTotal) const
{
  mprintf("TIME:");
  for (int i = 0; i < indents; i++)
    mprintf("\t");
  mprintf("%s %.4f s", header, total_);
  if (FracTotal > 0.0)
    mprintf(" (%.2f%%)", (total_ / FracTotal) * 100.0);
  mprintf("\n");
}
