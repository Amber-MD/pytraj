#include <cstdio>
#include <cstdarg>
#ifdef MPI
#  include "MpiRoutines.h"
#endif

static bool worldsilent = false; // If true suppress all mprintf output.
// mflush()
/** Call flush on STDOUT only if this is the master thread */
void mflush() {
#ifdef MPI
  if (worldrank!=0) return;
#endif
  fflush(stdout);
}

/** Print message to STDOUT even if worldsilent */
void loudPrintf(const char* format, ...) {
  va_list args;
#ifdef MPI
  if (worldrank!=0) return;
#endif
  va_start(args,format);
  vfprintf(stdout,format,args);
  va_end(args);
}

// mprintf()
/** Print message to STDOUT only if this is the master thread */
void mprintf(const char *format, ...) {
  va_list args;
  if (worldsilent) return;
#ifdef MPI
  if (worldrank!=0) return;
#endif
  va_start(args,format);
  vfprintf(stdout,format,args);
  va_end(args);
}

// mprinterr()
/** Print message to STDERR only if this is the master thread */
void mprinterr(const char *format, ...) {
  va_list args;

#ifdef MPI
  if (worldrank!=0) return; 
#endif
  va_start(args,format);
  vfprintf(stderr,format,args);
  va_end(args);
}

// rprintf()
/** Print message to STDOUT for this worldrank */
void rprintf(const char *format, ...) {
  va_list args;
  if (worldsilent) return;
  va_start(args, format);
#ifdef MPI
  char buffer[1024];
  int nc = sprintf(buffer, "[%i]\t", worldrank);
  nc += vsprintf(buffer + nc, format, args);
  fwrite(buffer, 1, nc, stdout);
# else
  vfprintf(stdout,format,args);
# endif
  va_end(args);
  return;
}

// rprinterr()
/** Print message to STDERR for this worldrank */
void rprinterr(const char *format, ...) {
  va_list args;

#ifdef MPI
  fprintf(stderr,"[%i]\t",worldrank);
#endif
  va_start(args,format);
  vfprintf(stderr,format,args);
  va_end(args);
  return;
}

void SetWorldSilent(bool silentIn) {
  worldsilent = silentIn;
}

// printerr()
/** Print error message along with calling routine.  */
/*void printerr(const char *ROUTINE, const char *format, ...) {
  va_list args;

  va_start(args,format);
  fprintf(stdout,"Error: ");
  if (ROUTINE!=0)
    fprintf(stdout, "%s: ",ROUTINE);
  vfprintf(stdout,format,args);
  va_end(args);
  fprintf(stdout,"\n");
  return;
}*/

// printwar()
/** Print warning message along with calling routine.  */
/*void printwar(const char *ROUTINE, const char *format, ...) {
  va_list args;

  va_start(args,format);
  fprintf(stdout, "Warning: ");
  if (ROUTINE!=0)
    fprintf(stdout,"%s: ",ROUTINE);
  vfprintf(stdout,format,args);
  va_end(args);
  fprintf(stdout,"\n");
  return;
}*/
