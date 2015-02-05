/*
 * CPPTRAJ: Rewrite of PTRAJ in C++
 * Copyright (c) 2010-2014 Daniel R. Roe
 * For license information see the LICENSE file.
 * For a full list of contributing authors see the README file.
 */
#include "Cpptraj.h"
#include "MpiRoutines.h"
// ----------========== CPPTRAJ MAIN ROUTINE ==========----------
/// Main routine.
int main(int argc, char **argv) {
  Cpptraj Program;
  if (parallel_init(argc,argv) != 0) return 1;
  int err = Program.RunCpptraj(argc, argv);
  parallel_end();
  return err;
}
