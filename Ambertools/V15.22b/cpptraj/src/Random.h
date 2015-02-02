#ifndef INC_RANDOM_H
#define INC_RANDOM_H
// Class: Random_Number
/// Marsaglias random number generator as implemented in Amber 3.0 Rev A
/** \author George Seibel
  * 
  * Adapted from fortran code in $AMBERHOME/src/pmemd/src/random.F90
  *
  * This random number generator originally appeared in "Toward a Universal
  * Random Number Generator" by George Marsaglia and Arif Zaman.  Florida
  * State University Report: FSU-SCRI-87-50 (1987)
  *
  * It was later modified by F. James and published in "A Review of Pseudo-
  * random Number Generators"
  *
  * This is claimed to be the best known random number generator available.
  * It passes ALL of the tests for random number generators and has a
  * period of 2^144, is completely portable (gives bit identical results on
  * all machines with at least 24-bit mantissas in the floating point
  * representation).
  *
  * The algorithm is a combination of a Fibonacci sequence (with lags of 97
  * and 33, and operation "subtraction plus one, modulo one") and an
  * "arithmetic sequence" (using subtraction).
  */
class Random_Number {
  public:
    Random_Number();
    /// Initialize the random number generator with the given seed
    void rn_set(int);
    /// Initialize RN generator with 71277 (Amber default)
    void rn_set();
    /// Generate a random number between 0.0 and 1.0
    double rn_gen();
    /// Generate a pseudo-random Gaussian sequence.
    double rn_gauss(double,double);
    /// \return true if RN generator has been set up.
    bool IsSet() const { return RN_generator.iseed != -1; }
  private:
    /// Variables necessary for Marsaglia random number stream.
    /** This is placed in a struct in case the state of the random number
      * generator ever needs to be stored.
      */
    struct random_state {
      // Real variables in Marsaglia algorithm
      double u[97];
      double c;
      double cd;
      double cm;
      // Pointers into u() in Marsaglia algorithm
      int i97;
      int j97;
      // Random seed; if -1 random generator has not been set
      int iseed;
    };
    /// Hold the state of the random number generator
    random_state RN_generator;
};
#endif
