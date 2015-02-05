#ifndef INC_BYTEROUTINES_H
#define INC_BYTEROUTINES_H
/*! /file ByteRoutines.h
    /brief Collection of routines for manipulating byte "endianness".
 */
/// Union of 8 bytes that can be converted to 2 ints or 1 double.
union byte8 {
  unsigned char c[8];
  int i[2];
  double d;
};
/// Perform byte swaps on 4-byte segments.
void endian_swap(void *, long);
/// Perform byte swaps on 8-byte segments.
void endian_swap8(void *, long);
#endif
