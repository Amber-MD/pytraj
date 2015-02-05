/** \param x starting memory address of consecutive 4-byte segments to swap
  * \param nswap number of segments to swap
  */
void endian_swap(void *x, long nswap) {
  int *Xi;

  int *addr = (int *) x;
  for (long i=0; i<nswap; i++) {
    Xi = addr + i;
    *Xi=(
         ((*Xi>>24) & 0xFF)   |
         ((*Xi&0xFF)   << 24) |
         ((*Xi>>8)  & 0xFF00) |
         ((*Xi&0xFF00) << 8)
        );
  }
}

/** Break up the total swap into two 4-byte swaps.
  * \param x starting memory address of consecutive 8-byte segments to swap
  * \param nswap number of segments to swap
  */
void endian_swap8(void *x, long nswap) {
  int *Xi;
  int x0, x1;

  int *addr = (int *) x;
  for (long i=0; i<nswap; i++) {
    Xi = addr + (i<<1);
    // Perform swap on first 4 bytes
    x0 = Xi[0];
    x0=(
        ((x0>>24) & 0xFF)   |
        ((x0&0xFF)   << 24) |
        ((x0>>8)  & 0xFF00) |
        ((x0&0xFF00) << 8)
       );
    // Perform swap on second 4 bytes
    x1 = Xi[1];
    x1=(
        ((x1>>24) & 0xFF)   |
        ((x1&0xFF)   << 24) |
        ((x1>>8)  & 0xFF00) |
        ((x1&0xFF00) << 8)
       );
    // Re-assemble swapped bytes
    Xi[0] = x1;
    Xi[1] = x0;
  }
}

