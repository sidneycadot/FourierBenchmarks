
//////////////
// fft_v0.h //
//////////////

#ifndef fft_v0_h
#define fft_v0_h

#include <complex.h>

// This is not, in fact, an FFT, but rather a DFT with O(size **2) time complexity.

void fft_v0(complex double * z, unsigned size);

#endif // fft_v0_h
