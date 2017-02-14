
//////////////
// fft_v1.h //
//////////////

#ifndef fft_v1_h
#define fft_v1_h

#include <complex.h>

// A recursive FFT without any optimization.

void fft_v1(complex double * z, unsigned size);

#endif // fft_v1_h
