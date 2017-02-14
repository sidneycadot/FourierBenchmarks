
///////////
// czt.h //
///////////

#ifndef czt_h
#define czt_h

#include <complex.h>

void czt(complex double * z, unsigned n, complex double * ztrans, unsigned m, complex double w, complex double a, void (*fftfunc)(complex double *, unsigned));

void czt_fft(complex double * z, unsigned size, void (*fftfunc)(complex double *, unsigned));

#endif // czt_h
