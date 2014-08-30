
///////////////////////////////
// ReferenceImplementation.h //
///////////////////////////////

#ifndef FastFourierTransform_ReferenceImplementation_h
#define FastFourierTransform_ReferenceImplementation_h

// SUMMARY
// -------
//
// This code implements an arbitrary-precision FFT.
//
// AUTHOR
// ------
//
// Sidney Cadot <sidney@jigsaw.nl>
//
// REMARKS
// -------
//
// The FFT implemented has the following properties:
//
//  - both "forward" and "inverse" transforms are implemented;
//  - 1D complex input, 1D complex output;
//  - in-place (output overwrites input);
//  - stride may be specified;
//  - number of points n can be any integer (not restricted to powers of 2);
//  - arbitrary precision (using MPC / MPFR / GMP).
//
// This FFT code is optimized for accuracy, not performance:
//
// - For all n, the algorithm is O(n log n), so this is indeed a "Fast" FFT.
//   However, it is still pretty slow since we use the CZT and arbitrary
//   precision numbers.
//
// - The FFT is actually calculated using the Chirp Z-transform (CZT).
//
// - The CZT makes it possible to express any n-point FFT in terms of 3 larger
//   FFTs (2 forward, 1 inverse). The nice thing is that the larger FFTs can
//   be chosen to be of a power-of-2 size, making it possible to use a
//   straightforward Cooley-Tukey implementation.
//
// So, our code consists of 3 functions: a "power-of-2" FFT, the CZT transform,
// and a generic FFT routine that expresses a generic n-point FFT as a CZT
// transform.
//
// The idea to use the CZT transform to support arbitrary numbers of points was
// found on StackExchange; see (1), first answer.
//
// That answer refers to publications describing the CZT: (2) and (3).
//
// In order to implement the CZT, I consulted both the Octave implementation (4) as well
// as the Matlab implementation. They are nearly identical.
//
// Wikipedia has a nice article on the Chirp Z-tranform as well (5).
//
// REFERENCES
// ----------
//
// (1) http://math.stackexchange.com/questions/77118/non-power-of-2-ffts
// (2) http://www.alcatel-lucent.com/bstj/vol48-1969/articles/bstj48-5-1249.pdf
// (3) http://dx.doi.org/10.1109/TAU.1969.1162034
// (4) http://sourceforge.net/p/octave/signal/ci/default/tree/inst/czt.m
// (5) http://en.wikipedia.org/wiki/Bluestein's_FFT_algorithm

#include <mpc.h>

enum class FourierTransformDirection
{
    Forward,
    Inverse
};

void pow2_fft(const FourierTransformDirection direction, mpc_t * z, const unsigned n, const unsigned stride, const mpfr_prec_t precision);

void czt(const mpc_t * x, const unsigned n, const unsigned x_stride, mpc_t * y, const unsigned m, const mpc_t & w, const mpc_t & a, const mpfr_prec_t precision);

void generic_fft(const FourierTransformDirection direction, mpc_t * z, const unsigned n, const unsigned stride, const mpfr_prec_t precision);

#endif // FastFourierTransform_ReferenceImplementation_h
