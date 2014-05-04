
////////////////////////////
// TestSignalGenerator.cc //
////////////////////////////

#include <cstdlib>
#include <vector>
#include <complex>

#include <mpc.h>

#include "SignalGenerator.h"

int main()
{
    using namespace std;

    const unsigned SIZE = 32;

    const unsigned n = SIZE * SIZE;

    const mpfr_prec_t precision = 1024;

    vector<complex<double>> arr(n);

    vector<unsigned> dims    = {SIZE, SIZE};
    vector<signed  > strides = {SIZE, 1};

    GaussianNoiseSignalFast noise("", precision);

    sample(noise, arr.data(), strides, dims, precision);

    // Free MPFR cache.
    mpfr_free_cache();

    return EXIT_SUCCESS;
}
