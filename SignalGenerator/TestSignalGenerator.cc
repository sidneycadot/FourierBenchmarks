
////////////////////////////
// TestSignalGenerator.cc //
////////////////////////////

#include "SignalGenerator.h"

#include <iostream>
#include <vector>
#include <complex>

#include <mpc.h>

using namespace std;

int main()
{
    const unsigned SIZE = 1024;

    unsigned n = SIZE * SIZE;

    const mpfr_prec_t precision = 1024;

    vector<complex<float>> arr(n);

    for (unsigned i = 0; i < n; ++i)
    {
        arr[i] = -999.0;
    }

    std::vector<unsigned> dims    = {SIZE, SIZE};
    std::vector<signed  > strides = {SIZE, 1};

    GaussianNoiseSignalFast noise("", precision);
    //ZeroSignal zero;

    sample(noise, arr.data(), strides, dims, precision);

    if (0)
    for (unsigned i = 0; i < n; ++i)
    {
        cout << i << " " << arr[i] << endl;
    }

    // Free cache of MPFR, used for PI.
    mpfr_free_cache();

    return 0;
}
