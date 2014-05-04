
////////////////
// TestMKL.cc //
////////////////

#include "MklDftiUtils.h"

#include <iostream>
#include <vector>
#include <complex>
#include <cassert>
#include <chrono>

#include <mkl_dfti.h>
#include <mkl_service.h>

using namespace std;

template <typename real_type>
struct MklTraits
{
};

template <>
struct MklTraits<float>
{
    static const enum DFTI_CONFIG_VALUE precision = DFTI_SINGLE;
};

template <>
struct MklTraits<double>
{
    static const enum DFTI_CONFIG_VALUE precision = DFTI_DOUBLE;
};

template <typename real_type>
static double test_complex_1d(const unsigned & n, const unsigned & repeats, const bool & noisy)
{
    vector<complex<real_type>> tvec(n);

    DFTI_DESCRIPTOR_HANDLE descriptor;

    MKL_LONG status;

    status = DftiCreateDescriptor(&descriptor, MklTraits<real_type>::precision, DFTI_COMPLEX, 1, n);
    assert(status == DFTI_NO_ERROR);

    status = DftiSetValue(descriptor, DFTI_PLACEMENT, DFTI_NOT_INPLACE);
    assert(status == DFTI_NO_ERROR);

    status = DftiCommitDescriptor(descriptor);
    assert(status == DFTI_NO_ERROR);

    if (noisy)
    {
        print_mkl_dfti_descriptor_info(cout, descriptor);
    }

    vector<complex<real_type>> fvec(n);

    double min_duration = numeric_limits<double>::infinity();

    for (unsigned rep = 0; rep < repeats; ++rep)
    {
        std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();

        status = DftiComputeForward(descriptor, tvec.data(), fvec.data());
        assert(status == DFTI_NO_ERROR);

        std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();

        const double duration = std::chrono::duration_cast<std::chrono::nanoseconds>(t2 - t1).count() / 1e9;

        min_duration = std::min(min_duration, duration);
    }

    status = DftiFreeDescriptor(&descriptor);
    assert(status == DFTI_NO_ERROR);

    return min_duration;
}

int main()
{
    bool noisy = true;

    unsigned n = 0;

    for (n = 16; n <= 16; ++n)
    {
        double durationFloat = test_complex_1d<float>(n, 100, noisy);

        double durationDouble = test_complex_1d<double>(n, 100, noisy);

        cout << n << " " << durationFloat << " " << durationDouble << endl;
    }

    mkl_free_buffers();

    return EXIT_SUCCESS;
}
