
#include <cmath>
#include <iostream>
#include <vector>
#include <complex>
#include <cassert>
#include <chrono>
#include <fftw3.h>

double test_fftw(const unsigned n, const unsigned repeats)
{
    vector<complex<double>> v = mk_complex_test_vector(n);

    fftw_complex * x = fftw_alloc_complex(n);
    assert(x != nullptr);

    fftw_complex * y = fftw_alloc_complex(n);
    assert(y != nullptr);

    fftw_plan plan = fftw_plan_dft_1d(n, x, y, FFTW_FORWARD, FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
    assert(plan != nullptr);

    // Initialize the input

    // copy test data
    for (unsigned i = 0; i < n; ++i)
    {
        x[i][0] = real(v[i]);
        x[i][1] = imag(v[i]);
    }

    auto t1 = chrono::high_resolution_clock::now();

    for (unsigned rep = 0; rep < repeats; ++rep)
    {
        fftw_execute(plan);
    }

    auto t2 = chrono::high_resolution_clock::now();

    double duration = chrono::duration_cast<chrono::nanoseconds>(t2 - t1).count() / 1e9;

    fftw_destroy_plan(plan);

    fftw_free(y);
    fftw_free(x);

    fftw_cleanup();

    return duration / repeats;
}
