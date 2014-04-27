
/////////////////
// fftw_fft.h //
////////////////

#include <vector>
#include <complex>
#include <cassert>
#include <fftw3.h>

void run_fftw(const std::vector<std::complex<double>> & tvec, std::vector<std::complex<double>> & fvec)
{
    const unsigned n = tvec.size();

    assert(fvec.size() == n);

    fftw_complex * x = fftw_alloc_complex(n);
    assert(x != nullptr);

    fftw_complex * y = fftw_alloc_complex(n);
    assert(y != nullptr);

    // Plan before copy (because the planner may destroy the input)

    fftw_plan plan = fftw_plan_dft_1d(n, x, y, FFTW_FORWARD, FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
    assert(plan != nullptr);

    // copy input data
    for (unsigned i = 0; i < n; ++i)
    {
        x[i][0] = std::real(tvec[i]);
        x[i][1] = std::imag(tvec[i]);
    }

    // Execute FFT
    fftw_execute(plan);

    // copy output data
    for (unsigned i = 0; i < n; ++i)
    {
        fvec[i] = std::complex<double>(y[i][0], y[i][1]);
    }

    fftw_destroy_plan(plan);

    fftw_free(y);
    fftw_free(x);

    fftw_cleanup();
}
