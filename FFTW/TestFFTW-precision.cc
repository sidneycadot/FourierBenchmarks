
///////////////////////////
// TestFFTW-precision.cc //
///////////////////////////

#include <vector>
#include <cassert>

#include <iostream>
#include <iomanip>

#include "FftwUtils.h"
#include "SignalGenerator/SignalGenerator.h"
#include "ReferenceImplementation/ReferenceImplementation.h"

using namespace std;

void calc_errors(const mpc_t * approx, const mpc_t * perfect, const unsigned n, mpfr_t & rms_err, mpfr_t & max_err, const mpfr_prec_t precision)
{
    // Determine max error and rms error

    mpc_t diff;
    mpfr_t err;

    mpc_init2(diff, precision);
    mpfr_init2(err, precision);

    mpfr_set_zero(rms_err, +1);
    mpfr_set_zero(max_err, +1);

    for (unsigned i = 0; i < n; ++i)
    {
        mpc_sub(diff, approx[i], perfect[i], DEFAULT_MPC_ROUNDINGMODE);

        mpc_abs(err, diff, DEFAULT_MPFR_ROUNDINGMODE);
        mpfr_max(max_err, max_err, err, DEFAULT_MPFR_ROUNDINGMODE);

        mpc_norm(err, diff, DEFAULT_MPFR_ROUNDINGMODE);
        mpfr_add(rms_err, rms_err, err, DEFAULT_MPFR_ROUNDINGMODE);
    }

    mpfr_div_ui(rms_err, rms_err, n, DEFAULT_MPFR_ROUNDINGMODE);
    mpfr_sqrt(rms_err, rms_err, DEFAULT_MPFR_ROUNDINGMODE);

    mpfr_clear(err);
    mpc_clear(diff);
}

template <typename fftw_traits>
void execute_tests_r2c_1d(const unsigned & n_in, const unsigned & repeats)
{
    const mpfr_prec_t NOISE_PRECISION = 256; // We fix this to ensure reproducibility.

    const mpfr_prec_t precision = 256;

    const unsigned n_out = n_in / 2 + 1;

    typename fftw_traits::real_type * x = fftw_traits::alloc_real(n_in);
    assert(x != nullptr);

    typename fftw_traits::complex_type * y = fftw_traits::alloc_complex(n_out);
    assert(y != nullptr);

    // Plan FFT

    const unsigned flags = FFTW_PRESERVE_INPUT | FFTW_ESTIMATE;

    typename fftw_traits::plan plan = fftw_traits::plan_dft_r2c_1d(n_in, x, y, flags);

    assert(plan != nullptr);

    for (unsigned rep = 1; rep <= repeats; ++rep)
    {
        GaussianNoiseSignal noise(to_string(rep), NOISE_PRECISION);
        //ZeroSignal noise;

        // Initialize x with signal

        const vector<unsigned> dims = { n_in };
        const vector<int> strides = { 1 };

        sample(noise, x, strides, dims, precision);

        // Perform FFTW fft

        fftw_traits::execute(plan);

        // Copy result array to zy

        mpc_t * zy = new mpc_t[n_out];

        for (unsigned i = 0; i < n_out; ++i)
        {
            mpc_init2(zy[i], precision);
        }

        for (unsigned i = 0; i < n_out; ++i)
        {
            mpc_set_fp_fp(zy[i], y[i][0], y[i][1], DEFAULT_MPC_ROUNDINGMODE);
        }

        // Do reference implementation.

        // Prepare z (complex array used for both input and output).

        mpc_t * z = new mpc_t[n_in];

        for (unsigned i = 0; i < n_in; ++i)
        {
            mpc_init2(z[i], precision);
        }

        for (unsigned i = 0; i < n_in; ++i)
        {
            mpc_set_fp(z[i], x[i], DEFAULT_MPC_ROUNDINGMODE);
        }

        generic_fft(FourierTransformDirection::Forward, z, n_in, 1, precision);

        // Determine max error and rms error
        {
            mpfr_t max_err1;
            mpfr_t rms_err1;
            mpfr_t max_err2;
            mpfr_t rms_err2;

            mpfr_init2(max_err1, precision);
            mpfr_init2(rms_err1, precision);
            mpfr_init2(max_err2, precision);
            mpfr_init2(rms_err2, precision);

            calc_errors(zy, z, n_out, rms_err1, max_err1, precision);

            // Reduce 'z' to the precision of the underlying variable, giving the
            // 'best possible' value.

            for (unsigned i = 0; i < n_out; ++i)
            {
                std::complex<typename fftw_traits::real_type> rounded = mpc_get_complex_fp<typename fftw_traits::real_type>(z[i], DEFAULT_MPFR_ROUNDINGMODE);

                mpc_set_fp_fp(z[i], std::real(rounded), std::imag(rounded), DEFAULT_MPC_ROUNDINGMODE);
            }

            calc_errors(zy, z, n_out, rms_err2, max_err2, precision);

            // rms_err1/max_err1 show the error relative to the mathematical thruth;
            // rms_err2/max_err2 show the error relative to the best possible approximation, given the precision of the type where the FFT result is stored.

            cout << "precision"  " " << setw( 8) << precision           << " "
                    "rms_error1" " " << setw(20) << to_string(rms_err1) << " "
                    "max_error1" " " << setw(20) << to_string(max_err1) << " "
                    "rms_error2" " " << setw(20) << to_string(rms_err2) << " "
                    "max_error2" " " << setw(20) << to_string(max_err2) << endl;

            mpfr_clear(rms_err2);
            mpfr_clear(max_err2);
            mpfr_clear(rms_err1);
            mpfr_clear(max_err1);
        }

        for (unsigned i = 0; i < n_in; ++i)
        {
            mpc_clear(z[i]);
        }

        for (unsigned i = 0; i < n_out; ++i)
        {
            mpc_clear(zy[i]);
        }

        delete [] z;
        delete [] zy;
    }

    fftw_traits::destroy_plan(plan);
    fftw_traits::free(y);
    fftw_traits::free(x);

    fftw_traits::cleanup();
}

int main()
{
    const unsigned n = 1024;
    const unsigned repeats = 2;

    execute_tests_r2c_1d<fftw_traits<float>>(n, repeats);

    // Free MPFR cache.
    mpfr_free_cache();

    return EXIT_SUCCESS;
}
