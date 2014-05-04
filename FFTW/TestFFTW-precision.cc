
///////////////////////////
// TestFFTW-precision.cc //
///////////////////////////

#include <vector>
#include <cassert>

#include "FftwUtils.h"
#include "SignalGenerator.h"
#include "ReferenceImplementation.h"

using namespace std;

template <typename traits>
void execute_tests_r2c_1d(const unsigned & n_in, const unsigned & repeats)
{
    const mpfr_prec_t precision = 1024;

    GaussianNoiseSignal noise("", precision);

    const unsigned n_out = n_in / 2 + 1;

    typename traits::real_type * x = traits::alloc_real(n_in);
    assert(x != nullptr);

    typename traits::complex_type * y = traits::alloc_complex(n_out);
    assert(y != nullptr);

    // Plan FFT

    const unsigned flags = FFTW_PRESERVE_INPUT | FFTW_ESTIMATE;

    typename traits::plan plan = traits::plan_dft_r2c_1d(n_in, x, y, flags);

    assert(plan != nullptr);

    if (true)
    {
        // Initialize x with signal

        const vector<unsigned> dims = { n_in };
        const vector<int> strides = { 1 };

        sample(noise, x, strides, dims, precision);

        // Perform FFTW fft

        traits::execute(plan);

        // Do reference implementation.

        // Prepare Z
        mpc_t * z = new mpc_t[n_in];

        for (unsigned i = 0; i < n_in; ++i)
        {
            mpc_init2(z[i], precision);
        }

        for (unsigned i = 0; i < n_in; ++i)
        {
            mpc_set_d(z[i], x[i], DEFAULT_MPC_ROUNDINGMODE);
        }

        generic_fft(FourierTransformDirection::Forward, z, n_in, 1, precision);

#if 0
        // Determine max error and rms error
        {
            mpfr_set_zero(max_err, +1);
            mpfr_set_zero(rms_err, +1);

            for (unsigned i = 0; i < n_in; ++i)
            {
                mpc_sub(diff, zy[i], z[i], DEFAULT_MPC_ROUNDINGMODE);

                mpc_abs(err, diff, DEFAULT_MPFR_ROUNDINGMODE);
                mpfr_max(max_err, max_err, err, DEFAULT_MPFR_ROUNDINGMODE);

                mpc_norm(err, diff, DEFAULT_MPFR_ROUNDINGMODE);
                mpfr_add(rms_err, rms_err, err, DEFAULT_MPFR_ROUNDINGMODE);
            }

            mpfr_div_ui(rms_err, rms_err, n_in, DEFAULT_MPFR_ROUNDINGMODE);
            mpfr_sqrt(rms_err, rms_err, DEFAULT_MPFR_ROUNDINGMODE);

            // (7) Present result for this trial.

            cout << "precision"  << setw( 8) << precision          << "    "
                    << "num_points" << setw( 8) << num_points         << "    "
                    << "repeat"     << setw( 6) << repeat             << "    "
                    << "forward"    << setw(12) << duration_forward   << "    "
                    << "inverse"    << setw(12) << duration_inverse   << "    "
                    << "rms_error"  << setw(20) << to_string(rms_err) << "    "
                    << "max_error"  << setw(20) << to_string(max_err) << endl;
        }
#endif

        for (unsigned i = 0; i < n_in; ++i)
        {
            mpc_clear(z[i]);
        }

        delete [] z;
    }

    traits::destroy_plan(plan);
    traits::free(y);
    traits::free(x);

    traits::cleanup();
}

int main()
{
    const unsigned n = 16;
    const unsigned repeats = 100;

    execute_tests_r2c_1d<FFTW_Traits<float>>(n, repeats);

    // Free MPFR cache.
    mpfr_free_cache();

    return EXIT_SUCCESS;
}
