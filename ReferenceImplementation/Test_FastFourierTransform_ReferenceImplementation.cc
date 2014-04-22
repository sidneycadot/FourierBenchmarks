
//////////////////////////////////////////////////////////
// Test_FastFourierTransform_ReferenceImplementation.cc //
//////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>
#include <chrono>
#include <cassert>
#include <cstdlib>
#include <mpc.h>

#include "FastFourierTransform_ReferenceImplementation.h"

static bool is_power_of_two(unsigned n)
{
    if (n == 0)
    {
        return false;
    }

    while (n % 2 == 0)
    {
        n /= 2;
    }

    return (n == 1);
}

// We assume a 64-bit unsigned long type, for the gmp_randseed_ui() call.
static_assert(sizeof(unsigned long) == 8, "unsigned long is too small.");

int main()
{
    using namespace std;

    bool print = false;

    size_t print_precision = 40;

    const unsigned NUM_REPEATS = 10;

    const unsigned TIME_LIMIT_REPEATS_US = 20000000; // 20 seconds max. for all repeats.

    // Prepare random number generator of GMP.

    gmp_randstate_t rnd_state;
    gmp_randinit_default(rnd_state);

    // Outer trial loop.

    for (mpfr_prec_t precision = 16; precision <= 1048576; precision = is_power_of_two(precision) ? (precision * 3 / 2) : (precision * 4 / 3))
    {
        mpfr_t rnd_re;
        mpfr_t rnd_im;

        mpc_t  diff;
        mpfr_t err;
        mpfr_t max_err;
        mpfr_t rms_err;

        mpfr_init2(rnd_re, precision);
        mpfr_init2(rnd_im, precision);

        mpc_init2 (diff    , precision + 1024); // give extra precision to the numbers used to calculate the error.
        mpfr_init2(err     , precision + 1024);
        mpfr_init2(max_err , precision + 1024);
        mpfr_init2(rms_err , precision + 1024);

        for (unsigned num_points = 1;; ++num_points)
        {
            mpc_t * x = new mpc_t[num_points];
            mpc_t * y = new mpc_t[num_points];

            for (unsigned i = 0; i < num_points; ++i)
            {
                mpc_init2(x[i], precision);
                mpc_init2(y[i], precision);
            }

            unsigned total_duration = 0;

            // Initialize the seed depending on precision and num_points, to ensure reproducible runs.

            gmp_randseed_ui(rnd_state, precision * 0x100000000UL + num_points);

            for (unsigned repeat = 1; repeat <= NUM_REPEATS; ++repeat)
            {
                // Prepare a single trial.

                // (1) Fill x[] array with gaussian noise.

                for (unsigned i = 0; i < num_points; ++i)
                {
                    mpfr_grandom(rnd_re, rnd_im, rnd_state, DEFAULT_MPFR_ROUNDINGMODE);

                    mpc_set_fr_fr(x[i], rnd_re, rnd_im, DEFAULT_MPC_ROUNDINGMODE);
                }

                // (2) Copy x to y.
                //     We do this because our FFT is in-place, and we need the original data to determine our error.

                for (unsigned i = 0; i < num_points; ++i)
                {
                    mpc_set(y[i], x[i], DEFAULT_MPC_ROUNDINGMODE);
                }

                // (3) Execute forward FFT of y[].

                unsigned duration_forward;

                {
                    std::chrono::time_point<std::chrono::high_resolution_clock> t1 = std::chrono::high_resolution_clock::now();

                    generic_fft(FourierTransformDirection::Forward, y, num_points, 1, precision);

                    std::chrono::time_point<std::chrono::high_resolution_clock> t2 = std::chrono::high_resolution_clock::now();

                    duration_forward = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
                }

                // (4) Optional print of y[] after FFT.

                if (print)
                {
                    for (unsigned i = 0; i < num_points; ++i)
                    {
                        char * y_str = mpc_get_str(10, print_precision, y[i], DEFAULT_MPC_ROUNDINGMODE);

                        cout << "y[" << i << "] = " << y_str << endl;

                        mpc_free_str(y_str);
                    }
                }

                // (5) Execute inverse FFT of y[].

                unsigned duration_inverse;

                {
                    std::chrono::time_point<std::chrono::high_resolution_clock> t1 = std::chrono::high_resolution_clock::now();

                    generic_fft(FourierTransformDirection::Inverse, y, num_points, 1, precision);

                    std::chrono::time_point<std::chrono::high_resolution_clock> t2 = std::chrono::high_resolution_clock::now();

                    duration_inverse = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
                }

                // (6) Calculate max_err = max(abs(y[i] - x[i]).
                //           and rms_err = sqrt(mean(norm(y[i] - x[i])))

                mpfr_set_ui(max_err, 0, DEFAULT_MPFR_ROUNDINGMODE);
                mpfr_set_ui(rms_err, 0, DEFAULT_MPFR_ROUNDINGMODE);

                for (unsigned i = 0; i < num_points; ++i)
                {
                    mpc_sub(diff, y[i], x[i], DEFAULT_MPC_ROUNDINGMODE);
                    mpc_abs(err, diff, DEFAULT_MPFR_ROUNDINGMODE);
                    mpfr_max(max_err, max_err, err, DEFAULT_MPFR_ROUNDINGMODE);

                    mpc_norm(err, diff, DEFAULT_MPFR_ROUNDINGMODE);
                    mpfr_add(rms_err, rms_err, err, DEFAULT_MPFR_ROUNDINGMODE);
                }

                mpfr_div_ui(rms_err, rms_err, num_points, DEFAULT_MPFR_ROUNDINGMODE);
                mpfr_sqrt(rms_err, rms_err, DEFAULT_MPFR_ROUNDINGMODE);

                // (7) Present result for this trial.

                {
                    int mpfr_asprintf_result;

                    char * max_err_str;
                    mpfr_asprintf_result = mpfr_asprintf(&max_err_str, "%.6Re", max_err);
                    assert(mpfr_asprintf_result > 0);

                    char * rms_err_str;
                    mpfr_asprintf_result = mpfr_asprintf(&rms_err_str, "%.6Re", rms_err);
                    assert(mpfr_asprintf_result > 0);

                    cout << "precision"  << setw( 8) << precision                << "    "
                         << "num_points" << setw( 8) << num_points               << "    "
                         << "repeat"     << setw( 6) << repeat                   << "    "
                         << "forward"    << setw(12) << duration_forward         << "    "
                         << "inverse"    << setw(12) << duration_inverse         << "    "
                         << "rms_error"  << setw(20) << rms_err_str              << "    "
                         << "max_error"  << setw(20) << max_err_str              << endl;

                    mpfr_free_str(max_err_str);
                    mpfr_free_str(rms_err_str);
                }

                // Trial done. Note the time taken.

                total_duration += (duration_forward + duration_inverse);

            } // repeat loop

            for (unsigned i = 0; i < num_points; ++i)
            {
                mpc_clear(x[i]);
                mpc_clear(y[i]);
            }

            delete [] y;
            delete [] x;

            // If the total time for the repeats exceeds 'TIME_LIMIT_REPEATS_US',
            // we are done with this precision; we will start with the next one.

            if (total_duration > TIME_LIMIT_REPEATS_US)
            {
                break;
            }

        } // num_points loop

        mpfr_clear(rnd_re);
        mpfr_clear(rnd_im);
        mpfr_clear(rms_err);
        mpfr_clear(max_err);
        mpfr_clear(err);
        mpc_clear(diff);

    } // precision loop

    // Free cache of MPFR, used for PI.
    mpfr_free_cache();

    // Free the GMP random state.
    gmp_randclear(rnd_state);

    return EXIT_SUCCESS;
}
