
//////////////////////////////////////////////////////////
// Test_FastFourierTransform_ReferenceImplementation.cc //
//////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>
#include <chrono>
#include <mpc.h>

#include "FastFourierTransform_ReferenceImplementation.h"

using namespace std;

bool is_power_of_two(unsigned n)
{
    while (n != 0 && n % 2 == 0)
    {
        n >>= 1;
    }

    return (n == 1);
}

int main()
{
    bool print = false;

    size_t print_precision = 40;

    for (mpfr_prec_t precision = 16; precision <= 1048576; precision = is_power_of_two(precision) ? (precision * 3 / 2) : (precision * 4 / 3))
    {
        for (unsigned NUM_POINTS = 1;; ++NUM_POINTS)
        {
            mpc_t * z = new mpc_t[NUM_POINTS];

            for (unsigned i = 0; i < NUM_POINTS; ++i)
            {
                mpc_init2(z[i], precision);
                mpc_set_ui_ui(z[i], 10 + i, 20 + i * i, DEFAULT_MPC_ROUNDINGMODE);
            }

            std::chrono::time_point<std::chrono::high_resolution_clock> t1 = std::chrono::high_resolution_clock::now();

            generic_fft(FourierTransformDirection::Forward, z, NUM_POINTS, precision);

            std::chrono::time_point<std::chrono::high_resolution_clock> t2 = std::chrono::high_resolution_clock::now();

            unsigned duration_us = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();

            cout << setw(10) << precision << setw(10) << NUM_POINTS << " " << setw(20) << duration_us << " us" << endl;

            if (print)
            {
                for (unsigned i = 0; i < NUM_POINTS; ++i)
                {
                    char * z_str = mpc_get_str(10, print_precision, z[i], DEFAULT_MPC_ROUNDINGMODE);

                    cout << "z[" << i << "] = " << z_str << endl;

                    mpc_free_str(z_str);
                }
            }

            for (unsigned i = 0; i < NUM_POINTS; ++i)
            {
                mpc_clear(z[i]);
            }

            delete [] z;

            if (duration_us >= 1000000)
            {
                break;
            }
        }
    }

    mpfr_free_cache();

    return 0;
}
