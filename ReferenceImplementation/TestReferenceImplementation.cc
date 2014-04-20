
////////////////////////////////////
// TestReferenceImplementation.cc //
////////////////////////////////////

#include <iostream>
#include <mpc.h>

#include "ReferenceImplementation.h"

using namespace std;

int main()
{
    mpfr_prec_t precision = 1024;

    unsigned NUM_POINTS = 16;

    bool print = true;

    size_t print_precision = 0;

    mpc_t * z = new mpc_t[NUM_POINTS];

    for (unsigned i = 0; i < NUM_POINTS; ++i)
    {
        mpc_init2(z[i], precision);
        mpc_set_ui_ui(z[i], 10 + i, 20 + i * i, DEFAULT_MPC_ROUNDINGMODE);
    }

    pow2_recursive_fft<FourierTransformDirection::Forward>(z, NUM_POINTS, 1, precision);
    pow2_recursive_fft<FourierTransformDirection::Inverse>(z, NUM_POINTS, 1, precision);

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

    mpfr_free_cache();

    return 0;
}
