
//////////////
// fft_v0.c //
//////////////

#include <string.h>

#include "fft_v0.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif

void fft_v0(complex double * z, unsigned size)
{
    complex double z_orig[size];

    memcpy(z_orig, z, sizeof(z_orig));

    for (unsigned i = 0; i < size; ++i)
    {
        z[i] = 0;

        for (unsigned j = 0; j < size; ++j)
        {
            z[i] += cexp(-2.0 * M_PI * I * i * j / size) * z_orig[j];
        }
    }
}
