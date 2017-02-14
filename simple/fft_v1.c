
//////////////
// fft_v1.c //
//////////////

#include <assert.h>
#include <string.h>
#include <stdio.h>

#include "fft_v1.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif

// This implementation is still very naive, but it actually *is* a true FFT.
// It uses recursion to decompose an even-sized FFT onto two half-size
// FFTs, and combining their result.
//
// This version uses a copy for each generation 

static void fft_v1_recursive(complex double * z, unsigned size, unsigned stride)
{
    if (size <= 1)
    {
        return;
    }

    // Fail on FFT sizes that are not decomposable in factors 2.

    assert(size % 2 == 0);

    // Do sub-fft's

    fft_v1_recursive(z         , size / 2, stride * 2);
    fft_v1_recursive(z + stride, size / 2, stride * 2);

    // Copy the sub-fft results to a local, temporary array.

    complex double zsub[size];

    for (unsigned k = 0; k < size; ++k)
    {
        zsub[k] = z[stride * k];
    }

    // Combine the sub-FFTs naively.

    for (unsigned k = 0; k < size; ++k)
    {
        const unsigned left  = (k * 2    ) % size;
        const unsigned right = (k * 2 + 1) % size;

        const complex double w = cexp(-2 * M_PI * I * k / size);

        z[stride * k] = zsub[left] + w * zsub[right];
    }
}

void fft_v1(complex double * z, unsigned size)
{
    fft_v1_recursive(z, size, 1);
}
