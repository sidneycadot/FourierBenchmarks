
///////////
// czt.c //
///////////

#include "czt.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846264338327
#endif

void czt(complex double * z, unsigned n, complex double * ztrans, unsigned m, complex double w, complex double a, void (*fftfunc)(complex double *, unsigned))
{
    // Determine next-biggest power-of-two that fits the (n + m - 1) entries we need.

    unsigned fft_size = 1;
    while (fft_size < n + m - 1)
    {
        fft_size *= 2;
    }

    complex double zz[fft_size];

    // Initialize zz.

    for (unsigned k = 0; k < fft_size; ++k)
    {
        if (k < n)
        {
            const complex double w1 = cpow(w, 0.5 * k * k) / cpow(a, k);
            zz[k] = w1 * z[k];
        }
        else
        {
            zz[k] = 0;
        }
    }

    // Do forward FFT of zz.

    fftfunc(zz, fft_size);

    // Allocate and initialize w2 that we will convolve with.

    complex double w2[fft_size];

    for (unsigned k = 0; k < fft_size; ++k)
    {
        if (k < n + m - 1)
        {
            const int kshift = k - (n - 1);

            w2[k] = cpow(w, -0.5 * kshift * kshift);
        }
        else
        {
            w2[k] = 0;
        }
    }

    // Do forward FFT of w2.

    fftfunc(w2, fft_size);

    // Do convolution: zz[i] = zz[i] * w2[i]

    for (unsigned k = 0; k < fft_size; ++k)
    {
        zz[k] *= w2[k];
    }

    // Do inverse FFT of (zz * w2), put result in zz.

    fftfunc(zz, fft_size); // forward FFT

    // Make an inverse FFT from the forward FFT.
    // - scale all elements by 1 / fft_size;
    // - reverse elements 1 .. (fft_size - 1).

    for (unsigned k = 0; k < fft_size; ++k)
    {
        zz[k] /= fft_size;
    }

    for (unsigned k = 1; k < fft_size - k; ++k)
    {
        const unsigned kswap = fft_size - k;

        const complex double temp = zz[k];
        zz[k]     = zz[kswap];
        zz[kswap] = temp;
    }

    // Extract output:

    for (unsigned k = 0; k < m; ++k)
    {
        const complex double w3 = cpow(w, (0.5 * k * k));
        ztrans[k] = w3 * zz[n - 1 + k];
    }
}

void czt_fft(complex double * z, unsigned size, void (*fftfunc)(complex double *, unsigned))
{
    if (size == 0)
    {
        return;
    }

    // Check if size is a power of two.

    unsigned sz = size;
    while (sz % 2 == 0)
    {
        sz /= 2;
    }

    if (sz == 1)
    {
        // Size is a power of two. Defer to the given radix-2 fftfunc.
        fftfunc(z, size);
    }
    else
    {
        // Size is not a power of two.
        // Calculate the FFT via the Chirp-Z transform.

        const complex double w = cexp(-2.0 * M_PI * I / size);
        const complex double a = 1;

        czt(z, size, z, size, w, a, fftfunc);
    }
}
