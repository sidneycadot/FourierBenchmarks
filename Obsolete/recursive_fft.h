
/////////////////////
// recursive_fft.h //
/////////////////////

#ifndef recursive_fft_h
#define recursive_fft_h

#include <complex>

template <typename real_type>
void recursive_fft(std::complex<real_type> * z, const unsigned n, const unsigned stride)
{
    if (n == 1)
    {
        return;
    }

    assert (n % 2 == 0);

    recursive_fft(z         , n / 2, 2 * stride);
    recursive_fft(z + stride, n / 2, 2 * stride);

    // Now, we need to combine the values of the sub-FFTs

    std::complex<real_type> temp[n];

    for (unsigned i = 0; i < n / 2; ++i)
    {
        const real_type turn = -2.0 * M_PI * i / n;

        const std::complex<real_type> coeff = std::polar(1.0, turn);

        temp[2 * i + 0] = z[stride * 2 * i] + coeff * z[stride * (2 * i + 1)];
        temp[2 * i + 1] = z[stride * 2 * i] - coeff * z[stride * (2 * i + 1)];
    }

    // And re-arrange them ...

    for (unsigned i = 0; i < n; ++i)
    {
        const unsigned j = (i / 2) + (i % 2) * (n / 2);

        assert(j < n);

        z[j * stride] = temp[i];
    }
}

template <typename real_type>
void recursive_inverse_fft(std::complex<real_type> * z, const unsigned n, const unsigned stride)
{
    if (n == 1)
    {
        return;
    }

    assert (n % 2 == 0);

    recursive_inverse_fft(z         , n / 2, 2 * stride);
    recursive_inverse_fft(z + stride, n / 2, 2 * stride);

    // Now, we need to combine the values of the sub-FFTs

    std::complex<real_type> temp[n];

    for (unsigned i = 0; i < n / 2; ++i)
    {
        const real_type turn = +2.0 * M_PI * i / n;

        const std::complex<real_type> coeff = std::polar(1.0, turn);

        temp[2 * i + 0] = 0.5 * (z[stride * 2 * i] + coeff * z[stride * (2 * i + 1)]);
        temp[2 * i + 1] = 0.5 * (z[stride * 2 * i] - coeff * z[stride * (2 * i + 1)]);
    }

    // And re-arrange them ...

    for (unsigned i = 0; i < n; ++i)
    {
        const unsigned j = (i / 2) + (i % 2) * (n / 2);

        assert(j < n);

        z[j * stride] = temp[i];
    }
}

#endif // recursive_fft_h
