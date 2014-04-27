
///////////////////////////////
// naive_fourier_transform.h //
///////////////////////////////

#ifndef naive_fourier_transform_h
#define naive_fourier_transform_h

#include <cassert>
#include <vector>
#include <complex>

template <typename real_type>
void naive_fourier_transform(const std::vector<std::complex<real_type>> & tvec, std::vector<std::complex<real_type>> & fvec)
{
    const unsigned n = tvec.size();

    // Make sure the destination vector has the correct size.
    assert(fvec.size() == n);

    for (unsigned fi = 0; fi < n; ++fi)
    {
        fvec[fi] = 0;

        for (unsigned ti = 0; ti < n; ++ti)
        {
            fvec[fi] += std::polar(1.0, -2.0 * M_PI * fi * ti / n) * tvec[ti];
        }
    }
}

template <typename real_type>
void naive_inverse_fourier_transform(const std::vector<std::complex<real_type>> & fvec, std::vector<std::complex<real_type>> & tvec)
{
    const unsigned n = fvec.size();

    // Make sure the destination vector has the correct size.
    assert(tvec.size() == n);

    for (unsigned ti = 0; ti < n; ++ti)
    {
        tvec[ti] = 0;

        for (unsigned fi = 0; fi < n; ++fi)
        {
            tvec[ti] += std::polar(1.0 / n, +2.0 * M_PI * fi * ti / n) * fvec[fi];
        }
    }
}

#endif // naive_fourier_transform_h
