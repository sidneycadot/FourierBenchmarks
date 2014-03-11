
#include <cstdio>
#include <iostream>

#include <vector>
#include <complex>

#include <mpc.h>
#include <gmp.h>

using namespace std;

template <typename real_type>
inline real_type to_native(const mpfr_t & x, const mpfr_rnd_t & rnd);

template <>
inline float to_native<float>(const mpfr_t & x, const mpfr_rnd_t & rnd)
{
    return mpfr_get_flt(x, rnd);
}

template <>
inline double to_native<double>(const mpfr_t & x, const mpfr_rnd_t & rnd)
{
    return mpfr_get_d(x, rnd);
}

template <>
inline long double to_native<long double>(const mpfr_t & x, const mpfr_rnd_t & rnd)
{
    return mpfr_get_ld(x, rnd);
}

template <typename real_type>
vector<complex<real_type>> mk_complex_gauss_vector(gmp_randstate_t & randstate, const unsigned & n)
{
    vector<complex<real_type>> v;
    v.reserve(n);

    mpfr_t x, y;

    mpfr_init2(x, 1000);
    mpfr_init2(y, 1000);

    while (v.size() < n)
    {
        mpfr_grandom(x, y, randstate, MPFR_RNDN);

        const real_type xnative = to_native<real_type>(x, MPFR_RNDN);
        const real_type ynative = to_native<real_type>(y, MPFR_RNDN);

        v.push_back(complex<real_type>(xnative, ynative));
    }

    mpfr_clear(y);
    mpfr_clear(x);

    return v;
}

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

    recursive_fft(z         , n / 2, 2 * stride);
    recursive_fft(z + stride, n / 2, 2 * stride);

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

template <typename Collection>
ostream & print_vector(ostream & os, const Collection & collection)
{
    bool first = true;

    os << "{";

    for (const auto & element : collection)
    {
        if (first)
        {
            first = false;
        }
        else
        {
            os << ", ";
        }

        os << element;
    }

    os << "}";

    return os;
}

int main()
{
    gmp_randstate_t randstate;

    gmp_randinit_default(randstate);

    vector<complex<float>>  vf = mk_complex_gauss_vector<float>(randstate, 10);

    cout << "vf: ";
    print_vector(cout, vf);
    cout << endl;

    gmp_randclear(randstate);

    gmp_randinit_default(randstate);

    vector<complex<double>> vd = mk_complex_gauss_vector<double>(randstate, 10);

    cout << "vd: ";
    print_vector(cout, vd);
    cout << endl;

    gmp_randclear(randstate);

    gmp_randinit_default(randstate);

    vector<complex<long double>> vld = mk_complex_gauss_vector<long double>(randstate, 10);

    gmp_randclear(randstate);

    cout << "vld: ";
    print_vector(cout, vld);
    cout << endl;

    // Get gaussian random number.

    mpfr_free_cache();

    return 0;
}
