
#include <cmath>
#include <iostream>
#include <vector>
#include <complex>
#include <cassert>
#include <chrono>
#include <fftw3.h>

using namespace std;

// Variations:
//
//     Number of transforms
//     Forward vs inverse
//     Scaling
//     Real vs complex
//     Basic data type (float, double, long double)
//
//     Dimensionality of transform
//     Memory layout
//     In-place or not
//
//     Multithreading
//     Optimization (planning) level (FFTW)
//     Preserve input or not

// Implementations:
//
// CPU:
//
//   MKL
//   IPP
//   FFTW
//   FFTS
//   KissFFT
//
//   Slow
//   Basic
//
// GPU:
//
//   CUDA
//   OpenCL?
//
// Precision of the implementations?
//
// Different HW plaforms? (Intel/AMD, AMD/nVidia)
//

bool is_prime(unsigned n)
{
    if (n < 2)
    {
        return false;
    }

    for (unsigned d = 2; d * d <= n; ++d)
    {
        if (n % d == 0)
        {
            return false;
        }
    }

    return true;
}

vector<double> mk_real_test_vector(unsigned n)
{
    vector<double> v;
    v.reserve(n);

    unsigned p = 0;
    while (v.size() < n)
    {
        while (!is_prime(p))
        {
            ++p;
        }

        v.push_back(sqrt(p));

        ++p;
    }

    return v;
}

vector<complex<double>> mk_complex_test_vector(unsigned n)
{
    vector<complex<double>> v;
    v.reserve(n);

    unsigned p = 0;
    while (v.size() < n)
    {
        while (!is_prime(p))
        {
            ++p;
        }

        v.push_back(complex<double>(sqrt(p), cbrt(p)));

        ++p;
    }

    return v;
}

void slow_fourier_1d(const vector<complex<double>> & v, vector<complex<double>> & fv)
{
    const unsigned n = v.size();

    // Make sure the destination vector has the correct size.
    assert(fv.size() == n);

    for (unsigned fi = 0; fi < n; ++fi)
    {
        fv[fi] = 0;

        for (unsigned ti = 0; ti < n; ++ti)
        {
            fv[fi] += polar(1.0, -2.0 * M_PI * fi * ti / n) * v[ti];
        }
    }
}

void slow_inverse_fourier_1d(const vector<complex<double>> & fv, vector<complex<double>> & ifv)
{
    const unsigned n = fv.size();

    // Make sure the destination vector has the correct size.
    assert(ifv.size() == n);

    for (unsigned ti = 0; ti < n; ++ti)
    {
        ifv[ti] = 0;

        for (unsigned fi = 0; fi < n; ++fi)
        {
            ifv[ti] += polar(1.0, +2.0 * M_PI * fi * ti / n) * fv[fi] / static_cast<double>(n);
        }
    }
}

template <typename T>
ostream & print_vector(ostream & os, const vector<T> & v)
{
    os << "{";
    for (unsigned i = 0; i < v.size(); ++i)
    {
        if (i != 0)
        {
            os << ", ";
        }
        os << v[i];
    }
    os << "}";

    return os;
}

double test_naive(const unsigned n, const unsigned repeats)
{
    vector<complex<double>> v = mk_complex_test_vector(n);

    vector<complex<double>> fv(n);

    auto t1 = chrono::high_resolution_clock::now();

    for (unsigned rep = 0; rep < repeats; ++rep)
    {
        slow_fourier_1d(v, fv);
    }

    auto t2 = chrono::high_resolution_clock::now();

    double duration = chrono::duration_cast<chrono::nanoseconds>(t2 - t1).count() / 1e9;

    return duration / repeats;
}

double test_fftw(const unsigned n, const unsigned repeats)
{
    vector<complex<double>> v = mk_complex_test_vector(n);

    fftw_complex * x = fftw_alloc_complex(n);
    assert(x != nullptr);

    fftw_complex * y = fftw_alloc_complex(n);
    assert(y != nullptr);

    fftw_plan plan = fftw_plan_dft_1d(n, x, y, FFTW_FORWARD, FFTW_ESTIMATE | FFTW_PRESERVE_INPUT);
    assert(plan != nullptr);

    // Initialize the input

    // copy test data
    for (unsigned i = 0; i < n; ++i)
    {
        x[i][0] = real(v[i]);
        x[i][1] = imag(v[i]);
    }

    auto t1 = chrono::high_resolution_clock::now();

    for (unsigned rep = 0; rep < repeats; ++rep)
    {
        fftw_execute(plan);
    }

    auto t2 = chrono::high_resolution_clock::now();

    double duration = chrono::duration_cast<chrono::nanoseconds>(t2 - t1).count() / 1e9;

    fftw_destroy_plan(plan);

    fftw_free(y);
    fftw_free(x);

    fftw_cleanup();

    return duration / repeats;
}

void show()
{
    const unsigned n = 8;

    vector<complex<double>> v = mk_complex_test_vector(n);

    print_vector(cout, v); cout << endl;

    vector<complex<double>> fv(n);

    slow_fourier_1d(v, fv);

    print_vector(cout, fv); cout << endl;
}

int main()
{
    show();
    return 0;

    unsigned rep = 1000;

    for (unsigned n = 1; n < 10; ++n)
    {
        double duration = test_fftw(n, rep);

        cout << n << "\t" << rep << "\t" << duration << endl;

        rep = min(1000u, static_cast<unsigned>(0.25 / duration));
    }

    return 0;
}
